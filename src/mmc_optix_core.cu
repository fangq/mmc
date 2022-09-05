#include <iostream>
#include <limits>
#include <math.h>
#include <stdint.h>
#include <optix.h>
#include <optix_device.h>
#include <optix_types.h>
#include <sutil/vec_math.h>

#include "mmc_optix_ray.h"
#include "mmc_optix_launchparam.h"

constexpr float R_C0 = 3.335640951981520e-12f; //1/C0 in s/mm
constexpr float MAX_ACCUM = 1000.0f;

// simulation configuration and medium optical properties
extern "C" {
    __constant__ MMCParam gcfg;
}

/**
 * @brief Init RNG seed for each thread
 */
__device__ __forceinline__ void initRNGSeed(optixray &r, const int &idx) {
    r.random = mcx::Random(((uint4*)gcfg.seedbuffer)[idx]);
}

/**
 * @brief Launch a new photon
 */
__device__ __forceinline__ void launchPhoton(optixray &r) {
    r.p0 = gcfg.srcpos;
    r.dir = gcfg.srcdir;
    r.slen = r.random.rand_next_scatlen();
    r.weight = 1.0f;
    r.photontimer = 0.0f;
    r.mediumid = gcfg.mediumid0;
}

/**
 * @brief Move a photon one step forward
 */
__device__ __forceinline__ void movePhoton(optixray &r) {
    optixTrace(gcfg.gashandle, r.p0, r.dir, 0.00001f, std::numeric_limits<float>::max(),
        0.0f, OptixVisibilityMask(255), OptixRayFlags::OPTIX_RAY_FLAG_NONE, 0, 1, 0,
        *(uint32_t*)&(r.p0.x), *(uint32_t*)&(r.p0.y), *(uint32_t*)&(r.p0.z),
        *(uint32_t*)&(r.dir.x), *(uint32_t*)&(r.dir.y), *(uint32_t*)&(r.dir.z),
        *(uint32_t*)&(r.slen), *(uint32_t*)&(r.weight), *(uint32_t*)&(r.photontimer),
        r.random.intSeed.x, r.random.intSeed.y, r.random.intSeed.z, r.random.intSeed.w,
        r.mediumid);
}

/**
 * @brief Rotate a vector with given azimuth and zenith angles
 */
__device__ __forceinline__ float3 rotateVector(const float3 &vec, const float2 &zen, 
    const float2 &azi) {
    if (vec.z > -1.0f + std::numeric_limits<float>::epsilon() && 
        vec.z < 1.0f - std::numeric_limits<float>::epsilon()) {
        float tmp0 = 1.0f - vec.z * vec.z;
        float tmp1 = zen.x * rsqrtf(tmp0);
        return tmp1 * (azi.y * make_float3(vec.x, vec.y, -tmp0) * 
            make_float3(vec.z, vec.z, 1.f) + azi.x * make_float3(-vec.y, vec.x, 0.0f))
             + zen.y * vec;
    } else {
        return make_float3(zen.x * azi.y, zen.x * azi.x, (vec.z > 0.0f) ? zen.y : -zen.y);
    }
}

/**
 * @brief Returns the sine and cosine from the Henyey-Greenstein distribution
 */
__device__ __forceinline__ float2 henyeyGreenstein(const float &g, mcx::Random& rand) {
    float ctheta;

    if (fabs(g) > std::numeric_limits<float>::epsilon()) {
        ctheta = (1.0f - g * g) / (1.0f - g + 2.0f * g * rand.uniform(0.0f, 1.0f));
        ctheta *= ctheta;
        ctheta = (1.0f + g * g - ctheta) / (2.0f * g);
        ctheta = fmax(-1.0f, fmin(1.0f, ctheta));
    } else {
        ctheta = 2.0f * rand.uniform(0.0f, 1.0f) - 1.0f;
    }

    return make_float2(sinf(acosf(ctheta)), ctheta);
}

/**
 * @brief Update ray direction after a scattering event
 */
__device__ __forceinline__ float3 selectScatteringDirection(const float3 &dir, 
    const float &g, mcx::Random& rand) {
    float2 aziScat;
    sincosf(rand.uniform(0.0f, 2.0f * M_PIf), &aziScat.x, &aziScat.y);

    float2 zenScat = henyeyGreenstein(g, rand);

    return rotateVector(dir, zenScat, aziScat);
}

/**
 * @brief Accumualte output quantities to a 3D grid
 */
__device__ __forceinline__ void accumulateOutput(const float3 &p0, const float3 &dir, 
    const float &currweight, const float &lmove, const Medium &prop, const float &tof) {
    // update photon weight
    float nextweight = currweight * expf(-prop.mua *lmove);
    optixSetPayload_7(__float_as_uint(nextweight));

    // divide path into segments of equal length
    int segcount = ((int)(lmove * gcfg.dstep) + 1) << 1;
    float seglen = lmove / segcount;
    float segdecay = expf(-prop.mua * seglen);
    float segloss = currweight * (1.0f - segdecay);
    float3 step = seglen * dir;

    // deposit weight loss of each segment to the corresponding grid
    float3 segmid = p0 - gcfg.nmin + 0.5f * step; // segment midpoint
    float currtof = tof + seglen * R_C0 * prop.n; // current time of flight
    for (int i = 0; i < segcount; ++i) {
        // find the index of the grid to store the absorbed weight
        int3 idx = make_int3(__float2int_rd(segmid.x * gcfg.dstep),
                             __float2int_rd(segmid.y * gcfg.dstep),
                             __float2int_rd(segmid.z * gcfg.dstep));
        int tshift = min(((int)((currtof - gcfg.tstart) * gcfg.Rtstep)),
            gcfg.maxgate - 1) * gcfg.crop0.z;
        int idx1d = (idx.z * gcfg.crop0.y + idx.y * gcfg.crop0.x + idx.x) + tshift;

        // to minimize numerical error, use the same trick as MCX
        // becomes much slower when using atomicAdd(*double, double)
        float oldval = atomicAdd(&((float*)gcfg.outputbuffer)[idx1d], segloss);
        if (oldval > MAX_ACCUM) {
            if (atomicAdd(&((float*)gcfg.outputbuffer)[idx1d], -oldval) < 0.0f) {
                atomicAdd(&((float*)gcfg.outputbuffer)[idx1d], oldval);
            } else {
                atomicAdd(&((float*)gcfg.outputbuffer)[idx1d + gcfg.crop0.w], oldval);
            }
        }

        segloss *= segdecay;
        segmid += step;
        currtof += seglen * R_C0 * prop.n;
    }
}

/**
 * @brief Set photon position
 */
__device__ __forceinline__ void setPosition(const float3 &p) {
    optixSetPayload_0(__float_as_uint(p.x));
    optixSetPayload_1(__float_as_uint(p.y));
    optixSetPayload_2(__float_as_uint(p.z));
}

/**
 * @brief Set ray direction
 */
__device__ __forceinline__ void setDirection(const float3 &v) {
    optixSetPayload_3(__float_as_uint(v.x));
    optixSetPayload_4(__float_as_uint(v.y));
    optixSetPayload_5(__float_as_uint(v.z));
}

/**
 * @brief Update time of flight for a photon
 */
__device__ __forceinline__ void setPhotonTimer(const float &tof,
    const float &lmove, const Medium &prop) {
    optixSetPayload_8(__float_as_uint(tof + lmove * R_C0 * prop.n));
}

/**
 * @brief Update RNG seed
 */
__device__ __forceinline__ void setRNGSeed(const mcx::Random &random) {
    optixSetPayload_9(random.intSeed.x);
    optixSetPayload_10(random.intSeed.y);
    optixSetPayload_11(random.intSeed.z);
    optixSetPayload_12(random.intSeed.w);
}

/**
 * @brief Set medium id
 */
__device__ __forceinline__ void setMediumID(const unsigned int &id) {
    optixSetPayload_13(id);
}

/**
 * @brief print ray information
 */
__device__ __forceinline__ void printRay(const optixray &r) {
    printf("pos:[%f %f %f], dir:[%f %f %f], slen:%f, weight:%f, tof:%fns, type:%u\n",
        r.p0.x, r.p0.y, r.p0.z,
        r.dir.x, r.dir.y, r.dir.z,
        r.slen, r.weight, r.photontimer * 1e9, r.mediumid);
}

/**
 * @brief Launch photon and trace ray iteratively
 */
extern "C" __global__ void __raygen__rg() {
    uint3 launchindex = optixGetLaunchIndex();

    optixray r;
    initRNGSeed(r, launchindex.x);
    launchPhoton(r);

    int ndone = 0;  // number of simulated photons
    while (ndone < (gcfg.threadphoton + (launchindex.x < gcfg.oddphoton))) {
        movePhoton(r);

        // when a photon escapes or tof reaches the upper limit
        if (!(r.mediumid && r.photontimer < gcfg.tend)) {
            launchPhoton(r);
            ++ndone;
        }
    }
}

/**
 * @brief when a photon hits a triangle
 */
extern "C" __global__ void __closesthit__ch() {
    // get photon and ray information from payload
    const float3 p0 = make_float3(__uint_as_float(optixGetPayload_0()), 
        __uint_as_float(optixGetPayload_1()), __uint_as_float(optixGetPayload_2()));
    const float3 dir = make_float3(__uint_as_float(optixGetPayload_3()), 
        __uint_as_float(optixGetPayload_4()), __uint_as_float(optixGetPayload_5()));
    const float slen = __uint_as_float(optixGetPayload_6());
    const float weight = __uint_as_float(optixGetPayload_7());
    const float tof = __uint_as_float(optixGetPayload_8());
    const Medium currprop = gcfg.medium[optixGetPayload_13()];

    // distance to intersection
    const float hitlen = optixGetRayTmax();
    
    // determine path length and save output
    const float lmove = (slen > hitlen * currprop.mus) ? hitlen : slen / currprop.mus;
    accumulateOutput(p0, dir, weight, lmove, currprop, tof);

    // update remaining scattering length
    optixSetPayload_6(__float_as_uint(slen - lmove * currprop.mus));

    // next photon position
    float3 p1 = p0 + dir * lmove;

    // update photon direction if needed
    if (slen > hitlen * currprop.mus) {
        // after hitting a boundary
        const TriangleMeshSBTData &sbtData = 
            *(const TriangleMeshSBTData*)optixGetSbtDataPointer();
        const int primid = optixGetPrimitiveIndex();

        // get triangle information
        const uint4 index = sbtData.face[primid];
        const float3 &v0 = sbtData.node[index.x];
        const float3 &v1 = sbtData.node[index.y];
        const float3 &v2 = sbtData.node[index.z];

        // get intersection (barycentric coordinate)
        const float2 bary = optixGetTriangleBarycentrics();
        p1 = (1 - bary.x - bary.y) * v0 + bary.x * v1 + bary.y * v2;

        // update medium id (assume matched boundary)
        if (optixIsFrontFaceHit()) {
            setMediumID(index.w & 0xFF); // back medium
        } else {
            setMediumID(index.w >> 16);  // front medium
        }
        // todo: update ray direction at a mismatched boundary
    } else {
        // after a scattering event
        mcx::Random random = mcx::Random(make_uint4(optixGetPayload_9(), 
                                                    optixGetPayload_10(),
                                                    optixGetPayload_11(),
                                                    optixGetPayload_12()));

        // update direction and scattering length
        setDirection(selectScatteringDirection(dir, currprop.g, random));
        optixSetPayload_6(__float_as_uint(random.rand_next_scatlen()));

        // update RNG seed
        setRNGSeed(random);
    }

    // update photon timer
    setPhotonTimer(tof, lmove, currprop);

    // update photon position and ray direction
    setPosition(p1);
}

extern "C" __global__ void __miss__ms() {
    optixSetPayload_13(0);
}