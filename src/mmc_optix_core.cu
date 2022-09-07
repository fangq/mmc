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
__device__ __forceinline__ void initRNGSeed(mcx::Random &rng, const int &idx) {
    rng = mcx::Random(((uint4*)gcfg.seedbuffer)[idx]);
}

/**
 * @brief Launch a new photon
 */
__device__ __forceinline__ void launchPhoton(optixray &r, mcx::Random &rng) {
    r.p0 = gcfg.srcpos;
    r.dir = gcfg.srcdir;
    r.slen = rng.rand_next_scatlen();
    r.weight = 1.0f;
    r.photontimer = 0.0f;
    r.mediumid = gcfg.mediumid0;
}

/**
 * @brief Move a photon one step forward
 */
__device__ __forceinline__ void movePhoton(optixray &r, mcx::Random &rng) {
    optixTrace(gcfg.gashandle, r.p0, r.dir, 0.0001f, std::numeric_limits<float>::max(),
        0.0f, OptixVisibilityMask(255), OptixRayFlags::OPTIX_RAY_FLAG_NONE, 0, 1, 0,
        *(uint32_t*)&(r.p0.x), *(uint32_t*)&(r.p0.y), *(uint32_t*)&(r.p0.z),
        *(uint32_t*)&(r.dir.x), *(uint32_t*)&(r.dir.y), *(uint32_t*)&(r.dir.z),
        *(uint32_t*)&(r.slen), *(uint32_t*)&(r.weight), *(uint32_t*)&(r.photontimer),
        r.mediumid,
        rng.intSeed.x, rng.intSeed.y, rng.intSeed.z, rng.intSeed.w);
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
__device__ __forceinline__ void accumulateOutput(const optixray &r, const Medium &prop,
    const float &lmove) {
    // divide path into segments of equal length
    int segcount = ((int)(lmove * gcfg.dstep) + 1) << 1;
    float seglen = lmove / segcount;
    float segdecay = expf(-prop.mua * seglen);
    float segloss = r.weight * (1.0f - segdecay);

    // deposit weight loss of each segment to the corresponding grid
    float3 step = seglen * r.dir;
    float3 segmid = r.p0 - gcfg.nmin + 0.5f * step; // segment midpoint
    float currtof = r.photontimer + seglen * R_C0 * prop.n; // current time of flight
    for (int i = 0; i < segcount; ++i) {
        // find the index of the grid to store the absorbed weight
        int3 idx = make_int3(segmid.x > 0.0f ? __float2int_rd(segmid.x * gcfg.dstep) : 0,
                             segmid.y > 0.0f ? __float2int_rd(segmid.y * gcfg.dstep) : 0,
                             segmid.z > 0.0f ? __float2int_rd(segmid.z * gcfg.dstep) : 0);
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
 * @brief Launch photon and trace ray iteratively
 */
extern "C" __global__ void __raygen__rg() {
    uint3 launchindex = optixGetLaunchIndex();

    // init RNG seed for each thread
    mcx::Random rng;
    initRNGSeed(rng, launchindex.x);

    // init a ray
    optixray r;
    launchPhoton(r, rng);

    int ndone = 0;  // number of simulated photons
    while (ndone < (gcfg.threadphoton + (launchindex.x < gcfg.oddphoton))) {
        movePhoton(r, rng);

        // when a photon escapes or tof reaches the upper limit
        if (!(r.mediumid && r.photontimer < gcfg.tend)) {
            launchPhoton(r, rng);
            ++ndone;
        }
    }
}

/**
 * @brief when a photon hits a triangle
 */
extern "C" __global__ void __closesthit__ch() {
    // get photon and ray information from payload
    optixray r = getRay();

    // get rng
    mcx::Random rng = getRNG();

    // get medium properties
    const Medium currprop = gcfg.medium[r.mediumid];

    // distance to intersection
    const float hitlen = optixGetRayTmax();
    
    // determine path length
    const float lmove = (r.slen > hitlen * currprop.mus) ?
        hitlen : r.slen / currprop.mus;

    // save output
    accumulateOutput(r, currprop, lmove);

    // update photon position
    r.p0 += r.dir * lmove;

    // update photon weight
    r.weight *= expf(-currprop.mua * lmove);

    // update photon timer
    r.photontimer += lmove * R_C0 * currprop.n;

    // update photon direction if needed
    if (r.slen > hitlen * currprop.mus) {
        // after hitting a boundary, update remaining scattering length
        r.slen -= lmove * currprop.mus;

        // get triangle information
        const TriangleMeshSBTData &sbtData = 
            *(const TriangleMeshSBTData*)optixGetSbtDataPointer();
        const int primid = optixGetPrimitiveIndex();
        const uint4 index = sbtData.face[primid];
        const float3 &v0 = sbtData.node[index.x];
        const float3 &v1 = sbtData.node[index.y];
        const float3 &v2 = sbtData.node[index.z];

        // get intersection (barycentric coordinate)
        const float2 bary = optixGetTriangleBarycentrics();
        r.p0 = (1.0f - bary.x - bary.y) * v0 + bary.x * v1 + bary.y * v2;

        // update medium id (assume matched boundary)
        if (optixIsFrontFaceHit()) {
            r.mediumid = (index.w & 0xFF); // back medium
        } else {
            r.mediumid = (index.w >> 16);  // front medium
        }

        // todo: update ray direction at a mismatched boundary
    } else {
        // after a scattering event, new direction and scattering length
        r.dir = selectScatteringDirection(r.dir, currprop.g, rng);
        r.slen = rng.rand_next_scatlen();
    }

    // update rng
    setRNG(rng);

    // update ray
    setRay(r);
}

extern "C" __global__ void __miss__ms() {
    // concave case needs further investigation
    setMediumID(0);
}