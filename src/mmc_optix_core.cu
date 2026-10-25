#include <iostream>
#include <limits>
#include <math.h>

/* M_PIf is new in glibc 2.35 (the CentOS 7 release builder has 2.17), and MSVC
   defines neither it nor M_PI without _USE_MATH_DEFINES */
#ifndef M_PIf
    #define M_PIf 3.14159265358979323846f
#endif
#include <stdint.h>
#include <optix.h>
#include <optix_device.h>
#include <optix_types.h>
#include "mmc_vector_ops.h"

#include "mmc_optix_ray.h"
#include "mmc_optix_launchparam.h"

enum TOutputType {otFlux, otFluence, otEnergy, otJacobian, otWL, otWP,
                  otRF, otRFmus,
                  otAdjoint, otAdjointDcoeff, otAdjointMus, otAdjointMusp, otAdjointMuaD, otAdjointMuaMusp
                 };

constexpr float R_C0 = 3.335640951981520e-12f; // 1/C0 in s/mm
constexpr float MAX_ACCUM = 1000.0f;
constexpr float SAFETY_DISTANCE = 0.0001f; // heuristic to ensure ray cut through triangle
constexpr float DOUBLE_SAFETY_DISTANCE = SAFETY_DISTANCE * 2.0f;
constexpr float TWO_PI = 6.28318530717959f;       //2*pi
#ifndef FLT_EPSILON   // a macro of <float.h>, which MSVC's headers pull in
constexpr float FLT_EPSILON = 1.19209290E-07F;    //round-off limit
#endif
constexpr float EPS = FLT_EPSILON;                //round-off limit
// simulation configuration and medium optical properties
extern "C" {
    __constant__ MMCParam gcfg;
}

/**
 * @brief Cross product of two vectors
 */
__device__ __forceinline__ float3 crossVec(const float3& a, const float3& b) {
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/**
 * @brief Normalize a vector
 */
__device__ __forceinline__ float3 unitVec(const float3& v) {
    return v * rsqrtf(dot(v, v));
}

__device__ __forceinline__ float3 rotateVector(const float3& vec, const float2& zen, const float2& azi);

/**
 * @brief Init RNG seed for each thread
 */
__device__ __forceinline__ void initRNGSeed(mcx::Random& rng, const int& idx) {
    rng = mcx::Random(((uint4*)gcfg.seedbuffer)[idx]);
}

/**
 * @brief Determine the initial medium ID where the photon is launched
 */
__device__ __forceinline__ void getInitialMediumId(optixray& r, mcx::Random& rng) {
    optixTrace(gcfg.gashandle, r.p0, r.dir, 0.0f, std::numeric_limits<float>::max(),
               0.0f, OptixVisibilityMask(255), OptixRayFlags::OPTIX_RAY_FLAG_NONE, 0, 1, 0,
               *(uint32_t*) & (r.p0.x), *(uint32_t*) & (r.p0.y), *(uint32_t*) & (r.p0.z),
               *(uint32_t*) & (r.dir.x), *(uint32_t*) & (r.dir.y), *(uint32_t*) & (r.dir.z),
               *(uint32_t*) & (r.slen), *(uint32_t*) & (r.weight), *(uint32_t*) & (r.photontimer),
               r.mediumid,
               rng.intSeed.x, rng.intSeed.y, rng.intSeed.z, rng.intSeed.w,
               *(uint32_t*) & (r.weight_im));
}

/**
 * @brief Launch a new photon
 */
__device__ __forceinline__ void launchPhoton(optixray& r, mcx::Random& rng) {

    r.p0 = gcfg.srcpos;
    r.dir = gcfg.srcdir;
    r.weight = 1.0f;
    r.weight_im = 0.0f;  /**< initialize imaginary weight to zero */

    /* multi-source mode, as in the OpenCL/CUDA kernel: srcid < 0 picks a source slot at
       random for each photon, whose deposits go to the field of that slot; srcid > 0
       launches all photons from slot srcid-1 into a single field */
    if (gcfg.extrasrclen > 0 && (gcfg.srcid < 0 || (gcfg.srcid > 0 && gcfg.srcid <= gcfg.extrasrclen))) {
        int srcslot;

        if (gcfg.srcid < 0) {
            srcslot = min((int)(rng.uniform(0.0f, 1.0f) * gcfg.extrasrclen), gcfg.extrasrclen - 1);
            ((uint*)gcfg.slotbuffer)[optixGetLaunchIndex().x] = srcslot;
        } else {
            srcslot = gcfg.srcid - 1;
        }

        ExtraSrc src = ((ExtraSrc*)gcfg.srcdatabuffer)[srcslot];
        r.p0 = make_float3(src.srcpos.x, src.srcpos.y, src.srcpos.z);
        r.dir = make_float3(src.srcdir.x, src.srcdir.y, src.srcdir.z);
        r.weight = src.srcpos.w;  /**< importance weight from srcpos.w */

        /* disk source sampling for adjoint detector-as-source */
        float radius = src.srcparam1.x;

        if (radius > 0.f) {
            float sphi, cphi;
            float phi = TWO_PI * rng.uniform(0.0f, 1.0f);
            sincosf(phi, &sphi, &cphi);
            float r0 = sqrtf(rng.uniform(0.0f, 1.0f)) * radius;
            float3 d = r.dir;

            if (d.z > -1.f + EPS && d.z < 1.f - EPS) {
                float tmp0 = 1.f - d.z * d.z;
                float tmp1 = r0 * rsqrt(tmp0);
                r.p0.x += tmp1 * (d.x * d.z * cphi - d.y * sphi);
                r.p0.y += tmp1 * (d.y * d.z * cphi + d.x * sphi);
                r.p0.z -= tmp1 * tmp0 * cphi;
            } else {
                r.p0.x += r0 * cphi;
                r.p0.y += r0 * sphi;
            }
        }

        /* focal point focusing: srcdir.w stores focal length */
        if (src.srcdir.w != 0.f) {
            float3 focal = make_float3(r.p0.x + src.srcdir.w * r.dir.x,
                                       r.p0.y + src.srcdir.w * r.dir.y,
                                       r.p0.z + src.srcdir.w * r.dir.z);
            float3 newdir = make_float3(focal.x - r.p0.x,
                                        focal.y - r.p0.y,
                                        focal.z - r.p0.z);
            float len = rsqrtf(dot(newdir, newdir));
            r.dir = newdir * len;
        }

        r.slen = rng.rand_next_scatlen();
        r.photontimer = 0.0f;

        /* the host stores the medium enclosing each source in srcparam2.w, < 0 if unknown */
        r.mediumid = (src.srcparam2.w >= 0.f) ? (unsigned int)src.srcparam2.w : INITIAL_MEDIUM_UNKNOWN;

        if (r.mediumid == INITIAL_MEDIUM_UNKNOWN) {
            getInitialMediumId(r, rng);
        }

        return;
    }

    float3 origin = r.p0;  /* the focal reference point of a wide-field source */
    bool canfocus = true;
    uint posidx = 0;       /* the pattern pixel of the photon, used when srcnum > 1 */

    if (gcfg.srctype == MCX_SRC_PLANAR || gcfg.srctype == MCX_SRC_PATTERN || gcfg.srctype == MCX_SRC_FOURIER
            || gcfg.srctype == MCX_SRC_PENCILARRAY || gcfg.srctype == MCX_SRC_PATTERN3D) {
        float rx = rng.uniform(0.0f, 1.0f);
        float ry = rng.uniform(0.0f, 1.0f);

        if (gcfg.srctype == MCX_SRC_PATTERN3D) {
            float rz = rng.uniform(0.0f, 1.0f);
            int nx = (int)gcfg.srcparam1.x, ny = (int)gcfg.srcparam1.y, nz = (int)gcfg.srcparam1.z;
            r.p0 = gcfg.srcpos + make_float3(rx * gcfg.srcparam1.x, ry * gcfg.srcparam1.y, rz * gcfg.srcparam1.z);
            posidx = min((int)(rz * nz), nz - 1) * ny * nx + min((int)(ry * ny), ny - 1) * nx + min((int)(rx * nx), nx - 1);
            r.weight = (gcfg.srcnum > 1) ? 1.f : ((float*)gcfg.srcpattern)[posidx];
            origin += 0.5f * make_float3(gcfg.srcparam1.x, gcfg.srcparam1.y, gcfg.srcparam1.z);
        } else {
            float3 v1 = make_float3(gcfg.srcparam1.x, gcfg.srcparam1.y, gcfg.srcparam1.z);
            float3 v2 = make_float3(gcfg.srcparam2.x, gcfg.srcparam2.y, gcfg.srcparam2.z);

            if (gcfg.srctype == MCX_SRC_PENCILARRAY) {
                r.p0 = gcfg.srcpos + floorf(rx * gcfg.srcparam1.w) / (gcfg.srcparam1.w - 1.f) * v1
                       + floorf(ry * gcfg.srcparam2.w) / (gcfg.srcparam2.w - 1.f) * v2;
            } else {
                r.p0 = gcfg.srcpos + rx * v1 + ry * v2;
            }

            if (gcfg.srctype == MCX_SRC_PATTERN) {
                int xsize = (int)gcfg.srcparam1.w, ysize = (int)gcfg.srcparam2.w;
                posidx = min((int)(ry * ysize), ysize - 1) * xsize + min((int)(rx * xsize), xsize - 1);
                r.weight = (gcfg.srcnum > 1) ? 1.f : ((float*)gcfg.srcpattern)[posidx];
            } else if (gcfg.srctype == MCX_SRC_FOURIER) {
                r.weight = (cosf((floorf(gcfg.srcparam1.w) * rx + floorf(gcfg.srcparam2.w) * ry + gcfg.srcparam1.w - floorf(gcfg.srcparam1.w)) * TWO_PI)
                            * (1.f - gcfg.srcparam2.w + floorf(gcfg.srcparam2.w)) + 1.f) * 0.5f;
            }

            origin += 0.5f * (v1 + v2);
        }
    } else if (gcfg.srctype == MCX_SRC_FOURIERX || gcfg.srctype == MCX_SRC_FOURIERX2D) {
        /* srcparam1: [v1, |v2|], v2 = |v2| * unit(srcdir x v1); srcparam2: [kx, ky, phi0, M] or [kx, ky, phix, phiy] */
        float rx = rng.uniform(0.0f, 1.0f);
        float ry = rng.uniform(0.0f, 1.0f);
        float3 v1 = make_float3(gcfg.srcparam1.x, gcfg.srcparam1.y, gcfg.srcparam1.z);
        float3 v2 = gcfg.srcparam1.w * rsqrtf(dot(v1, v1)) * crossVec(gcfg.srcdir, v1);
        r.p0 = gcfg.srcpos + rx * v1 + ry * v2;

        if (gcfg.srctype == MCX_SRC_FOURIERX2D) {
            r.weight = (sinf((gcfg.srcparam2.x * rx + gcfg.srcparam2.z) * TWO_PI) * sinf((gcfg.srcparam2.y * ry + gcfg.srcparam2.w) * TWO_PI) + 1.f) * 0.5f;
        } else {
            r.weight = (cosf((gcfg.srcparam2.x * rx + gcfg.srcparam2.y * ry + gcfg.srcparam2.z) * TWO_PI) * (1.f - gcfg.srcparam2.w) + 1.f) * 0.5f;
        }

        origin += 0.5f * (v1 + v2);
    } else if (gcfg.srctype == MCX_SRC_DISK || gcfg.srctype == MCX_SRC_RING || gcfg.srctype == MCX_SRC_GAUSSIAN) {
        float sphi, cphi, phi;

        if (gcfg.srctype == MCX_SRC_RING && (gcfg.srcparam1.z > 0.f || gcfg.srcparam1.w > 0.f)) {
            phi = fabsf(gcfg.srcparam1.z - gcfg.srcparam1.w) * rng.uniform(0.0f, 1.0f) + fminf(gcfg.srcparam1.z, gcfg.srcparam1.w);
        } else {
            phi = TWO_PI * rng.uniform(0.0f, 1.0f);
        }

        sincosf(phi, &sphi, &cphi);

        float r0;

        if (gcfg.srctype != MCX_SRC_GAUSSIAN) {
            /* a disk of radius srcparam1.x, or a ring between the radii srcparam1.x and srcparam1.y */
            r0 = sqrtf(rng.uniform(0.0f, 1.0f) * fabsf(gcfg.srcparam1.x * gcfg.srcparam1.x - gcfg.srcparam1.y * gcfg.srcparam1.y)
                       + gcfg.srcparam1.y * gcfg.srcparam1.y);
        } else if (fabsf(gcfg.focus) < 1e-5f || fabsf(gcfg.srcparam1.y) < 1e-5f) {
            r0 = sqrtf(-logf((rng.uniform(0.0f, 1.0f)))) * gcfg.srcparam1.x;
        } else {
            float z0 = gcfg.srcparam1.x * gcfg.srcparam1.x * M_PIf / gcfg.srcparam1.y; //Rayleigh range
            r0 = sqrtf(-logf(rng.uniform(0.0f, 1.0f)) * (1.f + (gcfg.focus * gcfg.focus / (z0 * z0)))) * gcfg.srcparam1.x;
        }

        if (gcfg.srcdir.z > -1.f + EPS && gcfg.srcdir.z < 1.f - EPS) {
            float tmp0 = 1.f - gcfg.srcdir.z * gcfg.srcdir.z;
            float tmp1 = r0 * rsqrt(tmp0);
            r.p0.x = gcfg.srcpos.x + tmp1 * (gcfg.srcdir.x * gcfg.srcdir.z * cphi - gcfg.srcdir.y * sphi);
            r.p0.y = gcfg.srcpos.y + tmp1 * (gcfg.srcdir.y * gcfg.srcdir.z * cphi + gcfg.srcdir.x * sphi);
            r.p0.z = gcfg.srcpos.z - tmp1 * tmp0 * cphi;
        } else {
            r.p0.x += r0 * cphi;
            r.p0.y += r0 * sphi;
        }
    }  else if (gcfg.srctype == MCX_SRC_CONE || gcfg.srctype == MCX_SRC_ISOTROPIC || gcfg.srctype == MCX_SRC_ARCSINE) {
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rng.uniform(0.0f, 1.0f); //next arimuth angle
        sincosf(ang, &sphi, &cphi);

        if (gcfg.srctype == MCX_SRC_CONE) {
            do {
                ang = (gcfg.srcparam1.y > 0) ? TWO_PI * rng.uniform(0.0f, 1.0f) : acos(2.f * rng.uniform(0.0f, 1.0f) - 1.f); //sine distribution
            } while (ang > gcfg.srcparam1.x);
        } else {
            if (gcfg.srctype == MCX_SRC_ISOTROPIC) { // uniform sphere
                ang = acosf(2.f * rng.uniform(0.0f, 1.0f) - 1.f);    //sine distribution
            } else {
                ang = M_PIf * rng.uniform(0.0f, 1.0f);    //uniform distribution in zenith angle, arcsine (single precision, as the rest of the kernel)
            }
        }

        sincosf(ang, &stheta, &ctheta);
        r.dir.x = stheta * cphi;
        r.dir.y = stheta * sphi;
        r.dir.z = ctheta;
        canfocus = false;
    }  else if (gcfg.srctype == MCX_SRC_ZGAUSSIAN) {
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rng.uniform(0.0f, 1.0f); //next arimuth angle
        sincosf(ang, &sphi, &cphi);
        ang = sqrtf(-2.f * logf(rng.uniform(0.0f, 1.0f))) * (1.f - 2.f * rng.uniform(0.0f, 1.0f)) * gcfg.srcparam1.x;
        sincosf(ang, &stheta, &ctheta);
        r.dir.x = stheta * cphi;
        r.dir.y = stheta * sphi;
        r.dir.z = ctheta;
        canfocus = false;
    } else if (gcfg.srctype == MCX_SRC_HYPERBOLOID) {
        /* a hyperboloid Gaussian beam: srcparam1.x is the waist radius, srcparam1.y the
           distance from srcpos to the waist and srcparam1.z the Rayleigh range */
        float sphi, cphi;
        sincosf(TWO_PI * rng.uniform(0.0f, 1.0f), &sphi, &cphi);
        float r0 = sqrtf(0.5f * rng.rand_next_scatlen()) * gcfg.srcparam1.x;

        /* position offset and direction in the beam frame */
        float tilt = (gcfg.srcparam1.z > 0.f) ? (-gcfg.srcparam1.y / gcfg.srcparam1.z) : 0.f;
        float3 bdir = unitVec(make_float3(-r0 * sphi, r0 * cphi, gcfg.srcparam1.z));
        float dp_x = r0 * (cphi - tilt * sphi);
        float dp_y = r0 * (sphi + tilt * cphi);

        /* rotate them to the srcdir frame */
        if (gcfg.srcdir.z > -1.f + EPS && gcfg.srcdir.z < 1.f - EPS) {
            float sin_theta = sqrtf(1.f - gcfg.srcdir.z * gcfg.srcdir.z);
            float cphi_beam = gcfg.srcdir.x / sin_theta;
            float sphi_beam = gcfg.srcdir.y / sin_theta;

            r.p0.x = gcfg.srcpos.x + dp_x * cphi_beam * gcfg.srcdir.z - dp_y * sphi_beam;
            r.p0.y = gcfg.srcpos.y + dp_x * sphi_beam * gcfg.srcdir.z + dp_y * cphi_beam;
            r.p0.z = gcfg.srcpos.z - dp_x * sin_theta;

            r.dir.x = bdir.x * cphi_beam * gcfg.srcdir.z - bdir.y * sphi_beam + bdir.z * cphi_beam * sin_theta;
            r.dir.y = bdir.x * sphi_beam * gcfg.srcdir.z + bdir.y * cphi_beam + bdir.z * sphi_beam * sin_theta;
            r.dir.z = -bdir.x * sin_theta + bdir.z * gcfg.srcdir.z;
        } else {
            r.p0 = gcfg.srcpos + make_float3(dp_x, dp_y, 0.f);
            r.dir = make_float3(bdir.x, bdir.y, (gcfg.srcdir.z > 0.f) ? bdir.z : -bdir.z);
        }

        r.dir = unitVec(r.dir);
        canfocus = false;
    } else if (gcfg.srctype == MCX_SRC_LINE || gcfg.srctype == MCX_SRC_SLIT) {
        float t = rng.uniform(0.0f, 1.0f);
        r.p0.x += t * gcfg.srcparam1.x;
        r.p0.y += t * gcfg.srcparam1.y;
        r.p0.z += t * gcfg.srcparam1.z;

        if (gcfg.srctype == MCX_SRC_LINE) {
            float s, p;
            t = 1.f - 2.f * rng.uniform(0.0f, 1.0f);
            s = 1.f - 2.f * rng.uniform(0.0f, 1.0f);
            p = sqrtf(1.f - r.dir.x * r.dir.x - r.dir.y * r.dir.y) * (rng.uniform(0.0f, 1.0f) > 0.5f ? 1.f : -1.f);
            float3 vv;
            vv.x = r.dir.y * p - r.dir.z * s;
            vv.y = r.dir.z * t - r.dir.x * p;
            vv.z = r.dir.x * s - r.dir.y * t;
            r.dir = vv;
            //*((float3*)&(r.dir))=(float3)(r.dir.y*p-r.dir.z*s,r.vec.z*t-r.vec.x*p,r.vec.x*s-r.vec.y*t);
        }

        origin += 0.5f * make_float3(gcfg.srcparam1.x, gcfg.srcparam1.y, gcfg.srcparam1.z);
        canfocus = (gcfg.srctype == MCX_SRC_SLIT);
    }

    /* the incident angles of a focused, isotropic or Lambertian beam, as the other backends do */
    if (canfocus && gcfg.focus != 0.f) {
        if (isnan(gcfg.focus)) {
            float ang, stheta, ctheta, sphi, cphi;
            sincosf(TWO_PI * rng.uniform(0.0f, 1.0f), &sphi, &cphi);
            ang = acosf(2.f * rng.uniform(0.0f, 1.0f) - 1.f); //sine distribution
            sincosf(ang, &stheta, &ctheta);
            r.dir = rotateVector(r.dir, make_float2(stheta, ctheta), make_float2(sphi, cphi));
        } else if (isinf(gcfg.focus) && gcfg.focus < 0.f) { // Lambertian (cosine distribution)
            float stheta, sphi, cphi;
            sincosf(TWO_PI * rng.uniform(0.0f, 1.0f), &sphi, &cphi);
            stheta = sqrtf(rng.uniform(0.0f, 1.0f));
            r.dir = rotateVector(r.dir, make_float2(stheta, sqrtf(1.f - stheta * stheta)), make_float2(sphi, cphi));
        } else {
            origin += gcfg.focus * r.dir;
            r.dir = unitVec((gcfg.focus < 0.f) ? (r.p0 - origin) : (origin - r.p0)); // diverging or converging beam
        }
    }

    if (gcfg.srcnum > 1) {
        ((uint*)gcfg.slotbuffer)[optixGetLaunchIndex().x] = posidx;
    }

    r.slen = rng.rand_next_scatlen();
    r.photontimer = 0.0f;
    r.mediumid = gcfg.mediumid0;

    // if the initial medium ID is unknown, determine it at runtime
    if (r.mediumid == INITIAL_MEDIUM_UNKNOWN) {
        getInitialMediumId(r, rng);
    }
}

/**
 * @brief Move a photon one step forward
 */
__device__ __forceinline__ void movePhoton(optixray& r, mcx::Random& rng) {
    optixTrace(gcfg.gashandle, r.p0, r.dir, 0.0f,
               gcfg.medium[r.mediumid].mus ? r.slen / gcfg.medium[r.mediumid].mus : std::numeric_limits<float>::max(),
               0.0f, OptixVisibilityMask(255), OptixRayFlags::OPTIX_RAY_FLAG_NONE, 0, 1, 0,
               *(uint32_t*) & (r.p0.x), *(uint32_t*) & (r.p0.y), *(uint32_t*) & (r.p0.z),
               *(uint32_t*) & (r.dir.x), *(uint32_t*) & (r.dir.y), *(uint32_t*) & (r.dir.z),
               *(uint32_t*) & (r.slen), *(uint32_t*) & (r.weight), *(uint32_t*) & (r.photontimer),
               r.mediumid,
               rng.intSeed.x, rng.intSeed.y, rng.intSeed.z, rng.intSeed.w,
               *(uint32_t*) & (r.weight_im));
}

/**
 * @brief Rotate a vector with given azimuth and zenith angles
 */
__device__ __forceinline__ float3 rotateVector(const float3& vec, const float2& zen,
        const float2& azi) {
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
__device__ __forceinline__ float2 henyeyGreenstein(const float& g, mcx::Random& rand) {
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
__device__ __forceinline__ float3 selectScatteringDirection(const float3& dir,
        const float& g, mcx::Random& rand) {
    float2 aziScat;
    sincosf(rand.uniform(0.0f, 2.0f * M_PIf), &aziScat.x, &aziScat.y);

    float2 zenScat = henyeyGreenstein(g, rand);

    return rotateVector(dir, zenScat, aziScat);
}

/**
 * @brief Convert 3D indices to 1D
 */
__device__ __forceinline__ uint subToInd(const uint3& idx3d) {
    return idx3d.z * gcfg.crop0.y + idx3d.y * gcfg.crop0.x + idx3d.x;
}

/**
 * @brief Get the index of the voxel that encloses p
 */
__device__ __forceinline__ uint getVoxelIdx(const float3& p) {
    return subToInd(make_uint3(p.x > 0.0f ? __float2int_rd(min(p.x, gcfg.nmax.x) * gcfg.dstep) : 0,
                               p.y > 0.0f ? __float2int_rd(min(p.y, gcfg.nmax.y) * gcfg.dstep) : 0,
                               p.z > 0.0f ? __float2int_rd(min(p.z, gcfg.nmax.z) * gcfg.dstep) : 0));
}

/**
 * @brief Get the offset of the current time frame
 */
__device__ __forceinline__ uint getTimeFrame(const float& tof) {
    return min(((int)((tof - gcfg.tstart) * gcfg.Rtstep)),
               gcfg.maxgate - 1) * gcfg.crop0.z;
}

/**
 * @brief Save output to a buffer (real part or CW)
 */
__device__ __forceinline__ void saveToBuffer(const uint& eid, const float& w) {
    // to minimize numerical error, use the same trick as MCX
    // becomes much slower when using atomicAdd(*double, double)
    float accum = atomicAdd(&((float*)gcfg.outputbuffer)[eid], w);

    if (accum > MAX_ACCUM) {
        if (atomicAdd(&((float*)gcfg.outputbuffer)[eid], -accum) < 0.0f) {
            atomicAdd(&((float*)gcfg.outputbuffer)[eid], accum);
        } else {
            atomicAdd(&((float*)gcfg.outputbuffer)[eid + gcfg.crop0.w], accum);
        }
    }
}

/**
 * @brief Save imaginary component to a buffer (RF forward simulation)
 * Imaginary buffer is located at outputbuffer + 2*crop0.w
 */
__device__ __forceinline__ void saveToBufferIm(const uint& eid, const float& w) {
    uint imeid = eid + 2 * gcfg.crop0.w;
    float accum = atomicAdd(&((float*)gcfg.outputbuffer)[imeid], w);

    if (accum > MAX_ACCUM) {
        if (atomicAdd(&((float*)gcfg.outputbuffer)[imeid], -accum) < 0.0f) {
            atomicAdd(&((float*)gcfg.outputbuffer)[imeid], accum);
        } else {
            atomicAdd(&((float*)gcfg.outputbuffer)[imeid + gcfg.crop0.w], accum);
        }
    }
}

/**
 * @brief The pattern pixel of the current photon of a pattern source shared by srcnum
 * patterns, or the source slot of the photon in the multi-source mode
 */
__device__ __forceinline__ uint threadSlot() {
    return ((uint*)gcfg.slotbuffer)[optixGetLaunchIndex().x];
}

/**
 * @brief Save a deposit to the output; a pattern source with srcnum > 1 patterns saves one
 * deposit per pattern, scaled by the pattern intensity at the launch pixel of the photon,
 * in the layout of the other backends: [(gate*nvox+voxel)*srcnum+pattern]
 */
__device__ __forceinline__ void saveOutput(uint eid, const float& w, const float& w_im, const bool& isrf) {
    if (gcfg.srcnum > 1) {
        const float* pat = (float*)gcfg.srcpattern + threadSlot() * gcfg.srcnum;
        eid *= gcfg.srcnum;

        for (int k = 0; k < gcfg.srcnum; k++) {
            saveToBuffer(eid + k, w * pat[k]);

            if (isrf) {
                saveToBufferIm(eid + k, w_im * pat[k]);
            }
        }

        return;
    }

    /* multi-source mode: the field of each source slot follows that of the previous slot */
    if (gcfg.srcid < 0 && gcfg.extrasrclen > 0) {
        eid += threadSlot() * gcfg.crop0.z * gcfg.maxgate;
    }

    saveToBuffer(eid, w);

    if (isrf) {
        saveToBufferIm(eid, w_im);
    }
}

/**
 * @brief The per-thread detected photon scratch: nscat per medium, ppath per medium and the
 * initial weight (2*maxmedia+1 floats), the same as the ppath array of the OpenCL/CUDA kernels
 */
__device__ __forceinline__ float* threadPath() {
    return (float*)gcfg.ppathbuffer + optixGetLaunchIndex().x * (2 * gcfg.maxmedia + 1);
}

/**
 * @brief Clear the detected photon scratch of a newly launched photon
 */
__device__ __forceinline__ void clearPath(const float& w0) {
    float* ppath = threadPath();

    for (int i = 0; i < 2 * gcfg.maxmedia; i++) {
        ppath[i] = 0.f;
    }

    ppath[2 * gcfg.maxmedia] = w0;
}

/**
 * @brief Save an exiting photon if it leaves through a detector, in the record layout of the
 * other backends: detector id, nscat and ppath per medium, exit position and direction
 * (if issaveexit), initial weight
 */
__device__ __forceinline__ void saveDetPhoton(const optixray& r) {
    for (int i = 0; i < gcfg.detnum; i++) {
        float4 det = ((float4*)gcfg.detposbuffer)[i];
        float dx = det.x - r.p0.x, dy = det.y - r.p0.y, dz = det.z - r.p0.z;

        if (dx * dx + dy * dy + dz * dz < det.w * det.w) {
            unsigned int idx = atomicAdd((unsigned int*)gcfg.detbuffer, 1u);

            if (idx < (unsigned int)gcfg.maxdetphoton) {
                const int reclen = 2 * gcfg.maxmedia + 2 + gcfg.issaveexit * 6;
                float* rec = (float*)gcfg.detbuffer + 1 + (size_t)idx * reclen;
                float* ppath = threadPath();

                rec[0] = (float)(i + 1);

                for (int j = 0; j < 2 * gcfg.maxmedia; j++) {
                    rec[1 + j] = ppath[j];
                }

                if (gcfg.issaveexit) {
                    rec[2 * gcfg.maxmedia + 1] = r.p0.x;
                    rec[2 * gcfg.maxmedia + 2] = r.p0.y;
                    rec[2 * gcfg.maxmedia + 3] = r.p0.z;
                    rec[2 * gcfg.maxmedia + 4] = r.dir.x;
                    rec[2 * gcfg.maxmedia + 5] = r.dir.y;
                    rec[2 * gcfg.maxmedia + 6] = r.dir.z;
                }

                rec[reclen - 1] = ppath[2 * gcfg.maxmedia];
            }

            return;
        }
    }
}

/**
 * @brief Save a trajectory point (-D M) in the record layout of the other backends: photon
 * index, position, weight and, in place of the element index, the medium label
 */
__device__ __forceinline__ void saveTraj(const optixray& r, const uint& id) {
    uint pos = atomicAdd((unsigned int*)gcfg.trajbuffer, 1u);

    if (pos < (uint)gcfg.maxjumpdebug) {
        float* rec = (float*)gcfg.trajbuffer + 1 + (size_t)pos * 6;
        rec[0] = __uint_as_float(id);
        rec[1] = r.p0.x;
        rec[2] = r.p0.y;
        rec[3] = r.p0.z;
        rec[4] = r.weight;
        rec[5] = __uint_as_float((r.mediumid >= DEAD_MEDIUM_ID) ? 0u : r.mediumid);  /* 0 once the photon has ended */
    }
}

/**
 * @brief Accumulate output quantities to a 3D grid
 */
__device__ __forceinline__ void accumulateOutput(optixray& r, const Medium& prop,
        const float& lmove) {
    bool isrf = (gcfg.omega != 0.f);

    // divide path into segments of equal length
    int segcount = ((int)(lmove * gcfg.dstep) + 1) << 1;
    float seglen = lmove / segcount;
    float segdecay = expf(-prop.mua * seglen);

    // RF: per-segment complex rotation angle = omega * n/c0 * seglen
    float rfcos = 1.0f, rfsin = 0.0f;

    if (isrf) {
        float phase = gcfg.omega * prop.n * gcfg.oneoverc0 * seglen;
        sincosf(phase, &rfsin, &rfcos);
    }

    float segloss, segloss_im = 0.0f;
    float w_re = r.weight, w_im = r.weight_im;

    if (isrf) {
        // complex Beer-Lambert per segment: w_new = w * segdecay * exp(-i*phase), the sign
        // that the deposit (w0-w)/(mua + i*omega*n/c0) below and the OpenCL/CUDA kernels use
        float a_im = gcfg.omega * prop.n * gcfg.oneoverc0;
        float a_mag2 = prop.mua * prop.mua + a_im * a_im;
        // new weight after one segment
        float new_re = segdecay * (w_re * rfcos + w_im * rfsin);
        float new_im = segdecay * (w_im * rfcos - w_re * rfsin);
        float dw_re = w_re - new_re;
        float dw_im = w_im - new_im;
        segloss    = (a_mag2 > 0.f) ? (dw_re * prop.mua + dw_im * a_im) / a_mag2 : 0.f;
        segloss_im = (a_mag2 > 0.f) ? (dw_im * prop.mua - dw_re * a_im) / a_mag2 : 0.f;
        w_re = new_re;
        w_im = new_im;
    } else {
        segloss = (gcfg.outputtype == otEnergy) ? r.weight * (1.0f - segdecay) :
                  (prop.mua < EPS ? r.weight * seglen :
                   __fdividef(r.weight * (1.0f - segdecay), prop.mua));
    }

    // deposit weight loss of each segment to the corresponding grid
    float3 step = seglen * r.dir;
    float3 segmid = r.p0 - gcfg.nmin + 0.5f * step; // segment midpoint
    float currtof = r.photontimer + seglen * R_C0 * prop.n; // current time of flight

    // find the information of the first segment
    uint oldeid = getTimeFrame(currtof) + getVoxelIdx(segmid);

    /* adjoint J_mus/J_musp: the medium label of each voxel the photon deposits in, for the
       chain-rule factor applied per voxel by the adjoint post-processing */
    if (gcfg.savevoxtype) {
        ((int*)gcfg.voxtypebuffer)[getVoxelIdx(segmid)] = r.mediumid;
    }

    float oldweight = segloss;
    float oldweight_im = segloss_im;

    // iterate over the rest of the segments
    for (int i = 1; i < segcount; ++i) {
        if (isrf) {
            float new_re2 = segdecay * (w_re * rfcos + w_im * rfsin);
            float new_im2 = segdecay * (w_im * rfcos - w_re * rfsin);
            float dw_re2 = w_re - new_re2;
            float dw_im2 = w_im - new_im2;
            float a_im = gcfg.omega * prop.n * gcfg.oneoverc0;
            float a_mag2 = prop.mua * prop.mua + a_im * a_im;
            segloss    = (a_mag2 > 0.f) ? (dw_re2 * prop.mua + dw_im2 * a_im) / a_mag2 : 0.f;
            segloss_im = (a_mag2 > 0.f) ? (dw_im2 * prop.mua - dw_re2 * a_im) / a_mag2 : 0.f;
            w_re = new_re2;
            w_im = new_im2;
        } else {
            segloss *= segdecay;
        }

        segmid += step;
        currtof += seglen * R_C0 * prop.n;
        uint neweid = getTimeFrame(currtof) + getVoxelIdx(segmid);

        // save when entering a new element or during the last segment
        if (neweid != oldeid) {
            saveOutput(oldeid, oldweight, oldweight_im, isrf);

            if (gcfg.savevoxtype) {
                ((int*)gcfg.voxtypebuffer)[getVoxelIdx(segmid)] = r.mediumid;
            }

            // reset oldeid and weight bucket
            oldeid = neweid;
            oldweight = 0.0f;
            oldweight_im = 0.0f;
        }

        oldweight += segloss;
        oldweight_im += segloss_im;
    }

    // save the weight loss of the last segment
    saveOutput(oldeid, oldweight, oldweight_im, isrf);

    // update both real and imaginary weight in ray after all segments
    if (isrf) {
        r.weight    = w_re;
        r.weight_im = w_im;
    }
}

/**
 * @brief reflect or refract a ray at mismatched boundary
 */
__device__ __forceinline__ bool reflectray(const float3& norm, const float& n1,
        const float& n2, mcx::Random& rng, optixray& r) {
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2;

    Icos = fabs(dot(r.dir, norm));
    tmp0 = n1 * n1;
    tmp1 = n2 * n2;
    tmp2 = 1.f - tmp0 / tmp1 * (1.f - Icos * Icos);

    if (tmp2 > 0.0f) {
        // if no total internal reflection
        Re = tmp0 * Icos * Icos + tmp1 * tmp2; /*transmission angle*/
        tmp2 = sqrtf(tmp2); /*to save one sqrt*/
        Im = 2.f * n1 * n2 * Icos * tmp2;
        Rtotal = (Re - Im) / (Re + Im); /*Rp*/
        Re = tmp1 * Icos * Icos + tmp0 * tmp2 * tmp2;
        Rtotal = (Rtotal + (Re - Im) / (Re + Im)) * 0.5f; /*(Rp+Rs)/2*/

        if (rng.uniform(0.0f, 1.0f) <= Rtotal) {
            // do reflection, correct photon position
            r.p0 -= r.dir * DOUBLE_SAFETY_DISTANCE;
            r.dir += -2.0f * Icos * norm;
        } else {
            // do transmission
            r.dir += -Icos * norm;
            r.dir = tmp2 * norm + n1 / n2 * r.dir;
            return false;
        }
    } else {
        // total internel reflection, correct photon position
        r.p0 -= r.dir * DOUBLE_SAFETY_DISTANCE;
        r.dir += -2.0f * Icos * norm;
    }

    // normalize new direction
    tmp0 = rsqrtf(dot(r.dir, r.dir));
    r.dir *= tmp0;
    return true;
}

/**
 * @brief Add a photon weight to a float accumulator pair with the MCX dual single-precision
 * buffer: once the primary sum exceeds MAX_ACCUM, it is moved to the secondary one
 */
__device__ __forceinline__ void accumEnergy(float& acc, float& acc2, const float& w) {
    acc += w;

    if (acc > MAX_ACCUM) {
        acc2 += acc;
        acc = 0.f;
    }
}

/**
 * @brief Add a thread's energy sums to the 4 floats (per pattern) after the output grids (2
 * or, for RF, 4 blocks of crop0.w): [tot, esc] and their secondary sums at [2, 3], using the
 * same dual single-precision buffer as saveToBuffer
 */
__device__ __forceinline__ void saveEnergy(const int& idx, const float& acc, const float& acc2) {
    float* buf = (float*)gcfg.outputbuffer + gcfg.crop0.w * ((gcfg.omega != 0.f) ? 4 : 2);
    float accum = atomicAdd(buf + idx, acc);

    if (accum > MAX_ACCUM) {
        if (atomicAdd(buf + idx, -accum) < 0.0f) {
            atomicAdd(buf + idx, accum);
        } else {
            atomicAdd(buf + idx + 2, accum);
        }
    }

    atomicAdd(buf + idx + 2, acc2);
}

/**
 * @brief Add the weight of a launched (isesc=0) or escaped (isesc=1) photon of a pattern
 * source shared by srcnum patterns to the per-thread energy sums of each pattern, kept in
 * global memory as [tot, esc, tot2, esc2] per pattern
 */
__device__ __forceinline__ void accumPatternEnergy(const int& isesc, const float& w) {
    float* e = (float*)gcfg.energybuffer + (size_t)optixGetLaunchIndex().x * 4 * gcfg.srcnum;
    const float* pat = (float*)gcfg.srcpattern + threadSlot() * gcfg.srcnum;

    for (int k = 0; k < gcfg.srcnum; k++) {
        accumEnergy(e[k * 4 + isesc], e[k * 4 + isesc + 2], w * pat[k]);
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

    // number of photons of this thread; a thread without photons must not launch one,
    // which would also count its weight in the launched energy
    const int nphoton = gcfg.threadphoton + (launchindex.x < gcfg.oddphoton);

    if (nphoton == 0) {
        return;
    }

    optixray r;
    int ndone = 0;  // number of simulated photons

    /* launched and escaped energy of this thread, as the other backends tally them: the
       weight a photon still carries when it ends (escape or time limit) has escaped;
       |w| for RF, which the phase rotation leaves unchanged */
    float energytot = 0.f, energytot2 = 0.f, energyesc = 0.f, energyesc2 = 0.f;

    /* photons are launched at one place in the loop, so the source code is inlined once */
    r.mediumid = DEAD_MEDIUM_ID;
    r.photontimer = 0.f;
    bool islaunch = true;

    while (true) {
        if (islaunch) {
            launchPhoton(r, rng);
            islaunch = false;

            if (gcfg.savetraj) {
                /* the photon index of the OpenCL/CUDA kernels: consecutive in each thread */
                uint id = launchindex.x * gcfg.threadphoton + min(launchindex.x, (uint)gcfg.oddphoton) + ndone;
                ((uint*)gcfg.photonidbuffer)[launchindex.x] = id;
                saveTraj(r, id);
            }

            if (gcfg.srcnum > 1) {
                accumPatternEnergy(0, r.weight);
            } else {
                accumEnergy(energytot, energytot2, r.weight);
            }

            if (gcfg.issavedet) {
                clearPath((gcfg.srcnum > 1) ? __uint_as_float(threadSlot()) : r.weight);  /* the pattern pixel, as the other backends save it */
            }
        }

        if (r.mediumid != DEAD_MEDIUM_ID) {
            movePhoton(r, rng);
        }

        // when a photon escapes or tof reaches the upper limit
        if (!(r.mediumid != DEAD_MEDIUM_ID && r.photontimer < gcfg.tend)) {
            float wesc = (gcfg.omega != 0.f) ? sqrtf(r.weight * r.weight + r.weight_im * r.weight_im) : r.weight;

            if (gcfg.savetraj) {
                saveTraj(r, ((uint*)gcfg.photonidbuffer)[launchindex.x]);
            }

            if (gcfg.srcnum > 1) {
                accumPatternEnergy(1, wesc);
            } else {
                accumEnergy(energyesc, energyesc2, wesc);
            }

            if (++ndone >= nphoton) {
                break;
            }

            islaunch = true;
        }
    }

    if (gcfg.srcnum > 1) {
        const float* e = (float*)gcfg.energybuffer + (size_t)launchindex.x * 4 * gcfg.srcnum;

        for (int k = 0; k < gcfg.srcnum; k++) {
            saveEnergy(k * 4, e[k * 4], e[k * 4 + 2]);
            saveEnergy(k * 4 + 1, e[k * 4 + 1], e[k * 4 + 3]);
        }
    } else {
        saveEnergy(0, energytot, energytot2);
        saveEnergy(1, energyesc, energyesc2);
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

    // distance to the intersection
    const float hitlen = optixGetRayTmax();

    // intersected triangle id
    const int primid = optixGetPrimitiveIndex();

    // hit front or back face
    const bool isfronthit = optixIsFrontFaceHit();

    // get info of triangle, including face normal, back and front media
    const TriangleMeshSBTData& sbtData =
        *(const TriangleMeshSBTData*)optixGetSbtDataPointer();
    float4 fnorm = sbtData.fnorm[primid];
    const uint media = __float_as_uint(fnorm.w);

    // immediately return once the initial medium ID unknown -> known
    if (r.mediumid == INITIAL_MEDIUM_UNKNOWN) {
        r.mediumid = isfronthit ? (media >> 16) : (media & 0xFFFF);
        setRay(r);
        return;
    }

    // current medium id
    r.mediumid = isfronthit ? (media >> 16) : (media & 0xFFFF);

    // get medium properties
    const Medium currprop = gcfg.medium[r.mediumid];

    // resolve some edge cases where hitlen is sligthly larger than slen
    const float lmove = min(hitlen, currprop.mus ? r.slen / currprop.mus : std::numeric_limits<float>::max());

    // save output; for RF mode, also updates r.weight and r.weight_im inside.
    // nothing is deposited in the ambient medium outside the mesh, e.g. the gap a
    // wide-field source photon crosses from its launch plane to the tissue surface
    if (r.mediumid != AMBIENT_MEDIUM_ID) {
        accumulateOutput(r, currprop, lmove);

        if (gcfg.issavedet) {
            threadPath()[gcfg.maxmedia + r.mediumid - 1] += lmove;  /* the partial path per medium */
        }
    }

    // update photon position
    r.p0 += r.dir * (lmove + SAFETY_DISTANCE);

    // for CW mode, update photon weight (RF weights already updated inside accumulateOutput)
    if (gcfg.omega == 0.f) {
        r.weight *= expf(-currprop.mua * lmove);
    }

    // update photon timer
    r.photontimer += lmove * R_C0 * currprop.n;

    // after hitting a boundary, update remaining scattering length
    r.slen -= lmove * currprop.mus;

    // update medium id (assume matched boundary)
    r.mediumid = isfronthit ? (media & 0xFFFF) : (media >> 16);

    // update ray direction at mismatched boundary
    if (gcfg.isreflect && currprop.n != gcfg.medium[r.mediumid].n) {
        fnorm = isfronthit ? -fnorm : fnorm;

        if (reflectray(*(float3*)&fnorm, currprop.n, gcfg.medium[r.mediumid].n,
                       rng, r)) {
            // if the ray is reflected
            r.mediumid = isfronthit ? (media >> 16) : (media & 0xFFFF);
        }

    }

    // terminate once a ray leaves the tissue for the ambient medium, also at a
    // matched boundary or without -b 1, the same as the other tracers do on exit
    if (r.mediumid == AMBIENT_MEDIUM_ID) {
        if (gcfg.issavedet) {
            saveDetPhoton(r);
        }

        r.mediumid = DEAD_MEDIUM_ID;
    }

    // update rng
    setRNG(rng);

    // update ray
    setRay(r);
}

extern "C" __global__ void __miss__ms() {
    // get photon and ray information from payload
    optixray r = getRay();

    // terminate if a ray's initial medium is unknown and
    // misses the entire tissue volume
    if (r.mediumid == INITIAL_MEDIUM_UNKNOWN) {
        r.mediumid = DEAD_MEDIUM_ID;
        setRay(r);
        return;
    }

    // get rng
    mcx::Random rng = getRNG();

    // get medium properties
    const Medium currprop = gcfg.medium[r.mediumid];

    // a ray that hits no triangle in a non-scattering medium (the void outside the
    // mesh, mus=0) has left the domain for good: slen/mus is inf, and depositing
    // that segment wrote NaN into the 8 corners of the output grid
    if (currprop.mus == 0.f) {
        r.mediumid = DEAD_MEDIUM_ID;
        setRay(r);
        return;
    }

    // distance to the scattering event site
    float lmove = r.slen / currprop.mus;

    // save output; for RF mode, also updates r.weight and r.weight_im inside
    accumulateOutput(r, currprop, lmove);

    // the partial path and the scattering count of this medium
    if (gcfg.issavedet) {
        float* ppath = threadPath();
        ppath[gcfg.maxmedia + r.mediumid - 1] += lmove;
        ppath[r.mediumid - 1] += 1.f;
    }

    // update photon position
    r.p0 += r.dir * lmove;

    // for CW mode, update photon weight (RF weights already updated inside accumulateOutput)
    if (gcfg.omega == 0.f) {
        r.weight *= expf(-currprop.mua * lmove);
    }

    // update photon timer
    r.photontimer += lmove * R_C0 * currprop.n;

    // scattering event
    r.dir = selectScatteringDirection(r.dir, currprop.g, rng);
    r.slen = rng.rand_next_scatlen();

    if (gcfg.savetraj) {
        saveTraj(r, ((uint*)gcfg.photonidbuffer)[optixGetLaunchIndex().x]);
    }

    // update rng
    setRNG(rng);

    // update ray
    setRay(r);
}
