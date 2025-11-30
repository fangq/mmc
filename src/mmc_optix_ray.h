#ifndef _MMC_OPTIX_RAY_H
#define _MMC_OPTIX_RAY_H

#include "mmc_optix_rand.cu"

/**
 * @brief struct for Ray information
 */
typedef struct __attribute__((aligned(16))) MMC_OptiX_Ray {
    float3 p0;                    /**< current photon position */
    float3 dir;                   /**< current photon direction vector */
    float slen;                   /**< the remaining unitless scattering length = length*mus */
    float weight;                 /**< photon current weight */
    float photontimer;            /**< the total time-of-fly of the photon */
    unsigned int mediumid;        /**< ID of current medium type */
} optixray;

/**
 * @brief Get photon position
 */
__device__ __forceinline__ float3 getPosition() {
    return make_float3(__uint_as_float(optixGetPayload_0()),
                       __uint_as_float(optixGetPayload_1()),
                       __uint_as_float(optixGetPayload_2()));
}

/**
 * @brief Set photon position
 */
__device__ __forceinline__ void setPosition(const float3& p) {
    optixSetPayload_0(__float_as_uint(p.x));
    optixSetPayload_1(__float_as_uint(p.y));
    optixSetPayload_2(__float_as_uint(p.z));
}

/**
 * @brief Get ray direction
 */
__device__ __forceinline__ float3 getDirection() {
    return make_float3(__uint_as_float(optixGetPayload_3()),
                       __uint_as_float(optixGetPayload_4()),
                       __uint_as_float(optixGetPayload_5()));
}

/**
 * @brief Set ray direction
 */
__device__ __forceinline__ void setDirection(const float3& v) {
    optixSetPayload_3(__float_as_uint(v.x));
    optixSetPayload_4(__float_as_uint(v.y));
    optixSetPayload_5(__float_as_uint(v.z));
}

/**
 * @brief Get remaining scattering length
 */
__device__ __forceinline__ float getSlen() {
    float slen = __uint_as_float(optixGetPayload_6());
    return slen;
}

/**
 * @brief Set remaining scattering length
 */
__device__ __forceinline__ void setSlen(const float& slen) {
    optixSetPayload_6(__float_as_uint(slen));
}

/**
 * @brief Get photon weight
 */
__device__ __forceinline__ float getWeight() {
    float w = __uint_as_float(optixGetPayload_7());
    return w;
}

/**
 * @brief Set photon weight
 */
__device__ __forceinline__ void setWeight(const float& w) {
    optixSetPayload_7(__float_as_uint(w));
}

/**
 * @brief Get time of flight for a photon
 */
__device__ __forceinline__ float getPhotonTimer() {
    float tof = __uint_as_float(optixGetPayload_8());
    return tof;
}

/**
 * @brief Update time of flight for a photon
 */
__device__ __forceinline__ void setPhotonTimer(const float& tof) {
    optixSetPayload_8(__float_as_uint(tof));
}

/**
 * @brief Get medium id
 */
__device__ __forceinline__ unsigned int getMediumID() {
    unsigned int id = optixGetPayload_9();
    return id;
}

/**
 * @brief Set medium id
 */
__device__ __forceinline__ void setMediumID(const unsigned int& id) {
    optixSetPayload_9(id);
}

/**
 * @brief Get RNG seed
 */
__device__ __forceinline__ uint4 getRNGSeed() {
    return make_uint4(optixGetPayload_10(),
                      optixGetPayload_11(),
                      optixGetPayload_12(),
                      optixGetPayload_13());
}

/**
 * @brief Set RNG seed
 */
__device__ __forceinline__ void setRNGSeed(const uint4& seed) {
    optixSetPayload_10(seed.y);
    optixSetPayload_11(seed.z);
    optixSetPayload_12(seed.w);
    optixSetPayload_13(seed.x);
}

/**
 * @brief Get ray info
 */
__device__ __forceinline__ optixray getRay() {
    optixray r;
    r.p0 = getPosition();
    r.dir = getDirection();
    r.slen = getSlen();
    r.weight = getWeight();
    r.photontimer = getPhotonTimer();
    r.mediumid = getMediumID();
    return r;
}

/**
 * @brief Set ray info
 */
__device__ __forceinline__ void setRay(const optixray& r) {
    setPosition(r.p0);
    setDirection(r.dir);
    setSlen(r.slen);
    setWeight(r.weight);
    setPhotonTimer(r.photontimer);
    setMediumID(r.mediumid);
}

/**
 * @brief Get RNG
 */
__device__ __forceinline__ mcx::Random getRNG() {
    return mcx::Random(getRNGSeed());
}

/**
 * @brief Set RNG
 */
__device__ __forceinline__ void setRNG(const mcx::Random& rng) {
    return setRNGSeed(rng.intSeed);
}

/**
 * @brief print ray information
 */
__device__ __forceinline__ void printRay(const optixray& r) {
    printf("pos:[%f %f %f], dir:[%f %f %f], slen:%f, weight:%f, tof:%fns, type:%u\n",
           r.p0.x, r.p0.y, r.p0.z,
           r.dir.x, r.dir.y, r.dir.z,
           r.slen, r.weight, r.photontimer * 1e9, r.mediumid);
}

#endif
