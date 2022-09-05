#ifndef _MMC_OPTIX_RAY_H
#define _MMC_OPTIX_RAY_H

#include "random.cu"

/**
 * @brief struct for Ray information
 */
typedef struct __attribute__((aligned(16))) MMC_OptiX_Ray {
    mcx::Random random;
    float3 p0;                    /**< current photon position */
    float3 dir;                   /**< current photon direction vector */
    float slen;                   /**< the remaining unitless scattering length = length*mus */
    float weight;                 /**< photon current weight */
    float photontimer;            /**< the total time-of-fly of the photon */
    unsigned int mediumid;        /**< ID of current medium type */
} optixray;

#endif