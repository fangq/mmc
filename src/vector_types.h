/***************************************************************************//**
\file    vector_types.h

\brief   Definitions of the basic short vector data structures
*******************************************************************************/

#ifndef _MMC_VECTOR_H
#define _MMC_VECTOR_H

/**
 \struct MMC_float4 vector_types.h
 \brief  floating-point quadraplet {x,y,z,w}

 the data structure is 16byte aligned to facilitate SSE operations
*/

typedef struct MMC_float4{
    float x,y,z,w;
} float4 __attribute__ ((aligned(16)));

/**
 \struct MMC_float3 vector_types.h
 \brief  floating-point triplet {x,y,z}

 if SSE is enabled, float3 is identical to float4
*/

#ifdef MMC_USE_SSE
 typedef struct MMC_float4 float3;
#else
 typedef struct MMC_float3{
    float x,y,z;
 } float3;
#endif

/**
 \struct MMC_int2 vector_types.h
 \brief  integer pair {ix,iy}
*/

typedef struct MMC_int2{
    int x,y;
} int2;

/**
 \struct MMC_int3 vector_types.h
 \brief  integer triplet {ix,iy,iz}
*/

typedef struct MMC_int3{
    int x,y,z;
} int3;

/**
 \struct MMC_int4 vector_types.h
 \brief  unsigned integer quadraplet {ix,iy,iz,iw}
*/
typedef struct MMC_int4{
    int x,y,z,w;
} int4 __attribute__ ((aligned(16)));

/**
 \struct MMC_uint3 vector_types.h
 \brief  unsigned integer triplet {ix,iy,iz}
*/
typedef struct MMC_uint3{
    unsigned int x,y,z;
} uint3;

/**
 \struct MMC_uint3 vector_types.h
 \brief  unsigned integer triplet {ix,iy,iz}
*/
typedef struct MMC_uint4{
    unsigned int x,y,z,w;
} uint4 __attribute__ ((aligned(16)));

/**
 \struct MMC_uint2 vector_types.h
 \brief  unsigned integer pair {ix,iy}
*/

typedef struct MMC_uint2{
    unsigned int x,y;
} uint2;


typedef unsigned int uint;

#endif
