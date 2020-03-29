/***************************************************************************//**
\file    vector_types.h

\brief   Definitions of the basic short vector data structures
*******************************************************************************/

#ifndef _MMC_VECTOR_H
#define _MMC_VECTOR_H

#ifdef _MSC_VER
    #define PRE_ALIGN(x) __declspec(align(x))
    #define POST_ALIGN(x)
#else
    #define PRE_ALIGN(x)
    #define POST_ALIGN(x) __attribute__ ((aligned(x)))
#endif

/**
 \struct MMC_float4 vector_types.h
 \brief  floating-point quadraplet {x,y,z,w}

 the data structure is 16byte aligned to facilitate SSE operations
*/

typedef struct PRE_ALIGN(16) MMC_float4{
    float x,y,z,w;
} float4 POST_ALIGN(16);

/**
 \struct MMC_float3 vector_types.h
 \brief  floating-point triplet {x,y,z}

 if SSE is enabled, float3 is identical to float4
*/

#if defined(MMC_USE_SSE) || defined(USE_OPENCL)
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
typedef struct PRE_ALIGN(16) MMC_int4{
    int x,y,z,w;
} int4 POST_ALIGN(16);

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
typedef struct PRE_ALIGN(16) MMC_uint4{
    unsigned int x,y,z,w;
} uint4 POST_ALIGN(16);

/**
 \struct MMC_uint2 vector_types.h
 \brief  unsigned integer pair {ix,iy}
*/

typedef struct MMC_uint2{
    unsigned int x,y;
} uint2;


typedef unsigned int uint;

#endif
