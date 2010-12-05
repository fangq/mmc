/***************************************************************************//**
\file    vector_types.h

\brief   Definitions of the basic short vector data structures
*******************************************************************************/

#ifndef _MMC_VECTOR_H
#define _MMC_VECTOR_H

typedef struct MMC_float4{
    float x,y,z,w;
} float4 __attribute__ ((aligned(16)));

#ifdef MMC_USE_SSE
 typedef struct MMC_float4 float3;
#else
 typedef struct MMC_float3{
    float x,y,z;
 } float3;
#endif

typedef struct MMC_int2{
    int x,y;
} int2;

typedef struct MMC_int3{
    int x,y,z;
} int3;

typedef struct MMC_int4{
    int x,y,z,w;
} int4;

typedef struct MMC_uint3{
    unsigned int x,y,z;
} uint3;

typedef struct MMC_uint2{
    unsigned int x,y;
} uint2;

#endif
