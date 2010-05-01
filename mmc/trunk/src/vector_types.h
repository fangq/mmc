#ifndef _OMP_CUDA_VECTOR_H
#define _OMP_CUDA_VECTOR_H

typedef struct GPU_float4{
    float x,y,z,w;
} float4 __attribute__ ((aligned(16)));

#ifdef MMC_USE_SSE
 typedef struct GPU_float4 float3;
#else
 typedef struct GPU_float3{
    float x,y,z;
 } float3;
#endif

typedef struct GPU_int2{
    int x,y;
} int2;

typedef struct GPU_int3{
    int x,y,z;
} int3;

typedef struct GPU_int4{
    int x,y,z,w;
} int4;

typedef struct GPU_uint3{
    unsigned int x,y,z;
} uint3;

typedef struct GPU_uint2{
    unsigned int x,y;
} uint2;

#endif
