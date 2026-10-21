#ifndef _MMC_VECTOR_OPS_H
#define _MMC_VECTOR_OPS_H

/**
 * @brief arithmetic operators and dot products for the CUDA short-vector types (float3/float4)
 *
 * CUDA does not define arithmetic on its vector types. These are shared by the CUDA build of
 * mmc_core.cl and by the OptiX kernel (mmc_optix_core.cu), so neither needs sutil/vec_math.h
 * from the OptiX SDK samples.
 */

inline __device__ __host__ float3 operator *(float3 a, float3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __device__ __host__ float3 operator *(float f, float3 v) {
    return make_float3(v.x * f, v.y * f, v.z * f);
}
inline __device__ __host__ float3 operator *(float3 v, float f) {
    return make_float3(v.x * f, v.y * f, v.z * f);
}
inline __device__ __host__ void operator *=(float3& b, float f) {
    b.x *= f;
    b.y *= f;
    b.z *= f;
}
inline __device__ __host__ float3 operator +(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __device__ __host__ void operator +=(float3& b, float3 a) {
    b.x += a.x;
    b.y += a.y;
    b.z += a.z;
}
inline __device__ __host__ float3 operator -(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __device__ __host__ void operator -=(float3& b, float3 a) {
    b.x -= a.x;
    b.y -= a.y;
    b.z -= a.z;
}
inline __device__ __host__ float3 operator /(float3 v, float f) {
    float inv = 1.0f / f;
    return v * inv;
}
inline __device__ __host__ void operator /=(float3& b, float f) {
    float inv = 1.0f / f;
    b.x *= inv;
    b.y *= inv;
    b.z *= inv;
}
inline __device__ __host__ void operator *=(float3& b, float3 f) {
    b.x *= f.x;
    b.y *= f.y;
    b.z *= f.z;
}
inline __device__ __host__ float4 operator -(float4 a) {
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}
inline __device__ __host__ float4 operator *(float4 a, float4 b) {
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}
inline __device__ __host__ float4 operator *(float f, float4 v) {
    return make_float4(v.x * f, v.y * f, v.z * f, v.w * f);
}
inline __device__ __host__ float4 operator *(float4 v, float f) {
    return make_float4(v.x * f, v.y * f, v.z * f, v.w * f);
}
inline __device__ __host__ void operator *=(float4& b, float f) {
    b.x *= f;
    b.y *= f;
    b.z *= f;
    b.w *= f;
}
inline __device__ __host__ float4 operator +(float4 a, float4 b) {
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
inline __device__ __host__ void operator +=(float4& b, float4 a) {
    b.x += a.x;
    b.y += a.y;
    b.z += a.z;
    b.w += a.w;
}
inline __device__ __host__ float4 operator -(float4 a, float4 b) {
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
inline __device__ __host__ void operator -=(float4& b, float4 a) {
    b.x -= a.x;
    b.y -= a.y;
    b.z -= a.z;
    b.w -= a.w;
}
inline __device__ __host__ float3 operator /(float3 a, float3 b) {
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __device__ __host__ float4 operator /(float4 v, float f) {
    float inv = 1.0f / f;
    return v * inv;
}
inline __device__ __host__ void operator /=(float4& b, float f) {
    float inv = 1.0f / f;
    b.x *= inv;
    b.y *= inv;
    b.z *= inv;
    b.w *= inv;
}
#ifdef __USE_FAST_MATH__
inline __device__ float4 operator /(float4 a, float4 b) {
    return make_float4(__fdividef(a.x, b.x), __fdividef(a.y, b.y), __fdividef(a.z, b.z), __fdividef(a.w, b.w));
}
#else
inline __device__ __host__ float4 operator /(float4 a, float4 b) {
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}
#endif

inline __device__ __host__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __device__ __host__ float dot(float4 a, float4 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

#endif
