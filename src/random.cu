#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <stdint.h>
#include <optix.h>
#include <optix_device.h>
#include <optix_types.h>
#include <sutil/vec_math.h>
#include <limits>

namespace mcx {
union Random {
public:
    uint4 intSeed;
    uint64_t seed[2];

    __forceinline__ __device__ float xorshift128p_nextf();

    __forceinline__ __device__ Random();
    __forceinline__ __device__ Random(uint64_t a, uint64_t b);
    __forceinline__ __device__ Random(uint4 se);

    __forceinline__ __device__ float uniform(float a, float b);
    __forceinline__ __device__ float rand_next_scatlen();
};

__forceinline__ __device__ Random::Random() {
    seed[0] = 0;
    seed[1] = 0;
}

__forceinline__ __device__ Random::Random(uint64_t a, uint64_t b) {
    this->seed[0] = a;
    this->seed[1] = b;
}

__forceinline__ __device__ Random::Random(uint4 se) {
    this->intSeed = se;
}

__forceinline__ __device__ float Random::uniform(float a, float b) {
    return (b - a) * this->xorshift128p_nextf() + a;
}

__forceinline__ __device__ float Random::rand_next_scatlen() {
    return -logf(this->xorshift128p_nextf() + std::numeric_limits<float>::epsilon());
}

__forceinline__ __device__ float Random::xorshift128p_nextf() {
    union {
        uint64_t  i;
        float f[2];
        unsigned int  u[2];
    } s1;

    const uint64_t s0 = this->seed[1];
    s1.i = this->seed[0];
    this->seed[0] = s0;
    s1.i ^= s1.i << 23; // a
    this->seed[1] = s1.i ^ s0 ^ (s1.i >> 18) ^ (s0 >> 5); // b, c
    s1.i = this->seed[1] + s0;
    s1.u[0] = 0x3F800000U | (s1.u[0] >> 9);

    return s1.f[0] - 1.0f;
}
}