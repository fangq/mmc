/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    posix_randr.c

\brief   A POSIX Random Number Generator for multi-threaded applications
*******************************************************************************/

#ifndef _MMC_POSIX_RAND_H
#define _MMC_POSIX_RAND_H

#include <math.h>
#include <stdio.h>
#include "xorshift128p_rand.h"
#include "fastmath.h"

#ifdef MMC_USE_SSE_MATH
#include "sse_math/sse_math.h"
#include <smmintrin.h>
#endif

#define LOG_RNG_MAX         22.1807097779182f
#define IEEE754_DOUBLE_BIAS     0x3FF0000000000000ul /* Added to exponent.  */

static float xorshift128p_nextf (RandType t[RAND_BUF_LEN]){
   union {
        ulong  i;
	float f[2];
	uint  u[2];
   } s1;
   const ulong s0 = t[1];
   s1.i = t[0];
   t[0] = s0;
   s1.i ^= s1.i << 23; // a
   t[1] = s1.i ^ s0 ^ (s1.i >> 18) ^ (s0 >> 5); // b, c
   s1.i = t[1] + s0;
   s1.u[0] = 0x3F800000U | (s1.u[0] >> 9);

   return s1.f[0] - 1.0f;
}

static void xorshift128p_seed (uint *seed,RandType t[RAND_BUF_LEN]){
    t[0] = (ulong)seed[0] << 32 | seed[1] ;
    t[1] = (ulong)seed[2] << 32 | seed[3];
}

// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    return xorshift128p_nextf(t);
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    xorshift128p_seed(n_seed+idx*RAND_SEED_WORD_LEN,t);
}
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tbuf[RAND_BUF_LEN]){
}

#include "rng_common.h"

#endif
