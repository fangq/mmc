/*****************************************************************//**
*  A Pseudo-RNG using the SIMD-oriented Fast Mersenne Twister (SFMT) *
*                                                                    *
*  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>              *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    sfmt_rand.h

\brief   An interface to use the SFMT-19937 random number generator
*******************************************************************************/

#ifndef _MCEXTREME_SFMT_RAND_H
#define _MCEXTREME_SFMT_RAND_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define __device__ static inline

typedef unsigned int RandType;
typedef unsigned int uint;

#define MCX_RNG_NAME       "SFMT19937 RNG"

#define RAND_BUF_LEN       1248    //buffer length
#define RAND_SEED_LEN      2       //
#define MEXP               19937

__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx);

__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tbuf[RAND_BUF_LEN]);
// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]);
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]);
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]);
__device__ float rand_next_reflect(RandType t[RAND_BUF_LEN]);
__device__ float rand_do_roulette(RandType t[RAND_BUF_LEN]);

#ifdef MMC_USE_SSE_MATH
__device__ void rand_next_aangle_sincos(RandType t[RAND_BUF_LEN],float *si, float *co);
__device__ float rand_next_scatlen_ps(RandType t[RAND_BUF_LEN]);
#endif

#endif
