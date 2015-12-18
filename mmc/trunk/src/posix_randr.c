/*********************************************************************
* A POSIX Random Number Generator for multi-threaded applications    *
*                                                                    *
*  Author: Qianqian Fang <q.fang at neu.edu>                         *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    posix_randr.c

\brief   A POSIX Random Number Generator for multi-threaded applications
*******************************************************************************/

#ifndef _MMC_POSIX_RAND_H
#define _MMC_POSIX_RAND_H

#include <math.h>
#include <stdio.h>
#include "posix_randr.h"
#include "fastmath.h"

#ifdef MMC_USE_SSE_MATH
#include "sse_math/sse_math.h"
#include <smmintrin.h>
#endif

#define LOG_RNG_MAX          22.1807097779182f
#define INIT_MULT            1812433253

// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    double ran;
    drand48_r(t,&ran);
    return (float)ran;
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    unsigned short s[3], *si=(unsigned short *)n_seed;
    s[0]=si[0]; s[1]=si[2]; s[2]=(INIT_MULT * (n_seed[0] ^ (n_seed[0] >> 30)) + idx) & 0xFFFF;
    seed48_r(s, t);
}
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tbuf[RAND_BUF_LEN]){
}

#include "rng_common.h"

#endif
