/*****************************************************************//**
*  A Pseudo-RNG using the SIMD-oriented Fast Mersenne Twister (SFMT) *
*                                                                    *
*  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>              *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    logistic_rand.c

\brief   A Random Number Generator based on coupled chaotic Logistic lattice
*******************************************************************************/

#ifndef _MMC_SFMT_RAND_H
#define _MMC_SFMT_RAND_H

#include "sfmt_rand.h"
#include "SFMT/SFMT.h"
#include "fastmath.h"

#ifdef MMC_USE_SSE_MATH
#include "sse_math/sse_math.h"
#include <smmintrin.h>
#endif

#define MAX_SFMT_RAND        4294967296           //2^32
#define R_MAX_SFMT_RAND      2.3283064365387e-10f //1/2^32
#define LOG_RNG_MAX          22.1807097779182f    //log(2^32)

// generate random number for the next zenith angle
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tnew[RAND_BUF_LEN]){
}

__device__ void sfmt_init(RandType *t,RandType *tnew,uint n_seed[],uint idx){
     uint32_t ini[2];
     ini[0]=n_seed[0];
     ini[1]=n_seed[1];
     init_by_array(ini, 2);
     fill_array32(t,RAND_BUF_LEN);
}
// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    static __thread unsigned int pos;
    if(pos>=RAND_BUF_LEN) {
	fill_array32(t,RAND_BUF_LEN);
	pos=0;
    }
    return t[pos++]*R_MAX_SFMT_RAND;
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    sfmt_init(t,tnew,n_seed,idx);
}

#include "rng_common.h"

#endif
