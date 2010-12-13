/*********************************************************************
*A Random Number Generator based on coupled chaotic Logistic lattice *
*                                                                    *
*  (both double and single precision random numbers are supported)   *
*                                                                    *
*  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>              *
*                                                                    *
*  History: 2009/03/02  CUDA version based on Neal Wagner 1993       *
*         http://www.cs.utsa.edu/~wagner/pubs/logistic/logistic.pdf  *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    logistic_rand.c

\brief   A Random Number Generator based on coupled chaotic Logistic lattice
*******************************************************************************/

#include "sfmt_rand.h"
#include "SFMT/SFMT.h"
#include "fastmath.h"

#define MAX_SFMT_RAND        4294967296           //2^32
#define R_MAX_SFMT_RAND      2.3283064365387e-10f //1/2^32
#define LOG_SFMT_MAX         22.1807097779182f    //log(2^32)

// generate random number for the next zenith angle
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tnew[RAND_BUF_LEN]){
}

__device__ void sfmt_init(RandType *t,RandType *tnew,uint n_seed[],uint idx){
     uint32_t ini[2];
     ini[0]=n_seed[0];
     ini[1]=n_seed[1];
     init_by_array(ini, 2);
}
// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    static unsigned int pos;
    if(pos>=RAND_BUF_LEN) {
	fill_array32(t,RAND_BUF_LEN);
	pos=0;
    }
    return t[pos++]*R_MAX_SFMT_RAND;
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    sfmt_init(t,tnew,n_seed,idx);
}
// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]){
    float ran=rand_uniform01(t);
    return ((ran==0.f)?LOG_SFMT_MAX:(-logf(ran)));
}
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_next_reflect(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_do_roulette(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}

