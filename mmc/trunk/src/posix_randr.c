/*********************************************************************
* A POSIX Random Number Generator for multi-threaded applications    *
*                                                                    *
*  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>              *
*                                                                    *
*********************************************************************/

#ifndef _MCEXTREME_POSIX_RAND_H
#define _MCEXTREME_POSIX_RAND_H

#include <math.h>
#include <stdio.h>
#include "posix_randr.h"

#define POSIX_R_RAND_MAX    (1.f/RAND_MAX)
#define LOG_MT_MAX          22.1807097779182f
#define INIT_MULT           1812433253u      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */

// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    return (uint)(rand_r(t))*POSIX_R_RAND_MAX;
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    *t=*n_seed + idx;
    //printf("seeding for %d: %d\n",idx,*t);
}
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tbuf[RAND_BUF_LEN]){
}
// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]){
    float ran=rand_uniform01(t);
    return ((ran==0.f)?LOG_MT_MAX:(-logf(ran)));
}
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
#endif
