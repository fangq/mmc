/*****************************************************************//**
*A Random Number Generator based on coupled chaotic Logistic lattice *
*                                                                    *
*  (both double and single precision random numbers are supported)   *
*                                                                    *
*  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>              *
*                                                                    *
*  History: 2009/03/02  CUDA version based on Neal Wagner 1993       *
*         http://www.cs.utsa.edu/~wagner/pubs/logistic/logistic.pdf  *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    logistic_rand.h

\brief   An interface to use a coupled chaotic Logistic lattice RNG
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

#define RAND_BUF_LEN       624     //buffer length
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
#endif
