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
#define R_MAX_SFMT_RAND      2.3283064359966e-10f //1/(2^32+1)
#define LOG_RNG_MAX          22.1807097779182f    //log(2^32)
#define INIT_MULT            1812433253

// generate random number for the next zenith angle
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tnew[RAND_BUF_LEN]){
}

__device__ void sfmt_init(RandType *t,RandType *tnew,uint n_seed[],uint idx){
     uint32_t ini[2];
     ini[0]=n_seed[0];
     ini[1]=(INIT_MULT * (n_seed[0] ^ (n_seed[0] >> 30)) + idx);
     init_by_array(ini, 2);
     #pragma omp critical
     {
        fill_array32(t,RAND_BUF_LEN);
     }
}
// transform into [0,1] random number
__device__ float rand_uniform01(RandType t[RAND_BUF_LEN]){
    static __thread unsigned int pos;
    if(pos>=RAND_BUF_LEN) {
        #pragma omp critical
	{
	    fill_array32(t,RAND_BUF_LEN);
	}
	pos=0;
    }
    return t[pos++]*R_MAX_SFMT_RAND;
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    sfmt_init(t,tnew,n_seed,idx);
}

#include "rng_common.h"

#endif
