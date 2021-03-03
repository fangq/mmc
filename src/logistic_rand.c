/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2021
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plucker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**  \li \c (\b Fang2019) Qianqian Fang and Shijie Yan, 
**          <a href="http://dx.doi.org/10.1117/1.JBO.24.11.115002">
**          "Graphics processing unit-accelerated mesh-based Monte Carlo photon transport 
**           simulations,"</a> J. of Biomedical Optics, 24(11), 115002 (2019)
**  \li \c (\b Yuan2021) Yaoshen Yuan, Shijie Yan, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-12-1-147">
**          "Light transport modeling in highly complex tissues using the implicit 
**           mesh-based Monte Carlo algorithm,"</a> Biomed. Optics Express, 12(1) 147-161 (2021)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    logistic_rand.c

\brief   A Random Number Generator based on coupled chaotic Logistic lattice
*******************************************************************************/

#include "logistic_rand.h"
#include "fastmath.h"

#define R_PI               0.318309886183791f
#define INIT_LOGISTIC      100
#define INIT_MULT          1812433253u      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */

#ifndef DOUBLE_PREC_LOGISTIC
  #define FUN(x)               (4.f*(x)*(1.f-(x)))
  #define NU 1e-8f
  #define NU2 (1.f-2.f*NU)
  #define MIN_INVERSE_LIMIT 1e-7f
  #define logistic_uniform(v)  (acosf(1.f-2.f*(v))*R_PI)
  #define R_MAX_C_RAND       (0.5f/RAND_MAX)  /*RAND_MAX is max signed int, need *2 for unsigned int*/
  #define LOG_MT_MAX         22.1807097779182f
#else
  #define FUN(x)               (4.0*(x)*(1.0-(x)))
  #define NU 1e-14
  #define NU2 (1.0-2.0*NU)
  #define MIN_INVERSE_LIMIT 1e-12
  #define logistic_uniform(v)  (acos(1.0-2.0*(v))*R_PI)
  #define R_MAX_C_RAND       (1./RAND_MAX)
  #define LOG_MT_MAX         22.1807097779182
#endif

#define RING_FUN(x,y,z)      (NU2*(x)+NU*((y)+(z)))


__device__ void logistic_step(RandType *t, RandType *tnew, int len_1){
/*
    int i;
    for(i=0;i<=len_1;i++)
       t[i]=FUN(t[i]);
    tnew[0]=RING_FUN(t[0],t[1],t[len_1]);
    for(i=1;i<len_1;i++)
       tnew[i]=RING_FUN(t[i],t[i-1],t[i+1]);
    tnew[len_1]=RING_FUN(t[len_1],t[0],t[len_1-1]);
*/
    RandType tmp;
    t[0]=FUN(t[0]);
    t[1]=FUN(t[1]);
    t[2]=FUN(t[2]);
    t[3]=FUN(t[3]);
    t[4]=FUN(t[4]);
    tnew[3]=RING_FUN(t[0],t[4],t[1]);   /* shuffle the results by separation of 2*/
    tnew[4]=RING_FUN(t[1],t[0],t[2]);
    tnew[0]=RING_FUN(t[2],t[1],t[3]);
    tnew[1]=RING_FUN(t[3],t[2],t[4]);
    tnew[2]=RING_FUN(t[4],t[3],t[0]);
    tmp =t[0];
    t[0]=t[2];
    t[2]=t[4];
    t[4]=t[1];
    t[1]=t[3];
    t[3]=tmp;
}
// generate random number for the next zenith angle
__device__ void rand_need_more(RandType t[RAND_BUF_LEN],RandType tbuf[RAND_BUF_LEN]){
    logistic_step(t,tbuf,RAND_BUF_LEN-1);
    logistic_step(tbuf,t,RAND_BUF_LEN-1);
}

__device__ void logistic_init(RandType *t,RandType *tnew,uint n_seed[],uint idx){
     int i;
     for(i=0;i<RAND_BUF_LEN;i++)
	t[i]=n_seed[idx*RAND_SEED_WORD_LEN+i]*R_MAX_C_RAND;

     for(i=0;i<INIT_LOGISTIC;i++)  /*initial randomization*/
           rand_need_more(t,tnew);
}
// transform into [0,1] random number
__device__ RandType rand_uniform01(RandType v){
    return logistic_uniform(v);
}
__device__ void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN],uint *n_seed,int idx){
    logistic_init(t,tnew,n_seed,idx);
}
// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]){
    RandType ran=rand_uniform01(t[0]);
    if(ran==0.f) ran=rand_uniform01(t[1]);
    return ((ran==0.f)?LOG_MT_MAX:(-logf(ran)));
}
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t[2]);
}
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t[4]);
}
// generate random number for the next zenith angle
__device__ float rand_next_reflect(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t[3]);
}
// generate random number for the next zenith angle
__device__ float rand_do_roulette(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t[1]);
}

