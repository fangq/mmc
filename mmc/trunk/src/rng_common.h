#ifndef _MMC_RNG_COMMON_H
#define _MMC_RNG_COMMON_H

#define TWO_PI     (M_PI*2.0)

// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen(RandType t[RAND_BUF_LEN]){
    float ran=rand_uniform01(t);
    return ((ran==0.f)?LOG_RNG_MAX:(-logf(ran)));
}
// generate [0,1] random number for the next arimuthal angle
__device__ float rand_next_aangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_next_zangle(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for reflection test
__device__ float rand_next_reflect(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
// generate random number for the next zenith angle
__device__ float rand_do_roulette(RandType t[RAND_BUF_LEN]){
    return rand_uniform01(t);
}
#ifdef MMC_USE_SSE_MATH

#define MATH_BLOCK 8

// generate [0,1] random number for the sin/cos of arimuthal angles
__device__ void rand_next_aangle_sincos(RandType t[RAND_BUF_LEN],float *si, float *co){
    static __thread V4SF sine[MATH_BLOCK], cosine[MATH_BLOCK];
    static __thread int pos=MATH_BLOCK+1;
    if(pos>=(MATH_BLOCK<<2)){
    	V4SF ran[MATH_BLOCK];
	int i,j;
	for(i=0;i<MATH_BLOCK;i++)
          for(j=0;j<4;j++)
	    ran[i].f[j]=TWO_PI*rand_uniform01(t);
        for(i=0;i<MATH_BLOCK;i++)
	    sincos_ps(ran[i].v,&(sine[i].v),&(cosine[i].v));
    	pos=0;
    }
    *si=sine[0].f[pos];
    *co=cosine[0].f[pos++];
}
// generate [0,1] random number for the next scattering length
__device__ float rand_next_scatlen_ps(RandType t[RAND_BUF_LEN]){
    static __thread V4SF logval[MATH_BLOCK];
    static __thread int pos=MATH_BLOCK+1;
    float res;
    if(pos>=(MATH_BLOCK<<2)){
    	V4SF ran[MATH_BLOCK];
	int i,j;	
	for(i=0;i<MATH_BLOCK;i++)
	  for(j=0;j<4;j++)
	    ran[i].f[j]=rand_uniform01(t);
        for(i=0;i<MATH_BLOCK;i++){
	    logval[i].v=log_ps(ran[i].v);
	}
    	pos=0;
    }
    res=((logval[0].f[pos]!=logval[0].f[pos])?LOG_RNG_MAX:(-logval[0].f[pos]));
    pos++;
    return res;
}
#endif

#endif
