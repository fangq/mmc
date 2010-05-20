#ifndef MCX_FEM_TRACING_H
#define MCX_FEM_TRACING_H

#include "simpmesh.h"
#include "mcx_utils.h"

#ifdef MMC_LOGISTIC
  #include "logistic_rand.h"
#else
  #include "posix_randr.h"
#endif

#define MAX_TRIAL          3
#define FIX_PHOTON         1e-3f

void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void fixphoton(float3 *p,float3 *nodes, int *ee);
float onephoton(int id,tetplucker *plucker,tetmesh *mesh,Config *cfg,float rtstep,RandType *ran,RandType *ran0);

#endif
