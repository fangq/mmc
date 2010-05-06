#ifndef MCX_FEM_TRACING_H
#define MCX_FEM_TRACING_H

#include "simpmesh.h"
#include "mcx_utils.h"

#define R_C0               3.335640951981520e-12f  //1/C0 in s/mm


void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
float trackpos(float3 *p0,float3 *pvec,tetplucker *plucker,int eid /*start from 1*/, 
              float3 *pout, float slen, int *faceid, float *weight, 
	      int *isend,float *photontimer,float rtstep, Config *cfg);
float onephoton(tetplucker *plucker,tetmesh *mesh,Config *cfg,float rtstep);
#endif
