#ifndef MCX_FEM_TRACING_H
#define MCX_FEM_TRACING_H

#include "simpmesh.h"

void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
float getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void trackpos(float3 *p0,float3 *p1,tetplucker *plucker,int eid /*start from 1*/, 
              float3 *pout, int *faceid, float *weight, int *isend);
#endif
