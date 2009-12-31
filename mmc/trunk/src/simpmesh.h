#ifndef MCX_FEM_RAY_H
#define MCX_FEM_RAY_H

#include <stdio.h>
#include <math.h>
#include "mcx_utils.h"


#define QLIMIT (3.40282347e+38F)
#define R_RAND_MAX  (1.f/RAND_MAX)
#define TWO_PI  (M_PI*2.0)
#define EPS  1e-9f
#define LOG_MT_MAX 22.1807097779182f


typedef struct MMCMedium{
	float mua,musp,n,g;
} medium;

typedef struct femmesh{
	int nn; // number of nodes
	int ne; // number of elems
	int prop;
	float3 *node;
	int4 *elem;
	int  *type;
	int4 *facenb;
	medium *med;
	float *weight;
} tetmesh;

typedef struct tplucker{
	tetmesh *mesh;
	float3 *d;
	float3 *m;
} tetplucker;

void vec_diff(float3 *a,float3 *b,float3 *res);
void vec_cross(float3 *a,float3 *b,float3 *res);
void vec_mult_add(float3 *a,float3 *b,float sa,float sb,float3 *res);
float vec_dot(float3 *a,float3 *b);
float pinner(float3 *Pd,float3 *Pm,float3 *Ad,float3 *Am);

void mesh_init(tetmesh *mesh);
void mesh_loadnode(tetmesh *mesh,Config *cfg);
void mesh_loadelem(tetmesh *mesh,Config *cfg);
void mesh_loadfaceneighbor(tetmesh *mesh,Config *cfg);
void mesh_loadmedia(tetmesh *mesh,Config *cfg);

void mesh_clear(tetmesh *mesh);
void mesh_build(tetmesh *mesh);
void mesh_error(char *msg);
void mesh_filenames(char *format,char *foutput,Config *cfg);
void mesh_saveweight(tetmesh *mesh,Config *cfg);

void plucker_init(tetplucker *plucker,tetmesh *mesh);
void plucker_build(tetplucker *plucker);
void plucker_clear(tetplucker *plucker);
float dist2(float3 *p0,float3 *p1);
float dist(float3 *p0,float3 *p1);
void mc_next_scatter(float g, float musp, float3 *pnext);
float rand01();


#endif
