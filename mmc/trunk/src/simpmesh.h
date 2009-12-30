#ifndef MCX_FEM_RAY_H
#define MCX_FEM_RAY_H

#include <stdio.h>
#include <math.h>

#define QLIMIT (3.40282347e+38F)

typedef struct float_3{
	float x,y,z;
} float3;

typedef struct int_4{
	int x,y,z,w;
} int4;

typedef struct femmesh{
	int nn; // number of nodes
	int ne; // number of elems
	float3 *node;
	int4 *elem;
	int4 *facenb;
} tetmesh;

typedef struct tplucker{
	tetmesh *mesh;
	float3 *d;
	float3 *m;
} tetplucker;

void vec_diff(float3 *a,float3 *b,float3 *res);
void vec_cross(float3 *a,float3 *b,float3 *res);
float vec_dot(float3 *a,float3 *b);
float pinner(float3 *Pd,float3 *Pm,float3 *Ad,float3 *Am);

void mesh_init(tetmesh *mesh);
void mesh_loadnode(tetmesh *mesh,char *node);
void mesh_loadelem(tetmesh *mesh,char *node);
void mesh_loadfaceneighbor(tetmesh *mesh,char *ffacenb);
void mesh_clear(tetmesh *mesh);
void mesh_build(tetmesh *mesh);

void plucker_init(tetplucker *plucker,tetmesh *mesh);
void plucker_build(tetplucker *plucker);
void plucker_clear(tetplucker *plucker);
float dist2(float3 *p0,float3 *p1);

#endif
