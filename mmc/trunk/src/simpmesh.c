#include <stdlib.h>
#include "simpmesh.h"

void vec_diff(float3 *a,float3 *b,float3 *res){
	res->x=b->x-a->x;
	res->y=b->y-a->y;
	res->z=b->z-a->z;
}
void vec_cross(float3 *a,float3 *b,float3 *res){
	res->x=a->y*b->z-a->z*b->y;
	res->y=a->z*b->x-a->x*b->z;
	res->z=a->x*b->y-a->y*b->x;
}
float vec_dot(float3 *a,float3 *b){
	return a->x*b->x+a->y*b->y+a->z*b->z;
}
float pinner(float3 *Pd,float3 *Pm,float3 *Ad,float3 *Am){
	return vec_dot(Pd,Am)+vec_dot(Pm,Ad);
}

void mesh_init(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->node=NULL;
	mesh->elem=NULL;
	mesh->facenb=NULL;
}
void mesh_error(char *msg){
	fprintf(stderr,"%s\n",msg);
	exit(1);
}
void mesh_loadnode(tetmesh *mesh,char *fnode){
	FILE *fp;
	int tmp,len,i;
	if((fp=fopen(fnode,"rt"))==NULL){
		mesh_error("can not open node file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->nn));
	if(len!=2 || mesh->nn<=0){
		mesh_error("mesh file has wrong format");
	}
	mesh->node=(float3 *)malloc(sizeof(float3)*mesh->nn);
	for(i=0;i<mesh->nn;i++){
		if(fscanf(fp,"%d %f %f %f",&tmp,&(mesh->node[i].x),&(mesh->node[i].y),&(mesh->node[i].z))!=4)
			mesh_error("mesh file has wrong format");
	}
	fclose(fp);
}
void mesh_loadelem(tetmesh *mesh,char *felem){
	FILE *fp;
	int tmp,len,i;
	int4 *pe;

	if((fp=fopen(felem,"rt"))==NULL){
		mesh_error("can not open element file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		mesh_error("mesh file has wrong format");
	}
	mesh->elem=(int4 *)malloc(sizeof(int4)*mesh->ne);
	for(i=0;i<mesh->ne;i++){
		pe=mesh->elem+i;
		if(fscanf(fp,"%d %d %d %d %d %d",&tmp,&(pe->x),&(pe->y),&(pe->z),&(pe->w),&tmp)!=6)
			mesh_error("mesh file has wrong format");
	}
	fclose(fp);
}
void mesh_loadfaceneighbor(tetmesh *mesh,char *ffacenb){
	FILE *fp;
	int tmp,len,i;
	int4 *pe;

	if((fp=fopen(ffacenb,"rt"))==NULL){
		mesh_error("can not open element file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		mesh_error("mesh file has wrong format");
	}
	mesh->facenb=(int4 *)malloc(sizeof(int4)*mesh->ne);
	for(i=0;i<mesh->ne;i++){
		pe=mesh->facenb+i;
		if(fscanf(fp,"%d %d %d %d",&(pe->x),&(pe->y),&(pe->z),&(pe->w))!=4)
			mesh_error("mesh file has wrong format");
	}
	fclose(fp);
}
void mesh_clear(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	if(mesh->node){
		free(mesh->node);
		mesh->node=NULL;
	}
	if(mesh->elem){
		free(mesh->elem);
		mesh->elem=NULL;
	}
	if(mesh->facenb){
		free(mesh->facenb);
		mesh->facenb=NULL;
	}
}

void plucker_init(tetplucker *plucker,tetmesh *pmesh){
	plucker->d=NULL;
	plucker->m=NULL;
	plucker->mesh=pmesh;
	plucker_build(plucker);
}
void plucker_clear(tetplucker *plucker){
	if(plucker->d) {
		free(plucker->d);
		plucker->d=NULL;
	}
	if(plucker->m) {
		free(plucker->m);
		plucker->m=NULL;
	}	
	plucker->mesh=NULL;
}
void plucker_build(tetplucker *plucker){
	int nn,ne,i,j;
	int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	float3 *nodes;
	int *elems,ebase;
	int e1,e0;
	
	if(plucker->d || plucker->m || plucker->mesh==NULL) return;
	nn=plucker->mesh->nn;
	ne=plucker->mesh->ne;
	nodes=plucker->mesh->node;
	elems=(int *)(plucker->mesh->elem); // convert int4* to int*
	plucker->d=(float3*)malloc(sizeof(float3)*ne*6); // 6 edges/elem
	plucker->m=(float3*)malloc(sizeof(float3)*ne*6); // 6 edges/elem
	for(i=0;i<ne;i++){
		ebase=i<<2;
		for(j=0;j<6;j++){
			e1=elems[ebase+pairs[j][1]]-1;
			e0=elems[ebase+pairs[j][0]]-1;
			vec_diff(&nodes[e0],&nodes[e1],plucker->d+i*6+j);
			vec_cross(&nodes[e0],&nodes[e1],plucker->m+i*6+j);
		}
	}
}

float dist2(float3 *p0,float3 *p1){
    return (p1->x-p0->x)*(p1->x-p0->x)+(p1->y-p0->y)*(p1->y-p0->y)+(p1->z-p0->z)*(p1->z-p0->z);
}
