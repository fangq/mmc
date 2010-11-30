/*******************************************************************************
**  Mesh-based Monte Carlo (MMC)
**
**  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates," Biomed. Opt. Express, 1(1) 165-175 (2010)
**
**  (Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, 17(22) 20178-20190 (2009)
**
**  simpmesh.c: basic vector math and mesh operations
**
**  License: GPL v3, see LICENSE.txt for details
**
*******************************************************************************/

#include <stdlib.h>
#include "simpmesh.h"

#ifdef MMC_USE_SSE
#include <smmintrin.h>
#endif

#define SINCOSF(theta,stheta,ctheta) {stheta=sinf(theta);ctheta=cosf(theta);}

#ifdef WIN32
         char pathsep='\\';
#else
         char pathsep='/';
#endif

void vec_add(float3 *a,float3 *b,float3 *res){
	res->x=a->x+b->x;
	res->y=a->y+b->y;
	res->z=a->z+b->z;
}
void vec_diff(float3 *a,float3 *b,float3 *res){
        res->x=b->x-a->x;
        res->y=b->y-a->y;
        res->z=b->z-a->z;
}
void vec_mult(float3 *a,float sa,float3 *res){
        res->x=sa*a->x;
        res->y=sa*a->y;
        res->z=sa*a->z;
}
void vec_mult_add(float3 *a,float3 *b,float sa,float sb,float3 *res){
	res->x=sb*b->x+sa*a->x;
	res->y=sb*b->y+sa*a->y;
	res->z=sb*b->z+sa*a->z;
}
void vec_cross(float3 *a,float3 *b,float3 *res){
	res->x=a->y*b->z-a->z*b->y;
	res->y=a->z*b->x-a->x*b->z;
	res->z=a->x*b->y-a->y*b->x;
}

#ifndef MMC_USE_SSE
inline float vec_dot(float3 *a,float3 *b){
        return a->x*b->x+a->y*b->y+a->z*b->z;
}
#else

#ifndef __SSE4_1__
inline float vec_dot(float3 *a,float3 *b){
        float dot;
        __m128 na,nb,res;
        na=_mm_load_ps((float*)a);
        nb=_mm_load_ps((float*)b);
        res=_mm_mul_ps(na,nb);
        res=_mm_hadd_ps(res,res);
        res=_mm_hadd_ps(res,res);
        _mm_store_ss(&dot,res);
        return dot;   
}
#else
inline float vec_dot(float3 *a,float3 *b){
        float dot;
        __m128 na,nb,res;
        na=_mm_load_ps((float*)a);
        nb=_mm_load_ps((float*)b);
        res=_mm_dp_ps(na,nb,0x7f);
        _mm_store_ss(&dot,res);
        return dot;
}
#endif
        
#endif
 
inline float pinner(float3 *Pd,float3 *Pm,float3 *Ad,float3 *Am){
        return vec_dot(Pd,Am)+vec_dot(Pm,Ad);
}

void mesh_init(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->prop=0;
	mesh->node=NULL;
	mesh->elem=NULL;
	mesh->facenb=NULL;
	mesh->type=NULL;
	mesh->med=NULL;
	mesh->atte=NULL;
	mesh->weight=NULL;
	mesh->evol=NULL;
	mesh->nvol=NULL;
}
void mesh_error(char *msg){
	fprintf(stderr,"%s\n",msg);
	exit(1);
}
void mesh_filenames(char *format,char *foutput,Config *cfg){
	char filename[MAX_PATH_LENGTH];
	sprintf(filename,format,cfg->meshtag);

	if(cfg->rootpath[0]) 
		sprintf(foutput,"%s%c%s",cfg->rootpath,pathsep,filename);
	else
		sprintf(foutput,"%s",filename);
}
void mesh_loadnode(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int tmp,len,i;
	char fnode[MAX_PATH_LENGTH];
	mesh_filenames("node_%s.dat",fnode,cfg);
	if((fp=fopen(fnode,"rt"))==NULL){
		fprintf(stdout,"nodefile=%s\n",fnode);
		mesh_error("can not open node file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->nn));
	if(len!=2 || mesh->nn<=0){
		mesh_error("node file has wrong format");
	}
	mesh->node=(float3 *)calloc(sizeof(float3),mesh->nn);
	mesh->weight=(float *)calloc(sizeof(float)*mesh->nn,cfg->maxgate);

	for(i=0;i<mesh->nn;i++){
		if(fscanf(fp,"%d %f %f %f",&tmp,&(mesh->node[i].x),&(mesh->node[i].y),&(mesh->node[i].z))!=4)
			mesh_error("node file has wrong format");
	}
	if(cfg->unitinmm!=1.f)
	  for(i=0;i<mesh->nn;i++){
		mesh->node[i].x*=cfg->unitinmm;
		mesh->node[i].y*=cfg->unitinmm;
		mesh->node[i].z*=cfg->unitinmm;
	  }
	fclose(fp);
}

void mesh_loadmedia(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int tmp,len,i;
	char fmed[MAX_PATH_LENGTH];
	mesh_filenames("prop_%s.dat",fmed,cfg);
	if((fp=fopen(fmed,"rt"))==NULL){
		mesh_error("can not open media property file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->prop));
	if(len!=2 || mesh->prop<=0){
		mesh_error("property file has wrong format");
	}
	mesh->med=(medium *)calloc(sizeof(medium),mesh->prop);
        mesh->atte=(float *)calloc(sizeof(float),mesh->prop);

	for(i=0;i<mesh->prop;i++){
		if(fscanf(fp,"%d %f %f %f %f",&tmp,&(mesh->med[i].mua),&(mesh->med[i].mus),
		                                   &(mesh->med[i].g),&(mesh->med[i].n))!=5)
			mesh_error("property file has wrong format");
		mesh->atte[i]=expf(-cfg->minstep*mesh->med[i].mua);
		/*user input musp, MMCM converts to mus
		mesh->med[i].mus=musp/(1.f-mesh->med[i].g); */
	}
	fclose(fp);
}
void mesh_loadelem(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int tmp,len,i;
	int4 *pe;
	char felem[MAX_PATH_LENGTH];
	mesh_filenames("elem_%s.dat",felem,cfg);
	if((fp=fopen(felem,"rt"))==NULL){
		mesh_error("can not open element file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		mesh_error("element file has wrong format");
	}
	mesh->elem=(int4 *)malloc(sizeof(int4)*mesh->ne);
	mesh->type=(int  *)malloc(sizeof(int )*mesh->ne);

	for(i=0;i<mesh->ne;i++){
		pe=mesh->elem+i;
		if(fscanf(fp,"%d %d %d %d %d %d",&tmp,&(pe->x),&(pe->y),&(pe->z),&(pe->w),mesh->type+i)!=6)
			mesh_error("element file has wrong format");
	}
	fclose(fp);
}
void mesh_loadelemvol(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int tmp,len,i,j,*ee;
	char fvelem[MAX_PATH_LENGTH];
	mesh_filenames("velem_%s.dat",fvelem,cfg);
	if((fp=fopen(fvelem,"rt"))==NULL){
		mesh_error("can not open element volume file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		mesh_error("mesh file has wrong format");
	}
        mesh->evol=(float *)malloc(sizeof(float)*mesh->ne);
	mesh->nvol=(float *)calloc(sizeof(float),mesh->nn);

	for(i=0;i<mesh->ne;i++){
		if(fscanf(fp,"%d %f",&tmp,mesh->evol+i)!=2)
			mesh_error("mesh file has wrong format");
		ee=(int *)(mesh->elem+i);
		for(j=0;j<4;j++)
			mesh->nvol[ee[j]-1]+=mesh->evol[i]*0.25f;
	}
	if(cfg->unitinmm!=1.f){
	  float unit3d=cfg->unitinmm*cfg->unitinmm*cfg->unitinmm;
	  for(i=0;i<mesh->ne;i++)
		mesh->evol[i]*=unit3d;
	  for(i=0;i<mesh->nn;i++)
		mesh->nvol[i]*=unit3d;
	}
	fclose(fp);
}
void mesh_loadfaceneighbor(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int tmp,len,i;
	int4 *pe;
	char ffacenb[MAX_PATH_LENGTH];
	mesh_filenames("facenb_%s.dat",ffacenb,cfg);

	if((fp=fopen(ffacenb,"rt"))==NULL){
		mesh_error("can not open face-neighbor list file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		mesh_error("mesh file has wrong format");
	}
	mesh->facenb=(int4 *)malloc(sizeof(int4)*mesh->ne);
	for(i=0;i<mesh->ne;i++){
		pe=mesh->facenb+i;
		if(fscanf(fp,"%d %d %d %d",&(pe->x),&(pe->y),&(pe->z),&(pe->w))!=4)
			mesh_error("face-neighbor list file has wrong format");
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
	if(mesh->type){
		free(mesh->type);
		mesh->type=NULL;
	}
	if(mesh->med){
		free(mesh->med);
		mesh->med=NULL;
	}
        if(mesh->atte){
                free(mesh->atte);
                mesh->atte=NULL;
        }
	if(mesh->weight){
		free(mesh->weight);
		mesh->weight=NULL;
	}
        if(mesh->evol){
                free(mesh->evol);
                mesh->evol=NULL;
        }
	if(mesh->nvol){
		free(mesh->nvol);
		mesh->nvol=NULL;
	}
}

void plucker_init(tetplucker *plucker,tetmesh *pmesh,int mode){
	plucker->d=NULL;
	plucker->m=NULL;
	plucker->mesh=pmesh;
	plucker->isplucker=mode;
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
	const int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

	float3 *nodes;
	int *elems,ebase;
	int e1,e0;

	if(plucker->d || plucker->m || plucker->mesh==NULL) return;
        if(plucker->mesh->node==NULL||plucker->mesh->elem==NULL||plucker->mesh->facenb==NULL||plucker->mesh->med==NULL)
                mesh_error("encountered error while loading mesh files");

	nn=plucker->mesh->nn;
	ne=plucker->mesh->ne;
	nodes=plucker->mesh->node;
	elems=(int *)(plucker->mesh->elem); // convert int4* to int*
	if(plucker->isplucker){
		plucker->d=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		plucker->m=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		for(i=0;i<ne;i++){
			ebase=i<<2;
			for(j=0;j<6;j++){
				e1=elems[ebase+pairs[j][1]]-1;
				e0=elems[ebase+pairs[j][0]]-1;
				vec_diff(&nodes[e0],&nodes[e1],plucker->d+i*6+j);
				vec_cross(&nodes[e0],&nodes[e1],plucker->m+i*6+j);
			}
		}
	}else{
		int ea,eb,ec;
		const int out[4][3]={{0,3,1},{3,2,1},{0,2,3},{0,1,2}};
		float3 vecAB={0.f},vecAC={0.f};

		plucker->d=NULL;
		plucker->m=(float3*)calloc(sizeof(float3),ne*12);
                for(i=0;i<ne;i++){
                        ebase=i<<2;
			for(j=0;j<4;j++){
				float3 *vecN=plucker->m+3*((i<<2)+j);
				float Rn2;

                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);

				vec_cross(&vecAB,&vecAC,vecN); /*N is defined as ACxAB in Jiri's code, but not the paper*/
                                vec_cross(&vecAC,vecN,vecN+1);
                                vec_cross(vecN,&vecAB,vecN+2);

				Rn2=1.f/(vecN->x*vecN->x+vecN->y*vecN->y+vecN->z*vecN->z);
				vec_mult(vecN+1,Rn2,vecN+1);
                                vec_mult(vecN+2,Rn2,vecN+2);
#ifdef MMC_USE_SSE
				vecN->w    = vec_dot(vecN,  &nodes[ea]);
				(vecN+1)->w=-vec_dot(vecN+1,&nodes[ea]);
                                (vecN+2)->w=-vec_dot(vecN+2,&nodes[ea]);
#endif
			}
                }
	}
}

inline float dist2(float3 *p0,float3 *p1){
    return (p1->x-p0->x)*(p1->x-p0->x)+(p1->y-p0->y)*(p1->y-p0->y)+(p1->z-p0->z)*(p1->z-p0->z);
}

inline float dist(float3 *p0,float3 *p1){
    return sqrt(dist2(p0,p1));
}

float mc_next_scatter(float g, float3 *dir,RandType *ran, RandType *ran0, Config *cfg){
    float nextslen;
    float sphi,cphi,tmp0,theta,stheta,ctheta,tmp1;
    float3 p;

    rand_need_more(ran,ran0);

    //random scattering length (normalized)
    nextslen=rand_next_scatlen(ran);

    //random arimuthal angle
    tmp0=TWO_PI*rand_next_aangle(ran); //next arimuth angle
    SINCOSF(tmp0,sphi,cphi);

    //Henyey-Greenstein Phase Function, "Handbook of Optical Biomedical Diagnostics",2002,Chap3,p234
    //see Boas2002

    if(g>EPS){  //if g is too small, the distribution of theta is bad
	tmp0=(1.f-g*g)/(1.f-g+2.f*g*rand_next_zangle(ran));
	tmp0*=tmp0;
	tmp0=(1.f+g*g-tmp0)/(2.f*g);

    	// when ran=1, CUDA will give me 1.000002 for tmp0 which produces nan later
    	if(tmp0> 1.f) tmp0=1.f;
        if(tmp0<-1.f) tmp0=-1.f;

	theta=acosf(tmp0);
	stheta=sinf(theta);
	ctheta=tmp0;
    }else{  //Wang1995 has acos(2*ran-1), rather than 2*pi*ran, need to check
	theta=M_PI*rand_next_zangle(ran);
    	SINCOSF(theta,stheta,ctheta);
    }

    if( dir->z>-1.f+EPS && dir->z<1.f-EPS ) {
	tmp0=1.f-dir->z*dir->z;   //reuse tmp to minimize registers
	tmp1=1.f/sqrtf(tmp0);
	tmp1=stheta*tmp1;

	p.x=tmp1*(dir->x*dir->z*cphi - dir->y*sphi) + dir->x*ctheta;
	p.y=tmp1*(dir->y*dir->z*cphi + dir->x*sphi) + dir->y*ctheta;
	p.z=-tmp1*tmp0*cphi			    + dir->z*ctheta;
    }else{
	p.x=stheta*cphi;
	p.y=stheta*sphi;
	p.z=(dir->z>0.f)?ctheta:-ctheta;
    }
    dir->x=p.x;
    dir->y=p.y;
    dir->z=p.z;
    return nextslen;
}

void mesh_saveweight(tetmesh *mesh,Config *cfg){
	FILE *fp;
	int i,j;
	float3 *pn;
	char fweight[MAX_PATH_LENGTH];
        if(cfg->rootpath[0])
                sprintf(fweight,"%s%c%s.dat",cfg->rootpath,pathsep,cfg->session);
        else
                sprintf(fweight,"%s.dat",cfg->session);

	if((fp=fopen(fweight,"wt"))==NULL){
		mesh_error("can not open weight file to write");
	}
	for(i=0;i<cfg->maxgate;i++)
	   for(j=0;j<mesh->nn;j++){
		pn=mesh->node+j;
		/*if(fprintf(fp,"%d %e %e %e %e\n",j+1,pn->x,pn->y,pn->z,mesh->weight[i*mesh->nn+j])==0)*/
		if(fprintf(fp,"%d\t%e\n",j+1,mesh->weight[i*mesh->nn+j])==0)
			mesh_error("can not write to weight file");
	   }
	fclose(fp);
}

/*see Eq (1) in Fang&Boas, Opt. Express, vol 17, No.22, pp. 20178-20190, Oct 2009*/
float mesh_normalize(tetmesh *mesh,Config *cfg, float Eabsorb, float Etotal){
        int i,j,k;
	float energydeposit=0.f, energyelem,normalizor;
	int *ee;

        for(i=0;i<cfg->maxgate;i++)
            for(j=0;j<mesh->nn;j++)
              mesh->weight[i*mesh->nn+j]/=mesh->nvol[j];

        for(i=0;i<mesh->ne;i++){
	   ee=(int *)(mesh->elem+i);
	   energyelem=0.f;
	   for(j=0;j<cfg->maxgate;j++)
	     for(k=0;k<4;k++)
		energyelem+=mesh->weight[j*mesh->nn+ee[k]-1]; /*1/4 factor is absorbed two lines below*/

	   energydeposit+=energyelem*mesh->evol[i]*mesh->med[mesh->type[i]-1].mua; /**mesh->med[mesh->type[i]-1].n;*/
	}

	energydeposit*=0.25f*cfg->tstep; /* unit conversions */
	normalizor=Eabsorb/(Etotal*energydeposit); /*scaling factor*/

        for(i=0;i<cfg->maxgate;i++)
           for(j=0;j<mesh->nn;j++)
	      mesh->weight[i*mesh->nn+j]*=normalizor;

	return normalizor;
}
