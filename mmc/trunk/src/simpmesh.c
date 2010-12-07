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
**  License: GPL v3, see LICENSE.txt for details
**
*******************************************************************************/

/***************************************************************************//**
\file    simpmesh.c

\brief   Basic vector math and mesh operations
*******************************************************************************/

#include <stdlib.h>
#include "simpmesh.h"

#ifdef WIN32
         char pathsep='\\';
#else
         char pathsep='/';
#endif


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
	if(cfg->basisorder) 
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
	if(!cfg->basisorder)
	   mesh->weight=(float *)calloc(sizeof(float)*mesh->ne,cfg->maxgate);

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

void tracer_init(raytracer *tracer,tetmesh *pmesh){
	tracer->d=NULL;
	tracer->m=NULL;
	tracer->mesh=pmesh;
	tracer_build(tracer);
}
void tracer_clear(raytracer *tracer){
	if(tracer->d) {
		free(tracer->d);
		tracer->d=NULL;
	}
	if(tracer->m) {
		free(tracer->m);
		tracer->m=NULL;
	}	
	tracer->mesh=NULL;
}
void tracer_build(raytracer *tracer){
	int nn,ne,i,j;
	const int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	float3 *nodes;
	int *elems,ebase;
	int e1,e0;
	
	if(tracer->d || tracer->m || tracer->mesh==NULL) return;
	nn=tracer->mesh->nn;
	ne=tracer->mesh->ne;
	nodes=tracer->mesh->node;
	elems=(int *)(tracer->mesh->elem); // convert int4* to int*
	tracer->d=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
	tracer->m=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
	for(i=0;i<ne;i++){
		ebase=i<<2;
		for(j=0;j<6;j++){
			e1=elems[ebase+pairs[j][1]]-1;
			e0=elems[ebase+pairs[j][0]]-1;
			vec_diff(&nodes[e0],&nodes[e1],tracer->d+i*6+j);
			vec_cross(&nodes[e0],&nodes[e1],tracer->m+i*6+j);
		}
	}
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
    mmc_sincosf(tmp0,&sphi,&cphi);

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
    	mmc_sincosf(theta,&stheta,&ctheta);
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
	if(cfg->basisorder)
	  for(i=0;i<cfg->maxgate;i++)
	   for(j=0;j<mesh->nn;j++){
		pn=mesh->node+j;
		/*if(fprintf(fp,"%d %e %e %e %e\n",j+1,pn->x,pn->y,pn->z,mesh->weight[i*mesh->nn+j])==0)*/
		if(fprintf(fp,"%d\t%e\n",j+1,mesh->weight[i*mesh->nn+j])==0)
			mesh_error("can not write to weight file");
	   }
	else
	  for(i=0;i<cfg->maxgate;i++)
	   for(j=0;j<mesh->ne;j++){
		if(fprintf(fp,"%d\t%e\n",j+1,mesh->weight[i*mesh->ne+j])==0)
			mesh_error("can not write to weight file");
	   }
	fclose(fp);
}

/*see Eq (1) in Fang&Boas, Opt. Express, vol 17, No.22, pp. 20178-20190, Oct 2009*/
float mesh_normalize(tetmesh *mesh,Config *cfg, float Eabsorb, float Etotal){
        int i,j,k;
	float energydeposit=0.f, energyelem,normalizor;
	int *ee;

	if(cfg->basisorder){
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
	    normalizor=Eabsorb/(Etotal*energydeposit*0.25f*cfg->tstep); /*scaling factor*/

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<mesh->nn;j++)
		  mesh->weight[i*mesh->nn+j]*=normalizor;
	}else{
            for(i=0;i<mesh->ne;i++)
	      for(j=0;j<cfg->maxgate;j++)
	         energydeposit+=mesh->weight[j*mesh->ne+i];

            for(i=0;i<mesh->ne;i++){
	      energyelem=mesh->evol[i]*mesh->med[mesh->type[i]-1].mua;
              for(j=0;j<cfg->maxgate;j++)
        	mesh->weight[j*mesh->ne+i]/=energyelem;
	    }
            normalizor=Eabsorb/(Etotal*energydeposit*cfg->tstep); /*scaling factor*/

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<mesh->ne;j++)
                  mesh->weight[i*mesh->ne+j]*=normalizor;
	}

	return normalizor;
}
