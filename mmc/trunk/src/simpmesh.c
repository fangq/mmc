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
#include <string.h>

#ifdef WIN32
         char pathsep='\\';
#else
         char pathsep='/';
#endif

const int out[4][3]={{0,3,1},{3,2,1},{0,2,3},{0,1,2}};
const int facemap[]={2,0,1,3};
const int ifacemap[]={1,2,0,3};

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

void mesh_init_from_cfg(tetmesh *mesh,mcconfig *cfg){
        mesh_init(mesh);
        mesh_loadnode(mesh,cfg);
        mesh_loadelem(mesh,cfg);
        mesh_loadfaceneighbor(mesh,cfg);
        mesh_loadmedia(mesh,cfg);
        mesh_loadelemvol(mesh,cfg);
	if(cfg->seed==SEED_FROM_FILE && cfg->seedfile[0]){
          mesh_loadseedfile(mesh,cfg);
        }
}

void mesh_error(char *msg){
#ifdef MCX_CONTAINER
        mmc_throw_exception(1,msg,__FILE__,__LINE__);
#else
	fprintf(stderr,"%s\n",msg);
        exit(1);
#endif
}
void mesh_filenames(char *format,char *foutput,mcconfig *cfg){
	char filename[MAX_PATH_LENGTH];
	sprintf(filename,format,cfg->meshtag);

	if(cfg->rootpath[0]) 
		sprintf(foutput,"%s%c%s",cfg->rootpath,pathsep,filename);
	else
		sprintf(foutput,"%s",filename);
}
void mesh_loadnode(tetmesh *mesh,mcconfig *cfg){
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
	   mesh->weight=(double *)calloc(sizeof(double)*mesh->nn,cfg->maxgate);

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

void mesh_loadmedia(tetmesh *mesh,mcconfig *cfg){
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
	mesh->med=(medium *)calloc(sizeof(medium),mesh->prop+1);
        mesh->atte=(float *)calloc(sizeof(float),mesh->prop+1);
	
	mesh->med[0].mua=0.f;
	mesh->med[0].mus=0.f;
	mesh->med[0].n=cfg->nout;
	mesh->med[0].g=1.f;

	for(i=1;i<=mesh->prop;i++){
		if(fscanf(fp,"%d %f %f %f %f",&tmp,&(mesh->med[i].mua),&(mesh->med[i].mus),
		                                   &(mesh->med[i].g),&(mesh->med[i].n))!=5)
			mesh_error("property file has wrong format");
		/*mesh->atte[i]=expf(-cfg->minstep*mesh->med[i].mua);*/
	}
	fclose(fp);
	cfg->his.maxmedia=mesh->prop; /*skip media 0*/
}
void mesh_loadelem(tetmesh *mesh,mcconfig *cfg){
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
	   mesh->weight=(double *)calloc(sizeof(double)*mesh->ne,cfg->maxgate);

	for(i=0;i<mesh->ne;i++){
		pe=mesh->elem+i;
		if(fscanf(fp,"%d %d %d %d %d %d",&tmp,&(pe->x),&(pe->y),&(pe->z),&(pe->w),mesh->type+i)!=6)
			mesh_error("element file has wrong format");
	}
	fclose(fp);
}
void mesh_loadelemvol(tetmesh *mesh,mcconfig *cfg){
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
void mesh_loadfaceneighbor(tetmesh *mesh,mcconfig *cfg){
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

void mesh_loadseedfile(tetmesh *mesh, mcconfig *cfg){
    history his;
    FILE *fp=fopen(cfg->seedfile,"rb");
    if(fp==NULL)
        mesh_error("can not open the specified history file");
    if(fread(&his,sizeof(history),1,fp)!=1)
        mesh_error("error when reading the history file");
    if(his.savedphoton==0 || his.seedbyte==0){
        fclose(fp);
        return;
    }
    if(his.maxmedia!=mesh->prop)
        mesh_error("the history file was generated with a different media setting");
    if(fseek(fp,his.savedphoton*his.colcount*sizeof(float),SEEK_CUR))
        mesh_error("illegal history file");
    cfg->photonseed=malloc(his.savedphoton*his.seedbyte);
    if(cfg->photonseed==NULL)
        mesh_error("can not allocate memory");
    if(fread(cfg->photonseed,his.seedbyte,his.savedphoton,fp)!=his.savedphoton)
        mesh_error("error when reading the seed data");
    cfg->seed=SEED_FROM_FILE;
    cfg->nphoton=his.savedphoton;

    if(cfg->outputtype==otJacobian || cfg->outputtype==otTaylor){ //cfg->replaydet>0
       int i,j;
       float *ppath=(float*)malloc(his.savedphoton*his.colcount*sizeof(float));
       cfg->replayweight=(float*)malloc(his.savedphoton*sizeof(float));
       fseek(fp,sizeof(his),SEEK_SET);
       if(fread(ppath,his.colcount*sizeof(float),his.savedphoton,fp)!=his.savedphoton)
           mesh_error("error when reading the seed data");

       cfg->nphoton=0;
       for(i=0;i<his.savedphoton;i++)
           if(cfg->replaydet==0 || cfg->replaydet==(int)(ppath[i*his.colcount])){
               memcpy((char *)(cfg->photonseed)+cfg->nphoton*his.seedbyte, (char *)(cfg->photonseed)+i*his.seedbyte, his.seedbyte);
               cfg->replayweight[cfg->nphoton]=1.f;
               for(j=2;j<his.maxmedia+2;j++)
                  cfg->replayweight[cfg->nphoton]*=expf(-mesh->med[j-1].mua*ppath[i*his.colcount+j]*his.unitinmm);
               cfg->nphoton++;
           }
	free(ppath);
        cfg->photonseed=realloc(cfg->photonseed, cfg->nphoton*his.seedbyte);
        cfg->replayweight=(float*)realloc(cfg->replayweight, cfg->nphoton*sizeof(float));
	cfg->minenergy=0.f;
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

void tracer_init(raytracer *tracer,tetmesh *pmesh,char methodid){
	tracer->d=NULL;
	tracer->m=NULL;
	tracer->n=NULL;
	tracer->mesh=pmesh;
	tracer->method=methodid;
	tracer_build(tracer);
}

void tracer_prep(raytracer *tracer,mcconfig *cfg){
	if(tracer->n==NULL && tracer->m==NULL && tracer->d==NULL){
	    if(tracer->mesh!=NULL)
		tracer_build(tracer);
	    else
	    	mesh_error("tracer is not associated with a mesh");
	}else{
            int eid=cfg->dim.x-1;
	    float3 vecS={0.f}, *nodes=tracer->mesh->node, vecAB, vecAC, vecN;
	    int i,ea,eb,ec;
	    float s=0.f, *bary=&(cfg->bary0.x);
	    int *elems=(int *)(tracer->mesh->elem+eid); // convert int4* to int*

	    for(i=0;i<4;i++){
            	ea=elems[out[i][0]]-1;
            	eb=elems[out[i][1]]-1;
	    	ec=elems[out[i][2]]-1;
            	vec_diff(&nodes[ea],&nodes[eb],&vecAB);
            	vec_diff(&nodes[ea],&nodes[ec],&vecAC);
	    	vec_diff(&nodes[ea],&(cfg->srcpos),&vecS);
            	vec_cross(&vecAB,&vecAC,&vecN);
	    	bary[facemap[i]]=-vec_dot(&vecS,&vecN);
	    }
	    if(cfg->debuglevel&dlWeight)
	       fprintf(cfg->flog,"initial bary-centric volumes [%e %e %e %e]\n",
	           bary[0]/6.,bary[1]/6.,bary[2]/6.,bary[3]/6.);
	    for(i=0;i<4;i++){
	        if(bary[i]<0.f)
		    mesh_error("initial element does not enclose the source!");
	        s+=bary[i];
	    }
	    for(i=0;i<4;i++){
	        bary[i]/=s;
		if(bary[i]<1e-5f)
		    cfg->dim.y=ifacemap[i]+1;
	    }
	}
}

void tracer_build(raytracer *tracer){
	int ne,i,j;
	const int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

	float3 *nodes;
	int *elems,ebase;
	int e1,e0;
	float Rn2;

	if(tracer->d || tracer->m || tracer->n || tracer->mesh==NULL) return;

        if(tracer->mesh->node==NULL||tracer->mesh->elem==NULL||
	   tracer->mesh->facenb==NULL||tracer->mesh->med==NULL)
                mesh_error("mesh is missing");

	ne=tracer->mesh->ne;
	nodes=tracer->mesh->node;
	elems=(int *)(tracer->mesh->elem); // convert int4* to int*
	if(tracer->method==rtPlucker){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f};

		tracer->d=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		tracer->m=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		tracer->n=(float3*)calloc(sizeof(float3),ne*4); // 4 face norms
		for(i=0;i<ne;i++){
			ebase=i<<2;
			for(j=0;j<6;j++){
				e1=elems[ebase+pairs[j][1]]-1;
				e0=elems[ebase+pairs[j][0]]-1;
				vec_diff(&nodes[e0],&nodes[e1],tracer->d+i*6+j);
				vec_cross(&nodes[e0],&nodes[e1],tracer->m+i*6+j);
			}
			for(j=0;j<4;j++){
                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);
				vec_cross(&vecAB,&vecAC,tracer->n+ebase+j);

				Rn2=1.f/sqrt(vec_dot(tracer->n+ebase+j,tracer->n+ebase+j));
				vec_mult(tracer->n+ebase+j,Rn2,tracer->n+ebase+j);
			}
		}
	}else if(tracer->method==rtHavel || tracer->method==rtBadouel){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f};

		tracer->d=NULL;
		tracer->m=(float3*)calloc(sizeof(float3),ne*12);
                for(i=0;i<ne;i++){
                        ebase=i<<2;
			for(j=0;j<4;j++){
				float3 *vecN=tracer->m+3*(ebase+j);

                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);

				vec_cross(&vecAB,&vecAC,vecN); /*N is defined as ACxAB in Jiri's code, but not the paper*/
                                vec_cross(&vecAC,vecN,vecN+1);
                                vec_cross(vecN,&vecAB,vecN+2);

				Rn2=1.f/sqrt(vec_dot(vecN,vecN));

				vec_mult(vecN,Rn2,vecN);
				
				Rn2*=Rn2;
				vec_mult(vecN+1,Rn2,vecN+1);
                                vec_mult(vecN+2,Rn2,vecN+2);
#ifdef MMC_USE_SSE
				vecN->w    = vec_dot(vecN,  &nodes[ea]);
				(vecN+1)->w=-vec_dot(vecN+1,&nodes[ea]);
                                (vecN+2)->w=-vec_dot(vecN+2,&nodes[ea]);
#endif
			}
                }
	}else if(tracer->method==rtBLBadouel){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f},vN={0.f};

		tracer->d=NULL;
		tracer->n=(float3*)calloc(sizeof(float3),ne*4);
                for(i=0;i<ne;i++){
                        ebase=i<<2;
			float *vecN=&(tracer->n[ebase].x);
			for(j=0;j<4;j++){
                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);

				vec_cross(&vecAB,&vecAC,&vN); /*N is defined as ACxAB in Jiri's code, but not the paper*/

				Rn2=1.f/sqrt(vec_dot(&vN,&vN));
				vec_mult(&vN,Rn2,&vN);

				vecN[j]=vN.x;
				vecN[j+4]=vN.y;
				vecN[j+8]=vN.z;
#ifdef MMC_USE_SSE
				vecN[j+12]    = vec_dot(&vN, &nodes[ea]);
#endif
			}
                }
	}
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
	if(tracer->n) {
		free(tracer->n);
		tracer->n=NULL;
	}	
	tracer->mesh=NULL;
}


float mc_next_scatter(float g, float3 *dir,RandType *ran, RandType *ran0, mcconfig *cfg, float *pmom){

    float nextslen;
    float sphi,cphi,tmp0,theta,stheta,ctheta,tmp1;
    float3 p;

    rand_need_more(ran,ran0);

    //random scattering length (normalized)
#ifdef MMC_USE_SSE_MATH
    nextslen=rand_next_scatlen_ps(ran);
#else
    nextslen=rand_next_scatlen(ran);
#endif

    //random arimuthal angle
#ifdef MMC_USE_SSE_MATH
    rand_next_aangle_sincos(ran,&sphi,&cphi);
#else
    tmp0=TWO_PI*rand_next_aangle(ran); //next arimuth angle
    mmc_sincosf(tmp0,&sphi,&cphi);
#endif

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
	stheta=sqrt(1.f-tmp0*tmp0);
	//stheta=sinf(theta);
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
    if (cfg->ismomentum)
        pmom[0]+=(1.f-ctheta);

    dir->x=p.x;
    dir->y=p.y;
    dir->z=p.z;
    return nextslen;
}

void mesh_saveweightat(tetmesh *mesh,mcconfig *cfg,int id){
	char sess[MAX_SESSION_LENGTH];
	int i,found=0;
	for(i=0;i<MAX_CHECKPOINT;i++){
           if(cfg->checkpt[i]==0)  return;
	   if(id==cfg->checkpt[i]) {
		found=1;
		break;
	   }
	}
        if(!found) return;
	memcpy(sess,cfg->session,MAX_SESSION_LENGTH);
	sprintf(cfg->session,"%s_%d",sess,id);
	mesh_saveweight(mesh,cfg);
	memcpy(cfg->session,sess,MAX_SESSION_LENGTH);
}

void mesh_saveweight(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int i,j;
	char fweight[MAX_PATH_LENGTH];
        if(cfg->rootpath[0])
                sprintf(fweight,"%s%c%s.dat",cfg->rootpath,pathsep,cfg->session);
        else
                sprintf(fweight,"%s.dat",cfg->session);

        if(cfg->outputformat==ofBin){
		if((fp=fopen(fweight,"wb"))==NULL)
         	        mesh_error("can not open weight file to write");
		if(fwrite((void*)mesh->weight,sizeof(mesh->weight[0]),mesh->nn*cfg->maxgate,fp)!=mesh->nn*cfg->maxgate)
			mesh_error("fail to write binary weight file");
		fclose(fp);
		return;
	}
	if((fp=fopen(fweight,"wt"))==NULL){
		mesh_error("can not open weight file to write");
	}
	if(cfg->basisorder)
	  for(i=0;i<cfg->maxgate;i++)
	   for(j=0;j<mesh->nn;j++){
		/*pn=mesh->node+j;
		if(fprintf(fp,"%d %e %e %e %e\n",j+1,pn->x,pn->y,pn->z,mesh->weight[i*mesh->nn+j])==0)*/
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

void mesh_savedetphoton(float *ppath, void *seeds, int count, int seedbyte, mcconfig *cfg){
	FILE *fp;
	char fhistory[MAX_PATH_LENGTH];
        if(cfg->rootpath[0])
                sprintf(fhistory,"%s%c%s.mch",cfg->rootpath,pathsep,cfg->session);
        else
                sprintf(fhistory,"%s.mch",cfg->session);

	if((fp=fopen(fhistory,"wb"))==NULL){
		mesh_error("can not open history file to write");
	}
	cfg->his.totalphoton=cfg->nphoton;
        cfg->his.unitinmm=cfg->unitinmm;
        cfg->his.detected=count;
	cfg->his.savedphoton=count;
	cfg->his.detnum=cfg->detnum;
	if(cfg->issaveseed && seeds!=NULL){
	   cfg->his.seedbyte=seedbyte;
        }
        cfg->his.colcount=(1+(cfg->ismomentum>0))*cfg->his.maxmedia+(cfg->issaveexit>0)*6+2; /*column count=maxmedia+2*/

	fwrite(&(cfg->his),sizeof(history),1,fp);
	fwrite(ppath,sizeof(float),count*cfg->his.colcount,fp);
	if(cfg->issaveseed && seeds!=NULL)
           fwrite(seeds,seedbyte,count,fp);
	fclose(fp);
}

/*see Eq (1) in Fang&Boas, Opt. Express, vol 17, No.22, pp. 20178-20190, Oct 2009*/
float mesh_normalize(tetmesh *mesh,mcconfig *cfg, float Eabsorb, float Etotal){
        int i,j,k;
	float energydeposit=0.f, energyelem,normalizor;
	int *ee;

	if(cfg->seed==SEED_FROM_FILE && (cfg->outputtype==otJacobian || cfg->outputtype==otTaylor)){
            int datalen=(cfg->basisorder) ? mesh->nn : mesh->ne;
            float normalizor=1.f/(DELTA_MUA*cfg->nphoton);
            if(cfg->outputtype==otTaylor)
               normalizor=1.f/cfg->nphoton; /*DELTA_MUA is not used in this mode*/

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<datalen;j++)
                  mesh->weight[i*datalen+j]*=normalizor;
           return normalizor;
        }
	if(cfg->outputtype==otEnergy){
            int datalen=(cfg->basisorder) ? mesh->nn : mesh->ne;
            normalizor=1.f/cfg->nphoton;
            
            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<datalen;j++)
                  mesh->weight[i*datalen+j]*=normalizor;
	    return normalizor;
        }
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

	      energydeposit+=energyelem*mesh->evol[i]*mesh->med[mesh->type[i]].mua; /**mesh->med[mesh->type[i]].n;*/
	    }
	    normalizor=Eabsorb/(Etotal*energydeposit*0.25f); /*scaling factor*/
	    if(cfg->outputtype==otFlux)
               normalizor/=cfg->tstep;

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<mesh->nn;j++)
		  mesh->weight[i*mesh->nn+j]*=normalizor;
	}else{
            for(i=0;i<mesh->ne;i++)
	      for(j=0;j<cfg->maxgate;j++)
	         energydeposit+=mesh->weight[j*mesh->ne+i];

            for(i=0;i<mesh->ne;i++){
	      energyelem=mesh->evol[i]*mesh->med[mesh->type[i]].mua;
              for(j=0;j<cfg->maxgate;j++)
        	mesh->weight[j*mesh->ne+i]/=energyelem;
	    }
            normalizor=Eabsorb/(Etotal*energydeposit); /*scaling factor*/
            if(cfg->outputtype==otFlux)
               normalizor/=cfg->tstep;

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<mesh->ne;j++)
                  mesh->weight[i*mesh->ne+j]*=normalizor;
	}
	return normalizor;
}
