/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mex.h"
#include "simpmesh.h"
#include "tictoc.h"
#include "tettracing.h"

#define GET_1ST_FIELD(x,y)  if(strcmp(name,#y)==0) {double *val=mxGetPr(item);x->y=val[0];printf("mmc.%s=%g;\n",#y,(float)(x->y));}
#define GET_ONE_FIELD(x,y)  else GET_1ST_FIELD(x,y)
#define GET_VEC3_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];\
                                 printf("mmc.%s=[%g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z));}
#define ABS(a)    ((a)<0?-(a):(a))
#define MAX(a,b)  ((a)>(b)?(a):(b))

void mmc_set_field(const mxArray *root,const mxArray *item,int idx, mcconfig *cfg, tetmesh *mesh);
void mmc_validate_config(mcconfig *cfg, tetmesh *mesh);
void mmclab_usage();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  mcconfig cfg;
  tetmesh mesh;
  raytracer tracer;
  visitor master={0.f,0.f,0,0,NULL};
  double Eabsorb=0.0;
  RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  unsigned int i;
  float raytri=0.f;
  unsigned int threadid=0,ncomplete=0,t0,dt;

  mxArray    *tmp;
  int        ifield, jstruct;
  int        ncfg, nfields;
  int        fielddim[4];
  const char       *outputtag[]={"data"};

  if (nrhs==0){
     mmclab_usage();
     return;
  }
  printf("Launching MMCLAB - Mesh-based Monte Carlo for MATLAB & GNU Octave ...\n");
  if (!mxIsStruct(prhs[0]))
     mexErrMsgTxt("Input must be a structure.");

  nfields = mxGetNumberOfFields(prhs[0]);
  ncfg = mxGetNumberOfElements(prhs[0]);

  if(nlhs>=1)
      plhs[0] = mxCreateStructMatrix(ncfg,1,1,outputtag);
  if(nlhs>=2)
      plhs[1] = mxCreateStructMatrix(ncfg,1,1,outputtag);

  for (jstruct = 0; jstruct < ncfg; jstruct++) {  /* how many configs */
    printf("Running simulations for configuration #%d ...\n", jstruct+1);

    t0=StartTimer();
    mcx_initcfg(&cfg);
    MMCDEBUG(&cfg,dlTime,(cfg.flog,"initializing ... "));
    mesh_init(&mesh);

    for (ifield = 0; ifield < nfields; ifield++) { /* how many input struct fields */
        tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield);
	if (tmp == NULL) {
		continue;
	}
	mmc_set_field(prhs[0],tmp,ifield,&cfg,&mesh);
    }
    mexEvalString("pause(.001);");
    cfg.issave2pt=(nlhs>=1);
    cfg.issavedet=(nlhs>=2);

    mmc_validate_config(&cfg,&mesh);

    tracer_init(&tracer,&mesh,cfg.method);
    tracer_prep(&tracer,&cfg);

    MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\nsimulating ... ",GetTimeMillis()-t0));

/** \subsection ssimu Parallel photon transport simulation */

#pragma omp parallel private(ran0,ran1,threadid)
{
    visitor visit={0.f,1.f/cfg.tstep,DET_PHOTON_BUF,0,NULL};
    int buflen=(1+((cfg.ismomentum)>0))*mesh.prop+2;
    if(cfg.issavedet) 
	visit.partialpath=(float*)calloc(visit.detcount*buflen,sizeof(float));
#ifdef _OPENMP
    threadid=omp_get_thread_num();	
#endif
    rng_init(ran0,ran1,(unsigned int *)&(cfg.seed),threadid);

    /*launch photons*/
    #pragma omp for reduction(+:Eabsorb) reduction(+:raytri)
    for(i=0;i<cfg.nphoton;i++){
	    visit.raytet=0.f;
	    Eabsorb+=onephoton(i,&tracer,&mesh,&cfg,ran0,ran1,&visit);
	    raytri+=visit.raytet;
	    #pragma omp atomic
	       ncomplete++;

	    if((cfg.debuglevel & dlProgress) && threadid==0)
		    mcx_progressbar(ncomplete,&cfg);
	    if(cfg.issave2pt && cfg.checkpt[0])
		    mesh_saveweightat(&mesh,&cfg,i+1);
    }
    if(cfg.issavedet){
	#pragma omp atomic
	    master.detcount+=visit.bufpos;
        #pragma omp barrier
	if(threadid==0){
	    if(nlhs>=2){
	        fielddim[0]=master.detcount*buflen; fielddim[1]=0; fielddim[2]=0; fielddim[3]=0; 
    		mxSetFieldByNumber(plhs[1],jstruct,0, mxCreateNumericArray(4,fielddim,mxSINGLE_CLASS,mxREAL));
		master.partialpath = (float*)mxGetPr(mxGetFieldByNumber(plhs[1],jstruct,0));
	    }else
	        master.partialpath=(float*)calloc(master.detcount*buflen,sizeof(float));
	}
        #pragma omp barrier
        #pragma omp critical
        {
	    memcpy(master.partialpath+master.bufpos*buflen,
		   visit.partialpath,visit.bufpos*buflen*sizeof(float));
	    master.bufpos+=visit.bufpos;
        }
        #pragma omp barrier
	free(visit.partialpath);
    }
}

    /** \subsection sreport Post simulation */

    if((cfg.debuglevel & dlProgress))
	    mcx_progressbar(cfg.nphoton,&cfg);

    dt=GetTimeMillis()-t0;
    MMCDEBUG(&cfg,dlProgress,(cfg.flog,"\n"));
    MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",dt));
    MMCDEBUG(&cfg,dlTime,(cfg.flog,"speed ...\t%.0f ray-tetrahedron tests\n",raytri));

    tracer_clear(&tracer);
    if(cfg.isnormalized && cfg.nphoton){
      fprintf(cfg.flog,"total simulated energy: %d\tabsorbed: %5.5f%%\tnormalizor=%g\n",
	    cfg.nphoton,100.f*Eabsorb/cfg.nphoton,mesh_normalize(&mesh,&cfg,Eabsorb,cfg.nphoton));
    }
    if(cfg.issavedet && master.partialpath)
	    free(master.partialpath);

    MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

    if(nlhs>=1){
    	if(!cfg.basisorder)
	    fielddim[0]=mesh.ne*cfg.maxgate; 
	else
	    fielddim[0]=mesh.nn*cfg.maxgate; 
	fielddim[1]=0; fielddim[2]=0; fielddim[3]=0; 
    	mxSetFieldByNumber(plhs[0],jstruct,0, mxCreateNumericArray(1,fielddim,mxDOUBLE_CLASS,mxREAL));
	double *output = (double*)mxGetPr(mxGetFieldByNumber(plhs[0],jstruct,0));
	memcpy(output,mesh.weight,fielddim[0]*sizeof(double));
    }
    if(nlhs>=2){
        fielddim[0]=(mesh.prop+1); fielddim[1]=cfg.his.savedphoton;
	fielddim[2]=0; fielddim[3]=0;
	if(cfg.his.savedphoton>0){
		mxSetFieldByNumber(plhs[1],jstruct,0, mxCreateNumericArray(2,fielddim,mxSINGLE_CLASS,mxREAL));
		memcpy((float*)mxGetPr(mxGetFieldByNumber(plhs[1],jstruct,0)),master.partialpath,
		     fielddim[0]*fielddim[1]*sizeof(float));
	}
    }

    /** \subsection sclean End the simulation */
    mesh_clear(&mesh);
    mcx_clearcfg(&cfg);
  }
  return;
}


void mmc_set_field(const mxArray *root,const mxArray *item,int idx, mcconfig *cfg, tetmesh *mesh){
    const char *name=mxGetFieldNameByNumber(root,idx);
    const int *arraydim;
    char *jsonshapes=NULL;
    int i,j;
    const char *srctypeid[]={"pencil","isotropic","cone","gaussian",""};

    cfg->flog=stderr;
    GET_1ST_FIELD(cfg,nphoton)
    GET_ONE_FIELD(cfg,seed)
    GET_ONE_FIELD(cfg,tstart)
    GET_ONE_FIELD(cfg,tstep)
    GET_ONE_FIELD(cfg,tend)
    GET_ONE_FIELD(cfg,isreflect)
    GET_ONE_FIELD(cfg,isspecular)
    GET_ONE_FIELD(cfg,outputtype)
    GET_ONE_FIELD(cfg,basisorder)
    GET_ONE_FIELD(cfg,outputformat)
    GET_ONE_FIELD(cfg,method)
    GET_ONE_FIELD(cfg,roulettesize)
    GET_ONE_FIELD(cfg,nout)
    GET_ONE_FIELD(cfg,isref3)
    GET_ONE_FIELD(cfg,isnormalized)
    GET_ONE_FIELD(cfg,minenergy)
    GET_ONE_FIELD(cfg,unitinmm)
    GET_VEC3_FIELD(cfg,srcpos)
    GET_VEC3_FIELD(cfg,srcdir)
    GET_VEC3_FIELD(cfg,steps)
    else if(strcmp(name,"e0")==0){
        double *val=mxGetPr(item);
	cfg->dim.x=val[0];
        printf("mmc.e0=%d;\n",cfg->dim.x);
    }else if(strcmp(name,"node")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=3)
            mexErrMsgTxt("the 'node' field must have 3 columns (x,y,z)");
        double *val=mxGetPr(item);
        mesh->nn=arraydim[0];
	if(mesh->node) free(mesh->node);
        mesh->node=(float3 *)calloc(sizeof(float3),mesh->nn);
        for(j=0;j<3;j++)
          for(i=0;i<mesh->nn;i++)
             ((float *)(&mesh->node[i]))[j]=val[j*mesh->nn+i];
        printf("mmc.nn=%d;\n",mesh->nn);
    }else if(strcmp(name,"elem")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=4)
            mexErrMsgTxt("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        double *val=mxGetPr(item);
        mesh->ne=arraydim[0];
	if(mesh->elem) free(mesh->elem);
        mesh->elem=(int4 *)calloc(sizeof(int4),mesh->ne);
        for(j=0;j<4;j++)
          for(i=0;i<mesh->ne;i++)
             ((int *)(&mesh->elem[i]))[j]=val[j*mesh->ne+i];
        printf("mmc.ne=%d;\n",mesh->ne);
    }else if(strcmp(name,"elemprop")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            mexErrMsgTxt("the 'elemprop' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ne=MAX(arraydim[0],arraydim[1]);
	if(mesh->type) free(mesh->type);
	mesh->type=(int  *)malloc(sizeof(int )*mesh->ne);
        for(i=0;i<mesh->ne;i++)
           mesh->type[i]=val[i];
        printf("mmc.ne=%d;\n",mesh->ne);
    }else if(strcmp(name,"facenb")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=4)
            mexErrMsgTxt("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        double *val=mxGetPr(item);
        mesh->ne=arraydim[0];
	if(mesh->facenb) free(mesh->facenb);
        mesh->facenb=(int4 *)malloc(sizeof(int4)*mesh->ne);
        for(j=0;j<4;j++)
          for(i=0;i<mesh->ne;i++)
             ((int *)(&mesh->facenb[i]))[j]=val[j*mesh->ne+i];
        printf("mmc.facenb=%d;\n",mesh->ne);
    }else if(strcmp(name,"evol")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            mexErrMsgTxt("the 'evol' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ne=MAX(arraydim[0],arraydim[1]);
	if(mesh->evol) free(mesh->evol);
        mesh->evol=(float *)malloc(sizeof(float)*mesh->ne);
        for(i=0;i<mesh->ne;i++)
           mesh->evol[i]=val[i];
        printf("mmc.evol=%d;\n",mesh->ne);
    }else if(strcmp(name,"detpos")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]>0 && arraydim[1]!=4)
            mexErrMsgTxt("the 'detpos' field must have 4 columns (x,y,z,radius)");
        double *val=mxGetPr(item);
        cfg->detnum=arraydim[0];
	if(cfg->detpos) free(cfg->detpos);
        cfg->detpos=(float4 *)malloc(cfg->detnum*sizeof(float4));
        for(j=0;j<4;j++)
          for(i=0;i<cfg->detnum;i++)
             ((float *)(&cfg->detpos[i]))[j]=val[j*cfg->detnum+i];
        printf("mmc.detnum=%d;\n",cfg->detnum);
    }else if(strcmp(name,"prop")==0){
        arraydim=mxGetDimensions(item);
        if(arraydim[0]>0 && arraydim[1]!=4)
            mexErrMsgTxt("the 'prop' field must have 4 columns (mua,mus,g,n)");
        double *val=mxGetPr(item);
        mesh->prop=arraydim[0]-1;
        if(mesh->med) free(mesh->med);
	if(mesh->atte) free(mesh->atte);
        mesh->med=(medium *)calloc(sizeof(medium),mesh->prop+1);
	mesh->atte=(float *)calloc(sizeof(float),mesh->prop+1);
        for(j=0;j<4;j++)
          for(i=0;i<=mesh->prop;i++)
             ((float *)(&mesh->med[i]))[j]=val[j*(mesh->prop+1)+i];
	/*for(i=0;i<=mesh->prop;i++)
             mesh->atte[i]=expf(-cfg->minstep*mesh->med[i].mua);*/
	cfg->his.maxmedia=mesh->prop;
        printf("mmc.prop=%d;\n",mesh->prop);
    }else if(strcmp(name,"debuglevel")==0){
        int len=mxGetNumberOfElements(item);
	char buf[MAX_SESSION_LENGTH];
        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'debuglevel' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'session' field is too long");
        int status = mxGetString(item, buf, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->debuglevel=mcx_parsedebugopt(buf);
	printf("mmc.session='%s';\n",buf);
    }else if(strcmp(name,"srctype")==0){
        int len=mxGetNumberOfElements(item);
	char buf[MAX_SESSION_LENGTH];
        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'srctype' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'srctype' field is too long");
        int status = mxGetString(item, buf, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->srctype=mcx_keylookup(buf,srctypeid);
	printf("mmc.srctype='%s';\n",buf);
    }else if(strcmp(name,"session")==0){
        int len=mxGetNumberOfElements(item);
        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'session' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'session' field is too long");
        int status = mxGetString(item, cfg->session, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");

	printf("mmc.session='%s';\n",cfg->session);
    }else if(strcmp(name,"shapes")==0){
        int len=mxGetNumberOfElements(item);
        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'shapes' field must be a non-empty string");

        jsonshapes=new char[len+1];
        mxGetString(item, jsonshapes, len+1);
        jsonshapes[len]='\0';
    }else{
        printf("WARNING: redundant field '%s'\n",name);
    }
}

void mmc_validate_config(mcconfig *cfg, tetmesh *mesh){
     int i,j,gates,idx1d,*ee;
     if(cfg->nphoton<=0){
         mexErrMsgTxt("cfg.nphoton must be defined as a positive number");
     }
     if(cfg->tstart>cfg->tend || cfg->tstep==0.f){
         mexErrMsgTxt("incorrect time gate settings");
     }
     if(ABS(cfg->srcdir.x*cfg->srcdir.x+cfg->srcdir.y*cfg->srcdir.y+cfg->srcdir.z*cfg->srcdir.z - 1.f)>1e-5)
         mexErrMsgTxt("field 'srcdir' must be a unitary vector");
     if(cfg->tend<=cfg->tstart)
         mexErrMsgTxt("field 'tend' must be greater than field 'tstart'");
     gates=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);
     if(mesh->prop==0)
        mexErrMsgTxt("you must define the 'prop' field in the input structure");
     if(mesh->nn==0||mesh->ne==0||mesh->evol==NULL || mesh->facenb==NULL)
        mexErrMsgTxt("a complete input mesh include 'node','elem','facenb' and 'evol'");

     mesh->nvol=(float *)calloc(sizeof(float),mesh->nn);
     for(i=0;i<mesh->ne;i++){
     	ee=(int *)(mesh->elem+i);
     	for(j=0;j<4;j++)
     	   	mesh->nvol[ee[j]-1]+=mesh->evol[i]*0.25f;
     }
     if(mesh->weight)
        free(mesh->weight);

     if(!cfg->basisorder)
        mesh->weight=(double*)calloc(mesh->ne*sizeof(double),cfg->maxgate);
     else
        mesh->weight=(double*)calloc(mesh->nn*sizeof(double),cfg->maxgate);

     if(cfg->unitinmm!=1.f){
	float unit3d=cfg->unitinmm*cfg->unitinmm*cfg->unitinmm;
	for(i=0;i<mesh->ne;i++)
	      mesh->evol[i]*=unit3d;
	for(i=0;i<mesh->nn;i++)
	      mesh->nvol[i]*=unit3d;
        for(i=1;i<mesh->prop;i++){
		mesh->med[i].mus*=cfg->unitinmm;
		mesh->med[i].mua*=cfg->unitinmm;
        }
	for(i=0;i<mesh->nn;i++){
	      mesh->node[i].x*=cfg->unitinmm;
	      mesh->node[i].y*=cfg->unitinmm;
	      mesh->node[i].z*=cfg->unitinmm;
	}
     }
     if(mesh->node==NULL || mesh->elem==NULL || mesh->prop==0){
	 mexErrMsgTxt("You must define 'mesh' and 'prop' fields.");
     }
     if(cfg->issavedet && cfg->detnum==0) 
      	cfg->issavedet=0;

     if(cfg->seed<0) cfg->seed=time(NULL);

     cfg->his.maxmedia=cfg->medianum-1; /*skip medium 0*/
     cfg->his.detnum=cfg->detnum;
     cfg->his.colcount=cfg->medianum+1; /*column count=maxmedia+2*/
}

extern "C" int mmc_throw_exception(const int id, const char *msg, const char *filename, const int linenum){
     printf("MMCLAB ERROR %d in unit %s:%d\n",id,filename,linenum);
     mexErrMsgTxt(msg);
     return id;
}

void mmclab_usage(){
     printf("Usage:\n    [flux,detphoton]=mmclab(cfg);\n\nPlease run 'help mmclab' for more details.\n");
}
