/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmclab.cpp

@brief   mex function for MMCLAB
*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <exception>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "mex.h"
#include "simpmesh.h"
#include "tictoc.h"
#include "tettracing.h"
#include "waitmex/waitmex.c"

//! Macro to read the 1st scalar cfg member
#define GET_1ST_FIELD(x,y)  if(strcmp(name,#y)==0) {double *val=mxGetPr(item);x->y=val[0];printf("mmc.%s=%g;\n",#y,(float)(x->y));}

//! Macro to read one scalar cfg member
#define GET_ONE_FIELD(x,y)  else GET_1ST_FIELD(x,y)

//! Macro to read one 3-element vector member of cfg
#define GET_VEC3_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];\
                                 printf("mmc.%s=[%g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z));}

//! Macro to read one 3- or 4-element vector member of cfg
#define GET_VEC34_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];if(mxGetNumberOfElements(item)==4) u->v.w=val[3];\
                                 printf("mmc.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w));}

//! Macro to read one 4-element vector member of cfg
#define GET_VEC4_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];u->v.w=val[3];\
                                 printf("mmc.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w));}
#define ABS(a)    ((a)<0?-(a):(a))                        //! Macro to calculate the absolute value
#define MAX(a,b)  ((a)>(b)?(a):(b))                       //! Macro to calculate the max of two floating points
#define MEXERROR(a)  mcx_error(999,a,__FILE__,__LINE__)   //! Macro to add unit name and line number in error printing

#if (! defined MX_API_VER) || (MX_API_VER < 0x07300000)
	typedef int dimtype;                              //! MATLAB before 2017 uses int as the dimension array
#else
	typedef size_t dimtype;                           //! MATLAB after 2017 uses size_t as the dimension array
#endif

void mmc_set_field(const mxArray *root,const mxArray *item,int idx, mcconfig *cfg, tetmesh *mesh);
void mmc_validate_config(mcconfig *cfg, tetmesh *mesh);
void mmclab_usage();

/** @brief Mex function for the MMC host function for MATLAB/Octave
 *  This is the master function to interface all MMC features inside MATLAB.
 *  In MMCLAB, all inputs are read from the cfg structure, which contains all
 *  simuation parameters and data.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  mcconfig cfg;
  tetmesh mesh;
  raytracer tracer={NULL,0,NULL,NULL,NULL};
  visitor master={0.f,0.f,0.f,0,0,0,NULL,NULL,0.f,0.f};
  RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
  unsigned int i;
  unsigned int threadid=0,t0,dt;
  float* detmap;

  mxArray    *tmp;
  int        ifield, jstruct;
  int        ncfg, nfields;
  dimtype     fielddim[4];
  int        usewaitbar=1;
  int        errorflag=0;
  const char       *outputtag[]={"data"};
#ifdef MATLAB_MEX_FILE
  waitbar    *hprop;
#endif

  /**
   * If no input is given for this function, it prints help information and return.
   */
  if (nrhs==0){
     mmclab_usage();
     return;
  }

  /**
   * If a structure is passed to this function, a simulation will be launched.
   */
  printf("Launching MMCLAB - Mesh-based Monte Carlo for MATLAB & GNU Octave ...\n");
  if (!mxIsStruct(prhs[0]))
     MEXERROR("Input must be a structure.");

  /**
   * Find out information about input and output.
   */
  nfields = mxGetNumberOfFields(prhs[0]);
  ncfg = mxGetNumberOfElements(prhs[0]);

  /**
   * The function can return 1-3 outputs (i.e. the LHS)
   */
  if(nlhs>=1)
      plhs[0] = mxCreateStructMatrix(ncfg,1,1,outputtag);
  if(nlhs>=2)
      plhs[1] = mxCreateStructMatrix(ncfg,1,1,outputtag);
  if(nlhs>=3)
      plhs[2] = mxCreateStructMatrix(ncfg,1,1,outputtag);

  if(mexEvalString("mmclab_waitbar_handle=figure('visible','off');")) // waitbar is not supported with nojvm after matlab R2013a
      usewaitbar=0;
  else
      mexEvalString("close(mmclab_waitbar_handle);");

  /**
   * Loop over each element of the struct if it is an array, each element is a simulation
   */
  for (jstruct = 0; jstruct < ncfg; jstruct++) {  /* how many configs */

    /** Enclose all simulation calls inside a try/catch construct for exception handling */
    try{
        printf("Running simulations for configuration #%d ...\n", jstruct+1);
        unsigned int ncomplete=0;
        double Eabsorb=0.0;
        float raytri=0.f,raytri0=0.f;;

        /** Initialize cfg with default values first */
	t0=StartTimer();
	mcx_initcfg(&cfg);
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"initializing ... "));
	mesh_init(&mesh);

        memset(&master,0,sizeof(visitor));

        /** Read each struct element from input and set value to the cfg configuration */
	for (ifield = 0; ifield < nfields; ifield++) { /* how many input struct fields */
            tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield);
	    if (tmp == NULL) {
		    continue;
	    }
	    mmc_set_field(prhs[0],tmp,ifield,&cfg,&mesh);
	}
	mexEvalString("pause(.001);");
	
        /** Overwite the output flags using the number of output present */
	cfg.issave2pt=(nlhs>=1);  /** save fluence rate to the 1st output if present */
	cfg.issavedet=(nlhs>=2);  /** save detected photon data to the 2nd output if present */
	cfg.issaveseed=(nlhs>=3); /** save detected photon seeds to the 3rd output if present */

#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
        cfg.issaveseed=0;
#endif
	mesh_srcdetelem(&mesh,&cfg);
	
	/** Validate all input fields, and warn incompatible inputs */
	mmc_validate_config(&cfg,&mesh);

	tracer_init(&tracer,&mesh,cfg.method);
	tracer_prep(&tracer,&cfg);

        dt=GetTimeMillis();
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\nsimulating ... ",dt-t0));

    /** \subsection ssimu Parallel photon transport simulation */

    /** Start multiple threads, one thread to run portion of the simulation on one CUDA GPU, all in parallel */
#pragma omp parallel private(ran0,ran1,threadid) shared(errorflag)
{
	visitor visit={0.f,0.f,1.f/cfg.tstep,DET_PHOTON_BUF,0,0,NULL,NULL,0.f,0.f};
	visit.reclen=(2+((cfg.ismomentum)>0))*mesh.prop+(cfg.issaveexit>0)*6+2;
	if(cfg.issavedet){
            if(cfg.issaveseed)
                visit.photonseed=calloc(visit.detcount,(sizeof(RandType)*RAND_BUF_LEN));
	    visit.partialpath=(float*)calloc(visit.detcount*visit.reclen,sizeof(float));
        }
    #ifdef _OPENMP
	threadid=omp_get_thread_num();	
    #endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg.seed),threadid);
    #ifdef MATLAB_MEX_FILE
        if((cfg.debuglevel & dlProgress) && threadid==0)
	  if(usewaitbar)
             hprop = waitbar_create (0, NULL);
    #endif

	/*launch photons*/
	#pragma omp for reduction(+:Eabsorb) reduction(+:raytri,raytri0)
	for(i=0;i<cfg.nphoton;i++){
            #pragma omp flush (errorflag)
	    if(errorflag){
		i=cfg.nphoton;
		continue;
	    }
            try{
		visit.raytet=0.f;
                visit.raytet0=0.f;
                if(cfg.seed==SEED_FROM_FILE)
                   Eabsorb+=onephoton(i,&tracer,&mesh,&cfg,((RandType *)cfg.photonseed)+i*RAND_BUF_LEN,ran1,&visit);
                else
                   Eabsorb+=onephoton(i,&tracer,&mesh,&cfg,ran0,ran1,&visit);
		raytri+=visit.raytet;
                raytri0+=visit.raytet0;
		#pragma omp atomic
		   ncomplete++;

		if((cfg.debuglevel & dlProgress) && threadid==0 && cfg.nphoton>0){
#ifdef MATLAB_MEX_FILE
                    int prog=ncomplete*100/cfg.nphoton;
                    static int oldprog=-1;
                    char percent[8]="";
                    sprintf(percent,"%d%%",prog);
                    if(prog!=oldprog){
		       if(usewaitbar)
                        waitbar_update (((double)ncomplete)/cfg.nphoton, hprop, percent);
		       else
		        mcx_progressbar(ncomplete,&cfg);
                    }
		    oldprog=prog;
#else
                    mcx_progressbar(ncomplete,&cfg);
#endif
                }
	    }catch(const char *err){
        	mexPrintf("Error from thread (%d): %s\n",threadid,err);
		errorflag++;
	    }catch(const std::exception &err){
        	mexPrintf("C++ Error from thread (%d): %s\n",threadid,err.what());
		errorflag++;
	    }catch(...){
        	mexPrintf("Unknown Exception from thread (%d)",threadid);
		errorflag++;
	    }
	}

        #pragma omp atomic
                master.totalweight += visit.totalweight;

        /** If no error is detected, generat all output data */
	if(cfg.issavedet && errorflag==0){
	    #pragma omp atomic
		master.detcount+=visit.bufpos;
            #pragma omp barrier
            if(master.detcount>0){
	      if(threadid==0){
		if(nlhs>=2){
		    if(cfg.issaveexit==2){
			fielddim[0]=cfg.detparam1.w; fielddim[1]=cfg.detparam2.w; fielddim[2]=cfg.maxgate; fielddim[3]=0;
			mxSetFieldByNumber(plhs[1],jstruct,0, mxCreateNumericArray(3,fielddim,mxSINGLE_CLASS,mxREAL));
			detmap = (float*)mxGetPr(mxGetFieldByNumber(plhs[1],jstruct,0));
			master.partialpath=(float*)calloc(master.detcount*visit.reclen,sizeof(float));
		    }
		    else{
	            	fielddim[0]=visit.reclen; fielddim[1]=master.detcount; fielddim[2]=0; fielddim[3]=0;
    		    	mxSetFieldByNumber(plhs[1],jstruct,0, mxCreateNumericArray(2,fielddim,mxSINGLE_CLASS,mxREAL));
		    	master.partialpath = (float*)mxGetPr(mxGetFieldByNumber(plhs[1],jstruct,0));
		    }
                    if(nlhs>=3 && cfg.issaveseed){
                        fielddim[0]=(sizeof(RandType)*RAND_BUF_LEN); fielddim[1]=master.detcount; fielddim[2]=0; fielddim[3]=0; 
                        mxSetFieldByNumber(plhs[2],jstruct,0, mxCreateNumericArray(2,fielddim,mxUINT8_CLASS,mxREAL));
                        master.photonseed = (unsigned char*)mxGetPr(mxGetFieldByNumber(plhs[2],jstruct,0));
		    }
		}else{
	            master.partialpath=(float*)calloc(master.detcount*visit.reclen,sizeof(float));
                    if(cfg.issaveseed)
                        master.photonseed=calloc(master.detcount,(sizeof(RandType)*RAND_BUF_LEN));
                }
	      }
              #pragma omp barrier
              #pragma omp critical
              {
		memcpy(master.partialpath+master.bufpos*visit.reclen,
		       visit.partialpath,visit.bufpos*visit.reclen*sizeof(float));
                if(nlhs>=3 && cfg.issaveseed)
                    memcpy((unsigned char*)master.photonseed+master.bufpos*(sizeof(RandType)*RAND_BUF_LEN),
                            visit.photonseed,visit.bufpos*(sizeof(RandType)*RAND_BUF_LEN));
		master.bufpos+=visit.bufpos;
              }
            }
            #pragma omp barrier
	    free(visit.partialpath);
            visit.partialpath=NULL;
	}
}

	/** \subsection sreport Post simulation */

#ifdef MATLAB_MEX_FILE
	if((cfg.debuglevel & dlProgress)){
	    if(usewaitbar)
                 waitbar_update (1.0, hprop, NULL);
	    else
		 mcx_progressbar(cfg.nphoton,&cfg);
	}
#endif
	dt=GetTimeMillis()-dt;
	MMCDEBUG(&cfg,dlProgress,(cfg.flog,"\n"));
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",dt));
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"speed ...\t%.2f photon/ms,%.0f ray-tetrahedron tests (%.0f overhead, %.2f test/ms)\n",(double)cfg.nphoton/dt,raytri,raytri0,raytri/dt));

    /** Clear up simulation data structures by calling the destructors */

#ifdef MATLAB_MEX_FILE
        if((cfg.debuglevel & dlProgress))
	   if(usewaitbar)
             waitbar_destroy (hprop) ;
#endif

	tracer_clear(&tracer);
	if(cfg.isnormalized && master.totalweight){
          cfg.his.normalizer=mesh_normalize(&mesh,&cfg,Eabsorb,master.totalweight);
	  printf("total simulated energy: %.0f\tabsorbed: %5.5f%%\tnormalizor=%g\n",
		master.totalweight,100.f*Eabsorb/master.totalweight,cfg.his.normalizer);
	}
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	if(nlhs>=1){
	    int datalen=(cfg.method==rtBLBadouelGrid) ? cfg.crop0.z : ( (cfg.basisorder) ? mesh.nn : mesh.ne);
            fielddim[0]=datalen;
	    fielddim[1]=cfg.maxgate; fielddim[2]=0; fielddim[3]=0; 
	    if(cfg.method==rtBLBadouelGrid){
		fielddim[0]=cfg.dim.x;
		fielddim[1]=cfg.dim.y; 
		fielddim[2]=cfg.dim.z; 
		fielddim[3]=cfg.maxgate;
	        mxSetFieldByNumber(plhs[0],jstruct,0, mxCreateNumericArray(4,fielddim,mxDOUBLE_CLASS,mxREAL));
	    }else{
    	        mxSetFieldByNumber(plhs[0],jstruct,0, mxCreateNumericArray(2,fielddim,mxDOUBLE_CLASS,mxREAL));
	    }
	    double *output = (double*)mxGetPr(mxGetFieldByNumber(plhs[0],jstruct,0));
	    memcpy(output,mesh.weight,datalen*cfg.maxgate*sizeof(double));
	}
	if(nlhs>=2 && cfg.issaveexit==2){
	    float *detimage=(float*)calloc(cfg.detparam1.w*cfg.detparam2.w*cfg.maxgate,sizeof(float));
	    mesh_getdetimage(detimage,master.partialpath,master.bufpos,&cfg,&mesh);
	    memcpy(detmap,detimage,cfg.detparam1.w*cfg.detparam2.w*cfg.maxgate*sizeof(float));
	    free(detimage);	detimage = NULL;
	}
        if(errorflag)
            mexErrMsgTxt("MMCLAB Terminated due to exception!");
    }catch(const char *err){
      mexPrintf("Error: %s\n",err);
    }catch(const std::exception &err){
      mexPrintf("C++ Error: %s\n",err.what());
    }catch(...){
      mexPrintf("Unknown Exception");
    }

    /** \subsection sclean End the simulation */
    mesh_clear(&mesh);
    mcx_clearcfg(&cfg);
  }
  return;
}

/** 
 * @brief Function to parse one subfield of the input structure
 *
 * This function reads in all necessary information from the cfg input structure.
 * it can handle single scalar inputs, short vectors (3-4 elem), strings and arrays.
 *
 * @param[in] root: the cfg input data structure
 * @param[in] item: the current element of the cfg input data structure
 * @param[in] idx: the index of the current element (starting from 0)
 * @param[out] cfg: the simulation configuration structure to store all input read from the parameters
 * @param[out] mesh: the mesh data structure
 */

void mmc_set_field(const mxArray *root,const mxArray *item,int idx, mcconfig *cfg, tetmesh *mesh){
    const char *name=mxGetFieldNameByNumber(root,idx);
    const dimtype *arraydim;
    char *jsonshapes=NULL;
    int i,j;
    dimtype k;

    if(strcmp(name,"nphoton")==0 && cfg->photonseed!=NULL)
	return;

    cfg->flog=stderr;
    GET_1ST_FIELD(cfg,nphoton)
    GET_ONE_FIELD(cfg,tstart)
    GET_ONE_FIELD(cfg,tstep)
    GET_ONE_FIELD(cfg,tend)
    GET_ONE_FIELD(cfg,isreflect)
    GET_ONE_FIELD(cfg,isspecular)
    GET_ONE_FIELD(cfg,ismomentum)
    GET_ONE_FIELD(cfg,issaveexit)
    GET_ONE_FIELD(cfg,issaveseed)
    GET_ONE_FIELD(cfg,isatomic)
    GET_ONE_FIELD(cfg,basisorder)
    GET_ONE_FIELD(cfg,outputformat)
    GET_ONE_FIELD(cfg,roulettesize)
    GET_ONE_FIELD(cfg,nout)
    GET_ONE_FIELD(cfg,isref3)
    GET_ONE_FIELD(cfg,isnormalized)
    GET_ONE_FIELD(cfg,debugphoton)
    GET_ONE_FIELD(cfg,minenergy)
    GET_ONE_FIELD(cfg,replaydet)
    GET_ONE_FIELD(cfg,unitinmm)
    GET_ONE_FIELD(cfg,voidtime)
    GET_ONE_FIELD(cfg,mcmethod)
    GET_VEC3_FIELD(cfg,srcpos)
    GET_VEC34_FIELD(cfg,srcdir)
    GET_VEC3_FIELD(cfg,steps)
    GET_VEC4_FIELD(cfg,srcparam1)
    GET_VEC4_FIELD(cfg,srcparam2)
    GET_VEC4_FIELD(cfg,detparam1)
    GET_VEC4_FIELD(cfg,detparam2)
    else if(strcmp(name,"e0")==0){
        double *val=mxGetPr(item);
	cfg->e0=val[0];
        printf("mmc.e0=%d;\n",cfg->e0);
    }else if(strcmp(name,"node")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]!=3)
            MEXERROR("the 'node' field must have 3 columns (x,y,z)");
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
	if(arraydim[0]<=0 || arraydim[1]<4)
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        double *val=mxGetPr(item);
        mesh->ne=arraydim[0];
	mesh->elemlen=arraydim[1];
	if(mesh->elem) free(mesh->elem);
        mesh->elem=(int *)calloc(sizeof(int)*arraydim[1],mesh->ne);
        for(j=0;j<mesh->elemlen;j++)
          for(i=0;i<mesh->ne;i++)
             mesh->elem[i*mesh->elemlen+j]=val[j*mesh->ne+i];
        printf("mmc.elem=[%d,%d];\n",mesh->ne,mesh->elemlen);
    }else if(strcmp(name,"elemprop")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'elemprop' field can not be empty");
        double *val=mxGetPr(item);
        mesh->ne=MAX(arraydim[0],arraydim[1]);
	if(mesh->type) free(mesh->type);
	mesh->type=(int  *)malloc(sizeof(int )*mesh->ne);
        for(i=0;i<mesh->ne;i++)
           mesh->type[i]=val[i];
        printf("mmc.ne=%d;\n",mesh->ne);
    }else if(strcmp(name,"facenb")==0){
        arraydim=mxGetDimensions(item);
	if(arraydim[0]<=0 || arraydim[1]<4)
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        double *val=mxGetPr(item);
        mesh->ne=arraydim[0];
	mesh->elemlen=arraydim[1];
	if(mesh->facenb) free(mesh->facenb);
        mesh->facenb=(int *)malloc(sizeof(int)*arraydim[1]*mesh->ne);
        for(j=0;j<arraydim[1];j++)
          for(i=0;i<mesh->ne;i++)
             mesh->facenb[i*arraydim[1]+j]=val[j*mesh->ne+i];
        printf("mmc.facenb=[%d,%d];\n",mesh->ne,mesh->elemlen);
    }else if(strcmp(name,"evol")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'evol' field can not be empty");
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
            MEXERROR("the 'detpos' field must have 4 columns (x,y,z,radius)");
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
            MEXERROR("the 'prop' field must have 4 columns (mua,mus,g,n)");
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
             MEXERROR("the 'debuglevel' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     MEXERROR("the 'debuglevel' field is too long");
        int status = mxGetString(item, buf, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->debuglevel=mcx_parsedebugopt(buf);
	printf("mmc.debuglevel='%s';\n",buf);
    }else if(strcmp(name,"srctype")==0){
        int len=mxGetNumberOfElements(item);
        const char *srctypeid[]={"pencil","isotropic","cone","gaussian","planar","pattern","fourier","arcsine","disk","fourierx","fourierx2d","zgaussian","line","slit",""};
        char strtypestr[MAX_SESSION_LENGTH]={'\0'};
        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'srctype' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'srctype' field is too long");
        int status = mxGetString(item, strtypestr, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->srctype=mcx_keylookup(strtypestr,srctypeid);
        if(cfg->srctype==-1)
             mexErrMsgTxt("the specified source type is not supported");
	printf("mmc.srctype='%s';\n",strtypestr);
    }else if(strcmp(name,"session")==0){
        int len=mxGetNumberOfElements(item);
        if(!mxIsChar(item) || len==0)
             MEXERROR("the 'session' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     MEXERROR("the 'session' field is too long");
        int status = mxGetString(item, cfg->session, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");

	printf("mmc.session='%s';\n",cfg->session);
    }else if(strcmp(name,"srcpattern")==0){
        arraydim=mxGetDimensions(item);
        double *val=mxGetPr(item);
        if(cfg->srcpattern) free(cfg->srcpattern);
        cfg->srcpattern=(float*)malloc(arraydim[0]*arraydim[1]*sizeof(float));
        for(k=0;k<arraydim[0]*arraydim[1];k++)
             cfg->srcpattern[k]=val[k];
        printf("mmc.srcpattern=[%d %d];\n",arraydim[0],arraydim[1]);
    }else if(strcmp(name,"method")==0){
        int len=mxGetNumberOfElements(item);
        const char *methods[]={"plucker","havel","badouel","elem","grid",""};
        char methodstr[MAX_SESSION_LENGTH]={'\0'};

        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'method' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'method' field is too long");
        int status = mxGetString(item, methodstr, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->method=mcx_keylookup(methodstr,methods);
        if(cfg->method==-1)
             mexErrMsgTxt("the specified method is not supported");
	printf("mmc.method='%s';\n",methodstr);
    }else if(strcmp(name,"outputtype")==0){
        int len=mxGetNumberOfElements(item);
        const char *outputtype[]={"flux","fluence","energy","jacobian","wl","wp",""};
        char outputstr[MAX_SESSION_LENGTH]={'\0'};

        if(!mxIsChar(item) || len==0)
             mexErrMsgTxt("the 'outputtype' field must be a non-empty string");
	if(len>MAX_SESSION_LENGTH)
	     mexErrMsgTxt("the 'outputtype' field is too long");
        int status = mxGetString(item, outputstr, MAX_SESSION_LENGTH);
        if (status != 0)
             mexWarnMsgTxt("not enough space. string is truncated.");
        cfg->outputtype=mcx_keylookup(outputstr,outputtype);
        if(cfg->outputtype==-1)
             mexErrMsgTxt("the specified output type is not supported");
	printf("mmc.outputtype='%s';\n",outputstr);
    }else if(strcmp(name,"shapes")==0){
        int len=mxGetNumberOfElements(item);
        if(!mxIsChar(item) || len==0)
             MEXERROR("the 'shapes' field must be a non-empty string");

        jsonshapes=new char[len+1];
        mxGetString(item, jsonshapes, len+1);
        jsonshapes[len]='\0';
    }else if(strcmp(name,"seed")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'seed' field can not be empty");
        if(!mxIsUint8(item)){
            double *val=mxGetPr(item);
            cfg->seed=val[0];
            printf("mmc.seed=%d;\n",cfg->seed);
        }else{
	    cfg->photonseed=malloc(arraydim[0]*arraydim[1]);
            if(arraydim[0]!=(sizeof(RandType)*RAND_BUF_LEN))
		MEXERROR("the row number of cfg.seed does not match RNG seed byte-length");
            memcpy(cfg->photonseed,mxGetData(item),arraydim[0]*arraydim[1]);
            cfg->seed=SEED_FROM_FILE;
            cfg->nphoton=arraydim[1];
            printf("mmc.nphoton=%d;\n",cfg->nphoton);
	}
    }else if(strcmp(name,"replayweight")==0){
        arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'replayweight' field can not be empty");
	cfg->his.detected=arraydim[0]*arraydim[1];
	cfg->replayweight=(float *)malloc(cfg->his.detected*sizeof(float));
        memcpy(cfg->replayweight,mxGetData(item),cfg->his.detected*sizeof(float));
        printf("mmc.replayweight=%d;\n",cfg->his.detected);
    }else if(strcmp(name,"replaytime")==0){
	arraydim=mxGetDimensions(item);
	if(MAX(arraydim[0],arraydim[1])==0)
            MEXERROR("the 'replaytime' field can not be empty");
	cfg->his.detected=arraydim[0]*arraydim[1];
	cfg->replaytime=(float *)malloc(cfg->his.detected*sizeof(float));
        memcpy(cfg->replaytime,mxGetData(item),cfg->his.detected*sizeof(float));
        printf("mmc.replaytime=%d;\n",cfg->his.detected);
    }else if(strcmp(name,"isreoriented")==0){
        /*internal flag, don't need to do anything*/
    }else{
        printf("WARNING: redundant field '%s'\n",name);
    }
}

/** 
 * @brief Validate all input fields, and warn incompatible inputs
 *
 * Perform self-checking and raise exceptions or warnings when input error is detected
 *
 * @param[in,out] cfg: the simulation configuration structure
 * @param[out] mesh: the mesh data structure
 */

void mmc_validate_config(mcconfig *cfg, tetmesh *mesh){
     int i,j,*ee,datalen;
     if(cfg->nphoton<=0){
         MEXERROR("cfg.nphoton must be a positive number");
     }
     if(cfg->tstart>cfg->tend || cfg->tstep==0.f){
         MEXERROR("incorrect time gate settings or missing tstart/tend/tstep fields");
     }
     if(cfg->tstep>cfg->tend-cfg->tstart){
         cfg->tstep=cfg->tend-cfg->tstart;
     }
     if(ABS(cfg->srcdir.x*cfg->srcdir.x+cfg->srcdir.y*cfg->srcdir.y+cfg->srcdir.z*cfg->srcdir.z - 1.f)>1e-5)
         MEXERROR("field 'srcdir' must be a unitary vector");
     if(cfg->tend<=cfg->tstart)
         MEXERROR("field 'tend' must be greater than field 'tstart'");
     cfg->maxgate=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);

     if(mesh->prop==0)
        MEXERROR("you must define the 'prop' field in the input structure");
     if(mesh->nn==0||mesh->ne==0||mesh->evol==NULL || mesh->facenb==NULL)
        MEXERROR("a complete input mesh include 'node','elem','facenb' and 'evol'");

     mesh->nvol=(float *)calloc(sizeof(float),mesh->nn);
     for(i=0;i<mesh->ne;i++){
        if(mesh->type[i]<=0)
		continue;
     	ee=(int *)(mesh->elem+i*mesh->elemlen);
     	for(j=0;j<4;j++)
     	   	mesh->nvol[ee[j]-1]+=mesh->evol[i]*0.25f;
     }
     if(mesh->weight)
        free(mesh->weight);

     if(cfg->method==rtBLBadouelGrid){
	mesh_createdualmesh(mesh,cfg);
	cfg->basisorder=0;
     }
     datalen=(cfg->method==rtBLBadouelGrid) ? cfg->crop0.z : ( (cfg->basisorder) ? mesh->nn : mesh->ne);
     mesh->weight=(double *)calloc(sizeof(double)*datalen,cfg->maxgate);

     if(cfg->srctype==stPattern && cfg->srcpattern==NULL)
        mexErrMsgTxt("the 'srcpattern' field can not be empty when your 'srctype' is 'pattern'");

     if(cfg->method!=rtBLBadouelGrid && cfg->unitinmm!=1.f){
        for(i=1;i<mesh->prop;i++){
		mesh->med[i].mus*=cfg->unitinmm;
		mesh->med[i].mua*=cfg->unitinmm;
        }
     }
     cfg->his.unitinmm=cfg->unitinmm;
     if(mesh->node==NULL || mesh->elem==NULL || mesh->prop==0){
	 MEXERROR("You must define 'mesh' and 'prop' fields.");
     }
     /*make medianum+1 the same as medium 0*/
     if(cfg->isextdet){
         mesh->med=(medium *)realloc(mesh->med, sizeof(medium)*(mesh->prop+2));
         memcpy(mesh->med+mesh->prop+1,mesh->med,sizeof(medium));
         for(i=0;i<mesh->ne;i++){
             if(mesh->type[i]==-2)
                   mesh->type[i]=mesh->prop+1;
         }
     }

     if(cfg->issavedet && cfg->detnum==0 && cfg->isextdet==0) 
      	cfg->issavedet=0;
     if(cfg->seed<0 && cfg->seed!=SEED_FROM_FILE) cfg->seed=time(NULL);
     if(cfg->issavedet==0){
        cfg->ismomentum=0;
        cfg->issaveexit=0;
     }
     if(cfg->seed==SEED_FROM_FILE && cfg->his.detected!=cfg->nphoton){
        cfg->his.detected=0;
	if(cfg->replayweight==NULL)
	    MEXERROR("You must define 'replayweight' when you specify 'seed'.");
	else if(cfg->replaytime==NULL)
	    MEXERROR("You must define 'replayweight' when you specify 'seed'.");
	else
	    MEXERROR("The dimension of the 'replayweight' OR 'replaytime' field does not match the column number of the 'seed' field.");
     }
     // cfg->his.maxmedia=cfg->medianum-1; /*skip medium 0*/
     cfg->his.detnum=cfg->detnum;
     cfg->his.colcount=(1+(cfg->ismomentum>0))*cfg->his.maxmedia+(cfg->issaveexit>0)*6+1;
}

/**
 * @brief Error reporting function in the mex function, equivallent to mcx_error in binary mode
 *
 * @param[in] id: a single integer for the types of the error
 * @param[in] msg: the error message string
 * @param[in] filename: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

extern "C" int mmc_throw_exception(const int id, const char *msg, const char *filename, const int linenum){
     printf("MMCLAB ERROR (%d): %s in unit %s:%d\n",id,msg,filename,linenum);
     throw(msg);
     return id;
}

/**
 * @brief Print a brief help information if nothing is provided
 */

void mmclab_usage(){
     printf("Usage:\n    [flux,detphoton]=mmclab(cfg);\n\nPlease run 'help mmclab' for more details.\n");
}

/**
 * @brief Force matlab refresh the command window to print all buffered messages
 */

extern "C" void mcx_matlab_flush(){
#ifndef MATLAB_MEX_FILE
	mexEvalString("fflush(stdout);");
#else
	mexEvalString("pause(.0001);");
#endif
}
