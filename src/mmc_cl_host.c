/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Pluker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmc_host.c

\brief   << Driver program of MMC >>
*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mmc_cl_host.h"
#include "mcx_const.h"
#include "tictoc.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the 
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

const char *sourceflag[]={"-DMCX_SRC_PENCIL","-DMCX_SRC_ISOTROPIC","-DMCX_SRC_CONE",
    "-DMCX_SRC_GAUSSIAN","-DMCX_SRC_PLANAR","-DMCX_SRC_PATTERN","-DMCX_SRC_FOURIER",
    "-DMCX_SRC_ARCSINE","-DMCX_SRC_DISK","-DMCX_SRC_FOURIERX","-DMCX_SRC_FOURIERX2D",
    "-DMCX_SRC_ZGAUSSIAN","-DMCX_SRC_LINE","-DMCX_SRC_SLIT","-DMCX_SRC_PENCILARRAY",
    "-DMCX_SRC_PATTERN3D"};

extern cl_event kernelevent;

/*
   master driver code to run MC simulations
*/

void mmc_run_cl(mcconfig *cfg,tetmesh *mesh, raytracer *tracer){

     cl_uint i,j,iter;
     cl_float t,twindow0,twindow1;
     cl_float fullload=0.f;
     cl_float *energy;
     cl_uint *progress=NULL;
     cl_uint detected=0,workdev;

     cl_uint tic,tic0,tic1,toc=0,fieldlen;

     cl_context mcxcontext;                 // compute mcxcontext
     cl_command_queue *mcxqueue;          // compute command queue
     cl_program mcxprogram;                 // compute mcxprogram
     cl_kernel *mcxkernel;                   // compute mcxkernel
     cl_int status = 0;
     cl_device_id devices[MAX_DEVICE];
     cl_event * waittoread;
     cl_platform_id platform = NULL;

     cl_uint  totalcucore;
     cl_uint  devid=0;
     cl_mem gnode,gelem,gtype,gfacenb,gsrcelem,gnormal,gproperty,gparam,gdetpos; /*read-only buffers*/
     cl_mem *gweight,*gdref,*gdetphoton,*gseed,*genergy,*greporter;          /*read-write buffers*/
     cl_mem *gprogress=NULL,*gdetected, *gsrcpattern;  /*read-write buffers*/

     cl_uint meshlen=((cfg->method==rtBLBadouelGrid) ? cfg->crop0.z : mesh->ne)<<cfg->nbuffer; // use 4 copies to reduce racing
     
     cl_float  *field,*dref=NULL;

     cl_uint   *Pseed;
     float  *Pdet;
     char opt[MAX_PATH_LENGTH]={'\0'};
     cl_uint detreclen=(2+((cfg->ismomentum)>0))*mesh->prop+(cfg->issaveexit>0)*6+1;
     cl_uint hostdetreclen=detreclen+1;
     GPUInfo *gpu=NULL;

     MCXParam param={{{cfg->srcpos.x,cfg->srcpos.y,cfg->srcpos.z}}, {{cfg->srcdir.x,cfg->srcdir.y,cfg->srcdir.z}},
		     cfg->tstart, cfg->tend, (uint)cfg->isreflect,(uint)cfg->issavedet,(uint)cfg->issaveexit,
		     (uint)cfg->ismomentum, (uint)cfg->isatomic, (uint)cfg->isspecular, 1.f/cfg->tstep, cfg->minenergy, 
		     cfg->maxdetphoton, mesh->prop, cfg->detnum, (uint)cfg->voidtime, (uint)cfg->srctype, 
		     {{cfg->srcparam1.x,cfg->srcparam1.y,cfg->srcparam1.z,cfg->srcparam1.w}},
		     {{cfg->srcparam2.x,cfg->srcparam2.y,cfg->srcparam2.z,cfg->srcparam2.w}},
		     cfg->issaveref,cfg->maxgate,(uint)cfg->debuglevel, detreclen, cfg->outputtype, mesh->elemlen, 
		     cfg->mcmethod, cfg->method, 1.f/cfg->unitinmm, cfg->srcpos.w, 
		     mesh->nn, mesh->ne, mesh->nf, {{mesh->nmin.x,mesh->nmin.y,mesh->nmin.z}}, cfg->nout,
		     cfg->roulettesize, cfg->srcnum, {{cfg->crop0.x,cfg->crop0.y,cfg->crop0.z}}, 
		     mesh->srcelemlen, {{cfg->bary0.x,cfg->bary0.y,cfg->bary0.z,cfg->bary0.w}}, 
		     cfg->e0, cfg->isextdet, meshlen, cfg->nbuffer, ((1 << cfg->nbuffer)-1)};

     MCXReporter reporter={0.f};
     platform=mcx_list_gpu(cfg,&workdev,devices,&gpu);

     if(workdev>MAX_DEVICE)
         workdev=MAX_DEVICE;
     if(devices == NULL || workdev==0)
         mcx_error(-99,(char*)("Unable to find devices!"),__FILE__,__LINE__);

     cl_context_properties cps[3]={CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};

     /* Use NULL for backward compatibility */
     cl_context_properties* cprops=(platform==NULL)?NULL:cps;
     OCL_ASSERT(((mcxcontext=clCreateContextFromType(cprops,CL_DEVICE_TYPE_ALL,NULL,NULL,&status),status)));

     mcxqueue= (cl_command_queue*)malloc(workdev*sizeof(cl_command_queue));
     waittoread=(cl_event *)malloc(workdev*sizeof(cl_event));

     gseed=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gweight=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdref=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdetphoton=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     genergy=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gprogress=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdetected=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gsrcpattern=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     greporter=(cl_mem *)malloc(workdev*sizeof(cl_mem));

     /* The block is to move the declaration of prop closer to its use */
     cl_command_queue_properties prop = CL_QUEUE_PROFILING_ENABLE;

     totalcucore=0;
     for(i=0;i<workdev;i++){
         OCL_ASSERT(((mcxqueue[i]=clCreateCommandQueue(mcxcontext,devices[i],prop,&status),status)));
         totalcucore+=gpu[i].core;
	 if(!cfg->autopilot){
	    gpu[i].autothread=cfg->nthread;
	    gpu[i].autoblock=cfg->nblocksize;
	    gpu[i].maxgate=cfg->maxgate;
	 }else{
             // persistent thread mode
             if (gpu[i].vendor == dvIntelGPU){ // Intel HD graphics GPU
                 gpu[i].autoblock  = 64;
                 gpu[i].autothread = gpu[i].autoblock * 7 * gpu[i].sm; // 7 thread x SIMD-16 per Exec Unit (EU)
	     }else if (gpu[i].vendor == dvAMD){ // AMD GPU 
		 gpu[i].autoblock  = 64;
		 gpu[i].autothread = 2560 * gpu[i].sm; // 40 wavefronts * 64 threads/wavefront
             }else if(gpu[i].vendor == dvNVIDIA){
	       if (gpu[i].major == 2 || gpu[i].major == 3) { // fermi 2.x, kepler 3.x : max 7 blks per SM, 8 works better
                 gpu[i].autoblock  = 128;
                 gpu[i].autothread = gpu[i].autoblock * 8 * gpu[i].sm;
               }else if (gpu[i].major == 5) { // maxwell 5.x
                 gpu[i].autoblock  = 64;
                 gpu[i].autothread = gpu[i].autoblock * 16 * gpu[i].sm;
               }else if (gpu[i].major >= 6) { // pascal 6.x : max 32 blks per SM
                 gpu[i].autoblock  = 64;
                 gpu[i].autothread = gpu[i].autoblock * 64 * gpu[i].sm;
	       }
             }
         }
	 if(gpu[i].autothread%gpu[i].autoblock)
     	    gpu[i].autothread=(gpu[i].autothread/gpu[i].autoblock)*gpu[i].autoblock;
         if(gpu[i].maxgate==0 && meshlen>0){
             int needmem=meshlen+gpu[i].autothread*sizeof(float4)*4+sizeof(float)*cfg->maxdetphoton*hostdetreclen+10*1024*1024; /*keep 10M for other things*/
             gpu[i].maxgate=(gpu[i].globalmem-needmem)/meshlen;
             gpu[i].maxgate=MIN(((cfg->tend-cfg->tstart)/cfg->tstep+0.5),gpu[i].maxgate);     
	 }
     }
     cfg->maxgate=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);
     param.maxgate=cfg->maxgate;
     cl_uint nflen=mesh->nf*cfg->maxgate;

     fullload=0.f;
     for(i=0;i<workdev;i++)
     	fullload+=cfg->workload[i];

     if(fullload<EPS){
	for(i=0;i<workdev;i++)
     	    cfg->workload[i]=gpu[i].core;
	fullload=totalcucore;
     }

     field=(cl_float *)calloc(sizeof(cl_float)*meshlen,cfg->maxgate);
     dref=(cl_float *)calloc(sizeof(cl_float)*mesh->nf,cfg->maxgate);
     Pdet=(float*)calloc(cfg->maxdetphoton*sizeof(float),hostdetreclen);

     fieldlen=meshlen*cfg->maxgate;

     if(cfg->seed>0)
     	srand(cfg->seed);
     else
        srand(time(0));

     OCL_ASSERT(((gnode=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float3)*(mesh->nn),mesh->node,&status),status)));
     OCL_ASSERT(((gelem=clCreateBuffer(mcxcontext,RO_MEM, sizeof(int4)*(mesh->ne),mesh->elem,&status),status)));
     OCL_ASSERT(((gtype=clCreateBuffer(mcxcontext,RO_MEM, sizeof(int)*(mesh->ne),mesh->type,&status),status)));
     OCL_ASSERT(((gfacenb=clCreateBuffer(mcxcontext,RO_MEM, sizeof(int4)*(mesh->ne),mesh->facenb,&status),status)));
     if(mesh->srcelemlen>0)
         OCL_ASSERT(((gsrcelem=clCreateBuffer(mcxcontext,RO_MEM, sizeof(int)*(mesh->srcelemlen),mesh->srcelem,&status),status)));
     else
         gsrcelem=NULL;
     OCL_ASSERT(((gnormal=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float4)*(mesh->ne)*4,tracer->n,&status),status)));
     if(cfg->detpos && cfg->detnum)
           OCL_ASSERT(((gdetpos=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
     else
         gdetpos=NULL;

     OCL_ASSERT(((gproperty=clCreateBuffer(mcxcontext,RO_MEM, (mesh->prop+1+cfg->isextdet)*sizeof(medium),mesh->med,&status),status)));
     OCL_ASSERT(((gparam=clCreateBuffer(mcxcontext,RO_MEM, sizeof(MCXParam),&param,&status),status)));
     cl_mem (*clCreateBufferNV)(cl_context,cl_mem_flags, cl_mem_flags_NV, size_t, void*, cl_int*) = (cl_mem (*)(cl_context,cl_mem_flags, cl_mem_flags_NV, size_t, void*, cl_int*)) clGetExtensionFunctionAddressForPlatform(platform, "clCreateBufferNV");
     if (clCreateBufferNV == NULL)
         OCL_ASSERT(((gprogress[0]=clCreateBuffer(mcxcontext,RW_PTR, sizeof(cl_uint),NULL,&status),status)));
     else
         OCL_ASSERT(((gprogress[0]=clCreateBufferNV(mcxcontext,CL_MEM_READ_WRITE, NV_PIN, sizeof(cl_uint),NULL,&status),status)));
     progress = (cl_uint *)clEnqueueMapBuffer(mcxqueue[0], gprogress[0], CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(cl_uint), 0, NULL, NULL, NULL);
     *progress=0;

     for(i=0;i<workdev;i++){
       Pseed=(cl_uint*)malloc(sizeof(cl_uint)*gpu[i].autothread*RAND_SEED_WORD_LEN);
       energy=(cl_float*)calloc(sizeof(cl_float),gpu[i].autothread<<1);
       for (j=0; j<gpu[i].autothread*RAND_SEED_WORD_LEN;j++)
	   Pseed[j]=rand();
       OCL_ASSERT(((gseed[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(cl_uint)*gpu[i].autothread*RAND_SEED_WORD_LEN,Pseed,&status),status)));
       OCL_ASSERT(((gweight[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(float)*fieldlen,field,&status),status)));
       OCL_ASSERT(((gdref[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(float)*nflen,dref,&status),status)));
       OCL_ASSERT(((gdetphoton[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(float)*cfg->maxdetphoton*hostdetreclen,Pdet,&status),status)));
       OCL_ASSERT(((genergy[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(float)*(gpu[i].autothread<<1),energy,&status),status)));
       OCL_ASSERT(((gdetected[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(cl_uint),&detected,&status),status)));
       OCL_ASSERT(((greporter[i]=clCreateBuffer(mcxcontext,RW_MEM, sizeof(MCXReporter),&reporter,&status),status)));
       if(cfg->srctype==MCX_SRC_PATTERN)
           OCL_ASSERT(((gsrcpattern[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float)*(int)(cfg->srcparam1.w*cfg->srcparam2.w),cfg->srcpattern,&status),status)));
       else if(cfg->srctype==MCX_SRC_PATTERN3D)
           OCL_ASSERT(((gsrcpattern[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float)*(int)(cfg->srcparam1.x*cfg->srcparam1.y*cfg->srcparam1.z),cfg->srcpattern,&status),status)));
       else
           gsrcpattern[i]=NULL;
       free(Pseed);
       free(energy);
     }

     mcx_printheader(cfg);

     tic=StartTimer();
     if(cfg->issavedet){
         MMC_FPRINTF(cfg->flog,"- variant name: [%s] compiled with OpenCL version [%d]\n",
             "MMC-OpenCL",CL_VERSION_1_0);
     }else{
         MMC_FPRINTF(cfg->flog,"- code name: [MMC-OpenCL] compiled with OpenCL version [%d]\n",
             CL_VERSION_1_0);
     }
     MMC_FPRINTF(cfg->flog,"- compiled with: [RNG] %s [Seed Length] %d\n",MCX_RNG_NAME,RAND_SEED_WORD_LEN);
     MMC_FPRINTF(cfg->flog,"initializing streams ...\t");

     MMC_FPRINTF(cfg->flog,"init complete : %d ms\n",GetTimeMillis()-tic);fflush(cfg->flog);

     OCL_ASSERT(((mcxprogram=clCreateProgramWithSource(mcxcontext, 1,(const char **)&(cfg->clsource), NULL, &status),status)));

     if(cfg->optlevel>=1)
         sprintf(opt,"%s ","-cl-mad-enable -DMCX_USE_NATIVE");
     if(cfg->optlevel>=3)
         sprintf(opt+strlen(opt),"%s ","-DMCX_SIMPLIFY_BRANCH -DMCX_VECTOR_INDEX");
     
     if((uint)cfg->srctype<sizeof(sourceflag)/sizeof(sourceflag[0]))
         sprintf(opt+strlen(opt),"%s ",sourceflag[(uint)cfg->srctype]);

     sprintf(opt+strlen(opt),"%s",cfg->compileropt);
     if(cfg->isatomic)
         sprintf(opt+strlen(opt)," -DUSE_ATOMIC");
     if(cfg->issave2pt==0)
         sprintf(opt+strlen(opt)," -DMCX_SKIP_VOLUME");
     if(cfg->issavedet)
         sprintf(opt+strlen(opt)," -DMCX_SAVE_DETECTORS");
     if(cfg->issaveref)
         sprintf(opt+strlen(opt)," -DMCX_SAVE_DREF");
     if(cfg->isreflect)
         sprintf(opt+strlen(opt)," -DMCX_DO_REFLECTION");
     if(cfg->method==rtBLBadouelGrid)
         sprintf(opt+strlen(opt)," -DUSE_DMMC");

     MMC_FPRINTF(cfg->flog,"Building kernel with option: %s\n",opt);
     status=clBuildProgram(mcxprogram, 0, NULL, opt, NULL, NULL);

     size_t len;
     // get the details on the error, and store it in buffer
     clGetProgramBuildInfo(mcxprogram,devices[devid],CL_PROGRAM_BUILD_LOG,0,NULL,&len);
     if(len>0){
         char *msg;
	 int i;
         msg=(char *)calloc(len,1);
         clGetProgramBuildInfo(mcxprogram,devices[devid],CL_PROGRAM_BUILD_LOG,len,msg,NULL);
         for(i=0;i<(int)len;i++)
             if(msg[i]<='z' && msg[i]>='A'){
                 MMC_FPRINTF(cfg->flog,"Kernel build log:\n%s\n", msg);
                 break;
             }
	 free(msg);
     }
     if(status!=CL_SUCCESS)
	 mcx_error(-(int)status,(char*)("Error: Failed to build program executable!"),__FILE__,__LINE__);

     MMC_FPRINTF(cfg->flog,"build program complete : %d ms\n",GetTimeMillis()-tic);fflush(cfg->flog);

     mcxkernel=(cl_kernel*)malloc(workdev*sizeof(cl_kernel));

     for(i=0;i<workdev;i++){
         cl_int threadphoton, oddphotons;

         threadphoton=(int)(cfg->nphoton*cfg->workload[i]/(fullload*gpu[i].autothread*cfg->respin));
         oddphotons=(int)(cfg->nphoton*cfg->workload[i]/(fullload*cfg->respin)-threadphoton*gpu[i].autothread);

         MMC_FPRINTF(cfg->flog,"- [device %d(%d): %s] threadph=%d oddphotons=%d np=%.1f nthread=%d nblock=%d repetition=%d\n",i, gpu[i].id, gpu[i].name,threadphoton,oddphotons,
               cfg->nphoton*cfg->workload[i]/fullload,(int)gpu[i].autothread,(int)gpu[i].autoblock,cfg->respin);

	 OCL_ASSERT(((mcxkernel[i] = clCreateKernel(mcxprogram, "mmc_main_loop", &status),status)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 0, sizeof(cl_uint),(void*)&threadphoton)));
         OCL_ASSERT((clSetKernelArg(mcxkernel[i], 1, sizeof(cl_uint),(void*)&oddphotons)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 2, sizeof(cl_mem), (void*)&gparam)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 3, cfg->issavedet? sizeof(cl_float)*((int)gpu[i].autoblock)*detreclen : sizeof(int), NULL)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 4, sizeof(cl_mem), (void*)&gnode)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 5, sizeof(cl_mem), (void*)&gelem)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 6, sizeof(cl_mem), (void*)(gweight+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 7, sizeof(cl_mem), (void*)(gdref+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 8, sizeof(cl_mem), (void*)&gtype)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 9, sizeof(cl_mem), (void*)&gfacenb)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],10, sizeof(cl_mem), (void*)&gsrcelem)));	 
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],11, sizeof(cl_mem), (void*)&gnormal)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],12, sizeof(cl_mem), (void*)&gproperty)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],13, sizeof(cl_mem), (void*)&gdetpos)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],14, sizeof(cl_mem), (void*)(gdetphoton+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],15, sizeof(cl_mem), (void*)(gdetected+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],16, sizeof(cl_mem), (void*)(gseed+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],17, sizeof(cl_mem), (void*)(gprogress))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],18, sizeof(cl_mem), (void*)(genergy+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],19, sizeof(cl_mem), (void*)(greporter+i))));
	// OCL_ASSERT((clSetKernelArg(mcxkernel[i], 8, sizeof(cl_mem), (void*)(gsrcpattern+i))));
     }
     MMC_FPRINTF(cfg->flog,"set kernel arguments complete : %d ms %d\n",GetTimeMillis()-tic, param.method);fflush(cfg->flog);

     if(cfg->exportfield==NULL)
         cfg->exportfield=mesh->weight;
     if(cfg->exportdetected==NULL)
         cfg->exportdetected=(float*)malloc(hostdetreclen*cfg->maxdetphoton*sizeof(float));

     cfg->energytot=0.f;
     cfg->energyesc=0.f;
     cfg->runtime=0;

     //simulate for all time-gates in maxgate groups per run

     tic0=GetTimeMillis();

     for(t=cfg->tstart;t<cfg->tend;t+=cfg->tstep*cfg->maxgate){
       twindow0=t;
       twindow1=t+cfg->tstep*cfg->maxgate;

       MMC_FPRINTF(cfg->flog,"lauching mcx_main_loop for time window [%.1fns %.1fns] ...\n"
           ,twindow0*1e9,twindow1*1e9);fflush(cfg->flog);

       //total number of repetition for the simulations, results will be accumulated to field
       for(iter=0;iter<cfg->respin;iter++){
           MMC_FPRINTF(cfg->flog,"simulation run#%2d ... \n",iter+1); fflush(cfg->flog);fflush(cfg->flog);
	   param.tstart=twindow0;
	   param.tend=twindow1;

           for(devid=0;devid<workdev;devid++){
               OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],gparam,CL_TRUE,0,sizeof(MCXParam),&param, 0, NULL, NULL)));
               OCL_ASSERT((clSetKernelArg(mcxkernel[devid],2, sizeof(cl_mem), (void*)&gparam)));
               // launch mcxkernel
#ifndef USE_OS_TIMER
               OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid],mcxkernel[devid],1,NULL,&gpu[devid].autothread,&gpu[devid].autoblock, 0, NULL, &kernelevent)));
#else
               OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid],mcxkernel[devid],1,NULL,&gpu[devid].autothread,&gpu[devid].autoblock, 0, NULL, NULL)));
#endif
               OCL_ASSERT((clFlush(mcxqueue[devid])));
           }
           if((cfg->debuglevel & MCX_DEBUG_PROGRESS)){
	     int p0 = 0, ndone=-1;

	     mcx_progressbar(-0.f,cfg);

	     do{
               ndone = *progress;

	       if (ndone > p0){
		  mcx_progressbar((float)ndone/gpu[0].autothread*cfg->nphoton,cfg);
		  p0 = ndone;
	       }
               sleep_ms(100);
	     }while (p0 < gpu[0].autothread);
             mcx_progressbar(cfg->nphoton,cfg);
             MMC_FPRINTF(cfg->flog,"\n");

           }
           clEnqueueUnmapMemObject(mcxqueue[0], gprogress[0], progress, 0, NULL, NULL);

           //clWaitForEvents(workdev,waittoread);
           for(devid=0;devid<workdev;devid++)
               OCL_ASSERT((clFinish(mcxqueue[devid])));

           tic1=GetTimeMillis();
	   toc+=tic1-tic0;
           MMC_FPRINTF(cfg->flog,"kernel complete:  \t%d ms\nretrieving flux ... \t",tic1-tic);fflush(cfg->flog);

           if(cfg->runtime<tic1-tic)
               cfg->runtime=tic1-tic;

           for(devid=0;devid<workdev;devid++){
	     MCXReporter rep;
	     OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],greporter[devid],CL_TRUE,0,sizeof(MCXReporter),
                                            &rep, 0, NULL, waittoread+devid)));
	     reporter.raytet+=rep.raytet;
             if(cfg->issavedet){
                OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gdetected[devid],CL_FALSE,0,sizeof(uint),
                                            &detected, 0, NULL, NULL)));
                OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gdetphoton[devid],CL_TRUE,0,sizeof(float)*cfg->maxdetphoton*hostdetreclen,
	                                        Pdet, 0, NULL, NULL)));
		if(detected>cfg->maxdetphoton){
			MMC_FPRINTF(cfg->flog,"WARNING: the detected photon (%d) \
is more than what your have specified (%d), please use the -H option to specify a greater number\t"
                           ,detected,cfg->maxdetphoton);
		}else{
			MMC_FPRINTF(cfg->flog,"detected %d photons, total: %d\t",detected,cfg->detectedcount+detected);
		}
                cfg->his.detected+=detected;
                detected=MIN(detected,cfg->maxdetphoton);
		if(cfg->exportdetected){
                        cfg->exportdetected=(float*)realloc(cfg->exportdetected,(cfg->detectedcount+detected)*hostdetreclen*sizeof(float));
	                memcpy(cfg->exportdetected+cfg->detectedcount*(hostdetreclen),Pdet,detected*(hostdetreclen)*sizeof(float));
                        cfg->detectedcount+=detected;
		}
	     }
	     if(cfg->issaveref){
	        float *rawdref=(float*)calloc(sizeof(float),nflen);
	        OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gdref[devid],CL_TRUE,0,sizeof(float)*nflen,
	                                        rawdref, 0, NULL, NULL)));
		for(i=0;i<nflen;i++)  //accumulate field, can be done in the GPU
	            dref[i]+=rawdref[i]; //+rawfield[i+fieldlen];
	        free(rawdref);			
	     }
	     //handling the 2pt distributions
             if(cfg->issave2pt){
                float *rawfield=(float*)malloc(sizeof(float)*fieldlen);

        	OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gweight[devid],CL_TRUE,0,sizeof(cl_float)*fieldlen,
	                                         rawfield, 0, NULL, NULL)));
        	MMC_FPRINTF(cfg->flog,"transfer complete:        %d ms\n",GetTimeMillis()-tic);  fflush(cfg->flog);

                for(i=0;i<fieldlen;i++)  //accumulate field, can be done in the GPU
	            field[(i>>cfg->nbuffer)]+=rawfield[i]; //+rawfield[i+fieldlen];

	        free(rawfield);

/*        	if(cfg->respin>1){
                    for(i=0;i<fieldlen;i++)  //accumulate field, can be done in the GPU
                       field[fieldlen+i]+=field[i];
        	}
        	if(iter+1==cfg->respin){ 
                    if(cfg->respin>1)  //copy the accumulated fields back
                	memcpy(field,field+fieldlen,sizeof(cl_float)*fieldlen);
        	}
*/
        	if(cfg->isnormalized){
                    energy=(cl_float*)calloc(sizeof(cl_float),gpu[devid].autothread<<1);
                    OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],genergy[devid],CL_TRUE,0,sizeof(cl_float)*(gpu[devid].autothread<<1),
	                                     energy, 0, NULL, NULL)));
                    for(i=0;i<gpu[devid].autothread;i++){
                	cfg->energyesc+=energy[(i<<1)];
       	       		cfg->energytot+=energy[(i<<1)+1];
                	//eabsorp+=Plen[i].z;  // the accumulative absorpted energy near the source
                    }
		    free(energy);
        	}
             }
	     if(cfg->respin>1 && RAND_SEED_WORD_LEN>1){
               Pseed=(cl_uint*)malloc(sizeof(cl_uint)*gpu[devid].autothread*RAND_SEED_WORD_LEN);
               for (i=0; i<gpu[devid].autothread*RAND_SEED_WORD_LEN; i++)
		   Pseed[i]=rand();
               OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],gseed[devid],CL_TRUE,0,sizeof(cl_uint)*gpu[devid].autothread*RAND_SEED_WORD_LEN,
	                                        Pseed, 0, NULL, NULL)));
	       OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 15, sizeof(cl_mem), (void*)(gseed+devid))));
	       free(Pseed);
	     }
             OCL_ASSERT((clFinish(mcxqueue[devid])));
           }// loop over work devices
       }// iteration
     }// time gates

     fieldlen=(fieldlen>>cfg->nbuffer);
     field=realloc(field,sizeof(field[0])*fieldlen);
     if(cfg->exportfield){
         if(cfg->basisorder==0 || cfg->method==rtBLBadouelGrid){
             for(i=0;i<fieldlen;i++)
	         cfg->exportfield[i]+=field[i];
	 }else{
             for(i=0;i<cfg->maxgate;i++)
	       for(j=0;j<mesh->ne;j++){
		 float ww=field[i*mesh->ne+j]*0.25f;
		 int k;
	         for(k=0;k<mesh->elemlen;k++)
	             cfg->exportfield[i*mesh->nn+mesh->elem[j*mesh->elemlen+k]-1]+=ww;
	       }
	 }
     }
     
     if(cfg->issaveref && mesh->dref){
        for(i=0;i<nflen;i++)
	   mesh->dref[i]+=dref[i];
     }

     if(cfg->isnormalized){
         MMC_FPRINTF(cfg->flog,"normalizing raw data ...\t");fflush(cfg->flog);

         cfg->energyabs=cfg->energytot-cfg->energyesc;
	 mesh_normalize(mesh,cfg,cfg->energyabs,cfg->energytot,0);
     }
     if(cfg->issave2pt && cfg->parentid==mpStandalone){
         MMC_FPRINTF(cfg->flog,"saving data to file ...\t");
         mesh_saveweight(mesh,cfg,0);
         MMC_FPRINTF(cfg->flog,"saving data complete : %d ms\n\n",GetTimeMillis()-tic);
         fflush(cfg->flog);
     }
     if(cfg->issavedet && cfg->parentid==mpStandalone && cfg->exportdetected){
         cfg->his.unitinmm=cfg->unitinmm;
         cfg->his.savedphoton=cfg->detectedcount;
         cfg->his.detected=cfg->detectedcount;
         mesh_savedetphoton(cfg->exportdetected,NULL,cfg->detectedcount,(sizeof(RandType)*RAND_BUF_LEN),cfg);
     }
     if(cfg->issaveref){
	MMC_FPRINTF(cfg->flog,"saving surface diffuse reflectance ...");
	mesh_saveweight(mesh,cfg,1);
     }
     // total energy here equals total simulated photons+unfinished photons for all threads
     MMC_FPRINTF(cfg->flog,"simulated %ld photons (%ld) with %d devices (ray-tet %.0f)\nMCX simulation speed: %.2f photon/ms\n",
             cfg->nphoton,cfg->nphoton,workdev, reporter.raytet,(double)cfg->nphoton/toc);
     MMC_FPRINTF(cfg->flog,"total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to initial specular reflection is excluded in the total)\n",
             cfg->energytot,(cfg->energytot-cfg->energyesc)/cfg->energytot*100.f);
     fflush(cfg->flog);

     clReleaseMemObject(gnode);
     clReleaseMemObject(gelem);
     clReleaseMemObject(gtype);
     clReleaseMemObject(gfacenb);
     clReleaseMemObject(gsrcelem);
     clReleaseMemObject(gnormal);
     clReleaseMemObject(gproperty);
     clReleaseMemObject(gparam);
     if(cfg->detpos) clReleaseMemObject(gdetpos);

     for(i=0;i<workdev;i++){
         clReleaseMemObject(gseed[i]);
         clReleaseMemObject(gdetphoton[i]);
         clReleaseMemObject(gweight[i]);
	 clReleaseMemObject(gdref[i]);
         clReleaseMemObject(genergy[i]);
         clReleaseMemObject(gprogress[i]);
         clReleaseMemObject(gdetected[i]);
         if(gsrcpattern[i]) clReleaseMemObject(gsrcpattern[i]);
         clReleaseMemObject(greporter[i]);
         clReleaseKernel(mcxkernel[i]);
     }
     free(gseed);
     free(gdetphoton);
     free(gweight);
     free(gdref);
     free(genergy);
     free(gprogress);
     free(gdetected);
     free(gsrcpattern);
     free(mcxkernel);

     free(waittoread);

     if(gpu)
        free(gpu);

     for(devid=0;devid<workdev;devid++)
        clReleaseCommandQueue(mcxqueue[devid]);

     free(mcxqueue);
     clReleaseProgram(mcxprogram);
     clReleaseContext(mcxcontext);
#ifndef USE_OS_TIMER
     clReleaseEvent(kernelevent);
#endif
     free(field);
     if(Pdet)free(Pdet);
     free(dref);
}
