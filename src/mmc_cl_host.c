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
#include "mmc_cl_host.h"
#include "tictoc.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the 
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

extern cl_event kernelevent;

int mmc_init_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer,int argc, char**argv){
        mcx_initcfg(cfg);
        mcx_parsecmd(argc,argv,cfg);
	MMCDEBUG(cfg,dlTime,(cfg->flog,"initializing from commands ... "));
	mesh_init_from_cfg(mesh,cfg);
        return 0;
}

int mmc_reset_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	mmc_cleanup(cfg,mesh,tracer);
	mcx_initcfg(cfg);
	mesh_init(mesh);
        return 0;
}
int mmc_cleanup_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	tracer_clear(tracer);
	mesh_clear(mesh);
        mcx_clearcfg(cfg);
        return 0;
}
int mmc_prep(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	tracer_init(tracer,mesh,cfg->method);
	tracer_prep(tracer,cfg);
        return 0;
}

int mmc_run_cl (mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	double Eabsorb=0.0;
	RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
	RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
	unsigned int i, j;
	float raytri=0.f,raytri0=0.f;
	unsigned int threadid=0,ncomplete=0,t0,dt;
	visitor master={0.f,0.f,0.f,0,0,0,NULL,NULL,NULL,NULL,NULL,NULL};
     cl_uint detected=0,workdev;
     cl_float fullload=0.f;
     cl_uint   *Pseed;
     float  *Pdet;

     cl_context mcxcontext;                 // compute mcxcontext
     cl_command_queue *mcxqueue;          // compute command queue
     cl_program mcxprogram;                 // compute mcxprogram
     cl_kernel *mcxkernel;                   // compute mcxkernel
     cl_int status = 0;
     cl_device_id devices[MAX_DEVICE];
     cl_event * waittoread;
     cl_platform_id platform = NULL;

     cl_uint *cucount,totalcucore;
     cl_uint  devid=0;

     size_t mcgrid[1], mcblock[1];
     
     char opt[MAX_PATH_LENGTH]={'\0'};
     cl_uint detreclen=cfg->medianum+1;

     gmcconfig gcfg;
     gtetmesh  gmesh;
     gmmcdata  gdata; /*organizer for gpu global mem data*/

     platform=mcx_list_gpu(cfg,&workdev,devices);

     if(workdev>MAX_DEVICE)
         workdev=MAX_DEVICE;

     if(devices == NULL){
         OCL_ASSERT(-1);
     }

     cl_context_properties cps[3]={CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};

     /* Use NULL for backward compatibility */
     cl_context_properties* cprops=(platform==NULL)?NULL:cps;
     OCL_ASSERT(((mcxcontext=clCreateContextFromType(cprops,CL_DEVICE_TYPE_ALL,NULL,NULL,&status),status)));

     mcxqueue= (cl_command_queue*)malloc(workdev*sizeof(cl_command_queue));
     waittoread=(cl_event *)malloc(workdev*sizeof(cl_event));
     cucount=(cl_uint *)calloc(workdev,sizeof(cl_uint));

     gdata.detpos   =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.srcpattern=(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.node     =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.elem     =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.elem2    =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.srcelem  =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.detelem  =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.type     =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.facenb   =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.weight   =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.evol     =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.nvol     =(cl_mem *)malloc(workdev*sizeof(cl_mem));
     gdata.gvisitor =(cl_mem *)malloc(workdev*sizeof(cl_mem));

     /* The block is to move the declaration of prop closer to its use */
     cl_command_queue_properties prop = CL_QUEUE_PROFILING_ENABLE;

     totalcucore=0;
     for(i=0;i<workdev;i++){
         char pbuf[100]={'\0'};
         OCL_ASSERT(((mcxqueue[i]=clCreateCommandQueue(mcxcontext,devices[i],prop,&status),status)));
         OCL_ASSERT((clGetDeviceInfo(devices[i],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),(void*)(cucount+i),NULL)));
         OCL_ASSERT((clGetDeviceInfo(devices[i],CL_DEVICE_NAME,100,(void*)&pbuf,NULL)));
         if(strstr(pbuf,"ATI")){
            cucount[i]*=(80/5); // an ati core typically has 80 SP, and 80/5=16 VLIW
	 }else if(strstr(pbuf,"GeForce") || strstr(pbuf,"Quadro") || strstr(pbuf,"Tesla")){
            cucount[i]*=8;  // an nvidia MP typically has 8 SP
         }
         totalcucore+=cucount[i];
     }
     fullload=0.f;
     for(i=0;i<workdev;i++)
     	fullload+=cfg->workload[i];

     if(fullload<EPS){
	for(i=0;i<workdev;i++)
     	    cfg->workload[i]=cucount[i];
	fullload=totalcucore;
     }

     if(cfg->nthread%cfg->nblocksize)
        cfg->nthread=(cfg->nthread/cfg->nblocksize)*cfg->nblocksize;
     
     mcgrid[0]=cfg->nthread;
     mcblock[0]=cfg->nblocksize;

     Pseed=(cl_uint*)malloc(sizeof(cl_uint)*cfg->nthread*RAND_SEED_LEN);

     if(cfg->seed>0)
     	srand(cfg->seed);
     else
        srand(time(0));

     OCL_ASSERT(((gdata.gcfg=weightclCreateBuffer(mcxcontext,RO_MEM, sizeof(gmcconfig),&gcfg,&status),status)));
     OCL_ASSERT(((gdata.gmesh=clCreateBuffer(mcxcontext,RO_MEM, sizeof(gmcconfig),&gmesh.gcfg,&status),status)));

     for(i=0;i<workdev;i++){
       for (j=0; j<cfg->nthread*RAND_SEED_LEN;j++)
	   Pseed[j]=rand();
       OCL_ASSERT(((gdata.detpos[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(cl_uint)*cfg->nthread*RAND_SEED_LEN,Pseed,&status),status)));
       OCL_ASSERT(((gdata.srcpattern[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(cl_float)*(dimxyz)*cfg->maxgate,field,&status),status)));
       OCL_ASSERT(((gdata.node[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float)*cfg->maxdetphoton*(cfg->medianum+1),Pdet,&status),status)));
       OCL_ASSERT(((gdata.elem[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(float)*(cfg->nthread<<1),energy,&status),status)));
       OCL_ASSERT(((gdata.elem2[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(cl_uint),&stopsign,&status),status)));
       OCL_ASSERT(((gdata.facenb[i]=clCreateBuffer(mcxcontext,RO_MEM, sizeof(cl_uint),&detected,&status),status)));
       OCL_ASSERT(((gdata.weight[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.srcelem[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.detelem[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.evol[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.nvol[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.type[i]=clCreateBuffer(mcxcontext,RO_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
       OCL_ASSERT(((gdata.gvisitor[i]=clCreateBuffer(mcxcontext,RW_MEM, cfg->detnum*sizeof(float4),cfg->detpos,&status),status)));
     }

     fprintf(cfg->flog,"\
===============================================================================\n\
=                Mesh-based Monte Carlo eXtreme (MMCX) -- OpenCL              =\n\
=              Copyright (c) 2016 Qianqian Fang <q.fang at neu.edu>           =\n\
=                                                                             =\n\
=                    Computational Imaging Laboratory (CIL)                   =\n\
=             Department of Bioengineering, Northeastern University           =\n\
===============================================================================\n\
$MCXCL$Rev::    $ Last Commit $Date::                     $ by $Author::      $\n\
===============================================================================\n");

     tic=StartTimer();

     if(cfg->issavedet)
         fprintf(cfg->flog,"- variant name: [%s] compiled with OpenCL version [%d]\n",
             "Detective MCXCL",CL_VERSION_1_0);
     else
         fprintf(cfg->flog,"- code name: [Vanilla MCXCL] compiled with OpenCL version [%d]\n",
             CL_VERSION_1_0);

     fprintf(cfg->flog,"- compiled with: [RNG] %s [Seed Length] %d\n",MCX_RNG_NAME,RAND_SEED_LEN);
     fprintf(cfg->flog,"initializing streams ...\t");
     fflush(cfg->flog);

     fprintf(cfg->flog,"init complete : %d ms\n",GetTimeMillis()-tic);

     OCL_ASSERT(((mcxprogram=clCreateProgramWithSource(mcxcontext, 1,(const char **)&(cfg->clsource), NULL, &status),status)));

     sprintf(opt,"-cl-mad-enable -cl-fast-relaxed-math %s",cfg->compileropt);
     if(cfg->issavedet)
         sprintf(opt+strlen(opt)," -D MCX_SAVE_DETECTORS");
     if(cfg->isreflect)
         sprintf(opt+strlen(opt)," -D MCX_DO_REFLECTION");
     sprintf(opt+strlen(opt)," %s",cfg->compileropt);

     status=clBuildProgram(mcxprogram, 0, NULL, opt, NULL, NULL);
     
     if(status!=CL_SUCCESS){
	 size_t len;
	 char *msg;
	 // get the details on the error, and store it in buffer
	 clGetProgramBuildInfo(mcxprogram,devices[devid],CL_PROGRAM_BUILD_LOG,0,NULL,&len); 
	 msg=new char[len];
	 clGetProgramBuildInfo(mcxprogram,devices[devid],CL_PROGRAM_BUILD_LOG,len,msg,NULL); 
	 fprintf(cfg->flog,"Kernel build error:\n%s\n", msg);
	 mcx_error(-(int)status,(char*)("Error: Failed to build program executable!"),__FILE__,__LINE__);
	 delete msg;
     }
     fprintf(cfg->flog,"build program complete : %d ms\n",GetTimeMillis()-tic);

     mcxkernel=(cl_kernel*)malloc(workdev*sizeof(cl_kernel));

     for(i=0;i<workdev;i++){
         cl_int threadphoton, oddphotons;

         threadphoton=(int)(cfg->nphoton*cfg->workload[i]/(fullload*cfg->nthread*cfg->respin));
         oddphotons=(int)(cfg->nphoton*cfg->workload[i]/(fullload*cfg->respin)-threadphoton*cfg->nthread);
         fprintf(cfg->flog,"- [device %d] threadph=%d oddphotons=%d np=%.1f nthread=%d repetition=%d\n",i,threadphoton,oddphotons,
               cfg->nphoton*cfg->workload[i]/fullload,cfg->nthread,cfg->respin);

	 OCL_ASSERT(((mcxkernel[i] = clCreateKernel(mcxprogram, "mcx_main_loop", &status),status)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 0, sizeof(cl_uint),(void*)&threadphoton)));
         OCL_ASSERT((clSetKernelArg(mcxkernel[i], 1, sizeof(cl_uint),(void*)&oddphotons)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 2, sizeof(cl_mem), (void*)&gmedia)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 3, sizeof(cl_mem), (void*)(gfield+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 4, sizeof(cl_mem), (void*)(genergy+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 5, sizeof(cl_mem), (void*)(gseed+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 6, sizeof(cl_mem), (void*)(gdetphoton+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 7, sizeof(cl_mem), (void*)&gproperty)));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 8, sizeof(cl_mem), (void*)(gdetpos+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i], 9, sizeof(cl_mem), (void*)(gstopsign+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],10, sizeof(cl_mem), (void*)(gdetected+i))));
	 OCL_ASSERT((clSetKernelArg(mcxkernel[i],11, cfg->issavedet? sizeof(cl_float)*cfg->nblocksize*param.maxmedia : 1, NULL)));
     }
     fprintf(cfg->flog,"set kernel arguments complete : %d ms\n",GetTimeMillis()-tic);

     if(cfg->exportfield==NULL)
         cfg->exportfield=(float *)calloc(sizeof(float)*cfg->dim.x*cfg->dim.y*cfg->dim.z,cfg->maxgate*2);
     if(cfg->exportdetected==NULL)
         cfg->exportdetected=(float*)malloc((cfg->medianum+1)*cfg->maxdetphoton*sizeof(float));

     cfg->energytot=0.f;
     cfg->energyesc=0.f;
     cfg->runtime=0;

     //simulate for all time-gates in maxgate groups per run

     cl_float Vvox;
     Vvox=cfg->steps.x*cfg->steps.y*cfg->steps.z;
     tic0=GetTimeMillis();

     for(t=cfg->tstart;t<cfg->tend;t+=cfg->tstep*cfg->maxgate){
       twindow0=t;
       twindow1=t+cfg->tstep*cfg->maxgate;

       fprintf(cfg->flog,"lauching mcx_main_loop for time window [%.1fns %.1fns] ...\n"
           ,twindow0*1e9,twindow1*1e9);

       //total number of repetition for the simulations, results will be accumulated to field
       for(iter=0;iter<cfg->respin;iter++){
           fprintf(cfg->flog,"simulation run#%2d ... \t",iter+1); fflush(cfg->flog);
	   param.twin0=twindow0;
	   param.twin1=twindow1;
           for(devid=0;devid<workdev;devid++){
               OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],gparam,CL_TRUE,0,sizeof(MCXParam),&param, 0, NULL, NULL)));
               OCL_ASSERT((clSetKernelArg(mcxkernel[devid],12, sizeof(cl_mem), (void*)&gparam)));
               // launch mcxkernel
#ifndef USE_OS_TIMER
               OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid],mcxkernel[devid],1,NULL,mcgrid,mcblock, 0, NULL, &kernelevent)));
#else
               OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid],mcxkernel[devid],1,NULL,mcgrid,mcblock, 0, NULL, NULL)));
#endif
               OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gdetected[devid],CL_FALSE,0,sizeof(uint),
                                            &detected, 0, NULL, waittoread+devid)));
           }
           clWaitForEvents(workdev,waittoread);
           tic1=GetTimeMillis();
	   toc+=tic1-tic0;
           fprintf(cfg->flog,"kernel complete:  \t%d ms\nretrieving flux ... \t",tic1-tic);

           for(devid=0;devid<workdev;devid++){
             if(cfg->issavedet){
                OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gdetphoton[devid],CL_TRUE,0,sizeof(float)*cfg->maxdetphoton*(cfg->medianum+1),
	                                        Pdet, 0, NULL, NULL)));
		if(detected>cfg->maxdetphoton){
			fprintf(cfg->flog,"WARNING: the detected photon (%d) \
is more than what your have specified (%d), please use the -H option to specify a greater number\t"
                           ,detected,cfg->maxdetphoton);
		}else{
			fprintf(cfg->flog,"detected %d photons, total: %d\t",detected,cfg->detectedcount+detected);
		}
                cfg->his.detected+=detected;
                detected=MIN(detected,cfg->maxdetphoton);
		if(cfg->exportdetected){
                        cfg->exportdetected=(float*)realloc(cfg->exportdetected,(cfg->detectedcount+detected)*detreclen*sizeof(float));
	                memcpy(cfg->exportdetected+cfg->detectedcount*(detreclen),Pdet,detected*(detreclen)*sizeof(float));
                        cfg->detectedcount+=detected;
		}
	     }
	     //handling the 2pt distributions
             if(cfg->issave2pt){
               OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],gfield[devid],CL_TRUE,0,sizeof(cl_float)*dimxyz*cfg->maxgate,
	                                        field, 0, NULL, NULL)));
               fprintf(cfg->flog,"transfer complete:\t%d ms\n",GetTimeMillis()-tic);  fflush(cfg->flog);

               if(cfg->respin>1){
                   for(i=0;i<fieldlen;i++)  //accumulate field, can be done in the GPU
                      field[fieldlen+i]+=field[i];
               }
               if(iter+1==cfg->respin){ 
                   if(cfg->respin>1)  //copy the accumulated fields back
                       memcpy(field,field+fieldlen,sizeof(cl_float)*fieldlen);
               }
                   if(cfg->isnormalized){

                       OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid],genergy[devid],CL_TRUE,0,sizeof(cl_float)*(cfg->nthread<<1),
	                                        energy, 0, NULL, NULL)));
                       for(i=0;i<cfg->nthread;i++){
                           cfg->energyesc+=energy[(i<<1)];
       	       	       	   cfg->energytot+=energy[(i<<1)+1];
                           //eabsorp+=Plen[i].z;  // the accumulative absorpted energy near the source
                       }
                   }
	           if(cfg->exportfield){
	               for(i=0;i<fieldlen;i++)
		           cfg->exportfield[i]+=field[i];
                   }
             }
	     //initialize the next simulation
	     if(twindow1<cfg->tend && iter<cfg->respin){
                  memset(field,0,sizeof(cl_float)*dimxyz*cfg->maxgate);
                  OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],gfield[devid],CL_TRUE,0,sizeof(cl_float)*dimxyz*cfg->maxgate,
                                                field, 0, NULL, NULL)));
		  OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 3, sizeof(cl_mem), (void*)(gfield+devid))));
	     }
	     if(cfg->respin>1 && RAND_SEED_LEN>1){
               for (i=0; i<cfg->nthread*RAND_SEED_LEN; i++)
		   Pseed[i]=rand();
               OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],gseed[devid],CL_TRUE,0,sizeof(cl_uint)*cfg->nthread*RAND_SEED_LEN,
	                                        Pseed, 0, NULL, NULL)));
	       OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 5, sizeof(cl_mem), (void*)(gseed+devid))));
	     }
             OCL_ASSERT((clFinish(mcxqueue[devid])));
           }// loop over work devices
       }// iteration
       if(twindow1<cfg->tend){
	    cl_float *tmpenergy=(cl_float*)calloc(sizeof(cl_float),cfg->nthread*3);
            OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid],genergy[devid],CL_TRUE,0,sizeof(cl_float)*(cfg->nthread<<1),
                                        tmpenergy, 0, NULL, NULL)));
	    OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 4, sizeof(cl_mem), (void*)(genergy+devid))));	
	    free(tmpenergy);
       }
     }// time gates

     if(cfg->isnormalized){
	   float scale=0.f;
           fprintf(cfg->flog,"normalizing raw data ...\t");

           if(cfg->outputtype==otFlux || cfg->outputtype==otFluence){
               scale=1.f/(cfg->energytot*Vvox*cfg->tstep);
	       if(cfg->unitinmm!=1.f)
		   scale*=cfg->unitinmm; /* Vvox (in mm^3 already) * (Tstep) * (Eabsorp/U) */

               if(cfg->outputtype==otFluence)
		   scale*=cfg->tstep;
	   }else if(cfg->outputtype==otEnergy || cfg->outputtype==otJacobian)
	       scale=1.f/cfg->energytot;

	 fprintf(cfg->flog,"normalization factor alpha=%f\n",scale);  fflush(cfg->flog);
         mcx_normalize(cfg->exportfield,scale,fieldlen);
     }
     if(cfg->issave2pt && cfg->parentid==mpStandalone){
         fprintf(cfg->flog,"saving data to file ... %d %d\t",fieldlen,cfg->maxgate);
         mcx_savedata(cfg->exportfield,fieldlen,0,"mc2",cfg);
         fprintf(cfg->flog,"saving data complete : %d ms\n\n",GetTimeMillis()-tic);
         fflush(cfg->flog);
     }
     if(cfg->issavedet && cfg->parentid==mpStandalone && cfg->exportdetected){
         cfg->his.unitinmm=cfg->unitinmm;
         cfg->his.savedphoton=cfg->detectedcount;
         cfg->his.detected=cfg->detectedcount;
         mcx_savedetphoton(cfg->exportdetected,cfg->seeddata,cfg->detectedcount,0,cfg);
     }

     // total energy here equals total simulated photons+unfinished photons for all threads
     fprintf(cfg->flog,"simulated %d photons (%d) with %d CUs with %d threads (repeat x%d)\nMCX simulation speed: %.2f photon/ms\n",
             cfg->nphoton,cfg->nphoton,workdev,cfg->nthread, cfg->respin,(double)cfg->nphoton/toc); fflush(cfg->flog);
     fprintf(cfg->flog,"total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to initial specular reflection is excluded in the total)\n",
             cfg->energytot,(cfg->energytot-cfg->energyesc)/cfg->energytot*100.f);fflush(cfg->flog);
     fflush(cfg->flog);

     clReleaseMemObject(gmedia);
     clReleaseMemObject(gproperty);
     clReleaseMemObject(gparam);

     for(i=0;i<workdev;i++){
         clReleaseMemObject(gdata.[i]);
         clReleaseMemObject(gseed[i]);
         clReleaseMemObject(genergy[i]);
         clReleaseMemObject(gstopsign[i]);
         clReleaseMemObject(gdetected[i]);
         clReleaseMemObject(gdetpos[i]);
         clReleaseKernel(mcxkernel[i]);
     }
     free(gfield);
     free(gseed);
     free(genergy);
     free(gstopsign);
     free(gdetected);
     free(gdetpos);
     free(mcxkernel);

     free(waittoread);

     for(devid=0;devid<workdev;devid++)
        clReleaseCommandQueue(mcxqueue[devid]);

     free(mcxqueue);
     clReleaseProgram(mcxprogram);
     clReleaseContext(mcxcontext);
#ifndef USE_OS_TIMER
     clReleaseEvent(kernelevent);
#endif
     free(Pseed);
     free(energy);
     free(field);

        dt=GetTimeMillis();
        MMCDEBUG(cfg,dlTime,(cfg->flog,"seed=%u\nsimulating ... ",cfg->seed));

	/***************************************************************************//**
	The master thread then spawn multiple work-threads depending on your
	OpenMP settings. By default, the total thread number (master + work) is 
	your total CPU core number. For example, if you have a dual-core CPU, 
	the total thread number is 2; if you have two quad-core CPUs, the total 
	thread number is 8. If you want to set the total thread number manually, 
	you need to set the OMP_NUM_THREADS environment variable. For example, 
	\c OMP_NUM_THREADS=3 sets the total thread number to 3.
	*******************************************************************************/

/** \subsection ssimu Parallel photon transport simulation */

#pragma omp parallel private(ran0,ran1,threadid)
{
	visitor visit={0.f,0.f,1.f/cfg->tstep,DET_PHOTON_BUF,0,0,NULL,NULL,0.f,0.f};
	visit.reclen=(2+((cfg->ismomentum)>0))*mesh->prop+(cfg->issaveexit>0)*6+2;
	if(cfg->issavedet){
	    if(cfg->issaveseed)
	        visit.photonseed=calloc(visit.detcount,(sizeof(RandType)*RAND_BUF_LEN));
	    visit.partialpath=(float*)calloc(visit.detcount*visit.reclen,sizeof(float));
	}
#ifdef _OPENMP
	threadid=omp_get_thread_num();
#endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg->seed),threadid);

	/*launch photons*/
        #pragma omp for reduction(+:Eabsorb) reduction(+:raytri,raytri0)
	for(i=0;i<cfg->nphoton;i++){
		visit.raytet=0.f;
		visit.raytet0=0.f;
		if(cfg->seed==SEED_FROM_FILE)
		    Eabsorb+=onephoton(i,tracer,mesh,cfg,((RandType *)cfg->photonseed)+i*RAND_BUF_LEN,ran1,&visit);
		else
		    Eabsorb+=onephoton(i,tracer,mesh,cfg,ran0,ran1,&visit);
		raytri+=visit.raytet;
		raytri0+=visit.raytet0;

		#pragma omp atomic
		   ncomplete++;

		if((cfg->debuglevel & dlProgress) && threadid==0)
			mcx_progressbar(ncomplete,cfg);
		if(cfg->issave2pt && cfg->checkpt[0])
			mesh_saveweightat(mesh,cfg,i+1);
	}

	#pragma omp atomic
		master.totalweight += visit.totalweight;

	if(cfg->issavedet){
	    #pragma omp atomic
		master.detcount+=visit.bufpos;
            #pragma omp barrier
	    if(threadid==0){
		master.partialpath=(float*)calloc(master.detcount*visit.reclen,sizeof(float));
	        if(cfg->issaveseed)
        	    master.photonseed=calloc(master.detcount,(sizeof(RandType)*RAND_BUF_LEN));
            }
            #pragma omp barrier
            #pragma omp critical
            {
		memcpy(master.partialpath+master.bufpos*visit.reclen,
		       visit.partialpath,visit.bufpos*visit.reclen*sizeof(float));
                if(cfg->issaveseed)
                    memcpy((unsigned char*)master.photonseed+master.bufpos*(sizeof(RandType)*RAND_BUF_LEN),
                            visit.photonseed,visit.bufpos*(sizeof(RandType)*RAND_BUF_LEN));
		master.bufpos+=visit.bufpos;
            }
            #pragma omp barrier
	    free(visit.partialpath);
            if(cfg->issaveseed && visit.photonseed)
                 free(visit.photonseed);
	}
}

        /** \subsection sreport Post simulation */

	if((cfg->debuglevel & dlProgress))
		mcx_progressbar(cfg->nphoton,cfg);

	dt=GetTimeMillis()-dt;
	MMCDEBUG(cfg,dlProgress,(cfg->flog,"\n"));
        MMCDEBUG(cfg,dlTime,(cfg->flog,"\tdone\t%d\n",dt));
        MMCDEBUG(cfg,dlTime,(cfg->flog,"speed ...\t%.2f photon/ms, %.0f ray-tetrahedron tests (%.0f were overhead)\n",(double)cfg->nphoton/dt,raytri,raytri0));
        if(cfg->issavedet)
           fprintf(cfg->flog,"detected %d photons\n",master.detcount);

	if(cfg->isnormalized){
          cfg->his.normalizer=mesh_normalize(mesh,cfg,Eabsorb,master.totalweight);
          fprintf(cfg->flog,"total simulated energy: %f\tabsorbed: %5.5f%%\tnormalizor=%g\n",
		master.totalweight,100.f*Eabsorb/master.totalweight,cfg->his.normalizer);
	}
	if(cfg->issave2pt){
		switch(cfg->outputtype){
		    case otFlux:   MMCDEBUG(cfg,dlTime,(cfg->flog,"saving flux ...")); break;
                    case otFluence:MMCDEBUG(cfg,dlTime,(cfg->flog,"saving fluence ...")); break;
                    case otEnergy: MMCDEBUG(cfg,dlTime,(cfg->flog,"saving energy deposit ...")); break;
		}
		mesh_saveweight(mesh,cfg);
	}
	if(cfg->issavedet){
		MMCDEBUG(cfg,dlTime,(cfg->flog,"saving detected photons ..."));
		mesh_savedetphoton(master.partialpath,master.photonseed,master.bufpos,(sizeof(RandType)*RAND_BUF_LEN),cfg);
		free(master.partialpath);
                if(cfg->issaveseed && master.photonseed)
                    free(master.photonseed);
	}
        MMCDEBUG(cfg,dlTime,(cfg->flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	return 0;
}

