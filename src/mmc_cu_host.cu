/**
 **  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
 **
 **  \author Qianqian Fang <q.fang at neu.edu>
 **
 **  \section sref Reference:
 **  \li \c (\b Fang2010) Qianqian Fang, <a
 *href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
 **          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing
 **          in Pluker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175
 *(2010).
 **  \li \c (\b Fang2009) Qianqian Fang and David A. Boas,
 **          <a
 *href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
 **          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media
 *Accelerated
 **          by Graphics Processing Units,"</a> Optics Express, 17(22)
 *20178-20190 (2009).
 **
 **  \section slicense License
 **          GPL v3, see LICENSE.txt for details
 *******************************************************************************/

/**
  \file    mmc_host.c
  \brief   << Driver program of MMC >>
*/

#define inlinefun __device__

#include "mcx_const.h"
#include "mmc_cu_host.h"
#include "tictoc.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "mcx_const.h"

#include "mmc_core.cu"

/************************************************************************** In
this unit, we first launch a master thread and initialize the necessary data
structures.This include the command line options(cfg), tetrahedral mesh(mesh)
and the ray tracer precomputed data (tracer).
******************************************************************************/
#define CUDA_ASSERT(a)                                                         \
  mcx_cu_assess((a), __FILE__, __LINE__) ///< macro to report CUDA error
int mcx_corecount(int v1, int v2) {
  int v = v1 * 10 + v2;
  if (v < 20)
    return 8;
  else if (v < 21)
    return 32;
  else if (v < 30)
    return 48;
  else if (v < 50)
    return 192;
  else
    return 128;
}

int mcx_smxblock(int v1, int v2) {
  int v = v1 * 10 + v2;
  if (v < 30)
    return 8;
  else if (v < 50)
    return 16;
  else
    return 32;
}

/**
  assert cuda memory allocation result
 */
void mcx_cu_assess(cudaError_t cuerr, const char *file, const int linenum) {
  if (cuerr != cudaSuccess) {
    mcx_error(-(int)cuerr, (char *)cudaGetErrorString(cuerr), file, linenum);
  }
}

/*
   master driver code to run MC simulations
*/
int mcx_list_cu_gpu(mcconfig *cfg, GPUInfo **info) {
#if __DEVICE_EMULATION__
  return 1;
#else
  int dev;
  int deviceCount, activedev = 0;

  CUDA_ASSERT(cudaGetDeviceCount(&deviceCount));
  if (deviceCount == 0) {
    MMC_FPRINTF(stderr,
                S_RED "ERROR: No CUDA-capable GPU device found\n" S_RESET);
    return 0;
  }
  *info = (GPUInfo *)calloc(deviceCount, sizeof(GPUInfo));
  if (cfg->gpuid && cfg->gpuid > (uint)deviceCount) {
    MMC_FPRINTF(stderr,
                S_RED "ERROR: Specified GPU ID is out of range\n" S_RESET);
    return 0;
  }

  // scan from the first device
  for (dev = 0; dev < deviceCount; dev++) {
    cudaDeviceProp dp;
    CUDA_ASSERT(cudaGetDeviceProperties(&dp, dev));

    if (cfg->isgpuinfo == 3)
      activedev++;
    else if (cfg->deviceid[dev] == '1') {
      cfg->deviceid[dev] = '\0';
      cfg->deviceid[activedev] = dev + 1;
      activedev++;
    }

    strncpy((*info)[dev].name, dp.name, MAX_SESSION_LENGTH);
    (*info)[dev].id = dev + 1;
    (*info)[dev].devcount = deviceCount;
    (*info)[dev].major = dp.major;
    (*info)[dev].minor = dp.minor;
    (*info)[dev].globalmem = dp.totalGlobalMem;
    (*info)[dev].constmem = dp.totalConstMem;
    (*info)[dev].sharedmem = dp.sharedMemPerBlock;
    (*info)[dev].regcount = dp.regsPerBlock;
    (*info)[dev].clock = dp.clockRate;
    (*info)[dev].sm = dp.multiProcessorCount;
    (*info)[dev].core =
        dp.multiProcessorCount * mcx_corecount(dp.major, dp.minor);
    (*info)[dev].maxmpthread = dp.maxThreadsPerMultiProcessor;
    (*info)[dev].maxgate = cfg->maxgate;
    (*info)[dev].autoblock =
        (*info)[dev].maxmpthread / mcx_smxblock(dp.major, dp.minor);
    (*info)[dev].autothread = (*info)[dev].autoblock *
                              mcx_smxblock(dp.major, dp.minor) *
                              (*info)[dev].sm;

    if (strncmp(dp.name, "Device Emulation", 16)) {
      if (cfg->isgpuinfo) {
        MMC_FPRINTF(stdout,
                    S_BLUE "=============================   GPU Infomation  ================================\n" S_RESET);
        MMC_FPRINTF(stdout, "Device %d of %d:\t\t%s\n", (*info)[dev].id,
                    (*info)[dev].devcount, (*info)[dev].name);
        MMC_FPRINTF(stdout, "Compute Capability:\t%u.%u\n", (*info)[dev].major,
                    (*info)[dev].minor);
        MMC_FPRINTF(stdout,
                    "Global Memory:\t\t%u B\nConstant Memory:\t%u B\n"
                    "Shared Memory:\t\t%u B\nRegisters:\t\t%u\nClock "
                    "Speed:\t\t%.2f GHz\n",
                    (unsigned int)(*info)[dev].globalmem,
                    (unsigned int)(*info)[dev].constmem,
                    (unsigned int)(*info)[dev].sharedmem,
                    (unsigned int)(*info)[dev].regcount,
                    (*info)[dev].clock * 1e-6f);
#if CUDART_VERSION >= 2000
        MMC_FPRINTF(stdout, "Number of MPs:\t\t%u\nNumber of Cores:\t%u\n",
                    (*info)[dev].sm, (*info)[dev].core);
#endif
        MMC_FPRINTF(stdout, "SMX count:\t\t%u\n", (*info)[dev].sm);
      }
    }
  }

  if (cfg->isgpuinfo == 2 &&
      cfg->parentid == mpStandalone) { // list GPU info only
    exit(0);
  }

  if (activedev < MAX_DEVICE) {
    cfg->deviceid[activedev] = '\0';
  }

  return activedev;
#endif
}

void mmc_run_simulation(mcconfig *cfg, tetmesh *mesh, raytracer *tracer,GPUInfo *gpu,void (*progressfun)(float, void *),void *handle){
  uint i, j;
  float t, twindow0, twindow1;
  float fullload = 0.f;
  float *energy;

  uint detected = 0;
  int gpuid, threadid = 0;
  uint tic, tic0, tic1, toc = 0, fieldlen;
  int threadphoton, oddphotons;
  dim3 mcgrid, mcblock;

  float3 *gnode;
  int4 *gelem, *gfacenb;
  float4 *gnormal;
  int *gtype, *gsrcelem;
  uint *gseed, *gdetected;
  volatile int *progress, *gprogress;
  float *gweight, *gdref, *gdetphoton, *genergy, *gsrcpattern;
  RandType *gphotonseed=NULL, *greplayseed=NULL;
  float  *greplayweight=NULL, *greplaytime=NULL;

  MCXReporter *greporter;
  uint meshlen = ((cfg->method == rtBLBadouelGrid) ? cfg->crop0.z : mesh->ne)
                 << cfg->nbuffer; // use 4 copies to reduce racing

  float *field, *dref = NULL;

  uint *Pseed=NULL;
  float *Pdet=NULL;
  RandType *Pphotonseed=NULL;

  uint detreclen = (2 + ((cfg->ismomentum) > 0)) * mesh->prop +
                   (cfg->issaveexit > 0) * 6 + 1;
  uint hostdetreclen = detreclen + 1;

  MCXParam param = {make_float3(cfg->srcpos.x, cfg->srcpos.y, cfg->srcpos.z),
                    make_float3(cfg->srcdir.x, cfg->srcdir.y, cfg->srcdir.z),
                    cfg->tstart,
                    cfg->tend,
                    (uint)cfg->isreflect,
                    (uint)cfg->issavedet,
                    (uint)cfg->issaveexit,
                    (uint)cfg->ismomentum,
                    (uint)cfg->isatomic,
                    (uint)cfg->isspecular,
                    1.f / cfg->tstep,
                    cfg->minenergy,
                    cfg->maxdetphoton,
                    (uint)mesh->prop,
                    (uint)cfg->detnum,
                    (int)cfg->voidtime,
                    (int)cfg->srctype,
                    cfg->srcparam1,
                    cfg->srcparam2,
                    (uint)cfg->issaveref,
                    (uint)cfg->maxgate,
                    (uint)cfg->debuglevel,
                    (int)detreclen,
                    cfg->outputtype,
                    mesh->elemlen,
                    cfg->mcmethod,
                    cfg->method,
                    1.f / cfg->unitinmm,
                    cfg->srcdir.w,
                    mesh->nn,
                    mesh->ne,
                    mesh->nf,
                    make_float3(mesh->nmin.x, mesh->nmin.y, mesh->nmin.z),
                    cfg->nout,
                    cfg->roulettesize,
                    cfg->srcnum,
                    cfg->crop0,
                    mesh->srcelemlen,
                    cfg->bary0,
                    cfg->e0,
                    cfg->isextdet,
                    (int)meshlen,
                    cfg->nbuffer,
                    (uint)(mesh->prop + 1 + cfg->isextdet) + cfg->detnum,
		    (uint)(MIN((MAX_PROP-(mesh->prop + 1 + cfg->isextdet) - cfg->detnum), ((mesh->ne)<<2)) >>2),  /*max count of elem normal data in const mem*/
		    cfg->issaveseed,
		    cfg->seed
		    };

  MCXReporter reporter = {0.f};

  if(progressfun==NULL)
       cfg->debuglevel=cfg->debuglevel & (~MCX_DEBUG_PROGRESS);

#ifdef _OPENMP
  threadid = omp_get_thread_num();
#endif
  if (threadid < MAX_DEVICE && cfg->deviceid[threadid] == '\0')
    return;

  gpuid = cfg->deviceid[threadid] - 1;
  if (gpuid < 0)
    mcx_error(-1, "GPU ID must be non-zero", __FILE__, __LINE__);
  CUDA_ASSERT(cudaSetDevice(gpuid));

  #pragma omp master
  {
    if (cfg->exportfield == NULL)
      cfg->exportfield = mesh->weight;
    if (cfg->exportdetected == NULL)
      cfg->exportdetected = (float *)malloc(hostdetreclen * cfg->maxdetphoton * sizeof(float));
    if(cfg->issaveseed && cfg->exportseed==NULL)
      cfg->exportseed=(unsigned char*)malloc(cfg->maxdetphoton*(sizeof(RandType)*RAND_BUF_LEN));

    cfg->energytot = 0.f;
    cfg->energyesc = 0.f;
    cfg->runtime = 0;
  }
  #pragma omp barrier

  if (gpu[gpuid].maxgate == 0 && meshlen > 0) {
    int needmem = meshlen + gpu[gpuid].autothread * sizeof(float4) * 4 +
                  sizeof(float) * cfg->maxdetphoton * hostdetreclen +
                  10 * 1024 * 1024; /*keep 10M for other things*/
    gpu[gpuid].maxgate = (gpu[gpuid].globalmem - needmem) / meshlen;
    gpu[gpuid].maxgate =
        MIN(((cfg->tend - cfg->tstart) / cfg->tstep + 0.5), gpu[gpuid].maxgate);
  }
  
  if (!cfg->autopilot) {
    uint gates = (uint)((cfg->tend - cfg->tstart) / cfg->tstep + 0.5);
    gpu[gpuid].autothread = cfg->nthread;
    gpu[gpuid].autoblock = cfg->nblocksize;
    if (cfg->maxgate == 0)
      cfg->maxgate = gates;
    else if ((uint)cfg->maxgate > gates)
      cfg->maxgate = gates;
    gpu[gpuid].maxgate = cfg->maxgate;
  }
  if (gpu[gpuid].autothread % gpu[gpuid].autoblock)
    gpu[gpuid].autothread =
        (gpu[gpuid].autothread / gpu[gpuid].autoblock) * gpu[gpuid].autoblock;

  param.maxgate = gpu[gpuid].maxgate;

  uint nflen = mesh->nf * cfg->maxgate;
  #pragma omp master
  fullload = 0.f;
  for (i = 0; cfg->deviceid[i]; i++)
    fullload += cfg->workload[i];

  if (fullload < EPS) {
    for (i = 0; cfg->deviceid[i]; i++)
      cfg->workload[i] = gpu[cfg->deviceid[i]-1].core;
    
  }
  #pragma omp barrier
  fullload=0.f;
  for(i=0;cfg->deviceid[i];i++)
     if(cfg->workload[i]>0.f)
         fullload+=cfg->workload[i];
     else
         mcx_error(-1,"workload was unspecified for an active device",__FILE__,__LINE__);

  
  threadphoton = (int)(cfg->nphoton * cfg->workload[gpuid] /
                       (fullload * gpu[gpuid].autothread * cfg->respin));
  oddphotons =
      (int)(cfg->nphoton * cfg->workload[gpuid] / (fullload * cfg->respin) -
            threadphoton * gpu[gpuid].autothread);
  field = (float *)calloc(sizeof(float) * meshlen, cfg->maxgate);
  dref = (float *)calloc(sizeof(float) * mesh->nf, cfg->maxgate);
  Pdet = (float *)calloc(cfg->maxdetphoton * sizeof(float), hostdetreclen);

  mcgrid.x = gpu[gpuid].autothread / gpu[gpuid].autoblock;
  mcblock.x = gpu[gpuid].autoblock;
  fieldlen = meshlen * cfg->maxgate;

  if (cfg->seed > 0)
    srand(cfg->seed);
  else
    srand(time(0));

  // create gpu pointer
  // gnode,gelem,gtype,gfacenb,gsrcelem,gnormal,gdetpos,gproperty and copy the
  // data from cpu to gpu
  CUDA_ASSERT(cudaMalloc((void **)&gnode, sizeof(float3) * (mesh->nn)));
  CUDA_ASSERT(cudaMemcpy(gnode, mesh->node, sizeof(float3) * (mesh->nn),
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gelem, sizeof(int4) * (mesh->ne)));
  CUDA_ASSERT(cudaMemcpy(gelem, mesh->elem, sizeof(int4) * (mesh->ne),
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gtype, sizeof(int) * (mesh->ne)));
  CUDA_ASSERT(cudaMemcpy(gtype, mesh->type, sizeof(int) * (mesh->ne),
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gfacenb, sizeof(int4) * (mesh->ne)));
  CUDA_ASSERT(cudaMemcpy(gfacenb, mesh->facenb, sizeof(int4) * (mesh->ne),
                         cudaMemcpyHostToDevice));

  if (mesh->srcelemlen > 0) {
    CUDA_ASSERT(cudaMalloc((void **)&gsrcelem, sizeof(int) * (mesh->srcelemlen)));
    CUDA_ASSERT(cudaMemcpy(gsrcelem, mesh->srcelem,
                           sizeof(int) * (mesh->srcelemlen),
                           cudaMemcpyHostToDevice));
  } else
    gsrcelem = NULL;

  CUDA_ASSERT(cudaMalloc((void **)&gnormal, sizeof(float4) * (mesh->ne) * 4));
  CUDA_ASSERT(cudaMemcpy(gnormal, tracer->n, sizeof(float4) * (mesh->ne) * 4,
                         cudaMemcpyHostToDevice));

  // gparam
  CUDA_ASSERT(cudaMemcpyToSymbol(gcfg, &param, sizeof(MCXParam), 0,cudaMemcpyHostToDevice));
  CUDA_ASSERT(cudaMemcpyToSymbol(gmed, mesh->med,
                         (mesh->prop + 1 + cfg->isextdet) * sizeof(Medium),0,
                         cudaMemcpyHostToDevice));

  if (cfg->detpos && cfg->detnum) {
    if((mesh->prop + 1 + cfg->isextdet)+cfg->detnum >= MAX_PROP)
        mcx_error(-5, "Total tissue type and detector count must be less than 2000", __FILE__, __LINE__);
    CUDA_ASSERT(cudaMemcpyToSymbol(gmed, cfg->detpos, 
                         sizeof(float4)*cfg->detnum, (mesh->prop + 1 + cfg->isextdet) * sizeof(Medium), 
			 cudaMemcpyHostToDevice));
  }
  CUDA_ASSERT(cudaMemcpyToSymbol(gmed, tracer->n,
                         (param.normbuf<<2)*(sizeof(float4)), sizeof(float4)*param.maxpropdet,
			 cudaMemcpyHostToDevice));

  // gprogress
  CUDA_ASSERT(
      cudaHostAlloc((void **)&progress, sizeof(int), cudaHostAllocMapped));
  CUDA_ASSERT(cudaHostGetDevicePointer((int **)&gprogress, (int *)progress, 0));
  *progress = 0;

  Pseed =(uint *)malloc(sizeof(uint) * gpu[gpuid].autothread * RAND_SEED_WORD_LEN);
  energy = (float *)calloc(sizeof(float), gpu[gpuid].autothread << 1);
  for (j = 0; j < gpu[gpuid].autothread * RAND_SEED_WORD_LEN; j++)
    Pseed[j] = rand();

  CUDA_ASSERT(cudaMalloc((void **)&gseed, sizeof(uint) * gpu[gpuid].autothread *
                                              RAND_SEED_WORD_LEN));
  CUDA_ASSERT(cudaMemcpy(
      gseed, Pseed, sizeof(uint) * gpu[gpuid].autothread * RAND_SEED_WORD_LEN,
      cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gweight, sizeof(float) * fieldlen));
  CUDA_ASSERT(cudaMemcpy(gweight, field, sizeof(float) * fieldlen,
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gdref, sizeof(float) * nflen));
  CUDA_ASSERT(
      cudaMemcpy(gdref, dref, sizeof(float) * nflen, cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gdetphoton,
                         sizeof(float) * cfg->maxdetphoton * hostdetreclen));
  CUDA_ASSERT(cudaMemcpy(gdetphoton, Pdet,
                         sizeof(float) * cfg->maxdetphoton * hostdetreclen,
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&genergy,
                         sizeof(float) * (gpu[gpuid].autothread << 1)));
  CUDA_ASSERT(cudaMemcpy(genergy, energy,
                         sizeof(float) * (gpu[gpuid].autothread << 1),
                         cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&gdetected, sizeof(uint)));
  CUDA_ASSERT(cudaMemcpy(gdetected, &detected, sizeof(uint), cudaMemcpyHostToDevice));

  CUDA_ASSERT(cudaMalloc((void **)&greporter, sizeof(MCXReporter)));
  CUDA_ASSERT(cudaMemcpy(greporter, &reporter, sizeof(MCXReporter),
                         cudaMemcpyHostToDevice));

  if (cfg->srctype == MCX_SRC_PATTERN) {
    CUDA_ASSERT(cudaMalloc((void **)&gsrcpattern,
                   sizeof(float) * (int)(cfg->srcparam1.w * cfg->srcparam2.w)));
    CUDA_ASSERT(cudaMemcpy(gsrcpattern, cfg->srcpattern,
                   sizeof(float) * (int)(cfg->srcparam1.w * cfg->srcparam2.w),
                   cudaMemcpyHostToDevice));
  } else if (cfg->srctype == MCX_SRC_PATTERN3D) {
    CUDA_ASSERT(cudaMalloc((void **)&gsrcpattern,
                   sizeof(float) * (int)(cfg->srcparam1.x * cfg->srcparam1.y *
                                         cfg->srcparam1.z)));
    CUDA_ASSERT(cudaMemcpy(gsrcpattern, cfg->srcpattern,
                   sizeof(float) * (int)(cfg->srcparam1.x * cfg->srcparam1.y *
                                         cfg->srcparam1.z),
                   cudaMemcpyHostToDevice));
  } else {
    gsrcpattern = NULL;
  }

  if(cfg->issaveseed){
    Pphotonseed=(RandType *)calloc(cfg->maxdetphoton,(sizeof(RandType)*RAND_BUF_LEN));
    CUDA_ASSERT(cudaMalloc((void **)&gphotonseed, cfg->maxdetphoton*(sizeof(RandType)*RAND_BUF_LEN)));
  }
  if(cfg->seed==SEED_FROM_FILE){
    CUDA_ASSERT(cudaMalloc((void **)&greplayweight, sizeof(float)*cfg->nphoton));
    CUDA_ASSERT(cudaMemcpy(greplayweight, cfg->replayweight, sizeof(float)*cfg->nphoton, cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void **)&greplaytime, sizeof(float)*cfg->nphoton));
    CUDA_ASSERT(cudaMemcpy(greplaytime, cfg->replaytime, sizeof(float)*cfg->nphoton, cudaMemcpyHostToDevice));
  
    CUDA_ASSERT(cudaMalloc((void **)&greplayseed, (sizeof(RandType)*RAND_BUF_LEN)*cfg->nphoton));
    CUDA_ASSERT(cudaMemcpy(greplayseed, cfg->photonseed, (sizeof(RandType)*RAND_BUF_LEN)*cfg->nphoton, cudaMemcpyHostToDevice));
  }

  free(Pseed);
  free(energy);
  tic = StartTimer();

  #pragma omp master
  {
    mcx_printheader(cfg);

#ifdef MCX_TARGET_NAME
    MMC_FPRINTF(
        cfg->flog, "- variant name: [%s] compiled by nvcc [%d.%d] with CUDA [%d]\n",
        "Fermi", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, CUDART_VERSION);
#else
    MMC_FPRINTF(
        cfg->flog, "- code name: [Vanilla MCX] compiled by nvcc [%d.%d] with CUDA [%d]\n",
        __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, CUDART_VERSION);
#endif
    MMC_FPRINTF(cfg->flog, "- compiled with: [RNG] %s [Seed Length] %d\n",
                MCX_RNG_NAME, RAND_SEED_WORD_LEN);
    fflush(cfg->flog);
  }
  #pragma omp barrier

  MMC_FPRINTF(cfg->flog,
              "- [device %d(%d): %s] threadph=%d oddphotons=%d np=%.1f "
              "nthread=%d nblock=%d repetition=%d\n",
              gpuid + 1, gpu[gpuid].id, gpu[gpuid].name, threadphoton,
              oddphotons, cfg->nphoton * cfg->workload[gpuid] / fullload,
              (int)gpu[gpuid].autothread, (int)gpu[gpuid].autoblock,
              cfg->respin);

  // simulate for all time-gates in maxgate groups per run

  tic0 = GetTimeMillis();

  for (t = cfg->tstart; t < cfg->tend; t += cfg->tstep * cfg->maxgate) {
    twindow0 = t;
    twindow1 = t + cfg->tstep * cfg->maxgate;

    MMC_FPRINTF(cfg->flog,
                "lauching mcx_main_loop for time window [%.1fns %.1fns] ...\n",
                twindow0 * 1e9, twindow1 * 1e9);
    fflush(cfg->flog);

    // total number of repetition for the simulations, results will be
    // accumulated to field
    for (int iter = 0; iter < cfg->respin; iter++) {
      MMC_FPRINTF(cfg->flog, "simulation run#%2d ... \n", iter + 1);
      fflush(cfg->flog);
      fflush(cfg->flog);
      param.tstart = twindow0;
      param.tend = twindow1;

      // launch mcxkernel
      size_t sharedMemSize = sizeof(int);
      if (cfg->issavedet) {
        sharedMemSize = sizeof(float) * ((int)gpu[gpuid].autoblock) * detreclen;
      }

      mmc_main_loop<<<mcgrid, mcblock, sharedMemSize>>>(
          threadphoton, oddphotons, gnode, (int *)gelem, gweight, gdref,
          gtype, (int *)gfacenb, gsrcelem, gnormal,
          gdetphoton, gdetected, gseed, (int *)gprogress, genergy, greporter,
	  gsrcpattern, greplayweight, greplaytime, greplayseed, gphotonseed);

      #pragma omp master
      {
        if ((cfg->debuglevel & MCX_DEBUG_PROGRESS)) {
          int p0 = 0, ndone = -1;

          progressfun(-0.f,handle);

          do {
            ndone = *progress;

            if (ndone > p0) {
              progressfun((float)ndone / gpu[0].autothread,
                              handle);
              p0 = ndone;
            }
            sleep_ms(100);
          } while (p0 < (int)gpu[0].autothread);
          progressfun(1.f, handle);
          MMC_FPRINTF(cfg->flog, "\n");
        }
      }
      CUDA_ASSERT(cudaDeviceSynchronize());
      tic1 = GetTimeMillis();
      toc += tic1 - tic0;
      MMC_FPRINTF(cfg->flog,
                  "kernel complete:  \t%d ms\nretrieving flux ... \t",
                  tic1 - tic);
      fflush(cfg->flog);
      #pragma omp critical
      if (cfg->runtime < tic1 - tic)
        cfg->runtime = tic1 - tic;

      MCXReporter rep;
      CUDA_ASSERT(cudaMemcpy(&rep, greporter, sizeof(MCXReporter),
                             cudaMemcpyDeviceToHost));
      reporter.raytet += rep.raytet;
      if (cfg->issavedet) {
        CUDA_ASSERT(cudaMemcpy(&detected, gdetected, sizeof(uint), cudaMemcpyDeviceToHost));

        CUDA_ASSERT(cudaMemcpy(Pdet, gdetphoton, sizeof(float) * cfg->maxdetphoton * hostdetreclen,
            cudaMemcpyDeviceToHost));

        if (cfg->issaveseed) {
            CUDA_ASSERT(cudaMemcpy(Pphotonseed, gphotonseed, cfg->maxdetphoton*(sizeof(RandType)*RAND_BUF_LEN),
                cudaMemcpyDeviceToHost));
	}

        if (detected > cfg->maxdetphoton) {
          MMC_FPRINTF(cfg->flog, "WARNING: the detected photon (%d) \
              is more than what your have specified (%d), please use the -H option to specify a greater number\t",
                      detected, cfg->maxdetphoton);
        } else {
          MMC_FPRINTF(cfg->flog, "detected %d photons, total: %d\t", detected,
                      cfg->detectedcount + detected);
        }
        #pragma omp atomic
        cfg->his.detected += detected;
        detected = MIN(detected, cfg->maxdetphoton);
        if (cfg->exportdetected) {
          #pragma omp critical
          {
            cfg->exportdetected = (float *)realloc(
                cfg->exportdetected, (cfg->detectedcount + detected) *
                                         hostdetreclen * sizeof(float));
            memcpy(cfg->exportdetected + cfg->detectedcount * (hostdetreclen),
                   Pdet, detected * (hostdetreclen) * sizeof(float));
            if(cfg->issaveseed){
                 cfg->exportseed=(unsigned char *)realloc(cfg->exportseed,(cfg->detectedcount+detected)*(sizeof(RandType)*RAND_BUF_LEN));
                 memcpy(cfg->exportseed+cfg->detectedcount*sizeof(RandType)*RAND_BUF_LEN,Pphotonseed,detected*(sizeof(RandType)*RAND_BUF_LEN));
            }
            cfg->detectedcount += detected;
          }
        }
      }
      if (cfg->issaveref) {
        float *rawdref = (float *)calloc(sizeof(float), nflen);

        CUDA_ASSERT(cudaMemcpy(rawdref, gdref, sizeof(float) * nflen,
                               cudaMemcpyDeviceToHost));
        for (i = 0; i < nflen; i++) // accumulate field, can be done in the GPU
          dref[i] += rawdref[i];    //+rawfield[i+fieldlen];
        free(rawdref);
      }
      // handling the 2pt distributions
      if (cfg->issave2pt) {
        float *rawfield = (float *)malloc(sizeof(float) * fieldlen);

        CUDA_ASSERT(cudaMemcpy(rawfield, gweight, sizeof(float) * fieldlen,
                               cudaMemcpyDeviceToHost));
        MMC_FPRINTF(cfg->flog, "transfer complete:        %d ms\n",
                    GetTimeMillis() - tic);
        fflush(cfg->flog);

        for (i = 0; i < fieldlen; i++) // accumulate field, can be done in the GPU
          field[(i >> cfg->nbuffer)] += rawfield[i]; //+rawfield[i+fieldlen];

        free(rawfield);

	energy = (float *)calloc(sizeof(float), gpu[gpuid].autothread << 1);

	CUDA_ASSERT(cudaMemcpy(energy, genergy,
	                       sizeof(float) * (gpu[gpuid].autothread << 1),
			       cudaMemcpyDeviceToHost));
	#pragma omp critical
	{
	  for (i = 0; i < gpu[gpuid].autothread; i++) {
	    cfg->energyesc += energy[(i << 1)];
	    cfg->energytot += energy[(i << 1) + 1];
	    // eabsorp+=Plen[i].z;  // the accumulative absorpted energy near
	    // the source
	  }
	}
	free(energy);
      }
      if (cfg->respin > 1 && RAND_SEED_WORD_LEN > 1) {
        Pseed = (uint *)malloc(sizeof(uint) * gpu[gpuid].autothread *
                               RAND_SEED_WORD_LEN);
        for (i = 0; i < gpu[gpuid].autothread * RAND_SEED_WORD_LEN; i++)
          Pseed[i] = rand();
        CUDA_ASSERT(cudaMemcpy(gseed, Pseed,
                               sizeof(uint) * gpu[gpuid].autothread *
                                   RAND_SEED_WORD_LEN,
                               cudaMemcpyHostToDevice));
        free(Pseed);
      }

      // loop over work devices
    } // iteration
  }   // time gates
  #pragma omp master
  {
    int i,j;
    fieldlen = (fieldlen >> cfg->nbuffer);
    field = (float *)realloc(field, sizeof(field[0]) * fieldlen);
    if (cfg->exportfield) {
      if (cfg->basisorder == 0 || cfg->method == rtBLBadouelGrid) {
        for (uint i = 0; i < fieldlen; i++)
          #pragma omp atomic
          cfg->exportfield[i] += field[i];
      } else {
        for (i = 0; i < cfg->maxgate; i++)
          for (j = 0; j < mesh->ne; j++) {
            float ww = field[i * mesh->ne + j] * 0.25f;
            for (int k = 0; k < mesh->elemlen; k++)
              cfg->exportfield[i * mesh->nn +
                               mesh->elem[j * mesh->elemlen + k] - 1] += ww;
          }
      }
    }

    if (cfg->issaveref && mesh->dref) {
      for (uint i = 0; i < nflen; i++)
        mesh->dref[i] += dref[i];
    }

    if (cfg->isnormalized) {
      MMC_FPRINTF(cfg->flog, "normalizing raw data ...\t");
      fflush(cfg->flog);

      cfg->energyabs = cfg->energytot - cfg->energyesc;
      mesh_normalize(mesh, cfg, cfg->energyabs, cfg->energytot, 0);
    }
    if (cfg->issave2pt && cfg->parentid == mpStandalone) {
      MMC_FPRINTF(cfg->flog, "saving data to file ...\t");
      mesh_saveweight(mesh, cfg, 0);
      MMC_FPRINTF(cfg->flog, "saving data complete : %d ms\n\n",
                  GetTimeMillis() - tic);
      fflush(cfg->flog);
    }
    if (cfg->issavedet && cfg->parentid == mpStandalone &&
        cfg->exportdetected) {
      cfg->his.unitinmm = cfg->unitinmm;
      cfg->his.savedphoton = cfg->detectedcount;
      cfg->his.detected = cfg->detectedcount;
      mesh_savedetphoton(cfg->exportdetected, (void*)(cfg->exportseed), cfg->detectedcount,
                         (sizeof(uint64_t) * RAND_BUF_LEN), cfg);
    }
    if (cfg->issaveref) {
      MMC_FPRINTF(cfg->flog, "saving surface diffuse reflectance ...");
      mesh_saveweight(mesh, cfg, 1);
    }
    // total energy here equals total simulated photons+unfinished photons for
    // all threads
    MMC_FPRINTF(cfg->flog,
                "simulated %ld photons (%ld) with devices (ray-tet "
                "%.0f)\nMCX simulation speed: %.2f photon/ms\n",
                cfg->nphoton, cfg->nphoton, reporter.raytet,
                (double)cfg->nphoton / toc);
    MMC_FPRINTF(cfg->flog,
                "total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to "
                "initial specular reflection is excluded in the total)\n",
                cfg->energytot,
                (cfg->energytot - cfg->energyesc) / cfg->energytot * 100.f);
    fflush(cfg->flog);
  }
  #pragma omp barrier
  CUDA_ASSERT(cudaFree(gnode));
  CUDA_ASSERT(cudaFree(gelem));
  CUDA_ASSERT(cudaFree(gtype));
  CUDA_ASSERT(cudaFree(gfacenb));
  CUDA_ASSERT(cudaFree(gsrcelem));
  CUDA_ASSERT(cudaFree(gnormal));
  CUDA_ASSERT(cudaFree(gseed));
  CUDA_ASSERT(cudaFree(gdetphoton));
  CUDA_ASSERT(cudaFree(gweight));
  CUDA_ASSERT(cudaFree(gdref));
  CUDA_ASSERT(cudaFree(genergy));
  CUDA_ASSERT(cudaFree(gdetected));
  if (gsrcpattern)
    CUDA_ASSERT(cudaFree(gsrcpattern));
  if (greplayweight)
    CUDA_ASSERT(cudaFree(greplayweight));
  if (greplayseed)
    CUDA_ASSERT(cudaFree(greplayseed));
  if (greplaytime)
    CUDA_ASSERT(cudaFree(greplaytime));
  if (gphotonseed)
    CUDA_ASSERT(cudaFree(gphotonseed));
  CUDA_ASSERT(cudaFree(greporter));

  #pragma omp master
  if (gpu)
    free(gpu);

  free(field);
  if (Pdet)
    free(Pdet);
  if (Pphotonseed)
    free(Pphotonseed);
  free(dref);
}

void mmc_run_cu(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, void (*progressfun)(float, void *),void *handle){
    GPUInfo *gpuinfo=NULL;        /** gpuinfo: structure to store GPU information */
    unsigned int activedev=0;     /** activedev: count of total active GPUs to be used */
    if(!(activedev=mcx_list_cu_gpu(cfg,&gpuinfo))){
         mcx_error(-1,"No GPU device found\n",__FILE__,__LINE__);
    }

#ifdef _OPENMP
    /** 
        Now we are ready to launch one thread for each involked GPU to run the simulation 
     */
    omp_set_num_threads(activedev);
    #pragma omp parallel
    {
#endif

    /** 
        This line runs the main MCX simulation for each GPU inside each thread 
     */
    mmc_run_simulation(cfg,mesh,tracer,gpuinfo,progressfun,handle);

#ifdef _OPENMP
    }
#endif
}
