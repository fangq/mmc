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

/***************************************************************************//**
\file    mmc_cu_host.c

\brief   CUDA host code for NVIDIA GPUs
*******************************************************************************/


#define inlinefun __device__

#include "mmc_const.h"
#include "mmc_cu_host.h"
#include "mmc_tictoc.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "mmc_const.h"

#include "mmc_core.cu"

/******************************************************************************
In this unit, we first launch a master thread and initialize the necessary data
structures.This include the command line options(cfg), tetrahedral mesh(mesh)
and the ray tracer precomputed data (tracer).
******************************************************************************/
#define CUDA_ASSERT(a)                                                         \
    mcx_cu_assess((a), __FILE__, __LINE__) ///< macro to report CUDA error

/**
 * @brief Utility function to calculate the GPU stream processors (cores) per SM
 *
 * Obtain GPU core number per MP, this replaces
 * ConvertSMVer2Cores() in libcudautils to avoid
 * extra dependency.
 *
 * @param[in] v1: the major version of an NVIDIA GPU
 * @param[in] v2: the minor version of an NVIDIA GPU
 */

int mcx_corecount(int v1, int v2) {
    int v = v1 * 10 + v2;

    if (v < 20) {
        return 8;
    } else if (v < 21) {
        return 32;
    } else if (v < 30) {
        return 48;
    } else if (v < 50) {
        return 192;
    } else if (v < 60 || v == 61 || v >= 89) {
        return 128;
    } else {
        return 64;
    }
}

int mcx_smxblock(int v1, int v2) {
    int v = v1 * 10 + v2;

    if (v < 30) {
        return 8;
    } else if (v < 50) {
        return 16;
    } else {
        return 32;
    }
}

/**
  assert cuda memory allocation result
 */
void mcx_cu_assess(cudaError_t cuerr, const char* file, const int linenum) {
    if (cuerr != cudaSuccess) {
        mcx_error(-(int)cuerr, (char*)cudaGetErrorString(cuerr), file, linenum);
    }
}

/*
   master driver code to run MC simulations
*/
int mcx_list_cu_gpu(mcconfig* cfg, GPUInfo** info) {
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

    *info = (GPUInfo*)calloc(deviceCount, sizeof(GPUInfo));

    if (cfg->gpuid && cfg->gpuid > (uint)deviceCount) {
        MMC_FPRINTF(stderr,
                    S_RED "ERROR: Specified GPU ID is out of range\n" S_RESET);
        return 0;
    }

    // scan from the first device
    for (dev = 0; dev < deviceCount; dev++) {
        cudaDeviceProp dp;
        CUDA_ASSERT(cudaGetDeviceProperties(&dp, dev));

        if (cfg->isgpuinfo == 3) {
            activedev++;
        } else if (cfg->deviceid[dev] == '1') {
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

void mmc_run_simulation(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, GPUInfo* gpu) {
    uint i, j;
    float t, twindow0, twindow1;
    float fullload = 0.f;
    float* energy;

    uint detected = 0;
    int gpuid, threadid = 0;
    uint tic, tic0, tic1, toc = 0, fieldlen, debuglen = MCX_DEBUG_REC_LEN;
    int threadphoton, oddphotons;
    dim3 mcgrid, mcblock;

    float3* gnode;
    int4* gelem, *gfacenb;
    float4* gnormal;
    int* gtype, *gsrcelem;
    /* per-node optical property arrays for DOT reconstruction (NULL when unused) */
    float* gnodemua_cu = NULL;
    float* gnodemusp_cu = NULL;
    uint* gseed, *gdetected;
    volatile int* progress, *gprogress;
    float* gweight, *gdref, *gdetphoton, *genergy, *gsrcpattern, *gdebugdata;
    RandType* gphotonseed = NULL, *greplayseed = NULL;
    float*  greplayweight = NULL, *greplaytime = NULL;

    MCXReporter* greporter;
    int isrfforward = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE);
    /* J_mus/J_musp (U, V and the second half of Q) multiply J_D by a per-medium chain-rule factor;
       in grid mode the kernel records the medium label of each voxel it deposits into */
    int ismusjacob = (cfg->outputtype == otAdjointMus || cfg->outputtype == otAdjointMusp || cfg->outputtype == otAdjointMuaMusp);
    int* gvoxtype = NULL;
    int* groioff = NULL;          /* iMMC: ROI lists (edge/node) or slabs (face), see mesh_roilist() */
    float4* groilist = NULL;
    float* gfaceroi = NULL;
    /* For adjoint/multi-source mode, expand meshlen to cover all source slots.
     * cfg->srcid > 0 selects a single slot from srcdata, so the field buffer collapses to one slot. */
    uint nsrcslots = (cfg->extrasrclen > 0 && cfg->srcid <= 0) ? (uint)cfg->extrasrclen : 1u;
    uint meshlen = ((cfg->method == rtBLBadouelGrid) ? cfg->crop0.z : mesh->ne) * cfg->srcnum * nsrcslots;
    cfg->crop0.w = meshlen * cfg->maxgate; // offset for the second (double) buffer

    float* field, *dref = NULL;

    uint* Pseed = NULL;
    float* Pdet = NULL;
    RandType* Pphotonseed = NULL;

    uint detreclen = (2 + ((cfg->ismomentum) > 0)) * mesh->prop +
                     (cfg->issaveexit > 0) * 6 + 1;
    uint hostdetreclen = detreclen + 1;
    // launch mcxkernel
    size_t sharedmemsize = 0;
    double energytot = 0.0, energyesc = 0.0;

    MCXParam param = {
        cfg->srcparam1,
        cfg->srcparam2,
        cfg->crop0,
        cfg->bary0,
        make_float3(cfg->srcpos.x, cfg->srcpos.y, cfg->srcpos.z),
        make_float3(cfg->srcdir.x, cfg->srcdir.y, cfg->srcdir.z),
        make_float3(mesh->nmin.x, mesh->nmin.y, mesh->nmin.z),
        cfg->tstart,
        cfg->tend,
        (uint)cfg->isreflect,
        (uint)cfg->issavedet,
        (uint)cfg->issaveexit,
        (uint)cfg->ismomentum,
        (uint)cfg->isspecular,
        1.f / cfg->tstep,
        cfg->minenergy,
        cfg->maxdetphoton,
        (uint)mesh->prop,
        (uint)cfg->detnum,
        (int)cfg->voidtime,
        (int)cfg->srctype,
        (uint)cfg->issaveref,
        (uint)cfg->maxgate,
        (uint)cfg->debuglevel,
        (int)detreclen,
        cfg->outputtype,
        mesh->elemlen,
        cfg->mcmethod,
        cfg->method,
        1.f / cfg->steps.x,
        cfg->srcdir.w,
        mesh->nn,
        mesh->ne,
        mesh->nf,
        cfg->nout,
        cfg->roulettesize,
        cfg->srcnum,
        mesh->srcelemlen,
        cfg->e0,
        cfg->isextdet,
        (uint)(meshlen / (cfg->srcnum * nsrcslots)),  /* framelen = datalen (per-gate stride); slot stride is datalen*maxgate */
        (uint)(mesh->prop + 1 + cfg->isextdet) + (uint)(cfg->extrasrclen * 4) + cfg->detnum,
        0u, /*no element normals in __constant memory: the lanes of a warp read those of different elements, which serializes*/
        cfg->issaveseed,
        cfg->seed,
        cfg->maxjumpdebug,
        cfg->omega,
        (float)(3.335640951981520e-12),   /* oneoverc0 = 1/C0 in s/mm */
        (cfg->extrasrclen > 0 && cfg->srcid == 0) ? -1 : cfg->srcid,  /* srcid<0 multi-slot, srcid>0 single slot, srcid==0 with srcdata: auto-promote to -1 */
        cfg->extrasrclen,
        (int)(mesh->prop + 1 + cfg->isextdet),  /* srcpropoffset: gmed index where extra sources start */
        (uint)cfg->isnodalmua,
        (uint)cfg->isnodalmusp
    };

    MCXReporter reporter = {0.f, 0};

    if (cfg->issavedet) {
        sharedmemsize = sizeof(float) * detreclen;
    }

    param.reclen = sharedmemsize / sizeof(float);  //< the shared memory buffer length associated with detected photon

    sharedmemsize += sizeof(float) * (cfg->srcnum << 1);   /**< store energyesc/energytot */

    if (cfg->srctype == stPattern && cfg->srcnum > 1) {
        sharedmemsize += sizeof(float) * cfg->srcnum;
    }

#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif

    if (threadid < MAX_DEVICE && cfg->deviceid[threadid] == '\0') {
        return;
    }

    /* each OpenMP thread drives the GPU at its own slot in deviceid */
    gpuid = cfg->deviceid[threadid] - 1;

    sharedmemsize *= ((int)gpu[gpuid].autoblock);

    if (gpuid < 0) {
        mcx_error(-1, "GPU ID must be non-zero", __FILE__, __LINE__);
    }

    CUDA_ASSERT(cudaSetDevice(gpuid));

    #pragma omp master
    {
        if (cfg->exportfield == NULL) {
            cfg->exportfield = mesh->weight;
        }

        /* Pre-allocate exportadjoint for RF imaginary part or adjoint Jacobian */
        if (cfg->exportadjoint == NULL && isrfforward) {
            /* RF forward: the imaginary fluence, laid out as the real one; as in the
               OpenCL host, sized per node when basis order 1 spreads the per-element
               field onto nodes (nn can exceed ne).  Element- and node-based output
               lost the imaginary part when this was allocated for the grid alone. */
            size_t nodelen = (cfg->method != rtBLBadouelGrid && cfg->basisorder) ?
                             (size_t)mesh->nn * cfg->maxgate * cfg->srcnum * nsrcslots : 0;
            cfg->exportadjoint = (float*)calloc(MAX((size_t)meshlen * cfg->maxgate, nodelen), sizeof(float));  /* fieldlen is set later */
        }

        if (cfg->exportdetected == NULL) {
            cfg->exportdetected = (float*)malloc(hostdetreclen * cfg->maxdetphoton * sizeof(float));
        }

        if (cfg->issaveseed && cfg->exportseed == NULL) {
            cfg->exportseed = (unsigned char*)malloc(cfg->maxdetphoton * (sizeof(RandType) * RAND_BUF_LEN));
        }

        cfg->energytot = (double*)calloc(cfg->srcnum, sizeof(double));
        cfg->energyesc = (double*)calloc(cfg->srcnum, sizeof(double));
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

        if (cfg->maxgate == 0) {
            cfg->maxgate = gates;
        } else if ((uint)cfg->maxgate > gates) {
            cfg->maxgate = gates;
        }

        gpu[gpuid].maxgate = cfg->maxgate;
    }

    if (gpu[gpuid].autothread % gpu[gpuid].autoblock) {
        gpu[gpuid].autothread =
            (gpu[gpuid].autothread / gpu[gpuid].autoblock) * gpu[gpuid].autoblock;
    }

    param.maxgate = gpu[gpuid].maxgate;

    uint nflen = mesh->nf * cfg->maxgate;
    #pragma omp master
    fullload = 0.f;

    for (i = 0; cfg->deviceid[i]; i++) {
        fullload += cfg->workload[i];
    }

    if (fullload < EPS) {
        for (i = 0; cfg->deviceid[i]; i++) {
            cfg->workload[i] = gpu[cfg->deviceid[i] - 1].core;
        }

    }

    #pragma omp barrier
    fullload = 0.f;

    for (i = 0; cfg->deviceid[i]; i++) {
        if (cfg->workload[i] > 0.f) {
            fullload += cfg->workload[i];
        } else {
            mcx_error(-1, "workload was unspecified for an active device", __FILE__, __LINE__);
        }
    }

    threadphoton = (int)(cfg->nphoton * cfg->workload[threadid] /
                         (fullload * gpu[gpuid].autothread * cfg->respin));
    oddphotons =
        (int)(cfg->nphoton * cfg->workload[threadid] / (fullload * cfg->respin) -
              threadphoton * gpu[gpuid].autothread);
    field = (float*)calloc(sizeof(float) * meshlen * (isrfforward ? 4 : 2), cfg->maxgate);
    dref = (float*)calloc(sizeof(float) * mesh->nf, cfg->maxgate);
    Pdet = (float*)calloc(cfg->maxdetphoton * sizeof(float), hostdetreclen);

    mcgrid.x = gpu[gpuid].autothread / gpu[gpuid].autoblock;
    mcblock.x = gpu[gpuid].autoblock;
    fieldlen = meshlen * cfg->maxgate;

    MMCHostRNG hostrng;  /* fills the per-thread RNG seeds */

    if (cfg->seed > 0) {
        mmc_hostrng_seed(&hostrng, cfg->seed);
    } else {
        mmc_hostrng_seed(&hostrng, time(0));
    }

    // create gpu pointer
    // gnode,gelem,gtype,gfacenb,gsrcelem,gnormal,gdetpos,gproperty and copy the
    // data from cpu to gpu
    CUDA_ASSERT(cudaMalloc((void**)&gnode, sizeof(float3) * (mesh->nn)));
    CUDA_ASSERT(cudaMemcpy(gnode, mesh->node, sizeof(float3) * (mesh->nn),
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gelem, sizeof(int4) * (mesh->ne)));
    CUDA_ASSERT(cudaMemcpy(gelem, mesh->elem, sizeof(int4) * (mesh->ne),
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gtype, sizeof(int) * (mesh->ne)));
    CUDA_ASSERT(cudaMemcpy(gtype, mesh->type, sizeof(int) * (mesh->ne),
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gfacenb, sizeof(int4) * (mesh->ne)));
    CUDA_ASSERT(cudaMemcpy(gfacenb, mesh->facenb, sizeof(int4) * (mesh->ne),
                           cudaMemcpyHostToDevice));

    if (mesh->srcelemlen > 0) {
        CUDA_ASSERT(cudaMalloc((void**)&gsrcelem, sizeof(int) * (mesh->srcelemlen)));
        CUDA_ASSERT(cudaMemcpy(gsrcelem, mesh->srcelem,
                               sizeof(int) * (mesh->srcelemlen),
                               cudaMemcpyHostToDevice));
    } else {
        gsrcelem = NULL;
    }

    CUDA_ASSERT(cudaMalloc((void**)&gnormal, sizeof(float4) * (mesh->ne) * 4));
    CUDA_ASSERT(cudaMemcpy(gnormal, tracer->n, sizeof(float4) * (mesh->ne) * 4,
                           cudaMemcpyHostToDevice));

    /* upload per-node mua/musp for DOT reconstruction modes (size mesh->nn) */
    if (cfg->isnodalmua && cfg->nodemua) {
        CUDA_ASSERT(cudaMalloc((void**)&gnodemua_cu, sizeof(float) * (mesh->nn)));
        CUDA_ASSERT(cudaMemcpy(gnodemua_cu, cfg->nodemua, sizeof(float) * (mesh->nn),
                               cudaMemcpyHostToDevice));
    }

    if (cfg->isnodalmusp && cfg->nodemusp) {
        CUDA_ASSERT(cudaMalloc((void**)&gnodemusp_cu, sizeof(float) * (mesh->nn)));
        CUDA_ASSERT(cudaMemcpy(gnodemusp_cu, cfg->nodemusp, sizeof(float) * (mesh->nn),
                               cudaMemcpyHostToDevice));
    }

    // gparam
    CUDA_ASSERT(cudaMemcpyToSymbol(gcfg, &param, sizeof(MCXParam), 0, cudaMemcpyHostToDevice));

    if (cfg->implicit) {
        float2 immcparam = make_float2(mesh->roieps, (float)(mesh->nslab1 > 0));
        CUDA_ASSERT(cudaMemcpyToSymbol(gimmc, &immcparam, sizeof(float2), 0, cudaMemcpyHostToDevice));
    }

    /* Upload media, extra sources, and detectors to gmed[] constant memory.
     * Layout: [media 0..prop] | [extra sources, each 4 Medium slots] | [detector positions] | [normals] */
    uint srcpropoffset = (uint)(mesh->prop + 1 + cfg->isextdet);
    uint detpropoffset = srcpropoffset + (uint)(cfg->extrasrclen * 4);

    if (detpropoffset + cfg->detnum + (param.normbuf << 2) > MAX_PROP) {
        mcx_error(-5, "Total tissue types, extra sources, detectors and normals must fit in MAX_PROP slots", __FILE__, __LINE__);
    }

    CUDA_ASSERT(cudaMemcpyToSymbol(gmed, mesh->med,
                                   srcpropoffset * sizeof(Medium), 0,
                                   cudaMemcpyHostToDevice));

    /* Pack extra sources into gmed[] right after media; each ExtraSrc = 4 Medium slots */
    if (cfg->extrasrclen > 0 && cfg->srcdata) {
        CUDA_ASSERT(cudaMemcpyToSymbol(gmed, cfg->srcdata,
                                       sizeof(ExtraSrc) * cfg->extrasrclen,
                                       srcpropoffset * sizeof(Medium), cudaMemcpyHostToDevice));
    }

    if (cfg->detpos && cfg->detnum) {
        CUDA_ASSERT(cudaMemcpyToSymbol(gmed, cfg->detpos,
                                       sizeof(float4) * cfg->detnum,
                                       detpropoffset * sizeof(Medium),
                                       cudaMemcpyHostToDevice));
    }

    CUDA_ASSERT(cudaMemcpyToSymbol(gmed, tracer->n,
                                   (param.normbuf << 2) * (sizeof(float4)), sizeof(float4)*param.maxpropdet,
                                   cudaMemcpyHostToDevice));

    // gprogress
    CUDA_ASSERT(
        cudaHostAlloc((void**)&progress, sizeof(int), cudaHostAllocMapped));
    CUDA_ASSERT(cudaHostGetDevicePointer((int**)&gprogress, (int*)progress, 0));
    *progress = 0;

    Pseed = (uint*)malloc(sizeof(uint) * gpu[gpuid].autothread * RAND_SEED_WORD_LEN);
    energy = (float*)calloc(sizeof(float) * cfg->srcnum, gpu[gpuid].autothread << 1);

    for (j = 0; j < gpu[gpuid].autothread * RAND_SEED_WORD_LEN; j++) {
        Pseed[j] = mmc_rand_uint32(&hostrng);
    }

    CUDA_ASSERT(cudaMalloc((void**)&gseed, sizeof(uint) * gpu[gpuid].autothread*
                           RAND_SEED_WORD_LEN));
    CUDA_ASSERT(cudaMemcpy(
                    gseed, Pseed, sizeof(uint) * gpu[gpuid].autothread * RAND_SEED_WORD_LEN,
                    cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gweight, sizeof(float) * fieldlen * (isrfforward ? 4 : 2)));
    CUDA_ASSERT(cudaMemcpy(gweight, field, sizeof(float) * fieldlen * (isrfforward ? 4 : 2),
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gdref, sizeof(float) * nflen));
    CUDA_ASSERT(
        cudaMemcpy(gdref, dref, sizeof(float) * nflen, cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gdetphoton,
                           sizeof(float) * cfg->maxdetphoton * hostdetreclen));
    CUDA_ASSERT(cudaMemcpy(gdetphoton, Pdet,
                           sizeof(float) * cfg->maxdetphoton * hostdetreclen,
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&genergy,
                           sizeof(float) * (gpu[gpuid].autothread << 1) * cfg->srcnum));
    CUDA_ASSERT(cudaMemcpy(genergy, energy,
                           sizeof(float) * (gpu[gpuid].autothread << 1) * cfg->srcnum,
                           cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&gdetected, sizeof(uint)));
    CUDA_ASSERT(cudaMemcpy(gdetected, &detected, sizeof(uint), cudaMemcpyHostToDevice));

    CUDA_ASSERT(cudaMalloc((void**)&greporter, sizeof(MCXReporter)));
    CUDA_ASSERT(cudaMemcpy(greporter, &reporter, sizeof(MCXReporter),
                           cudaMemcpyHostToDevice));

    if (cfg->srctype == MCX_SRC_PATTERN) {
        CUDA_ASSERT(cudaMalloc((void**)&gsrcpattern,
                               sizeof(float) * (int)(cfg->srcparam1.w * cfg->srcparam2.w * cfg->srcnum)));
        CUDA_ASSERT(cudaMemcpy(gsrcpattern, cfg->srcpattern,
                               sizeof(float) * (int)(cfg->srcparam1.w * cfg->srcparam2.w * cfg->srcnum),
                               cudaMemcpyHostToDevice));
    } else if (cfg->srctype == MCX_SRC_PATTERN3D) {
        CUDA_ASSERT(cudaMalloc((void**)&gsrcpattern,
                               sizeof(float) * (int)(cfg->srcparam1.x * cfg->srcparam1.y*
                                       cfg->srcparam1.z * cfg->srcnum)));
        CUDA_ASSERT(cudaMemcpy(gsrcpattern, cfg->srcpattern,
                               sizeof(float) * (int)(cfg->srcparam1.x * cfg->srcparam1.y*
                                       cfg->srcparam1.z * cfg->srcnum),
                               cudaMemcpyHostToDevice));
    } else {
        gsrcpattern = NULL;
    }

    if (cfg->issaveseed) {
        Pphotonseed = (RandType*)calloc(cfg->maxdetphoton, (sizeof(RandType) * RAND_BUF_LEN));
        CUDA_ASSERT(cudaMalloc((void**)&gphotonseed, cfg->maxdetphoton * (sizeof(RandType)*RAND_BUF_LEN)));
    }

    if (cfg->debuglevel & dlTraj) {
        CUDA_ASSERT(cudaMalloc((void**) &gdebugdata, sizeof(float) * (debuglen * cfg->maxjumpdebug)));
    }

    if (cfg->implicit == 1 && mesh->roioff) {
        CUDA_ASSERT(cudaMalloc((void**) &groioff, sizeof(int) * mesh->roioff[mesh->ne]));
        CUDA_ASSERT(cudaMemcpy(groioff, mesh->roioff, sizeof(int) * mesh->roioff[mesh->ne], cudaMemcpyHostToDevice));
        CUDA_ASSERT(cudaMalloc((void**) &groilist, sizeof(float) * 8 * MAX(mesh->nroiprim, 1)));
        CUDA_ASSERT(cudaMemcpy(groilist, mesh->roilist, sizeof(float) * 8 * mesh->nroiprim, cudaMemcpyHostToDevice));
    } else if (cfg->implicit == 2 && mesh->faceroi) {
        CUDA_ASSERT(cudaMalloc((void**) &gfaceroi, sizeof(float) * 4 * mesh->ne));
        CUDA_ASSERT(cudaMemcpy(gfaceroi, mesh->faceroi, sizeof(float) * 4 * mesh->ne, cudaMemcpyHostToDevice));
    }

    if (ismusjacob && cfg->issave2pt && cfg->method == rtBLBadouelGrid) {
        CUDA_ASSERT(cudaMalloc((void**) &gvoxtype, sizeof(int) * cfg->crop0.z));
        CUDA_ASSERT(cudaMemset(gvoxtype, 0, sizeof(int) * cfg->crop0.z));
    }

    if (cfg->seed == SEED_FROM_FILE) {
        CUDA_ASSERT(cudaMalloc((void**)&greplayweight, sizeof(float)*cfg->nphoton));
        CUDA_ASSERT(cudaMemcpy(greplayweight, cfg->replayweight, sizeof(float)*cfg->nphoton, cudaMemcpyHostToDevice));

        CUDA_ASSERT(cudaMalloc((void**)&greplaytime, sizeof(float)*cfg->nphoton));
        CUDA_ASSERT(cudaMemcpy(greplaytime, cfg->replaytime, sizeof(float)*cfg->nphoton, cudaMemcpyHostToDevice));

        CUDA_ASSERT(cudaMalloc((void**)&greplayseed, (sizeof(RandType)*RAND_BUF_LEN)*cfg->nphoton));
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
            cfg->flog, "- code name: [%s] compiled by nvcc [%d.%d] with CUDA [%d]\n",
            "MMC-Trinity", __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, CUDART_VERSION);
#else
        MMC_FPRINTF(
            cfg->flog, "- code name: [MMC-Trinity] compiled by nvcc [%d.%d] with CUDA [%d]\n",
            __CUDACC_VER_MAJOR__, __CUDACC_VER_MINOR__, CUDART_VERSION);
#endif
        MMC_FPRINTF(cfg->flog, "- compiled with: [RNG] %s [Seed Length] %d\n",
                    MCX_RNG_NAME, RAND_SEED_WORD_LEN);
        mcx_fflush(cfg->flog);
    }
    #pragma omp barrier

    MMC_FPRINTF(cfg->flog,
                "- [device %d(%d): %s] threadph=%d oddphotons=%d np=%.1f "
                "nthread=%d nblock=%d repetition=%d\n",
                gpuid + 1, gpu[gpuid].id, gpu[gpuid].name, threadphoton,
                oddphotons, cfg->nphoton * cfg->workload[threadid] / fullload,
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

        MMC_FPRINTF(cfg->flog, "requesting %ld bytes of shared memory\n", sharedmemsize);

        mcx_fflush(cfg->flog);

        // total number of repetition for the simulations, results will be
        // accumulated to field
        for (int iter = 0; iter < cfg->respin; iter++) {
            MMC_FPRINTF(cfg->flog, "simulation run#%2d ... \n", iter + 1);
            mcx_fflush(cfg->flog);
            mcx_fflush(cfg->flog);
            param.tstart = twindow0;
            param.tend = twindow1;


            /* Template-dispatch on persistent per-photon state: each flag gates
             * register-resident state across the photon lifetime, letting ptxas
             * dead-code-eliminate unused features. NODAL_USE_MUA/MUSP are now
             * runtime reads from gcfg (their lifetime is too short to affect
             * peak register pressure). 8 instantiations enumerated below.
             *
             *   IS_RF          - cfg.omega>0 and not replay
             *   IS_MULTISRC    - cfg.srcnum>1 or (extrasrclen>0 && srcid<=0)
             *   SAVE_DETPHOTON - cfg.issavedet
             *   IMMC_MODE      - cfg.implicit (0, 1 edge/node, 2 face), 3 x 8 instantiations */
            const int is_rf = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE) ? 1 : 0;
            /* IS_MULTISRC includes pattern source (uses r->posidx for srcpattern[] lookup)
             * alongside pattern-multi-source (srcnum>1) and adjoint multi-source. */
            const int is_ms = (cfg->srctype == stPattern || cfg->srctype == stPattern3D
                               || cfg->srcnum > 1
                               || (cfg->extrasrclen > 0 && cfg->srcid <= 0)) ? 1 : 0;
            const int save_dp = cfg->issavedet ? 1 : 0;
            const int variant = (is_rf << 2) | (is_ms << 1) | save_dp;

#define MMC_DISPATCH_IM(RF, MS, SD, IM) \
    mmc_main_loop<RF, MS, SD, IM> <<< mcgrid, mcblock, sharedmemsize>>>( \
            threadphoton, oddphotons, gnode, (int*)gelem, gweight, gdref, \
            gtype, (int*)gfacenb, gsrcelem, gnormal, \
            gnodemua_cu, gnodemusp_cu, \
            gdetphoton, gdetected, gseed, (int*)gprogress, genergy, greporter, \
            gsrcpattern, greplayweight, greplaytime, greplayseed, gphotonseed, gdebugdata, gvoxtype, \
            groioff, groilist, gfaceroi)
#define MMC_DISPATCH(RF, MS, SD) \
    switch (cfg->implicit) { \
        case 1: MMC_DISPATCH_IM(RF, MS, SD, 1); break; \
        case 2: MMC_DISPATCH_IM(RF, MS, SD, 2); break; \
        default: MMC_DISPATCH_IM(RF, MS, SD, 0); \
    }

            switch (variant) {
                case 0:
                    MMC_DISPATCH(0, 0, 0);
                    break;

                case 1:
                    MMC_DISPATCH(0, 0, 1);
                    break;

                case 2:
                    MMC_DISPATCH(0, 1, 0);
                    break;

                case 3:
                    MMC_DISPATCH(0, 1, 1);
                    break;

                case 4:
                    MMC_DISPATCH(1, 0, 0);
                    break;

                case 5:
                    MMC_DISPATCH(1, 0, 1);
                    break;

                case 6:
                    MMC_DISPATCH(1, 1, 0);
                    break;

                case 7:
                    MMC_DISPATCH(1, 1, 1);
                    break;
            }

#undef MMC_DISPATCH
#undef MMC_DISPATCH_IM

            #pragma omp master
            {
                if ((cfg->debuglevel & MCX_DEBUG_PROGRESS)) {
                    int p0 = 0, ndone = -1;

                    mcx_progressbar(-0.f);

                    do {
                        ndone = *progress;

                        if (ndone > p0) {
                            mcx_progressbar((float)ndone / gpu[0].autothread);
                            p0 = ndone;
                        }

                        sleep_ms(100);
                    } while (p0 < (int)gpu[0].autothread);

                    mcx_progressbar(1.f);
                    MMC_FPRINTF(cfg->flog, "\n");
                }
            }
            CUDA_ASSERT(cudaDeviceSynchronize());
            tic1 = GetTimeMillis();
            toc += tic1 - tic0;
            MMC_FPRINTF(cfg->flog,
                        "kernel complete:  \t%d ms\nretrieving flux ... \t",
                        tic1 - tic);
            mcx_fflush(cfg->flog);
            #pragma omp critical

            if (cfg->runtime < tic1 - tic) {
                cfg->runtime = tic1 - tic;
            }

            MCXReporter rep;
            CUDA_ASSERT(cudaMemcpy(&rep, greporter, sizeof(MCXReporter),
                                   cudaMemcpyDeviceToHost));
            reporter.raytet += rep.raytet;
            reporter.jumpdebug += rep.jumpdebug;

            energy = (float*)calloc(sizeof(float) * cfg->srcnum, gpu[gpuid].autothread << 1);

            CUDA_ASSERT(cudaMemcpy(energy, genergy,
                                   sizeof(float) * (gpu[gpuid].autothread << 1) * cfg->srcnum,
                                   cudaMemcpyDeviceToHost));
            #pragma omp critical
            {

                for (i = 0; i < gpu[gpuid].autothread; i++) {
                    for (j = 0; j < (uint)cfg->srcnum; j++) {
                        cfg->energyesc[j] += energy[(i << 1) * cfg->srcnum + j];
                        cfg->energytot[j] += energy[((i << 1) + 1) * cfg->srcnum + j];
                        energyesc += energy[(i << 1) * cfg->srcnum + j];
                        energytot += energy[((i << 1) + 1) * cfg->srcnum + j];
                    }
                }
            }

            free(energy);

            /**
             * If '-D M' is specified, we retrieve photon trajectory data and store those to \c cfg.exportdebugdata and \c cfg.debugdatalen
             */
            if (cfg->debuglevel & dlTraj) {
                uint debugrec = reporter.jumpdebug;

                #pragma omp critical
                {
                    if (debugrec > 0) {
                        if (debugrec > cfg->maxjumpdebug) {
                            MMC_FPRINTF(cfg->flog, S_RED "WARNING: the saved trajectory positions (%d) \
are more than what your have specified (%d), please use the --maxjumpdebug option to specify a greater number\n" S_RESET
                                        , debugrec, cfg->maxjumpdebug);
                        } else {
                            MMC_FPRINTF(cfg->flog, "saved %u trajectory positions, total: %d\t", debugrec, cfg->debugdatalen + debugrec);
                        }

                        debugrec = min(debugrec, cfg->maxjumpdebug);
                        cfg->exportdebugdata = (float*)realloc(cfg->exportdebugdata, (cfg->debugdatalen + debugrec) * debuglen * sizeof(float));
                        CUDA_ASSERT(cudaMemcpy(cfg->exportdebugdata + cfg->debugdatalen, gdebugdata, sizeof(float)*debuglen * debugrec, cudaMemcpyDeviceToHost));
                        cfg->debugdatalen += debugrec;
                    }
                }
            }

            if (cfg->issavedet) {
                CUDA_ASSERT(cudaMemcpy(&detected, gdetected, sizeof(uint), cudaMemcpyDeviceToHost));

                CUDA_ASSERT(cudaMemcpy(Pdet, gdetphoton, sizeof(float) * cfg->maxdetphoton * hostdetreclen,
                                       cudaMemcpyDeviceToHost));

                if (cfg->issaveseed) {
                    CUDA_ASSERT(cudaMemcpy(Pphotonseed, gphotonseed, cfg->maxdetphoton * (sizeof(RandType)*RAND_BUF_LEN),
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
                        cfg->exportdetected = (float*)realloc(
                            cfg->exportdetected, (cfg->detectedcount + detected) *
                            hostdetreclen * sizeof(float));
                        memcpy(cfg->exportdetected + cfg->detectedcount * (hostdetreclen),
                               Pdet, detected * (hostdetreclen) * sizeof(float));

                        if (cfg->issaveseed) {
                            cfg->exportseed = (unsigned char*)realloc(cfg->exportseed, (cfg->detectedcount + detected) * (sizeof(RandType) * RAND_BUF_LEN));
                            memcpy(cfg->exportseed + cfg->detectedcount * sizeof(RandType)*RAND_BUF_LEN, Pphotonseed, detected * (sizeof(RandType)*RAND_BUF_LEN));
                        }

                        cfg->detectedcount += detected;
                    }
                }
            }

            if (cfg->issaveref) {
                float* rawdref = (float*)calloc(sizeof(float), nflen);

                CUDA_ASSERT(cudaMemcpy(rawdref, gdref, sizeof(float) * nflen,
                                       cudaMemcpyDeviceToHost));

                for (i = 0; i < nflen; i++) { // accumulate field, can be done in the GPU
                    dref[i] += rawdref[i];    //+rawfield[i+fieldlen];
                }

                free(rawdref);
            }

            // handling the 2pt distributions
            if (cfg->issave2pt) {
                int rfmul = isrfforward ? 4 : 2;
                float* rawfield = (float*)malloc(sizeof(float) * fieldlen * rfmul);

                CUDA_ASSERT(cudaMemcpy(rawfield, gweight, sizeof(float) * fieldlen * rfmul,
                                       cudaMemcpyDeviceToHost));
                MMC_FPRINTF(cfg->flog, "transfer complete:        %d ms\n",
                            GetTimeMillis() - tic);
                mcx_fflush(cfg->flog);

                for (i = 0; i < fieldlen; i++) { // accumulate real field (double-buffer re1+re2)
                    field[i] += rawfield[i] + rawfield[i + fieldlen];
                }

                if (isrfforward) {
                    /* accumulate imaginary part from im1+im2 into field[fieldlen..] */
                    for (i = 0; i < fieldlen; i++) {
                        field[i + fieldlen] += rawfield[i + fieldlen * 2] + rawfield[i + fieldlen * 3];
                    }
                }

                free(rawfield);
            }

            if (cfg->respin > 1 && RAND_SEED_WORD_LEN > 1) {
                Pseed = (uint*)malloc(sizeof(uint) * gpu[gpuid].autothread *
                                      RAND_SEED_WORD_LEN);

                for (i = 0; i < gpu[gpuid].autothread * RAND_SEED_WORD_LEN; i++) {
                    Pseed[i] = mmc_rand_uint32(&hostrng);
                }

                CUDA_ASSERT(cudaMemcpy(gseed, Pseed,
                                       sizeof(uint) * gpu[gpuid].autothread*
                                       RAND_SEED_WORD_LEN,
                                       cudaMemcpyHostToDevice));
                free(Pseed);
            }

            // loop over work devices
        } // iteration
    }   // time gates

    /* with multiple GPUs, every thread adds the field of its own device to the output */
    #pragma omp critical
    {
        int i, j, srcid;

        if (cfg->exportfield) {
            if (cfg->basisorder == 0 || cfg->method == rtBLBadouelGrid) {
                for (uint i = 0; i < fieldlen; i++)
                    #pragma omp atomic
                    cfg->exportfield[i] += field[i];

                /* RF forward: also store imaginary part in exportadjoint */
                if (isrfforward && cfg->exportadjoint) {
                    for (uint i = 0; i < fieldlen; i++) {
                        #pragma omp atomic
                        cfg->exportadjoint[i] += field[i + fieldlen];
                    }
                }
            } else {
                /* basisorder=1 mesh mode: redistribute per-element fluence onto nodes.
                 * Kernel layout per slot (adjoint, srcnum=1, nsrcslots>1):
                 *   field[eid + gate*ne + slot*ne*maxgate]
                 * Pattern source layout (srcnum>1, nsrcslots=1):
                 *   field[(gate*ne + eid)*srcnum + pidx]
                 * Output layout (matches grid's slot_stride = datalen * maxgate convention):
                 *   exportfield[node + gate*nn + slot*nn*maxgate]   (adjoint, srcnum=1)
                 *   exportfield[(gate*nn + node)*srcnum + pidx]      (pattern,   nsrcslots=1)
                 */
                uint nslots = (uint)((cfg->extrasrclen > 0) ? cfg->extrasrclen : 1);

                if (nslots > 1u && cfg->srcnum == 1) {
                    for (uint slot = 0; slot < nslots; slot++) {
                        size_t slot_off_f = (size_t)slot * (size_t)mesh->ne * (size_t)cfg->maxgate;
                        size_t slot_off_e = (size_t)slot * (size_t)mesh->nn * (size_t)cfg->maxgate;

                        for (i = 0; i < cfg->maxgate; i++) {
                            size_t f_gate_off = (size_t)i * (size_t)mesh->ne;
                            size_t e_gate_off = (size_t)i * (size_t)mesh->nn;

                            for (j = 0; j < mesh->ne; j++) {
                                float ww = field[slot_off_f + f_gate_off + j] * 0.25f;
                                /* RF: the imaginary part of each slot goes onto the nodes too, or the
                                   RF mesh adjoint Jacobians see a zero imaginary fluence */
                                float ww_im = (isrfforward && cfg->exportadjoint)
                                              ? field[slot_off_f + f_gate_off + j + fieldlen] * 0.25f : 0.f;
                                int k;

                                for (k = 0; k < mesh->elemlen; k++) {
                                    size_t nidx = slot_off_e + e_gate_off + mesh->elem[j * mesh->elemlen + k] - 1;
                                    cfg->exportfield[nidx] += ww;

                                    if (isrfforward && cfg->exportadjoint) {
                                        cfg->exportadjoint[nidx] += ww_im;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    for (i = 0; i < cfg->maxgate; i++) {
                        for (j = 0; j < mesh->ne; j++) {
                            for (srcid = 0; srcid < cfg->srcnum; srcid++) {
                                float ww = field[(i * mesh->ne + j) * cfg->srcnum + srcid] * 0.25f;
                                /* RF forward: the imaginary part goes onto the nodes too, or
                                   nodal RF output carries the real part only */
                                float ww_im = (isrfforward && cfg->exportadjoint && cfg->srcnum == 1)
                                              ? field[i * mesh->ne + j + fieldlen] * 0.25f : 0.f;
                                int k;

                                for (k = 0; k < mesh->elemlen; k++) {
                                    cfg->exportfield[(i * mesh->nn + mesh->elem[j * mesh->elemlen + k] - 1) * cfg->srcnum + srcid] += ww;

                                    if (ww_im != 0.f) {
                                        cfg->exportadjoint[i * mesh->nn + mesh->elem[j * mesh->elemlen + k] - 1] += ww_im;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (cfg->issaveref && mesh->dref) {
            for (uint i = 0; i < nflen; i++) {
                mesh->dref[i] += dref[i];
            }
        }
    }

    /* normalize only after all devices have added their fields and energies */
    #pragma omp barrier

    #pragma omp master
    {
        if (cfg->isnormalized) {
            double cur_normalizer, sum_normalizer = 0.0, energyabs = 0.0;

            for (j = 0; j < (uint)cfg->srcnum; j++) {
                /* with srcid<0, each photon is launched from one of the extrasrclen
                   slots at random: the outputs of a slot are normalized by the
                   energy launched (and absorbed) in it, as in mcx */
                double nslot = (cfg->extrasrclen > 0 && cfg->srcid < 0) ? cfg->extrasrclen : 1;

                energyabs =  cfg->energytot[j] - cfg->energyesc[j];
                cur_normalizer = mesh_normalize(mesh, cfg, energyabs / nslot, cfg->energytot[j] / nslot, j);
                sum_normalizer += cur_normalizer;
                MMCDEBUG(cfg, dlTime, (cfg->flog, "source %d\ttotal simulated energy: %f\tabsorbed: " S_BOLD "" S_BLUE "%5.5f%%" S_RESET "\tnormalizor=%g\n",
                                       j + 1, cfg->energytot[j], 100.f * energyabs / cfg->energytot[j], cur_normalizer));
            }

            cfg->his.normalizer = sum_normalizer / cfg->srcnum;

            /* Broadcast normalizor from slot 0 to the remaining adjoint slots
             * (srcnum..extrasrclen-1). mesh_normalize only processes pair=0..srcnum-1;
             * the multi-source/adjoint slots get scale = avg_normalizor (uniform launch
             * weight means w0/wk = 1, so no per-slot weight ratio). Includes the per-node
             * nvol division for basisorder=1 mesh mode. */
            if (cfg->extrasrclen > cfg->srcnum && cfg->exportfield) {
                double avg_normalizor = sum_normalizer / cfg->srcnum;
                int mesh_basis1 = (cfg->method != rtBLBadouelGrid) && (cfg->basisorder != 0);
                size_t datalen = (cfg->method == rtBLBadouelGrid)
                                 ? (size_t)cfg->crop0.z
                                 : (size_t)((cfg->basisorder) ? mesh->nn : mesh->ne);
                size_t slot_stride = datalen * (size_t)cfg->maxgate;

                for (int slot = cfg->srcnum; slot < cfg->extrasrclen; slot++) {
                    if (mesh_basis1) {
                        for (int t = 0; t < cfg->maxgate; t++) {
                            for (int jj = 0; jj < (int)datalen; jj++) {
                                size_t idx = (size_t)slot * slot_stride + (size_t)t * datalen + jj;

                                if (mesh->nvol[jj] > 0.f) {
                                    cfg->exportfield[idx] /= mesh->nvol[jj];
                                }

                                cfg->exportfield[idx] *= avg_normalizor;
                            }
                        }
                    } else {
                        for (size_t ki = 0; ki < slot_stride; ki++) {
                            cfg->exportfield[(size_t)slot * slot_stride + ki] *= avg_normalizor;
                        }
                    }
                }
            }

            /* RF imag fluence broadcast: same as the real-fluence block above
             * but applied to cfg->exportadjoint for slots srcnum..extrasrclen-1. */
            if (isrfforward && cfg->exportadjoint && cfg->extrasrclen > cfg->srcnum) {
                float normalizer = (float)(sum_normalizer / cfg->srcnum);
                int mesh_basis1 = (cfg->method != rtBLBadouelGrid) && (cfg->basisorder != 0);
                size_t datalen = (cfg->method == rtBLBadouelGrid)
                                 ? (size_t)cfg->crop0.z
                                 : (size_t)((cfg->basisorder) ? mesh->nn : mesh->ne);
                size_t slot_stride = datalen * (size_t)cfg->maxgate;

                for (int slot = cfg->srcnum; slot < cfg->extrasrclen; slot++) {
                    if (mesh_basis1) {
                        for (int t = 0; t < cfg->maxgate; t++) {
                            for (int jj = 0; jj < (int)datalen; jj++) {
                                size_t idx = (size_t)slot * slot_stride + (size_t)t * datalen + jj;

                                if (mesh->nvol[jj] > 0.f) {
                                    cfg->exportadjoint[idx] /= mesh->nvol[jj];
                                }

                                cfg->exportadjoint[idx] *= normalizer;
                            }
                        }
                    } else {
                        for (size_t ki = 0; ki < slot_stride; ki++) {
                            cfg->exportadjoint[(size_t)slot * slot_stride + ki] *= normalizer;
                        }
                    }
                }
            }
        }

        /* Mesh-mode adjoint Jacobian post-processing: FEM/nodal formulas on tet mesh.
         * Requires basisorder=1 (nodal fluence) so phi values are well-defined per node.
         * J_mua: either full FEM (rb_femjacobian) or nodal approximation (rbjacmuafast).
         * J_D:   always full FEM (no nodal-approx form). */
        if (cfg->issave2pt && MCX_IS_ADJOINT_TYPE(cfg->outputtype) && cfg->method != rtBLBadouelGrid &&
                cfg->basisorder && cfg->extrasrclen > 0 && cfg->detdir != NULL && cfg->exportfield) {
            unsigned int Ns = (unsigned int)(cfg->extrasrclen - cfg->detnum);
            unsigned int nch = mcx_jacchannels(cfg);   /* Ns*Nd, or the rows of cfg.sdmap */
            int isdual  = MCX_IS_DUAL_ADJOINT_TYPE(cfg->outputtype);
            int isnodal_approx = (cfg->adjointmode == 1);    /* nodal-adjoint approximation */
            /* output domain: node-based for both basisorder=1 full FEM (after elem→node scatter)
             * and the nodal approximation; element-based not currently exposed in mesh mode */
            unsigned int datalen = (unsigned int)mesh->nn;
            size_t adjointlen = (size_t)datalen * nch;
            size_t single_exportlen = adjointlen * (isrfforward ? 2 : 1);
            size_t exportlen_adj    = single_exportlen * (isdual ? 2 : 1);

            /* Convert double exportfield to float for GPU processing.
             * exportfield (= mesh->weight) is node-major for basisorder=1, sized
             * nn * extrasrclen * maxgate -- NOT fieldlen, which is sized on
             * mesh->ne (the kernel field buffer stride). Reading fieldlen
             * doubles walks past the allocation. */
            size_t nodefield_len = (size_t)datalen * (size_t)cfg->extrasrclen * (size_t)cfg->maxgate;

            float* hfield_re = (float*)malloc(sizeof(float) * nodefield_len);
            float* hfield_im = isrfforward ? (float*)malloc(sizeof(float) * nodefield_len) : NULL;

            for (size_t k = 0; k < nodefield_len; k++) {
                hfield_re[k] = (float)cfg->exportfield[k];
            }

            if (isrfforward && cfg->exportadjoint) {
                for (size_t k = 0; k < nodefield_len; k++) {
                    hfield_im[k] = cfg->exportadjoint[k];
                }
            }

            float* gfield_re_cu = NULL, *gfield_im_cu = NULL;
            float* gjmua_cu = NULL, *gjd_cu = NULL;
            CUDA_ASSERT(cudaMalloc((void**)&gfield_re_cu, sizeof(float) * nodefield_len));
            CUDA_ASSERT(cudaMemcpy(gfield_re_cu, hfield_re, sizeof(float) * nodefield_len, cudaMemcpyHostToDevice));
            free(hfield_re);

            if (hfield_im) {
                CUDA_ASSERT(cudaMalloc((void**)&gfield_im_cu, sizeof(float) * nodefield_len));
                CUDA_ASSERT(cudaMemcpy(gfield_im_cu, hfield_im, sizeof(float) * nodefield_len, cudaMemcpyHostToDevice));
                free(hfield_im);
            }

            /* Upload mesh helpers (elem, evol, deldotdel, nvol) to device */
            int* gelem_jac = NULL;
            float* gevol_jac = NULL;
            float* gdeldotdel_jac = NULL;
            float* gnvol_jac = NULL;
            int compute_jmua = (cfg->outputtype == otAdjoint || isdual);
            /* J_D is also the source of the mus/musp Jacobians (U, V), scaled per element in the kernel */
            int compute_jd   = (cfg->outputtype == otAdjointDcoeff || cfg->outputtype == otAdjointMus ||
                                cfg->outputtype == otAdjointMusp || isdual);

            /* upload elem (int4 -> int*) */
            CUDA_ASSERT(cudaMalloc((void**)&gelem_jac, sizeof(int) * mesh->ne * mesh->elemlen));
            CUDA_ASSERT(cudaMemcpy(gelem_jac, mesh->elem, sizeof(int) * mesh->ne * mesh->elemlen, cudaMemcpyHostToDevice));

            /* upload evol */
            CUDA_ASSERT(cudaMalloc((void**)&gevol_jac, sizeof(float) * mesh->ne));
            CUDA_ASSERT(cudaMemcpy(gevol_jac, mesh->evol, sizeof(float) * mesh->ne, cudaMemcpyHostToDevice));

            if (compute_jd || !isnodal_approx) {
                /* deldotdel needed for J_D (always) and for full-FEM J_mua */
                if (mesh->deldotdel == NULL) {
                    mesh_deldotdel(mesh);
                }

                float* deldotdel_f = (float*)malloc(sizeof(float) * mesh->ne * 10);

                for (size_t k = 0; k < (size_t)mesh->ne * 10; k++) {
                    deldotdel_f[k] = (float)mesh->deldotdel[k];
                }

                CUDA_ASSERT(cudaMalloc((void**)&gdeldotdel_jac, sizeof(float) * mesh->ne * 10));
                CUDA_ASSERT(cudaMemcpy(gdeldotdel_jac, deldotdel_f, sizeof(float) * mesh->ne * 10, cudaMemcpyHostToDevice));
                free(deldotdel_f);
            }

            if (isnodal_approx && compute_jmua) {
                CUDA_ASSERT(cudaMalloc((void**)&gnvol_jac, sizeof(float) * mesh->nn));
                CUDA_ASSERT(cudaMemcpy(gnvol_jac, mesh->nvol, sizeof(float) * mesh->nn, cudaMemcpyHostToDevice));
            }

            if (compute_jmua) {
                CUDA_ASSERT(cudaMalloc((void**)&gjmua_cu, sizeof(float) * single_exportlen));
                CUDA_ASSERT(cudaMemset(gjmua_cu, 0, sizeof(float) * single_exportlen));
            }

            if (compute_jd) {
                CUDA_ASSERT(cudaMalloc((void**)&gjd_cu, sizeof(float) * single_exportlen));
                CUDA_ASSERT(cudaMemset(gjd_cu, 0, sizeof(float) * single_exportlen));
            }

            size_t adjblocksize = 256;

            /* the 0-based (source, detector) pair of each Jacobian channel */
            unsigned int* hsdmap = mcx_getsdmap(cfg), *gsdmap_cu = NULL;
            CUDA_ASSERT(cudaMalloc((void**)&gsdmap_cu, sizeof(unsigned int) * 2 * (nch > 0 ? nch : 1)));
            CUDA_ASSERT(cudaMemcpy(gsdmap_cu, hsdmap, sizeof(unsigned int) * 2 * nch, cudaMemcpyHostToDevice));
            free(hsdmap);

            if (isnodal_approx && compute_jmua) {
                /* nodal approximation J_mua kernel: one thread per node */
                size_t adjgridsize = ((size_t)mesh->nn + adjblocksize - 1) / adjblocksize;
                mmc_adjoint_mesh_nodal_kernel <<< (unsigned int)adjgridsize, (unsigned int)adjblocksize>>>(
                    gfield_re_cu, gfield_im_cu, gnvol_jac, gjmua_cu,
                    (unsigned int)mesh->nn, (unsigned int)cfg->maxgate, Ns, gsdmap_cu, nch);
                CUDA_ASSERT(cudaDeviceSynchronize());
            }

            if ((!isnodal_approx && compute_jmua) || compute_jd) {
                /* full FEM kernel: one thread per element; nodal output via atomic scatter */
                size_t adjgridsize = ((size_t)mesh->ne + adjblocksize - 1) / adjblocksize;
                float* gmedscale_cu = NULL;

                /* J_mus/J_musp: the kernel multiplies each element's J_D by the factor of its label */
                if (ismusjacob) {
                    float* medscale = mesh_getmusscale(mesh, cfg->outputtype != otAdjointMus);
                    CUDA_ASSERT(cudaMalloc((void**)&gmedscale_cu, sizeof(float) * (mesh->prop + 1)));
                    CUDA_ASSERT(cudaMemcpy(gmedscale_cu, medscale, sizeof(float) * (mesh->prop + 1), cudaMemcpyHostToDevice));
                    free(medscale);
                }

                mmc_adjoint_mesh_full_kernel <<< (unsigned int)adjgridsize, (unsigned int)adjblocksize>>>(
                    gfield_re_cu, gfield_im_cu,
                    gelem_jac, gevol_jac, gdeldotdel_jac,
                    (isnodal_approx ? NULL : gjmua_cu),  /* skip J_mua here if nodal-approx already ran */
                    gjd_cu,
                    (unsigned int)mesh->ne, (unsigned int)mesh->nn, (unsigned int)cfg->maxgate,
                    Ns, gsdmap_cu, nch, (unsigned int)mesh->elemlen,
                    /*isnodal*/ 1, gtype, gmedscale_cu);
                CUDA_ASSERT(cudaDeviceSynchronize());

                if (gmedscale_cu) {
                    CUDA_ASSERT(cudaFree(gmedscale_cu));
                }
            }

            CUDA_ASSERT(cudaFree(gfield_re_cu));

            if (gfield_im_cu) {
                CUDA_ASSERT(cudaFree(gfield_im_cu));
            }

            CUDA_ASSERT(cudaFree(gelem_jac));
            CUDA_ASSERT(cudaFree(gevol_jac));
            CUDA_ASSERT(cudaFree(gsdmap_cu));

            if (gdeldotdel_jac) {
                CUDA_ASSERT(cudaFree(gdeldotdel_jac));
            }

            if (gnvol_jac) {
                CUDA_ASSERT(cudaFree(gnvol_jac));
            }

            /* Allocate output buffer and copy back */
            if (cfg->exportjacob) {
                free(cfg->exportjacob);
            }

            cfg->exportjacob = (float*)malloc(sizeof(float) * exportlen_adj);
            memset(cfg->exportjacob, 0, sizeof(float) * exportlen_adj);

            float* hmua = NULL, *hjd = NULL;

            if (compute_jmua) {
                hmua = (float*)malloc(sizeof(float) * single_exportlen);
                CUDA_ASSERT(cudaMemcpy(hmua, gjmua_cu, sizeof(float) * single_exportlen, cudaMemcpyDeviceToHost));
                CUDA_ASSERT(cudaFree(gjmua_cu));
            }

            if (compute_jd) {
                hjd = (float*)malloc(sizeof(float) * single_exportlen);
                CUDA_ASSERT(cudaMemcpy(hjd, gjd_cu, sizeof(float) * single_exportlen, cudaMemcpyDeviceToHost));
                CUDA_ASSERT(cudaFree(gjd_cu));
            }

            /* Pack into cfg->exportjacob using same layout as grid path:
             *   CW non-dual:    [Re_J1]
             *   CW dual:        [Re_J1, Re_J2]
             *   RF non-dual:    [Re_J1, Im_J1]
             *   RF dual:        [Re_J1, Re_J2, Im_J1, Im_J2] */
            if (isdual) {
                if (!isrfforward) {
                    memcpy(cfg->exportjacob,              hmua, adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen, hjd,  adjointlen * sizeof(float));
                } else {
                    memcpy(cfg->exportjacob,                  hmua,              adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen,     hjd,               adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 2 * adjointlen, hmua + adjointlen, adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 3 * adjointlen, hjd + adjointlen,  adjointlen * sizeof(float));
                }
            } else {
                float* hsrc = compute_jmua ? hmua : hjd;
                memcpy(cfg->exportjacob, hsrc, single_exportlen * sizeof(float));
            }

            if (hmua) {
                free(hmua);
            }

            if (hjd) {
                free(hjd);
            }

            MMC_FPRINTF(cfg->flog, "mesh adjoint Jacobian computation complete (%s): %d ms\n",
                        isnodal_approx ? "nodal approx" : "full FEM", GetTimeMillis() - tic);

#ifndef MCX_CONTAINER

            if (cfg->issave2pt && cfg->parentid == mpStandalone && cfg->exportjacob) {
                MMC_FPRINTF(cfg->flog, "saving mesh adjoint Jacobian to file ...\t");
                mesh_savejacob(cfg, mesh, cfg->exportjacob, (int)nch, isrfforward, isdual);
                MMC_FPRINTF(cfg->flog, "saving Jacobian complete : %d ms\n\n", GetTimeMillis() - tic);
                mcx_fflush(cfg->flog);
            }

#endif
        }

        /* Adjoint Jacobian post-processing (grid mode only) */
        mmc_cu_adjoint_grid(cfg, mesh, cfg->exportfield, fieldlen, isrfforward, gvoxtype, tic);

#ifndef MCX_CONTAINER

        if (cfg->issave2pt && cfg->parentid == mpStandalone) {
            MMC_FPRINTF(cfg->flog, "saving data to file ...\t");
            mesh_saveweight(mesh, cfg, 0);
            MMC_FPRINTF(cfg->flog, "saving data complete : %d ms\n\n",
                        GetTimeMillis() - tic);
            mcx_fflush(cfg->flog);
        }

        if (cfg->issavedet && cfg->parentid == mpStandalone &&
                cfg->exportdetected) {
            MMC_FPRINTF(cfg->flog, "saving detected photon data to file ...\t");
            cfg->his.totalphoton = cfg->nphoton;
            cfg->his.unitinmm = cfg->unitinmm;
            cfg->his.savedphoton = cfg->detectedcount;
            cfg->his.detected = cfg->detectedcount;
            cfg->his.colcount = (2 + (cfg->ismomentum > 0)) * cfg->his.maxmedia + (cfg->issaveexit > 0) * 6 + 2; /*column count=maxmedia+3*/
            cfg->his.seedbyte = (cfg->exportseed) ? (sizeof(RandType) * RAND_BUF_LEN) : 0;
            mcx_savedetphoton(cfg->exportdetected, (void*)(cfg->exportseed), cfg->detectedcount, 0, cfg);
            MMC_FPRINTF(cfg->flog, "saving detected photon data complete : %d ms\n\n", GetTimeMillis() - tic);
        }

        /**
         * If not running as a mex file, we need to save photon trajectory data, if enabled, as
         * a file, either as a .mct file, or a .jdat/.jbat file
         */

        if ((cfg->debuglevel & dlTraj) && cfg->parentid == mpStandalone && cfg->exportdebugdata) {
            MMC_FPRINTF(cfg->flog, "saving trajectory data to file ...\t");
            cfg->his.colcount = MCX_DEBUG_REC_LEN;
            cfg->his.savedphoton = cfg->debugdatalen;
            cfg->his.totalphoton = cfg->nphoton;
            cfg->his.detected = 0;
            mcx_savedetphoton(cfg->exportdebugdata, NULL, cfg->debugdatalen, 0, cfg);
            MMC_FPRINTF(cfg->flog, "saving trajectory data complete : %d ms\n\n", GetTimeMillis() - tic);
        }

        if (cfg->issaveref) {
            MMC_FPRINTF(cfg->flog, "saving surface diffuse reflectance ...");
            mesh_saveweight(mesh, cfg, 1);
        }

#endif

        // total energy here equals total simulated photons+unfinished photons for
        // all threads
        energytot = energyesc = 0.0;

        for (j = 0; j < (uint)cfg->srcnum; j++) {
            energytot += cfg->energytot[j];
            energyesc += cfg->energyesc[j];
        }

        MMC_FPRINTF(cfg->flog,
                    "simulated %ld photons (%ld) with devices (ray-tet "
                    "%.0f)\nMCX simulation speed: %.2f photon/ms\n",
                    cfg->nphoton, cfg->nphoton, reporter.raytet,
                    (double)cfg->nphoton / toc);
        MMC_FPRINTF(cfg->flog,
                    "total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to "
                    "initial specular reflection is excluded in the total)\n",
                    energytot,
                    (energytot - energyesc) / energytot * 100.f);
        mcx_fflush(cfg->flog);
    }
    #pragma omp barrier
    CUDA_ASSERT(cudaFree(gnode));
    CUDA_ASSERT(cudaFree(gelem));
    CUDA_ASSERT(cudaFree(gtype));
    CUDA_ASSERT(cudaFree(gfacenb));
    CUDA_ASSERT(cudaFree(gsrcelem));
    CUDA_ASSERT(cudaFree(gnormal));

    if (gnodemua_cu) {
        CUDA_ASSERT(cudaFree(gnodemua_cu));
    }

    if (gnodemusp_cu) {
        CUDA_ASSERT(cudaFree(gnodemusp_cu));
    }

    CUDA_ASSERT(cudaFree(gseed));
    CUDA_ASSERT(cudaFree(gdetphoton));
    CUDA_ASSERT(cudaFree(gweight));
    CUDA_ASSERT(cudaFree(gdref));
    CUDA_ASSERT(cudaFree(genergy));
    CUDA_ASSERT(cudaFree(gdetected));

    if (gsrcpattern) {
        CUDA_ASSERT(cudaFree(gsrcpattern));
    }

    if (greplayweight) {
        CUDA_ASSERT(cudaFree(greplayweight));
    }

    if (greplayseed) {
        CUDA_ASSERT(cudaFree(greplayseed));
    }

    if (greplaytime) {
        CUDA_ASSERT(cudaFree(greplaytime));
    }

    if (gphotonseed) {
        CUDA_ASSERT(cudaFree(gphotonseed));
    }

    if (cfg->debuglevel & dlTraj) {
        CUDA_ASSERT(cudaFree(gdebugdata));
    }

    if (gvoxtype) {
        CUDA_ASSERT(cudaFree(gvoxtype));
    }

    if (groioff) {
        CUDA_ASSERT(cudaFree(groioff));
        CUDA_ASSERT(cudaFree(groilist));
    }

    if (gfaceroi) {
        CUDA_ASSERT(cudaFree(gfaceroi));
    }

    CUDA_ASSERT(cudaFree(greporter));

    #pragma omp master
    {
        if (gpu) {
            free(gpu);
        }

        free(cfg->energytot);
        free(cfg->energyesc);
        cfg->energytot = NULL;
        cfg->energyesc = NULL;
    }

    free(field);

    if (Pdet) {
        free(Pdet);
    }

    if (Pphotonseed) {
        free(Pphotonseed);
    }

    free(dref);
}

void mmc_run_cu(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {
    GPUInfo* gpuinfo = NULL;      /** gpuinfo: structure to store GPU information */
    unsigned int activedev = 0;   /** activedev: count of total active GPUs to be used */

    if (!(activedev = mcx_list_cu_gpu(cfg, &gpuinfo))) {
        mcx_error(-1, "No GPU device found\n", __FILE__, __LINE__);
    }

    mmc_error_reset();

#ifdef _OPENMP
    /**
        Now we are ready to launch one thread for each involked GPU to run the simulation
     */
    omp_set_num_threads(activedev);
    #pragma omp parallel
    {
        /* an exception can not leave the parallel region: the error of a GPU thread is recorded
           (see mmc_utils.h) and raised again once all threads have finished */
        try {
            mmc_run_simulation(cfg, mesh, tracer, gpuinfo);
        } catch (...) {
            mmc_defer_error(-1, "unknown exception in a GPU thread", __FILE__, __LINE__);
        }
    }

    mmc_raise_deferred();
#else
    /**
        This line runs the main MCX simulation for each GPU inside each thread
     */
    mmc_run_simulation(cfg, mesh, tracer, gpuinfo);
#endif
}

/**
 * @brief Grid-mode adjoint Jacobian post-processing, shared by the CUDA and OptiX hosts
 *
 * Computes J_mua (-O A), J_D/J_mus/J_musp (-O D/U/V) or both (-O W/Q) from the
 * per-slot fields of the sources and the detectors, using the CUDA kernels of mmc_core.cu,
 * and stores them in cfg->exportjacob (saved to file in the standalone binary).
 *
 * \param[in] field: the normalized fields of all source/detector slots (datalen*maxgate per slot)
 * \param[in] fieldlen: the length of field
 * \param[in] isrfforward: 1 if cfg->exportadjoint holds the RF imaginary fields
 * \param[in] gvoxtype: device array of the medium label of each voxel (J_mus/J_musp), or NULL
 * \param[in] tic: the start time for the timing printouts
 */
void mmc_cu_adjoint_grid(mcconfig* cfg, tetmesh* mesh, double* field, size_t fieldlen, int isrfforward, int* gvoxtype, unsigned int tic) {
    if (!(cfg->issave2pt && MCX_IS_ADJOINT_TYPE(cfg->outputtype) && cfg->method == rtBLBadouelGrid &&
            cfg->extrasrclen > 0 && cfg->detdir != NULL && field)) {
        return;
    }

    unsigned int Ns = (unsigned int)(cfg->extrasrclen - cfg->detnum);
    unsigned int nch = mcx_jacchannels(cfg);   /* Ns*Nd, or the rows of cfg.sdmap */
    unsigned int pure_voxels = (unsigned int)cfg->crop0.z;
    int isdual  = MCX_IS_DUAL_ADJOINT_TYPE(cfg->outputtype);

    size_t adjointlen = (size_t)pure_voxels * nch;
    size_t single_exportlen = adjointlen * (isrfforward ? 2 : 1);
    size_t exportlen_adj    = single_exportlen * (isdual ? 2 : 1);

    /* Convert double exportfield to float for GPU processing */
    float* hfield_re = (float*)malloc(sizeof(float) * fieldlen);
    float* hfield_im = isrfforward ? (float*)malloc(sizeof(float) * fieldlen) : NULL;

    for (size_t k = 0; k < fieldlen; k++) {
        hfield_re[k] = (float)field[k];
    }

    if (isrfforward && cfg->exportadjoint) {
        for (size_t k = 0; k < fieldlen; k++) {
            hfield_im[k] = cfg->exportadjoint[k];
        }
    }

    float* gfield_re_cu = NULL, *gfield_im_cu = NULL;
    float* gadjoint_mua_cu = NULL, *gadjoint_tmp_cu = NULL;
    CUDA_ASSERT(cudaMalloc((void**)&gfield_re_cu, sizeof(float) * fieldlen));
    CUDA_ASSERT(cudaMemcpy(gfield_re_cu, hfield_re, sizeof(float) * fieldlen, cudaMemcpyHostToDevice));
    free(hfield_re);

    if (hfield_im) {
        CUDA_ASSERT(cudaMalloc((void**)&gfield_im_cu, sizeof(float) * fieldlen));
        CUDA_ASSERT(cudaMemcpy(gfield_im_cu, hfield_im, sizeof(float) * fieldlen, cudaMemcpyHostToDevice));
        free(hfield_im);
    }

    if (isdual) {
        CUDA_ASSERT(cudaMalloc((void**)&gadjoint_mua_cu, sizeof(float) * single_exportlen));
        CUDA_ASSERT(cudaMemset(gadjoint_mua_cu, 0, sizeof(float) * single_exportlen));
    }

    CUDA_ASSERT(cudaMalloc((void**)&gadjoint_tmp_cu, sizeof(float) * single_exportlen));
    CUDA_ASSERT(cudaMemset(gadjoint_tmp_cu, 0, sizeof(float) * single_exportlen));

    size_t adjblocksize = 256;
    size_t adjgridsize  = (pure_voxels + (unsigned int)adjblocksize - 1) / adjblocksize;

    /* the 0-based (source, detector) pair of each Jacobian channel */
    unsigned int* hsdmap = mcx_getsdmap(cfg), *gsdmap_cu = NULL;
    CUDA_ASSERT(cudaMalloc((void**)&gsdmap_cu, sizeof(unsigned int) * 2 * (nch > 0 ? nch : 1)));
    CUDA_ASSERT(cudaMemcpy(gsdmap_cu, hsdmap, sizeof(unsigned int) * 2 * nch, cudaMemcpyHostToDevice));
    free(hsdmap);

    if (isdual || cfg->outputtype == otAdjoint) {
        mmc_adjoint_kernel <<< (unsigned int)adjgridsize, (unsigned int)adjblocksize>>>(
            gfield_re_cu, gfield_im_cu,
            isdual ? gadjoint_mua_cu : gadjoint_tmp_cu,
            pure_voxels, (unsigned int)cfg->maxgate, Ns, gsdmap_cu, nch);
        CUDA_ASSERT(cudaDeviceSynchronize());
    }

    if (isdual || cfg->outputtype != otAdjoint) {
        float* gmedscale_cu = NULL;

        /* J_mus/J_musp: the kernel multiplies each voxel's J_D by the factor of the medium
           label recorded there during the simulation */
        if (gvoxtype) {
            float* medscale = mesh_getmusscale(mesh, cfg->outputtype != otAdjointMus);
            CUDA_ASSERT(cudaMalloc((void**)&gmedscale_cu, sizeof(float) * (mesh->prop + 1)));
            CUDA_ASSERT(cudaMemcpy(gmedscale_cu, medscale, sizeof(float) * (mesh->prop + 1), cudaMemcpyHostToDevice));
            free(medscale);
        }

        mmc_adjoint_dcoeff_kernel <<< (unsigned int)adjgridsize, (unsigned int)adjblocksize>>>(
            gfield_re_cu, gfield_im_cu, gadjoint_tmp_cu,
            pure_voxels, (unsigned int)cfg->maxgate, Ns, gsdmap_cu, nch,
            (unsigned int)cfg->dim.x, (unsigned int)cfg->dim.y,
            gvoxtype, gmedscale_cu);
        CUDA_ASSERT(cudaDeviceSynchronize());

        if (gmedscale_cu) {
            CUDA_ASSERT(cudaFree(gmedscale_cu));
        }
    }

    CUDA_ASSERT(cudaFree(gfield_re_cu));
    CUDA_ASSERT(cudaFree(gsdmap_cu));

    if (gfield_im_cu) {
        CUDA_ASSERT(cudaFree(gfield_im_cu));
    }

    /* Allocate separate Jacobian buffer; exportadjoint keeps RF imaginary fluence */
    if (cfg->exportjacob) {
        free(cfg->exportjacob);
    }

    cfg->exportjacob = (float*)malloc(sizeof(float) * exportlen_adj);

    float Vvox = cfg->unitinmm * cfg->unitinmm * cfg->unitinmm;

    if (isdual) {
        float* hmua    = (float*)malloc(sizeof(float) * single_exportlen);
        float* hsecond = (float*)malloc(sizeof(float) * single_exportlen);
        CUDA_ASSERT(cudaMemcpy(hmua,    gadjoint_mua_cu, sizeof(float) * single_exportlen, cudaMemcpyDeviceToHost));
        CUDA_ASSERT(cudaMemcpy(hsecond, gadjoint_tmp_cu, sizeof(float) * single_exportlen, cudaMemcpyDeviceToHost));
        CUDA_ASSERT(cudaFree(gadjoint_mua_cu));
        CUDA_ASSERT(cudaFree(gadjoint_tmp_cu));

        for (size_t k = 0; k < single_exportlen; k++) {
            hmua[k]    *= -Vvox;
            hsecond[k] *= -cfg->unitinmm;
        }

        if (!isrfforward) {
            memcpy(cfg->exportjacob,              hmua,    adjointlen * sizeof(float));
            memcpy(cfg->exportjacob + adjointlen, hsecond, adjointlen * sizeof(float));
        } else {
            memcpy(cfg->exportjacob,                   hmua,                 adjointlen * sizeof(float));
            memcpy(cfg->exportjacob + adjointlen,      hsecond,              adjointlen * sizeof(float));
            memcpy(cfg->exportjacob + 2 * adjointlen,  hmua + adjointlen,    adjointlen * sizeof(float));
            memcpy(cfg->exportjacob + 3 * adjointlen,  hsecond + adjointlen, adjointlen * sizeof(float));
        }

        free(hmua);
        free(hsecond);
    } else {
        CUDA_ASSERT(cudaMemcpy(cfg->exportjacob, gadjoint_tmp_cu, sizeof(float) * single_exportlen, cudaMemcpyDeviceToHost));
        CUDA_ASSERT(cudaFree(gadjoint_tmp_cu));

        float adj_scale = (cfg->outputtype == otAdjoint) ? -Vvox : -cfg->unitinmm;

        for (size_t k = 0; k < single_exportlen; k++) {
            cfg->exportjacob[k] *= adj_scale;
        }
    }

    MMC_FPRINTF(cfg->flog, "adjoint Jacobian computation complete: %d ms\n", GetTimeMillis() - tic);

#ifndef MCX_CONTAINER

    if (cfg->issave2pt && cfg->parentid == mpStandalone && cfg->exportjacob) {
        MMC_FPRINTF(cfg->flog, "saving adjoint Jacobian to file ...\t");
        mesh_savejacob(cfg, mesh, cfg->exportjacob, (int)nch, isrfforward, isdual);
        MMC_FPRINTF(cfg->flog, "saving Jacobian complete : %d ms\n\n", GetTimeMillis() - tic);
        mcx_fflush(cfg->flog);
    }

#endif
}
