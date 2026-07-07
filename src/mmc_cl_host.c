/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2025
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing
**          in Plucker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli,
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a>
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**  \li \c (\b Fang2019) Qianqian Fang and Shijie Yan,
**          <a href="http://dx.doi.org/10.1117/1.JBO.24.11.115002">
**          "Graphics processing unit-accelerated mesh-based Monte Carlo photon transport
**           simulations,"</a> J. of Biomedical Optics, 24(11), 115002 (2019)
**  \li \c (\b Yuan2021) Yaoshen Yuan, Shijie Yan, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-12-1-147">
**          "Light transport modeling in highly complex tissues using the implicit
**           mesh-based Monte Carlo algorithm,"</a> Biomed. Optics Express, 12(1) 147-161 (2021)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmc_cl_host.c

\brief   OpenCL host code for OpenCL based MMC simulations
*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mmc_cl_host.h"
#include "mmc_const.h"
#include "mmc_tictoc.h"

#define MAX_JIT_OPT_LEN  (MAX_PATH_LENGTH << 1)

#define IPARAM_TO_MACRO(macro,a,b)  snprintf(macro+strlen(macro),MAX_JIT_OPT_LEN-strlen(macro)-1," -Dgcfg%s=%u ",  #b,(a.b))
#define SIPARAM_TO_MACRO(macro,a,b) snprintf(macro+strlen(macro),MAX_JIT_OPT_LEN-strlen(macro)-1," -Dgcfg%s=%d ",  #b,(a.b))
#define FPARAM_TO_MACRO(macro,a,b)  snprintf(macro+strlen(macro),MAX_JIT_OPT_LEN-strlen(macro)-1," -Dgcfg%s=%.10ef ",#b,(a.b))

/* Acquisition-phase error handling. On failure these set 'status' and break out of
 * the do{...}while(0) guard that wraps the run, so the single NULL-safe teardown that
 * follows runs for both the success and the error path; the deferred error is then
 * re-reported at the very end. This matters when MMC runs as a library (mmclab/pmmc),
 * where mcx_error throws a C++ exception instead of calling exit(): a clCreateBuffer/
 * clBuildProgram/clCreateKernel failure (e.g. GPU out-of-memory) would otherwise leak
 * the context and all partially-allocated device buffers.
 *
 * NOTE: these expand to a bare 'if(...) break;' (break cannot live inside a do/while
 * wrapper here), so only use them as standalone statements or inside braced blocks,
 * never as the unbraced body of an if/else. 'break' inside one of the acquisition
 * for-loops exits that loop; an 'if (status != CL_SUCCESS) break;' after the loop then
 * leaves the guard. */
#define OCL_TRY(x)           if ((status = (x)) != CL_SUCCESS) { break; }
#define MMC_ASSERT_ALLOC(p)  if ((p) == NULL) { status = CL_OUT_OF_HOST_MEMORY; break; }

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

const char* sourceflag[] = {"-DMCX_SRC_PENCIL", "-DMCX_SRC_ISOTROPIC", "-DMCX_SRC_CONE",
                            "-DMCX_SRC_GAUSSIAN", "-DMCX_SRC_PLANAR", "-DMCX_SRC_PATTERN", "-DMCX_SRC_FOURIER",
                            "-DMCX_SRC_ARCSINE", "-DMCX_SRC_DISK", "-DMCX_SRC_FOURIERX", "-DMCX_SRC_FOURIERX2D",
                            "-DMCX_SRC_ZGAUSSIAN", "-DMCX_SRC_LINE", "-DMCX_SRC_SLIT", "-DMCX_SRC_PENCILARRAY",
                            "-DMCX_SRC_PATTERN3D", "-DMCX_SRC_HYPERBOLOID_GAUSSIAN", "-DMCX_SRC_RING"
                           };

extern cl_event kernelevent;

/*
   master driver code to run MC simulations
*/

void mmc_run_cl(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {

    cl_uint i, j, iter;
    cl_float t, twindow0, twindow1;
    cl_float fullload = 0.f;
    cl_float* energy = NULL;
    cl_uint* progress = NULL;
    cl_uint detected = 0, workdev = 0;

    cl_uint tic, tic0, tic1, toc = 0, debuglen = MCX_DEBUG_REC_LEN;
    size_t fieldlen = 0;                   /* size_t to avoid 32-bit overflow on large grid*gate*srcslot products */

    cl_context mcxcontext = NULL;          // compute mcxcontext
    cl_command_queue* mcxqueue = NULL;   // compute command queue
    cl_program mcxprogram = NULL;          // compute mcxprogram
    cl_kernel* mcxkernel = NULL;            // compute mcxkernel
    cl_int status = 0;
    cl_device_id devices[MAX_DEVICE];
    cl_platform_id platform = NULL;

    cl_uint  totalcucore;
    cl_uint  devid = 0;
    cl_mem* gnode = NULL, *gelem = NULL, *gtype = NULL, *gfacenb = NULL, *gsrcelem = NULL, *gnormal = NULL;
    cl_mem* gnodemua_cl = NULL, *gnodemusp_cl = NULL;   /* per-node mua/musp for DOT recon (NULL otherwise) */
    cl_mem* gproperty = NULL, *gparam = NULL, *gsrcpattern = NULL, *greplayweight = NULL, *greplaytime = NULL, *greplayseed = NULL; /*read-only buffers*/
    cl_mem* gweight = NULL, *gdref = NULL, *gdetphoton = NULL, *gseed = NULL, *genergy = NULL, *greporter = NULL, *gdebugdata = NULL;     /*read-write buffers*/
    cl_mem* gprogress = NULL, *gdetected = NULL, *gphotonseed = NULL; /*write-only buffers*/

    int isrfforward = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE);
    int savevolume = (cfg->issave2pt != 0);
    /* cfg->srcid > 0 selects a single slot from srcdata, collapsing the field buffer to one slot. */
    cl_uint nsrcslots = (cfg->extrasrclen > 0 && cfg->srcid <= 0) ? (cl_uint)cfg->extrasrclen : 1u;
    /* Compute the output length in 64-bit first and bail out cleanly if it cannot
     * be addressed by the kernel's 32-bit indices, rather than silently wrapping
     * meshlen/crop0.w and under-allocating the field buffer. */
    size_t meshlen64 = (size_t)((cfg->method == rtBLBadouelGrid) ? cfg->crop0.z : mesh->ne) * (size_t)cfg->srcnum * (size_t)nsrcslots;

    if (savevolume && meshlen64 * (size_t)((cfg->maxgate > 0) ? cfg->maxgate : 1) > 0xFFFFFFFFULL) {
        mcx_error(-6, (char*)("output buffer length exceeds 32-bit addressing; reduce grid size, time gates, or source count"), __FILE__, __LINE__);
    }

    /* MCX_SKIP_VOLUME keeps the kernel from touching gweight, but the kernel
     * signature still needs a valid buffer argument. Use one dummy slot when
     * volume output is disabled. */
    cl_uint meshlen = savevolume ? (cl_uint)meshlen64 : 1u;    /**< total output data length (per gate), or a dummy slot when volume saving is off */
    cfg->crop0.w = savevolume ? meshlen * cfg->maxgate : 1u;    /**< total output data length, before double-buffer expansion */

    cl_float*  field = NULL, *dref = NULL;

    cl_uint*   Pseed = NULL;
    float*     Pdet = NULL;
    RandType*  Pphotonseed = NULL;
    char opt[MAX_JIT_OPT_LEN] = {'\0'};
    cl_uint detreclen = (2 + ((cfg->ismomentum) > 0)) * mesh->prop + (cfg->issaveexit > 0) * 6 + 1;
    cl_uint hostdetreclen = detreclen + 1;
    int sharedmemsize = 0;

    GPUInfo* gpu = NULL;
    float4* propdet = NULL;
    double energytot = 0.0, energyesc = 0.0;

    MCXParam param = {{{cfg->srcparam1.x, cfg->srcparam1.y, cfg->srcparam1.z, cfg->srcparam1.w}},
        {{cfg->srcparam2.x, cfg->srcparam2.y, cfg->srcparam2.z, cfg->srcparam2.w}},
        {{cfg->crop0.x, cfg->crop0.y, cfg->crop0.z, cfg->crop0.w}},
        {{cfg->bary0.x, cfg->bary0.y, cfg->bary0.z, cfg->bary0.w}},
        {{cfg->srcpos.x, cfg->srcpos.y, cfg->srcpos.z}},
        {{cfg->srcdir.x, cfg->srcdir.y, cfg->srcdir.z}},
        {{mesh->nmin.x, mesh->nmin.y, mesh->nmin.z}},
        cfg->tstart, cfg->tend, (uint)cfg->isreflect, (uint)cfg->issavedet, (uint)cfg->issaveexit,
        (uint)cfg->ismomentum, (uint)cfg->isatomic, (uint)cfg->isspecular, 1.f / cfg->tstep, cfg->minenergy,
        cfg->maxdetphoton, mesh->prop, cfg->detnum, (uint)cfg->voidtime, (uint)cfg->srctype,
        cfg->issaveref, cfg->maxgate, (uint)cfg->debuglevel, detreclen, cfg->outputtype, mesh->elemlen,
        cfg->mcmethod, cfg->method, 1.f / cfg->steps.x,
#if defined(MMC_USE_SSE) || defined(USE_OPENCL)
        cfg->srcdir.w,
#else
        0.f,
#endif
        mesh->nn, mesh->ne, mesh->nf, cfg->nout, cfg->roulettesize, cfg->srcnum, mesh->srcelemlen,
        cfg->e0, cfg->isextdet, (meshlen / (cfg->srcnum * nsrcslots)),  /* framelen = datalen (per-gate stride); slot stride is datalen*maxgate */
        (mesh->prop + 1 + cfg->isextdet) + (cfg->extrasrclen * 4) + cfg->detnum,
        (MIN((MAX_PROP - param.maxpropdet), ((mesh->ne) << 2)) >> 2), /*max count of elem normal data in const mem*/
        cfg->issaveseed, cfg->seed, cfg->maxjumpdebug,
        cfg->omega,
        (float)(3.335640951981520e-12),   /* oneoverc0 */
        (cfg->extrasrclen > 0 && cfg->srcid == 0) ? -1 : cfg->srcid,  /* srcid<0 multi-slot, srcid>0 single slot, srcid==0 with srcdata: auto-promote to -1 */
        cfg->extrasrclen,
        (int)(mesh->prop + 1 + cfg->isextdet),  /* srcpropoffset */
        (uint)cfg->isnodalmua,
        (uint)cfg->isnodalmusp
    };

    MCXReporter reporter = {0.f, 0};
    platform = mcx_list_cl_gpu(cfg, &workdev, devices, &gpu);

    if (workdev > MAX_DEVICE) {
        workdev = MAX_DEVICE;
    }

    if (workdev == 0) {
        mcx_error(-99, (char*)("Unable to find devices!"), __FILE__, __LINE__);
    }

    if (cfg->issavedet) {
        sharedmemsize = sizeof(cl_float) * detreclen;
    }

    param.reclen = sharedmemsize / sizeof(float);

    sharedmemsize += sizeof(cl_float) * (cfg->srcnum << 1);   /**< store energyesc/energytot */

    if (cfg->srctype == stPattern && cfg->srcnum > 1) {
        sharedmemsize += sizeof(cl_float) * cfg->srcnum;
    }

    cl_context_properties cps[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};

    /* Use NULL for backward compatibility */
    cl_context_properties* cprops = (platform == NULL) ? NULL : cps;

    /* Guard the whole resource-acquisition + run region: any acquisition failure
     * (OCL_TRY/MMC_ASSERT_ALLOC) breaks out to the single NULL-safe teardown below
     * instead of throwing past it. See the OCL_TRY comment near the top of the file. */
    do {
        OCL_TRY(((mcxcontext = clCreateContext(cprops, workdev, devices, NULL, NULL, &status), status)));

        /* calloc (not malloc) so every per-device handle starts NULL: the cleanup
         * path can then walk these arrays and release only what was created, even if
         * allocation aborted partway through. */
        mcxqueue = (cl_command_queue*)calloc(workdev, sizeof(cl_command_queue));

        gseed = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gweight = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gdref = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gdetphoton = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        genergy = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gdetected = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gphotonseed = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        greporter = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gdebugdata = (cl_mem*)calloc(workdev, sizeof(cl_mem));

        gnode = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gelem = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gtype = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gfacenb = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gsrcelem = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gnormal = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gnodemua_cl  = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gnodemusp_cl = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gproperty = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gparam = (cl_mem*)calloc(workdev, sizeof(cl_mem));

        gprogress = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        gsrcpattern = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        greplayweight = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        greplaytime = (cl_mem*)calloc(workdev, sizeof(cl_mem));
        greplayseed = (cl_mem*)calloc(workdev, sizeof(cl_mem));

        if (!mcxqueue || !gseed || !gweight || !gdref || !gdetphoton || !genergy || !gdetected ||
                !gphotonseed || !greporter || !gdebugdata || !gnode || !gelem || !gtype || !gfacenb ||
                !gsrcelem || !gnormal || !gnodemua_cl || !gnodemusp_cl || !gproperty || !gparam ||
                !gprogress || !gsrcpattern || !greplayweight || !greplaytime || !greplayseed) {
            status = CL_OUT_OF_HOST_MEMORY;
            break;
        }

        /* The block is to move the declaration of prop closer to its use */
        cl_command_queue_properties prop = CL_QUEUE_PROFILING_ENABLE;

        totalcucore = 0;

        for (i = 0; i < workdev; i++) {
            OCL_TRY(((mcxqueue[i] = clCreateCommandQueue(mcxcontext, devices[i], prop, &status), status)));
            totalcucore += gpu[i].core;

            if (!cfg->autopilot) {
                gpu[i].autothread = cfg->nthread;
                gpu[i].autoblock = cfg->nblocksize;
                gpu[i].maxgate = cfg->maxgate;
            } else {
                // persistent thread mode
                if (gpu[i].vendor == dvIntelGPU || gpu[i].vendor == dvIntel) { // Intel HD graphics GPU
                    gpu[i].autoblock  = (gpu[i].vendor == dvIntelGPU) ? 64 : 1;
                    gpu[i].autothread = gpu[i].autoblock * 7 * gpu[i].sm; // 7 thread x SIMD-16 per Exec Unit (EU)
                } else if (gpu[i].vendor == dvAMD) { // AMD GPU
                    gpu[i].autoblock  = 64;
                    gpu[i].autothread = 2560 * gpu[i].sm; // 40 wavefronts * 64 threads/wavefront
                } else if (gpu[i].vendor == dvNVIDIA) {
                    if (gpu[i].major == 2 || gpu[i].major == 3) { // fermi 2.x, kepler 3.x : max 7 blks per SM, 8 works better
                        gpu[i].autoblock  = 128;
                        gpu[i].autothread = gpu[i].autoblock * 8 * gpu[i].sm;
                    } else if (gpu[i].major == 5) { // maxwell 5.x
                        gpu[i].autoblock  = 64;
                        gpu[i].autothread = gpu[i].autoblock * 16 * gpu[i].sm;
                    } else if (gpu[i].major >= 6) { // pascal 6.x : max 32 blks per SM
                        gpu[i].autoblock  = 64;
                        gpu[i].autothread = gpu[i].autoblock * 64 * gpu[i].sm;
                    }
                }
            }

            if (gpu[i].autothread % gpu[i].autoblock) {
                gpu[i].autothread = (gpu[i].autothread / gpu[i].autoblock) * gpu[i].autoblock;
            }

            if (gpu[i].maxgate == 0 && meshlen > 0) {
                int needmem = meshlen + gpu[i].autothread * sizeof(float4) * 4 + sizeof(float) * cfg->maxdetphoton * hostdetreclen + 10 * 1024 * 1024; /*keep 10M for other things*/
                gpu[i].maxgate = (gpu[i].globalmem - needmem) / meshlen;
                gpu[i].maxgate = MIN(((cfg->tend - cfg->tstart) / cfg->tstep + 0.5), gpu[i].maxgate);
            }
        }

        if (status != CL_SUCCESS) {
            break;    /* a clCreateCommandQueue above failed; leave the guard */
        }

        cfg->maxgate = (int)((cfg->tend - cfg->tstart) / cfg->tstep + 0.5);
        param.maxgate = cfg->maxgate;
        cl_uint nflen = mesh->nf * cfg->maxgate;

        fullload = 0.f;

        for (i = 0; i < workdev; i++) {
            fullload += cfg->workload[i];
        }

        if (fullload < EPS) {
            for (i = 0; i < workdev; i++) {
                cfg->workload[i] = gpu[i].core;
            }

            fullload = totalcucore;
        }

        field = savevolume ? (cl_float*)calloc(sizeof(cl_float) * (size_t)meshlen * (isrfforward ? 4 : 2), cfg->maxgate)
                : (cl_float*)calloc((size_t)(isrfforward ? 4 : 2), sizeof(cl_float));
        dref = (cl_float*)calloc(sizeof(cl_float) * mesh->nf, cfg->maxgate);
        Pdet = (float*)calloc(cfg->maxdetphoton * sizeof(float), hostdetreclen);
        MMC_ASSERT_ALLOC(field);
        MMC_ASSERT_ALLOC(dref);
        MMC_ASSERT_ALLOC(Pdet);

        if (cfg->issaveseed) {
            Pphotonseed = (RandType*)calloc(cfg->maxdetphoton, (sizeof(RandType) * RAND_BUF_LEN));
            MMC_ASSERT_ALLOC(Pphotonseed);
        }

        fieldlen = cfg->crop0.w;  /**< total float counts of the output buffer, before double-buffer expansion (x2) for improving saving accuracy */

        if (cfg->seed > 0) {
            srand(cfg->seed);
        } else {
            srand(time(0));
        }

        if (cfg->debuglevel & dlTraj && cfg->exportdebugdata == NULL) {
            cfg->exportdebugdata = (float*)calloc(sizeof(float), (debuglen * cfg->maxjumpdebug));
        }

        cl_mem (*clCreateBufferNV)(cl_context, cl_mem_flags, cl_mem_flags_NV, size_t, void*, cl_int*) = (cl_mem (*)(cl_context, cl_mem_flags, cl_mem_flags_NV, size_t, void*, cl_int*)) clGetExtensionFunctionAddressForPlatform(platform, "clCreateBufferNV");

        if (clCreateBufferNV == NULL) {
            OCL_TRY(((gprogress[0] = clCreateBuffer(mcxcontext, RW_PTR, sizeof(cl_uint), NULL, &status), status)));
        } else {
            gprogress[0] = clCreateBufferNV(mcxcontext, CL_MEM_READ_WRITE, NV_PIN, sizeof(cl_uint), NULL, &status);

            if (status != CL_SUCCESS) {
                OCL_TRY(((gprogress[0] = clCreateBuffer(mcxcontext, RW_PTR, sizeof(cl_uint), NULL, &status), status)));
            }
        }

        progress = (cl_uint*)clEnqueueMapBuffer(mcxqueue[0], gprogress[0], CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(cl_uint), 0, NULL, NULL, &status);
        OCL_TRY(status);
        *progress = 0;

        propdet = (float4*)malloc(MAX_PROP * sizeof(float4));
        MMC_ASSERT_ALLOC(propdet);

        /* Layout: [media 0..prop] | [extra sources, each 4 float4 slots] | [detector positions] | [normals] */
        cl_uint srcpropoffset = (cl_uint)(mesh->prop + 1 + cfg->isextdet);
        cl_uint detpropoffset = srcpropoffset + (cl_uint)(cfg->extrasrclen * 4);
        memcpy(propdet, mesh->med, srcpropoffset * sizeof(medium));

        /* Pack extra sources into propdet after media; sizeof(ExtraSrc) = 4*sizeof(float4) */
        if (cfg->extrasrclen > 0 && cfg->srcdata) {
            memcpy(propdet + srcpropoffset, cfg->srcdata, cfg->extrasrclen * sizeof(ExtraSrc));
        }

        if (cfg->detpos && cfg->detnum) {
            memcpy(propdet + detpropoffset, cfg->detpos, cfg->detnum * sizeof(float4));
        }

        memcpy(propdet + param.maxpropdet, tracer->n, (param.normbuf << 2)*sizeof(float4));

        for (i = 0; i < workdev; i++) {
            OCL_TRY(((gnode[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(FLOAT3) * (mesh->nn), mesh->node, &status), status)));
            OCL_TRY(((gelem[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(int4) * (mesh->ne), mesh->elem, &status), status)));
            OCL_TRY(((gtype[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(int) * (mesh->ne), mesh->type, &status), status)));
            OCL_TRY(((gfacenb[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(int4) * (mesh->ne), mesh->facenb, &status), status)));

            if (mesh->srcelemlen > 0) {
                OCL_TRY(((gsrcelem[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(int) * (mesh->srcelemlen), mesh->srcelem, &status), status)));
            } else {
                gsrcelem[i] = NULL;
            }

            OCL_TRY(((gnormal[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float4) * (mesh->ne) * 4, tracer->n, &status), status)));

            if (cfg->isnodalmua && cfg->nodemua) {
                OCL_TRY(((gnodemua_cl[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(cl_float) * (mesh->nn), cfg->nodemua, &status), status)));
            } else {
                gnodemua_cl[i] = (cl_mem)0;
            }

            if (cfg->isnodalmusp && cfg->nodemusp) {
                OCL_TRY(((gnodemusp_cl[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(cl_float) * (mesh->nn), cfg->nodemusp, &status), status)));
            } else {
                gnodemusp_cl[i] = (cl_mem)0;
            }

            OCL_TRY(((gproperty[i] = clCreateBuffer(mcxcontext, RO_MEM, MAX_PROP * sizeof(float4), propdet, &status), status)));
            OCL_TRY(((gparam[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(MCXParam), &param, &status), status)));

            Pseed = (cl_uint*)malloc(sizeof(cl_uint) * gpu[i].autothread * RAND_SEED_WORD_LEN);
            energy = (cl_float*)calloc(sizeof(cl_float) * cfg->srcnum, gpu[i].autothread << 1);
            MMC_ASSERT_ALLOC(Pseed);
            MMC_ASSERT_ALLOC(energy);

            for (j = 0; j < gpu[i].autothread * RAND_SEED_WORD_LEN; j++) {
                Pseed[j] = rand();
            }

            OCL_TRY(((gseed[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(cl_uint) * gpu[i].autothread * RAND_SEED_WORD_LEN, Pseed, &status), status)));
            OCL_TRY(((gweight[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * fieldlen * (isrfforward ? 4 : 2), field, &status), status)));
            OCL_TRY(((gdref[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * nflen, dref, &status), status)));
            OCL_TRY(((gdetphoton[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * cfg->maxdetphoton * hostdetreclen, Pdet, &status), status)));

            if (cfg->issaveseed) {
                OCL_TRY(((gphotonseed[i] = clCreateBuffer(mcxcontext, RW_MEM, cfg->maxdetphoton * (sizeof(RandType) * RAND_BUF_LEN), Pphotonseed, &status), status)));
            } else {
                gphotonseed[i] = NULL;
            }

            if (cfg->debuglevel & dlTraj) {
                OCL_TRY(((gdebugdata[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * (debuglen * cfg->maxjumpdebug), cfg->exportdebugdata, &status), status)));
            }

            OCL_TRY(((genergy[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * (gpu[i].autothread << 1) * cfg->srcnum, energy, &status), status)));
            OCL_TRY(((gdetected[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(cl_uint), &detected, &status), status)));
            OCL_TRY(((greporter[i] = clCreateBuffer(mcxcontext, RW_MEM, sizeof(MCXReporter), &reporter, &status), status)));

            if (cfg->srctype == MCX_SRC_PATTERN) {
                OCL_TRY(((gsrcpattern[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * (int)(cfg->srcparam1.w * cfg->srcparam2.w * cfg->srcnum), cfg->srcpattern, &status), status)));
            } else if (cfg->srctype == MCX_SRC_PATTERN3D) {
                OCL_TRY(((gsrcpattern[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * (int)(cfg->srcparam1.x * cfg->srcparam1.y * cfg->srcparam1.z * cfg->srcnum), cfg->srcpattern, &status), status)));
            } else {
                gsrcpattern[i] = NULL;
            }

            if (cfg->seed == SEED_FROM_FILE) {
                OCL_TRY(((greplayweight[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * cfg->nphoton, cfg->replayweight, &status), status)));
                OCL_TRY(((greplaytime[i] = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * cfg->nphoton, cfg->replaytime, &status), status)));
                OCL_TRY(((greplayseed[i] = clCreateBuffer(mcxcontext, RO_MEM, (sizeof(RandType) * RAND_BUF_LEN) * cfg->nphoton, cfg->photonseed, &status), status)));
            } else {
                greplayweight[i] = NULL;
                greplaytime[i] = NULL;
                greplayseed[i] = NULL;
            }

            free(Pseed);
            free(energy);
            Pseed = NULL;
            energy = NULL;
        }

        if (status != CL_SUCCESS) {
            break;    /* a per-device clCreateBuffer above failed; leave the guard */
        }

        free(propdet);
        propdet = NULL;

        mcx_printheader(cfg);

        tic = StartTimer();

        if (cfg->issavedet) {
            MMC_FPRINTF(cfg->flog, "- variant name: [%s] compiled with OpenCL version [%d]\n",
                        "MMC-OpenCL", CL_VERSION_1_0);
        } else {
            MMC_FPRINTF(cfg->flog, "- code name: [MMC-OpenCL] compiled with OpenCL version [%d]\n",
                        CL_VERSION_1_0);
        }

        MMC_FPRINTF(cfg->flog, "- compiled with: [RNG] %s [Seed Length] %d\n", MCX_RNG_NAME, RAND_SEED_WORD_LEN);
        MMC_FPRINTF(cfg->flog, "initializing streams ...\t");

        MMC_FPRINTF(cfg->flog, "init complete : %d ms\n", GetTimeMillis() - tic);
        mcx_fflush(cfg->flog);

        OCL_ASSERT(((mcxprogram = clCreateProgramWithSource(mcxcontext, 1, (const char**) & (cfg->clsource), NULL, &status), status)));

        if (cfg->optlevel >= 1) {
            snprintf(opt, MAX_JIT_OPT_LEN - strlen(opt), "%s ", "-cl-mad-enable -DMCX_USE_NATIVE");
        }

        if (cfg->optlevel >= 2) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), "%s ", "-DMCX_SIMPLIFY_BRANCH -DMCX_VECTOR_INDEX");
        }

        if (cfg->optlevel >= 3) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), "%s ", "-DUSE_MACRO_CONST");
        }

        if ((uint)cfg->srctype < sizeof(sourceflag) / sizeof(sourceflag[0])) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), "%s ", sourceflag[(uint)cfg->srctype]);
        }

        snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), "%s ", cfg->compileropt);

        if (cfg->isatomic) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DUSE_ATOMIC");
        }

        if (cfg->issave2pt == 0) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_SKIP_VOLUME");
        }

        if (cfg->issavedet) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_SAVE_DETECTORS");
        }

        if (cfg->issaveref) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_SAVE_DREF");
        }

        if (cfg->issaveseed) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_SAVE_SEED");
        }

        if (cfg->isreflect) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_DO_REFLECTION");
        }

        if (cfg->method == rtBLBadouelGrid) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DUSE_DMMC");
        }

        if (cfg->method == rtBLBadouel) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DUSE_BLBADOUEL");
        }

        if (cfg->srctype == stPattern && cfg->srcnum > 1) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DUSE_PHOTON_SHARING");
        }

        if (gpu[0].vendor == dvNVIDIA) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DUSE_NVIDIA_GPU");
        }

        /* enable per-node optical-property reads in the kernel (DOT recon) */
        if (cfg->isnodalmua) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_NODAL_MUA");
        }

        if (cfg->isnodalmusp) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DMCX_NODAL_MUSP");
        }

        /* IS_RF / IS_MULTISRC / SAVE_DETPHOTON: persistent per-photon state gates.
         * Setting these as JIT defines lets the OpenCL compiler dead-code-eliminate
         * the unused branches, matching the CUDA template-param register savings. */
        if (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DIS_RF=1");
        }

        if (cfg->srctype == stPattern || cfg->srctype == stPattern3D
                || cfg->srcnum > 1
                || (cfg->extrasrclen > 0 && cfg->srcid <= 0)) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DIS_MULTISRC=1");
        }

        if (cfg->issavedet) {
            snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt), " -DSAVE_DETPHOTON=1");
        }

        if (strstr(opt, "USE_MACRO_CONST")) {
            IPARAM_TO_MACRO(opt, param, debuglevel);
            FPARAM_TO_MACRO(opt, param, dstep);
            IPARAM_TO_MACRO(opt, param, e0);
            IPARAM_TO_MACRO(opt, param, elemlen);
            IPARAM_TO_MACRO(opt, param, framelen);
            IPARAM_TO_MACRO(opt, param, isextdet);
            IPARAM_TO_MACRO(opt, param, ismomentum);
            IPARAM_TO_MACRO(opt, param, isreflect);
            IPARAM_TO_MACRO(opt, param, issavedet);
            IPARAM_TO_MACRO(opt, param, issaveexit);
            IPARAM_TO_MACRO(opt, param, issaveref);
            IPARAM_TO_MACRO(opt, param, isspecular);
            IPARAM_TO_MACRO(opt, param, maxdetphoton);
            IPARAM_TO_MACRO(opt, param, maxjumpdebug);
            IPARAM_TO_MACRO(opt, param, maxmedia);
            IPARAM_TO_MACRO(opt, param, maxgate);
            IPARAM_TO_MACRO(opt, param, maxpropdet);
            IPARAM_TO_MACRO(opt, param, method);
            FPARAM_TO_MACRO(opt, param, minenergy);
            IPARAM_TO_MACRO(opt, param, normbuf);
            FPARAM_TO_MACRO(opt, param, nout);
            IPARAM_TO_MACRO(opt, param, outputtype);
            IPARAM_TO_MACRO(opt, param, reclen);
            FPARAM_TO_MACRO(opt, param, roulettesize);
            FPARAM_TO_MACRO(opt, param, Rtstep);
            IPARAM_TO_MACRO(opt, param, srcelemlen);
            IPARAM_TO_MACRO(opt, param, srctype);
            IPARAM_TO_MACRO(opt, param, voidtime);
            SIPARAM_TO_MACRO(opt, param, seed);
            IPARAM_TO_MACRO(opt, param, srcnum);
            IPARAM_TO_MACRO(opt, param, detnum);
            IPARAM_TO_MACRO(opt, param, issaveseed);
            IPARAM_TO_MACRO(opt, param, nf);
            IPARAM_TO_MACRO(opt, param, ne);
            IPARAM_TO_MACRO(opt, param, nn);
            FPARAM_TO_MACRO(opt, param, omega);
            FPARAM_TO_MACRO(opt, param, oneoverc0);
            SIPARAM_TO_MACRO(opt, param, srcid);
            IPARAM_TO_MACRO(opt, param, extrasrclen);
            SIPARAM_TO_MACRO(opt, param, srcpropoffset);

            if (param.focus != param.focus) {
                snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt) - 1, " -Dgcfgfocus=NAN");
            } else if (param.focus == INFINITY) {
                snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt) - 1, " -Dgcfgfocus=INFINITY");
            } else if (param.focus == -INFINITY) {
                snprintf(opt + strlen(opt), MAX_JIT_OPT_LEN - strlen(opt) - 1, " -Dgcfgfocus=-INFINITY");
            } else {
                FPARAM_TO_MACRO(opt, param, focus);
            }
        }

        MMC_FPRINTF(cfg->flog, "Building kernel with option: %s\n", opt);
        status = clBuildProgram(mcxprogram, 0, NULL, opt, NULL, NULL);

        size_t len;
        // get the details on the error, and store it in buffer
        clGetProgramBuildInfo(mcxprogram, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);

        /* Print the OpenCL build log only when the build failed or when verbose
         * timing/debug output is requested (dlTime). NVIDIA's OpenCL compiler
         * emits informational "Function X is a kernel, so overriding noinline"
         * lines on every successful build; those would otherwise clutter every
         * mmclab run. Real build errors still surface via mcx_error below. */
        if (len > 0 && (status != CL_SUCCESS || (cfg->debuglevel & dlTime))) {
            char* msg;
            int i;
            msg = (char*)calloc(len, 1);
            clGetProgramBuildInfo(mcxprogram, devices[0], CL_PROGRAM_BUILD_LOG, len, msg, NULL);

            for (i = 0; i < (int)len; i++)
                if (msg[i] <= 'z' && msg[i] >= 'A') {
                    MMC_FPRINTF(cfg->flog, "Kernel build log:\n%s\n", msg);
                    break;
                }

            free(msg);
        }

        if (status != CL_SUCCESS) {
            MMC_FPRINTF(cfg->flog, S_RED "Error: Failed to build program executable!\n" S_RESET);
            mcx_fflush(cfg->flog);
            break;   /* status holds the build error; reported after the teardown below */
        }

        MMC_FPRINTF(cfg->flog, "build program complete : %d ms\n", GetTimeMillis() - tic);
        mcx_fflush(cfg->flog);

        mcxkernel = (cl_kernel*)calloc(workdev, sizeof(cl_kernel));
        MMC_ASSERT_ALLOC(mcxkernel);

        for (i = 0; i < workdev; i++) {
            cl_int threadphoton, oddphotons;

            threadphoton = (int)(cfg->nphoton * cfg->workload[i] / (fullload * gpu[i].autothread * cfg->respin));
            oddphotons = (int)(cfg->nphoton * cfg->workload[i] / (fullload * cfg->respin) - threadphoton * gpu[i].autothread);

            MMC_FPRINTF(cfg->flog, "- [device %d(%d): %s] threadph=%d oddphotons=%d np=%.1f nthread=%d nblock=%d repetition=%d\n", i, gpu[i].id, gpu[i].name, threadphoton, oddphotons,
                        cfg->nphoton * cfg->workload[i] / fullload, (int)gpu[i].autothread, (int)gpu[i].autoblock, cfg->respin);

            MMC_FPRINTF(cfg->flog, "requesting %d bytes of shared memory\n", sharedmemsize * (int)gpu[i].autoblock);

            OCL_TRY(((mcxkernel[i] = clCreateKernel(mcxprogram, "mmc_main_loop", &status), status)));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 0, sizeof(cl_uint), (void*)&threadphoton)));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 1, sizeof(cl_uint), (void*)&oddphotons)));
            //OCL_ASSERT((clSetKernelArg(mcxkernel[i], 2, sizeof(cl_mem), (void*)(gparam+i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 3, sharedmemsize * (int)gpu[i].autoblock, NULL)));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 4, sizeof(cl_mem), (void*)(gproperty + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 5, sizeof(cl_mem), (void*)(gnode + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 6, sizeof(cl_mem), (void*)(gelem + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 7, sizeof(cl_mem), (void*)(gweight + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 8, sizeof(cl_mem), (void*)(gdref + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 9, sizeof(cl_mem), (void*)(gtype + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 10, sizeof(cl_mem), (void*)(gfacenb + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 11, sizeof(cl_mem), (void*)(gsrcelem + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 12, sizeof(cl_mem), (void*)(gnormal + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 13, sizeof(cl_mem), (void*)(gnodemua_cl + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 14, sizeof(cl_mem), (void*)(gnodemusp_cl + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 15, sizeof(cl_mem), (void*)(gdetphoton + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 16, sizeof(cl_mem), (void*)(gdetected + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 17, sizeof(cl_mem), (void*)(gseed + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 18, sizeof(cl_mem), (i == 0) ? ((void*)(gprogress)) : NULL)));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 19, sizeof(cl_mem), (void*)(genergy + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 20, sizeof(cl_mem), (void*)(greporter + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 21, sizeof(cl_mem), (void*)(gsrcpattern + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 22, sizeof(cl_mem), (void*)(greplayweight + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 23, sizeof(cl_mem), (void*)(greplaytime + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 24, sizeof(cl_mem), (void*)(greplayseed + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 25, sizeof(cl_mem), (void*)(gphotonseed + i))));
            OCL_ASSERT((clSetKernelArg(mcxkernel[i], 26, sizeof(cl_mem), ((cfg->debuglevel & dlTraj) ? (void*)(gdebugdata + i) : NULL) )));
        }

        if (status != CL_SUCCESS) {
            break;    /* a clCreateKernel above failed; leave the guard */
        }

        MMC_FPRINTF(cfg->flog, "set kernel arguments complete : %d ms %d\n", GetTimeMillis() - tic, param.method);
        mcx_fflush(cfg->flog);

        if (savevolume && cfg->exportfield == NULL) {
            cfg->exportfield = mesh->weight;
        }

        if (savevolume && cfg->exportadjoint == NULL && isrfforward) {
            /* Grid and mesh paths both need an imag fluence buffer when RF is on so
             * the adjoint Jacobian post-processing can read phi_im(node). */
            cfg->exportadjoint = (float*)calloc(fieldlen, sizeof(float));
        }

        if (cfg->exportdetected == NULL) {
            cfg->exportdetected = (float*)malloc(hostdetreclen * cfg->maxdetphoton * sizeof(float));
        }

        if (cfg->issaveseed && cfg->exportseed == NULL) {
            cfg->exportseed = (unsigned char*)malloc(cfg->maxdetphoton * (sizeof(RandType) * RAND_BUF_LEN));
        }

        if (cfg->debuglevel & dlTraj && cfg->exportdebugdata == NULL) {
            cfg->exportdebugdata = (float*)calloc(sizeof(float), (debuglen * cfg->maxjumpdebug));
        }

        cfg->energytot = (double*)calloc(cfg->srcnum, sizeof(double));
        cfg->energyesc = (double*)calloc(cfg->srcnum, sizeof(double));
        energytot = 0.0;
        energyesc = 0.0;
        cfg->runtime = 0;

        //simulate for all time-gates in maxgate groups per run

        tic0 = GetTimeMillis();

        for (t = cfg->tstart; t < cfg->tend; t += cfg->tstep * cfg->maxgate) {
            twindow0 = t;
            twindow1 = t + cfg->tstep * cfg->maxgate;

            MMC_FPRINTF(cfg->flog, "lauching mcx_main_loop for time window [%.1fns %.1fns] ...\n"
                        , twindow0 * 1e9, twindow1 * 1e9);
            mcx_fflush(cfg->flog);

            //total number of repetition for the simulations, results will be accumulated to field
            for (iter = 0; iter < cfg->respin; iter++) {
                MMC_FPRINTF(cfg->flog, "simulation run#%2d ... \n", iter + 1);
                mcx_fflush(cfg->flog);
                mcx_fflush(cfg->flog);
                param.tstart = twindow0;
                param.tend = twindow1;

                for (devid = 0; devid < workdev; devid++) {
                    OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid], gparam[devid], CL_TRUE, 0, sizeof(MCXParam), &param, 0, NULL, NULL)));
                    OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 2, sizeof(cl_mem), (void*)(gparam + devid))));

                    // launch mcxkernel
#ifndef USE_OS_TIMER

                    /* kernelevent feeds the profiling timer (GetTimeMillis) and is read
                     * repeatedly through this iteration, so it cannot be released until
                     * the next launch overwrites it. Release the previous one here to
                     * avoid leaking one cl_event per gate/repetition iteration. */
                    if (kernelevent) {
                        clReleaseEvent(kernelevent);
                        kernelevent = NULL;
                    }

                    OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid], mcxkernel[devid], 1, NULL, &gpu[devid].autothread, &gpu[devid].autoblock, 0, NULL, &kernelevent)));
#else
                    /* pass NULL for the event: the launch is synchronized below via
                     * clFinish (clWaitForEvents is disabled), so a returned cl_event
                     * would never be consumed nor released. On NVIDIA's OpenCL runtime
                     * a live event keeps an implicit reference to the queue/context,
                     * so leaking it here prevents clReleaseContext from ever freeing
                     * the GPU allocations -> per-run device-memory leak. */
                    OCL_ASSERT((clEnqueueNDRangeKernel(mcxqueue[devid], mcxkernel[devid], 1, NULL, &gpu[devid].autothread, &gpu[devid].autoblock, 0, NULL, NULL)));
#endif
                    OCL_ASSERT((clFlush(mcxqueue[devid])));
                }

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
                    } while (p0 < gpu[0].autothread);

                    mcx_progressbar(cfg->nphoton);
                    MMC_FPRINTF(cfg->flog, "\n");
                }

                /* NOTE: gprogress is mapped once before the loops and unmapped once in
                 * cleanup. Unmapping here (inside the gate/repetition loop) was a latent
                 * use-after-unmap: the next iteration's progress poll reads *progress
                 * after the buffer had already been unmapped (only safe today because
                 * the default single-gate/single-repetition case unmaps exactly once). */
                for (devid = 0; devid < workdev; devid++) {
                    OCL_ASSERT((clFinish(mcxqueue[devid])));
                }

                tic1 = GetTimeMillis();
                toc += tic1 - tic0;
                MMC_FPRINTF(cfg->flog, "kernel complete:  \t%d ms\nretrieving flux ... \t", tic1 - tic);
                mcx_fflush(cfg->flog);

                if (cfg->runtime < tic1 - tic) {
                    cfg->runtime = tic1 - tic;
                }

                for (devid = 0; devid < workdev; devid++) {
                    MCXReporter rep;
                    OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], greporter[devid], CL_TRUE, 0, sizeof(MCXReporter),
                                                    &rep, 0, NULL, NULL)));
                    reporter.raytet += rep.raytet;
                    reporter.jumpdebug += rep.jumpdebug;

                    energy = (cl_float*)calloc(sizeof(cl_float) * cfg->srcnum, gpu[devid].autothread << 1);
                    OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], genergy[devid], CL_TRUE, 0, sizeof(cl_float) * (gpu[devid].autothread << 1) * cfg->srcnum,
                                                    energy, 0, NULL, NULL)));

                    for (i = 0; i < gpu[devid].autothread; i++) {
                        for (j = 0; j < (int)cfg->srcnum; j++) {
                            cfg->energyesc[j] += energy[(i << 1) * cfg->srcnum + j];
                            cfg->energytot[j] += energy[((i << 1) + 1) * cfg->srcnum + j];
                            energyesc += energy[(i << 1) * cfg->srcnum + j];
                            energytot += energy[((i << 1) + 1) * cfg->srcnum + j];
                        }
                    }

                    free(energy);
                    energy = NULL;

                    if (cfg->debuglevel & dlTraj) {
                        uint debugrec = rep.jumpdebug;

                        if (debugrec > 0) {
                            if (debugrec > cfg->maxjumpdebug) {
                                MMC_FPRINTF(cfg->flog, S_RED "WARNING: the saved trajectory positions (%d) \
  are more than what your have specified (%d), please use the --maxjumpdebug option to specify a greater number\n" S_RESET
                                            , debugrec, cfg->maxjumpdebug);
                            } else {
                                MMC_FPRINTF(cfg->flog, "saved %u trajectory positions, total: %d\t", debugrec, cfg->debugdatalen + debugrec);
                            }

                            debugrec = MIN(debugrec, cfg->maxjumpdebug);
                            cfg->exportdebugdata = (float*)realloc(cfg->exportdebugdata, (cfg->debugdatalen + debugrec) * debuglen * sizeof(float));
                            OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gdebugdata[devid], CL_FALSE, 0, sizeof(float)*debuglen * debugrec,
                                                            cfg->exportdebugdata + cfg->debugdatalen, 0, NULL, NULL)));
                            cfg->debugdatalen += debugrec;
                        }
                    }

                    if (cfg->issavedet) {
                        OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gdetected[devid], CL_FALSE, 0, sizeof(uint),
                                                        &detected, 0, NULL, NULL)));
                        OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gdetphoton[devid], CL_TRUE, 0, sizeof(float)*cfg->maxdetphoton * hostdetreclen,
                                                        Pdet, 0, NULL, NULL)));

                        if (cfg->issaveseed) {
                            OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gphotonseed[devid], CL_TRUE, 0, cfg->maxdetphoton * (sizeof(RandType)*RAND_BUF_LEN),
                                                            Pphotonseed, 0, NULL, NULL)));
                        }

                        if (detected > cfg->maxdetphoton) {
                            MMC_FPRINTF(cfg->flog, "WARNING: the detected photon (%d) \
is more than what your have specified (%d), please use the -H option to specify a greater number\t"
                                        , detected, cfg->maxdetphoton);
                        } else {
                            MMC_FPRINTF(cfg->flog, "detected %d photons, total: %d\t", detected, cfg->detectedcount + detected);
                        }

                        cfg->his.detected += detected;
                        detected = MIN(detected, cfg->maxdetphoton);

                        if (cfg->exportdetected) {
                            cfg->exportdetected = (float*)realloc(cfg->exportdetected, (cfg->detectedcount + detected) * hostdetreclen * sizeof(float));
                            memcpy(cfg->exportdetected + cfg->detectedcount * (hostdetreclen), Pdet, detected * (hostdetreclen)*sizeof(float));

                            if (cfg->issaveseed) {
                                cfg->exportseed = (unsigned char*)realloc(cfg->exportseed, (cfg->detectedcount + detected) * (sizeof(RandType) * RAND_BUF_LEN));
                                memcpy(cfg->exportseed + cfg->detectedcount * sizeof(RandType)*RAND_BUF_LEN, Pphotonseed, detected * (sizeof(RandType)*RAND_BUF_LEN));
                            }

                            cfg->detectedcount += detected;
                        }
                    }

                    if (cfg->issaveref) {
                        float* rawdref = (float*)calloc(sizeof(float), nflen);
                        OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gdref[devid], CL_TRUE, 0, sizeof(float)*nflen,
                                                        rawdref, 0, NULL, NULL)));

                        //TODO: saving dref has not yet adopting double-buffer
                        for (i = 0; i < nflen; i++) { //accumulate field, can be done in the GPU
                            dref[i] += rawdref[i];    //+rawfield[i+fieldlen];
                        }

                        free(rawdref);
                    }

                    //handling the 2pt distributions
                    if (cfg->issave2pt) {
                        int rfmul = isrfforward ? 4 : 2;
                        float* rawfield = (float*)malloc(sizeof(float) * fieldlen * rfmul);

                        OCL_ASSERT((clEnqueueReadBuffer(mcxqueue[devid], gweight[devid], CL_TRUE, 0, sizeof(cl_float)*fieldlen * rfmul,
                                                        rawfield, 0, NULL, NULL)));
                        MMC_FPRINTF(cfg->flog, "transfer complete:        %d ms\n", GetTimeMillis() - tic);
                        mcx_fflush(cfg->flog);

                        for (i = 0; i < fieldlen; i++) { //accumulate real field
                            field[i] += rawfield[i] + rawfield[i + fieldlen];
                        }

                        if (isrfforward) {
                            for (i = 0; i < fieldlen; i++) { //accumulate imaginary part
                                field[i + fieldlen] += rawfield[i + fieldlen * 2] + rawfield[i + fieldlen * 3];
                            }
                        }

                        free(rawfield);
                    }

                    if (cfg->respin > 1 && RAND_SEED_WORD_LEN > 1) {
                        Pseed = (cl_uint*)malloc(sizeof(cl_uint) * gpu[devid].autothread * RAND_SEED_WORD_LEN);

                        for (i = 0; i < gpu[devid].autothread * RAND_SEED_WORD_LEN; i++) {
                            Pseed[i] = rand();
                        }

                        OCL_ASSERT((clEnqueueWriteBuffer(mcxqueue[devid], gseed[devid], CL_TRUE, 0, sizeof(cl_uint)*gpu[devid].autothread * RAND_SEED_WORD_LEN,
                                                         Pseed, 0, NULL, NULL)));
                        OCL_ASSERT((clSetKernelArg(mcxkernel[devid], 15, sizeof(cl_mem), (void*)(gseed + devid))));
                        free(Pseed);
                        Pseed = NULL;
                    }

                    OCL_ASSERT((clFinish(mcxqueue[devid])));
                }// loop over work devices
            }// iteration
        }// time gates

        if (cfg->exportfield) {
            if (cfg->basisorder == 0 || cfg->method == rtBLBadouelGrid) {
                for (i = 0; i < fieldlen; i++) {
                    cfg->exportfield[i] += field[i];
                }

                /* RF forward: store imaginary part in exportadjoint */
                if (isrfforward && cfg->exportadjoint) {
                    for (i = 0; i < fieldlen; i++) {
                        cfg->exportadjoint[i] += field[i + fieldlen];
                    }
                }
            } else {
                int srcid;
                /* basisorder=1 mesh mode: redistribute per-element fluence onto nodes.
                 * Kernel layout per slot (adjoint, srcnum=1, nsrcslots>1):
                 *   field[eid + gate*ne + slot*ne*maxgate]
                 * Pattern source layout (srcnum>1, nsrcslots=1):
                 *   field[(gate*ne + eid)*srcnum + pidx]
                 * Output (exportfield) per-slot stride = nn*maxgate to match grid convention.
                 */
                cl_uint nslots = (cfg->extrasrclen > 0) ? (cl_uint)cfg->extrasrclen : 1u;

                if (nslots > 1u && cfg->srcnum == 1) {
                    for (cl_uint slot = 0; slot < nslots; slot++) {
                        size_t slot_off_f = (size_t)slot * (size_t)mesh->ne * (size_t)cfg->maxgate;
                        size_t slot_off_e = (size_t)slot * (size_t)mesh->nn * (size_t)cfg->maxgate;

                        for (i = 0; i < cfg->maxgate; i++) {
                            size_t f_gate_off = (size_t)i * (size_t)mesh->ne;
                            size_t e_gate_off = (size_t)i * (size_t)mesh->nn;

                            for (j = 0; j < mesh->ne; j++) {
                                float ww = field[slot_off_f + f_gate_off + j] * 0.25f;
                                float ww_im = (isrfforward && cfg->exportadjoint)
                                              ? field[slot_off_f + f_gate_off + j + fieldlen] * 0.25f : 0.f;
                                int k;

                                for (k = 0; k < mesh->elemlen; k++) {
                                    size_t nidx = slot_off_e + e_gate_off
                                                  + mesh->elem[j * mesh->elemlen + k] - 1;
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
                                int k;

                                for (k = 0; k < mesh->elemlen; k++) {
                                    cfg->exportfield[(i * mesh->nn + mesh->elem[j * mesh->elemlen + k] - 1) * cfg->srcnum + srcid] += ww;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (cfg->issaveref && mesh->dref) {
            for (i = 0; i < nflen; i++) {
                mesh->dref[i] += dref[i];
            }
        }

        if (cfg->isnormalized && savevolume) {
            double cur_normalizer, sum_normalizer = 0.0, energyabs = 0.0;

            for (j = 0; j < cfg->srcnum; j++) {
                energyabs =  cfg->energytot[j] - cfg->energyesc[j];
                cur_normalizer = mesh_normalize(mesh, cfg, energyabs, cfg->energytot[j], j);
                sum_normalizer += cur_normalizer;
                MMCDEBUG(cfg, dlTime, (cfg->flog, "source %d\ttotal simulated energy: %f\tabsorbed: "S_BOLD""S_BLUE"%5.5f%%"S_RESET"\tnormalizor=%g\n",
                                       j + 1, cfg->energytot[j], 100.f * energyabs / cfg->energytot[j], cur_normalizer));
            }

            cfg->his.normalizer = sum_normalizer / cfg->srcnum;

            /* Broadcast normalizor from slot 0 to the remaining adjoint slots
             * (srcnum..extrasrclen-1). mesh_normalize only processes pair=0..srcnum-1;
             * the multi-source/adjoint slots get scale = avg_normalizor (uniform launch
             * weight means w0/wk = 1, so no per-slot weight ratio). Includes the per-node
             * nvol division for basisorder=1 mesh mode that mesh_normalize would have
             * applied to slot 0. */
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

            /* RF imag fluence: mesh_normalize applies the per-node nvol division +
             * scalar normalizor to slot 0 of cfg->exportadjoint. Mirror that for the
             * adjoint detector slots here (srcnum..extrasrclen-1). For mesh basis1
             * mode the imag fluence layout matches the real (slot-major), so we use
             * the same per-slot stride. */
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

        /* Mesh-mode adjoint Jacobian post-processing (OpenCL): FEM/nodal formulas on a tet mesh. */
        if (cfg->issave2pt && MCX_IS_ADJOINT_TYPE(cfg->outputtype) && cfg->method != rtBLBadouelGrid &&
                cfg->basisorder && cfg->extrasrclen > 0 && cfg->detdir != NULL && cfg->exportfield) {
            unsigned int Ns          = (unsigned int)(cfg->extrasrclen - cfg->detnum);
            unsigned int Nd          = (unsigned int)cfg->detnum;
            int isdual               = MCX_IS_DUAL_ADJOINT_TYPE(cfg->outputtype);
            int isnodal_approx       = (cfg->adjointmode == 1);
            unsigned int datalen     = (unsigned int)mesh->nn;
            size_t adjointlen        = (size_t)datalen * Ns * Nd;
            size_t single_exportlen  = adjointlen * (isrfforward ? 2 : 1);
            size_t exportlen_adj     = single_exportlen * (isdual ? 2 : 1);
            int compute_jmua = (cfg->outputtype == otAdjoint || isdual);
            int compute_jd   = (cfg->outputtype == otAdjointDcoeff || isdual);

            /* exportfield (= mesh->weight) is node-major for basisorder=1, sized
             * nn * extrasrclen * maxgate doubles -- NOT fieldlen, which is sized
             * on mesh->ne (the raw kernel field buffer stride). Reading fieldlen
             * doubles walks past the end. */
            size_t nodefield_len = (size_t)datalen * (size_t)cfg->extrasrclen * (size_t)cfg->maxgate;

            /* Convert double exportfield to float */
            float* hfield_re = (float*)malloc(sizeof(float) * nodefield_len);

            for (size_t k = 0; k < nodefield_len; k++) {
                hfield_re[k] = (float)cfg->exportfield[k];
            }

            cl_mem gcl_field_re = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * nodefield_len, hfield_re, &status);
            OCL_ASSERT(status);
            free(hfield_re);

            cl_mem gcl_field_im = (cl_mem)0;

            if (isrfforward && cfg->exportadjoint) {
                gcl_field_im = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * nodefield_len, cfg->exportadjoint, &status);
                OCL_ASSERT(status);
            }

            /* Upload mesh helpers */
            cl_mem gcl_elem = clCreateBuffer(mcxcontext, RO_MEM, sizeof(int) * mesh->ne * mesh->elemlen, mesh->elem, &status);
            OCL_ASSERT(status);
            cl_mem gcl_evol = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * mesh->ne, mesh->evol, &status);
            OCL_ASSERT(status);

            cl_mem gcl_deldotdel = (cl_mem)0, gcl_nvol = (cl_mem)0;

            if (compute_jd || !isnodal_approx) {
                if (mesh->deldotdel == NULL) {
                    mesh_deldotdel(mesh);
                }

                float* deldotdel_f = (float*)malloc(sizeof(float) * mesh->ne * 10);

                for (size_t k = 0; k < (size_t)mesh->ne * 10; k++) {
                    deldotdel_f[k] = (float)mesh->deldotdel[k];
                }

                gcl_deldotdel = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * mesh->ne * 10, deldotdel_f, &status);
                OCL_ASSERT(status);
                free(deldotdel_f);
            }

            if (isnodal_approx && compute_jmua) {
                gcl_nvol = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * mesh->nn, mesh->nvol, &status);
                OCL_ASSERT(status);
            }

            /* Zeroed output buffers */
            float* hzero = (float*)calloc(single_exportlen, sizeof(float));
            cl_mem gcl_jmua = (cl_mem)0, gcl_jd = (cl_mem)0;

            if (compute_jmua) {
                gcl_jmua = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * single_exportlen, hzero, &status);
                OCL_ASSERT(status);
            }

            if (compute_jd) {
                gcl_jd = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * single_exportlen, hzero, &status);
                OCL_ASSERT(status);
            }

            free(hzero);

            size_t adj_local = 256;

            /* nodal-approximation J_mua kernel (one thread per node) */
            if (isnodal_approx && compute_jmua) {
                cl_kernel kern = clCreateKernel(mcxprogram, "mmc_adjoint_mesh_nodal_kernel", &status);
                OCL_ASSERT(status);
                cl_uint nn_arg = (cl_uint)mesh->nn, maxgate_arg = (cl_uint)cfg->maxgate;
                OCL_ASSERT(clSetKernelArg(kern, 0, sizeof(cl_mem), &gcl_field_re));
                OCL_ASSERT(clSetKernelArg(kern, 1, sizeof(cl_mem), gcl_field_im ? &gcl_field_im : NULL));
                OCL_ASSERT(clSetKernelArg(kern, 2, sizeof(cl_mem), &gcl_nvol));
                OCL_ASSERT(clSetKernelArg(kern, 3, sizeof(cl_mem), &gcl_jmua));
                OCL_ASSERT(clSetKernelArg(kern, 4, sizeof(cl_uint), &nn_arg));
                OCL_ASSERT(clSetKernelArg(kern, 5, sizeof(cl_uint), &maxgate_arg));
                OCL_ASSERT(clSetKernelArg(kern, 6, sizeof(cl_uint), &Ns));
                OCL_ASSERT(clSetKernelArg(kern, 7, sizeof(cl_uint), &Nd));
                size_t adj_global = ((mesh->nn + adj_local - 1) / adj_local) * adj_local;
                OCL_ASSERT(clEnqueueNDRangeKernel(mcxqueue[0], kern, 1, NULL, &adj_global, &adj_local, 0, NULL, NULL));
                OCL_ASSERT(clFinish(mcxqueue[0]));
                clReleaseKernel(kern);
            }

            /* Full FEM kernel (one thread per element) */
            if ((!isnodal_approx && compute_jmua) || compute_jd) {
                cl_kernel kern = clCreateKernel(mcxprogram, "mmc_adjoint_mesh_full_kernel", &status);
                OCL_ASSERT(status);
                cl_uint ne_arg = (cl_uint)mesh->ne, nn_arg = (cl_uint)mesh->nn;
                cl_uint maxgate_arg = (cl_uint)cfg->maxgate;
                cl_uint elemlen_arg = (cl_uint)mesh->elemlen;
                cl_int  isnodal_arg = 1; /* output is nodal */
                cl_mem nulljmua = (cl_mem)0;
                cl_mem* jmua_arg = isnodal_approx ? &nulljmua : &gcl_jmua;
                OCL_ASSERT(clSetKernelArg(kern,  0, sizeof(cl_mem), &gcl_field_re));
                OCL_ASSERT(clSetKernelArg(kern,  1, sizeof(cl_mem), gcl_field_im ? &gcl_field_im : NULL));
                OCL_ASSERT(clSetKernelArg(kern,  2, sizeof(cl_mem), &gcl_elem));
                OCL_ASSERT(clSetKernelArg(kern,  3, sizeof(cl_mem), &gcl_evol));
                OCL_ASSERT(clSetKernelArg(kern,  4, sizeof(cl_mem), &gcl_deldotdel));
                OCL_ASSERT(clSetKernelArg(kern,  5, sizeof(cl_mem), isnodal_approx ? NULL : jmua_arg));
                OCL_ASSERT(clSetKernelArg(kern,  6, sizeof(cl_mem), compute_jd ? &gcl_jd : NULL));
                OCL_ASSERT(clSetKernelArg(kern,  7, sizeof(cl_uint), &ne_arg));
                OCL_ASSERT(clSetKernelArg(kern,  8, sizeof(cl_uint), &nn_arg));
                OCL_ASSERT(clSetKernelArg(kern,  9, sizeof(cl_uint), &maxgate_arg));
                OCL_ASSERT(clSetKernelArg(kern, 10, sizeof(cl_uint), &Ns));
                OCL_ASSERT(clSetKernelArg(kern, 11, sizeof(cl_uint), &Nd));
                OCL_ASSERT(clSetKernelArg(kern, 12, sizeof(cl_uint), &elemlen_arg));
                OCL_ASSERT(clSetKernelArg(kern, 13, sizeof(cl_int),  &isnodal_arg));
                size_t adj_global = ((mesh->ne + adj_local - 1) / adj_local) * adj_local;
                OCL_ASSERT(clEnqueueNDRangeKernel(mcxqueue[0], kern, 1, NULL, &adj_global, &adj_local, 0, NULL, NULL));
                OCL_ASSERT(clFinish(mcxqueue[0]));
                clReleaseKernel(kern);
            }

            OCL_ASSERT(clReleaseMemObject(gcl_field_re));

            if (gcl_field_im) {
                OCL_ASSERT(clReleaseMemObject(gcl_field_im));
            }

            clReleaseMemObject(gcl_elem);
            clReleaseMemObject(gcl_evol);

            if (gcl_deldotdel) {
                clReleaseMemObject(gcl_deldotdel);
            }

            if (gcl_nvol) {
                clReleaseMemObject(gcl_nvol);
            }

            /* Allocate output buffer */
            if (cfg->exportjacob) {
                free(cfg->exportjacob);
            }

            cfg->exportjacob = (float*)malloc(sizeof(float) * exportlen_adj);
            memset(cfg->exportjacob, 0, sizeof(float) * exportlen_adj);

            float* hmua = NULL, *hjd = NULL;

            if (compute_jmua) {
                hmua = (float*)malloc(sizeof(float) * single_exportlen);
                OCL_ASSERT(clEnqueueReadBuffer(mcxqueue[0], gcl_jmua, CL_TRUE, 0, sizeof(float) * single_exportlen, hmua, 0, NULL, NULL));
                clReleaseMemObject(gcl_jmua);
            }

            if (compute_jd) {
                hjd = (float*)malloc(sizeof(float) * single_exportlen);
                OCL_ASSERT(clEnqueueReadBuffer(mcxqueue[0], gcl_jd, CL_TRUE, 0, sizeof(float) * single_exportlen, hjd, 0, NULL, NULL));
                clReleaseMemObject(gcl_jd);
            }

            if (isdual) {
                if (!isrfforward) {
                    memcpy(cfg->exportjacob,              hmua, adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen, hjd,  adjointlen * sizeof(float));
                } else {
                    memcpy(cfg->exportjacob,                  hmua,              adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen,     hjd,               adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 2 * adjointlen, hmua + adjointlen, adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 3 * adjointlen, hjd  + adjointlen, adjointlen * sizeof(float));
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
                mesh_savejacob(cfg, mesh, cfg->exportjacob, (int)Ns, (int)Nd, isrfforward, isdual);
                MMC_FPRINTF(cfg->flog, "saving Jacobian complete : %d ms\n\n", GetTimeMillis() - tic);
                mcx_fflush(cfg->flog);
            }

#endif
        }

        /* Adjoint Jacobian post-processing: OpenCL equivalent of mmc_cu_host.cu lines 890-1007.
         * After simulation, exportfield has normalized multi-slot real fluence and exportadjoint
         * has normalized imaginary fluence (RF only). Launch adjoint/dcoeff kernels on GPU to
         * compute J_mua = phi_src * phi_det and/or J_D = grad(phi_src) . grad(phi_det). */
        if (cfg->issave2pt && MCX_IS_ADJOINT_TYPE(cfg->outputtype) && cfg->method == rtBLBadouelGrid &&
                cfg->extrasrclen > 0 && cfg->detdir != NULL && cfg->exportfield) {
            unsigned int Ns         = (unsigned int)(cfg->extrasrclen - cfg->detnum);
            unsigned int Nd         = (unsigned int)cfg->detnum;
            unsigned int pure_voxels = (unsigned int)cfg->crop0.z;
            int isdual              = MCX_IS_DUAL_ADJOINT_TYPE(cfg->outputtype);

            size_t adjointlen       = (size_t)pure_voxels * Ns * Nd;
            size_t single_exportlen = adjointlen * (isrfforward ? 2 : 1);
            size_t exportlen_adj    = single_exportlen * (isdual ? 2 : 1);

            /* Build float host arrays from exportfield (double) */
            float* hfield_re = (float*)malloc(sizeof(float) * fieldlen);

            for (size_t k = 0; k < fieldlen; k++) {
                hfield_re[k] = (float)cfg->exportfield[k];
            }

            /* Upload real field to GPU */
            cl_mem gcl_field_re = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * fieldlen, hfield_re, &status);
            OCL_ASSERT(status);
            free(hfield_re);

            /* Upload imaginary field to GPU (RF only) */
            cl_mem gcl_field_im = (cl_mem)0;

            if (isrfforward && cfg->exportadjoint) {
                gcl_field_im = clCreateBuffer(mcxcontext, RO_MEM, sizeof(float) * fieldlen, cfg->exportadjoint, &status);
                OCL_ASSERT(status);
            }

            /* Allocate zeroed output buffers for Jacobian */
            float* hzero = (float*)calloc(single_exportlen, sizeof(float));

            cl_mem gcl_adj_mua = (cl_mem)0;

            if (isdual) {
                gcl_adj_mua = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * single_exportlen, hzero, &status);
                OCL_ASSERT(status);
            }

            cl_mem gcl_adj_tmp = clCreateBuffer(mcxcontext, RW_MEM, sizeof(float) * single_exportlen, hzero, &status);
            OCL_ASSERT(status);
            free(hzero);

            size_t adj_local  = 256;
            size_t adj_global = ((pure_voxels + (unsigned int)adj_local - 1) / (unsigned int)adj_local) * adj_local;

            /* Launch mmc_adjoint_kernel (J_mua = phi_src * phi_det) for otAdjoint or dual modes */
            if (isdual || cfg->outputtype == otAdjoint) {
                cl_kernel adj_kern = clCreateKernel(mcxprogram, "mmc_adjoint_kernel", &status);
                OCL_ASSERT(status);
                cl_mem dest = isdual ? gcl_adj_mua : gcl_adj_tmp;
                OCL_ASSERT(clSetKernelArg(adj_kern, 0, sizeof(cl_mem), &gcl_field_re));
                OCL_ASSERT(clSetKernelArg(adj_kern, 1, sizeof(cl_mem), gcl_field_im ? &gcl_field_im : NULL));
                OCL_ASSERT(clSetKernelArg(adj_kern, 2, sizeof(cl_mem), &dest));
                OCL_ASSERT(clSetKernelArg(adj_kern, 3, sizeof(cl_uint), &pure_voxels));
                OCL_ASSERT(clSetKernelArg(adj_kern, 4, sizeof(cl_uint), (cl_uint*)&cfg->maxgate));
                OCL_ASSERT(clSetKernelArg(adj_kern, 5, sizeof(cl_uint), &Ns));
                OCL_ASSERT(clSetKernelArg(adj_kern, 6, sizeof(cl_uint), &Nd));
                OCL_ASSERT(clEnqueueNDRangeKernel(mcxqueue[0], adj_kern, 1, NULL, &adj_global, &adj_local, 0, NULL, NULL));
                OCL_ASSERT(clFinish(mcxqueue[0]));
                clReleaseKernel(adj_kern);
            }

            /* Launch mmc_adjoint_dcoeff_kernel (J_D = -grad(phi_src).grad(phi_det)) for dcoeff or dual */
            if (isdual || cfg->outputtype != otAdjoint) {
                cl_kernel dcoeff_kern = clCreateKernel(mcxprogram, "mmc_adjoint_dcoeff_kernel", &status);
                OCL_ASSERT(status);
                cl_uint Nx = (cl_uint)cfg->dim.x, Ny = (cl_uint)cfg->dim.y;
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 0, sizeof(cl_mem), &gcl_field_re));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 1, sizeof(cl_mem), gcl_field_im ? &gcl_field_im : NULL));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 2, sizeof(cl_mem), &gcl_adj_tmp));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 3, sizeof(cl_uint), &pure_voxels));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 4, sizeof(cl_uint), (cl_uint*)&cfg->maxgate));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 5, sizeof(cl_uint), &Ns));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 6, sizeof(cl_uint), &Nd));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 7, sizeof(cl_uint), &Nx));
                OCL_ASSERT(clSetKernelArg(dcoeff_kern, 8, sizeof(cl_uint), &Ny));
                OCL_ASSERT(clEnqueueNDRangeKernel(mcxqueue[0], dcoeff_kern, 1, NULL, &adj_global, &adj_local, 0, NULL, NULL));
                OCL_ASSERT(clFinish(mcxqueue[0]));
                clReleaseKernel(dcoeff_kern);
            }

            OCL_ASSERT(clReleaseMemObject(gcl_field_re));

            if (gcl_field_im) {
                OCL_ASSERT(clReleaseMemObject(gcl_field_im));
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
                OCL_ASSERT(clEnqueueReadBuffer(mcxqueue[0], gcl_adj_mua, CL_TRUE, 0, sizeof(float) * single_exportlen, hmua,    0, NULL, NULL));
                OCL_ASSERT(clEnqueueReadBuffer(mcxqueue[0], gcl_adj_tmp, CL_TRUE, 0, sizeof(float) * single_exportlen, hsecond, 0, NULL, NULL));
                clReleaseMemObject(gcl_adj_mua);
                clReleaseMemObject(gcl_adj_tmp);

                for (size_t k = 0; k < single_exportlen; k++) {
                    hmua[k]    *= -Vvox;
                    hsecond[k] *= -(cfg->unitinmm);
                }

                if (!isrfforward) {
                    memcpy(cfg->exportjacob,              hmua,    adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen, hsecond, adjointlen * sizeof(float));
                } else {
                    memcpy(cfg->exportjacob,                   hmua,                 adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + adjointlen,      hsecond,              adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 2 * adjointlen,  hmua    + adjointlen, adjointlen * sizeof(float));
                    memcpy(cfg->exportjacob + 3 * adjointlen,  hsecond + adjointlen, adjointlen * sizeof(float));
                }

                free(hmua);
                free(hsecond);
            } else {
                OCL_ASSERT(clEnqueueReadBuffer(mcxqueue[0], gcl_adj_tmp, CL_TRUE, 0, sizeof(float) * single_exportlen, cfg->exportjacob, 0, NULL, NULL));
                clReleaseMemObject(gcl_adj_tmp);

                float adj_scale = (cfg->outputtype == otAdjoint) ? -Vvox : -(cfg->unitinmm);

                for (size_t k = 0; k < single_exportlen; k++) {
                    cfg->exportjacob[k] *= adj_scale;
                }
            }

            MMC_FPRINTF(cfg->flog, "adjoint Jacobian computation complete: %d ms\n", GetTimeMillis() - tic);

#ifndef MCX_CONTAINER

            if (cfg->issave2pt && cfg->parentid == mpStandalone && cfg->exportjacob) {
                MMC_FPRINTF(cfg->flog, "saving adjoint Jacobian to file ...\t");
                mesh_savejacob(cfg, mesh, cfg->exportjacob, (int)Ns, (int)Nd, isrfforward, isdual);
                MMC_FPRINTF(cfg->flog, "saving Jacobian complete : %d ms\n\n", GetTimeMillis() - tic);
                mcx_fflush(cfg->flog);
            }

#endif
        }

#ifndef MCX_CONTAINER

        if (cfg->issave2pt && cfg->parentid == mpStandalone) {
            MMC_FPRINTF(cfg->flog, "saving data to file ...\t");
            mesh_saveweight(mesh, cfg, 0);
            MMC_FPRINTF(cfg->flog, "saving data complete : %d ms\n\n", GetTimeMillis() - tic);
            mcx_fflush(cfg->flog);
        }

        if (cfg->issavedet && cfg->parentid == mpStandalone && cfg->exportdetected) {
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

        if ((cfg->debuglevel & dlTraj) && cfg->parentid == mpStandalone && cfg->exportdebugdata) {
            MMC_FPRINTF(cfg->flog, "saving trajectory data to file ...\t");
            cfg->his.colcount = MCX_DEBUG_REC_LEN;
            cfg->his.savedphoton = cfg->debugdatalen;
            cfg->his.totalphoton = cfg->nphoton;
            cfg->his.detected = 0;
            mcx_savedetphoton(cfg->exportdebugdata, NULL, cfg->debugdatalen, 0, cfg);
            MMC_FPRINTF(cfg->flog, "saving trajectory data complete : %d ms\n\n", GetTimeMillis() - tic);
        }

        if (cfg->issaveref && cfg->parentid == mpStandalone) {
            MMC_FPRINTF(cfg->flog, "saving surface diffuse reflectance ...");
            mesh_saveweight(mesh, cfg, 1);
        }

#endif

        // total energy here equals total simulated photons+unfinished photons for all threads
        MMC_FPRINTF(cfg->flog, "simulated %zu photons (%zu) with %d devices (ray-tet %.0f)\nMCX simulation speed: %.2f photon/ms\n",
                    cfg->nphoton, cfg->nphoton, workdev, reporter.raytet, (double)cfg->nphoton / toc);
        MMC_FPRINTF(cfg->flog, "total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to initial specular reflection is excluded in the total)\n",
                    energytot, (energytot - energyesc) / energytot * 100.f);
        mcx_fflush(cfg->flog);

        status = CL_SUCCESS;   /* reached the end of the run with no error */
    } while (0);

    /* ---- single NULL-safe teardown, reached on both the success and error paths ----
     * Releases only what was actually created (handle arrays are calloc'd to NULL,
     * scalar handles are NULL-initialized), uses plain clRelease* (not OCL_ASSERT) so a
     * stray release error cannot mask the original failure, and re-reports any deferred
     * acquisition error at the very end. */
    if (progress && mcxqueue && mcxqueue[0] && gprogress && gprogress[0]) {
        /* gprogress is mapped once before the run and unmapped exactly once here */
        clEnqueueUnmapMemObject(mcxqueue[0], gprogress[0], progress, 0, NULL, NULL);
        clFinish(mcxqueue[0]);
    }

    if (gprogress && gprogress[0]) {
        clReleaseMemObject(gprogress[0]);
    }

    for (i = 0; i < workdev; i++) {
        if (gseed && gseed[i])               {
            clReleaseMemObject(gseed[i]);
        }

        if (gdetphoton && gdetphoton[i])     {
            clReleaseMemObject(gdetphoton[i]);
        }

        if (gweight && gweight[i])           {
            clReleaseMemObject(gweight[i]);
        }

        if (gdref && gdref[i])               {
            clReleaseMemObject(gdref[i]);
        }

        if (genergy && genergy[i])           {
            clReleaseMemObject(genergy[i]);
        }

        if (gdetected && gdetected[i])       {
            clReleaseMemObject(gdetected[i]);
        }

        if (gphotonseed && gphotonseed[i])   {
            clReleaseMemObject(gphotonseed[i]);
        }

        if (greporter && greporter[i])       {
            clReleaseMemObject(greporter[i]);
        }

        if (gnode && gnode[i])               {
            clReleaseMemObject(gnode[i]);
        }

        if (gelem && gelem[i])               {
            clReleaseMemObject(gelem[i]);
        }

        if (gtype && gtype[i])               {
            clReleaseMemObject(gtype[i]);
        }

        if (gfacenb && gfacenb[i])           {
            clReleaseMemObject(gfacenb[i]);
        }

        if (gsrcelem && gsrcelem[i])         {
            clReleaseMemObject(gsrcelem[i]);
        }

        if (gnormal && gnormal[i])           {
            clReleaseMemObject(gnormal[i]);
        }

        if (gnodemua_cl && gnodemua_cl[i])   {
            clReleaseMemObject(gnodemua_cl[i]);
        }

        if (gnodemusp_cl && gnodemusp_cl[i]) {
            clReleaseMemObject(gnodemusp_cl[i]);
        }

        if (gproperty && gproperty[i])       {
            clReleaseMemObject(gproperty[i]);
        }

        if (gparam && gparam[i])             {
            clReleaseMemObject(gparam[i]);
        }

        if (gsrcpattern && gsrcpattern[i])   {
            clReleaseMemObject(gsrcpattern[i]);
        }

        if (gdebugdata && gdebugdata[i])     {
            clReleaseMemObject(gdebugdata[i]);
        }

        if (greplayweight && greplayweight[i]) {
            clReleaseMemObject(greplayweight[i]);
        }

        if (greplayseed && greplayseed[i])   {
            clReleaseMemObject(greplayseed[i]);
        }

        if (greplaytime && greplaytime[i])   {
            clReleaseMemObject(greplaytime[i]);
        }

        if (mcxkernel && mcxkernel[i])       {
            clReleaseKernel(mcxkernel[i]);
        }
    }

    free(gseed);
    free(gdetphoton);
    free(gweight);
    free(gdref);
    free(genergy);
    free(gprogress);
    free(gdetected);
    free(gphotonseed);
    free(gdebugdata);
    free(greporter);

    free(gnode);
    free(gelem);
    free(gtype);
    free(gfacenb);
    free(gsrcelem);
    free(gnormal);
    free(gnodemua_cl);
    free(gnodemusp_cl);
    free(gproperty);
    free(gparam);

    free(gsrcpattern);
    free(greplayweight);
    free(greplayseed);
    free(greplaytime);
    free(mcxkernel);

    if (cfg->energytot) {
        free(cfg->energytot);
        cfg->energytot = NULL;
    }

    if (cfg->energyesc) {
        free(cfg->energyesc);
        cfg->energyesc = NULL;
    }

    if (gpu) {
        free(gpu);
    }

    for (devid = 0; devid < workdev; devid++) {
        if (mcxqueue && mcxqueue[devid]) {
            clFinish(mcxqueue[devid]);
            clReleaseCommandQueue(mcxqueue[devid]);
        }
    }

    free(mcxqueue);

    if (mcxprogram) {
        clReleaseProgram(mcxprogram);
    }

    if (mcxcontext) {
        clReleaseContext(mcxcontext);
    }

#ifndef USE_OS_TIMER

    /* free the last profiling event and NULL it so the next mmc_run_cl call does not
     * re-release a stale handle (kernelevent is a global reused across runs) */
    if (kernelevent) {
        clReleaseEvent(kernelevent);
        kernelevent = NULL;
    }

#endif
    free(field);
    free(Pdet);
    free(Pphotonseed);
    free(propdet);
    free(Pseed);
    free(energy);
    free(dref);

    if (status != CL_SUCCESS) {
        ocl_assess(status, __FILE__, __LINE__);   /* re-raise the deferred acquisition error after full cleanup */
    }
}
