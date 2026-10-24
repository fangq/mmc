#include <cstdlib>
#include <iostream>
#include <time.h>
#include <cstring>
#include <optix_function_table_definition.h>
#include <iomanip>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "mmc_optix_utils.h"
#include "mmc_cuda_query_gpu.h"
#include "mmc_tictoc.h"
#include "mmc_cu_host.h"
#include "incbin.h"

INCTXT(mmcShaderPtx, mmcShaderPtxSize, "built/mmc_optix_core.ptx");

void optix_run_simulation(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, GPUInfo* gpu) {
    uint tic0 = StartTimer();
    // ==================================================================
    // prepare optix pipeline
    // ==================================================================
    OptixParams optixcfg;

    initOptix();

    createContext(cfg, &optixcfg);
    MMC_FPRINTF(cfg->flog, "optix init complete:  \t%d ms\n", GetTimeMillis() - tic0);
    fflush(cfg->flog);

    prepBoundParams(cfg, mesh, &optixcfg);

    std::string ptxcodestr = std::string(mmcShaderPtx);
    createModule(cfg, &optixcfg, ptxcodestr);
    MMC_FPRINTF(cfg->flog, "optix module complete:  \t%d ms\n", GetTimeMillis() - tic0);
    fflush(cfg->flog);

    createRaygenPrograms(&optixcfg);
    createMissPrograms(&optixcfg);
    createHitgroupPrograms(&optixcfg);
    MMC_FPRINTF(cfg->flog, "optix device programs complete:  \t%d ms\n",
                GetTimeMillis() - tic0);
    fflush(cfg->flog);

    optixcfg.launchParams.gashandle = buildAccel(mesh, &optixcfg);
    MMC_FPRINTF(cfg->flog, "optix acceleration structure complete:  \t%d ms\n",
                GetTimeMillis() - tic0);
    fflush(cfg->flog);

    createPipeline(&optixcfg);
    MMC_FPRINTF(cfg->flog, "optix pipeline complete:  \t%d ms\n",
                GetTimeMillis() - tic0);
    fflush(cfg->flog);

    buildSBT(mesh, &optixcfg);
    MMC_FPRINTF(cfg->flog, "optix shader binding table complete:  \t%d ms\n",
                GetTimeMillis() - tic0);
    fflush(cfg->flog);

    // ==================================================================
    // prepare launch parameters
    // ==================================================================
    prepLaunchParams(cfg, mesh, gpu, &optixcfg);
    CUDA_ASSERT(cudaDeviceSynchronize());
    MMC_FPRINTF(cfg->flog, "optix launch parameters complete:  \t%d ms\n",
                GetTimeMillis() - tic0);
    fflush(cfg->flog);

    // ==================================================================
    // Launch simulation
    // ==================================================================
    MMC_FPRINTF(cfg->flog, "lauching OptiX for time window [%.1fns %.1fns] ...\n",
                cfg->tstart * 1e9, cfg->tend * 1e9);
    fflush(cfg->flog);
    uint kernel_start_time = GetTimeMillis();
    OPTIX_CHECK(optixLaunch(/*! pipeline we're launching launch: */
                    optixcfg.pipeline,
                    optixcfg.stream,
                    /*! parameters and SBT */
                    optixcfg.launchParamsBuffer.d_pointer(),
                    optixcfg.launchParamsBuffer.sizeInBytes,
                    &optixcfg.sbt,
                    /*! dimensions of the launch: */
                    optixcfg.launchWidth, 1, 1));
    CUDA_ASSERT(cudaDeviceSynchronize());
    uint kernel_runtime = GetTimeMillis() - kernel_start_time;
    MMC_FPRINTF(cfg->flog, "kernel complete:  \t%d ms\nretrieving flux ... \t",
                GetTimeMillis() - tic0);
    MMC_FPRINTF(cfg->flog, "simulated %zu photons (%zu) with 1 devices \nMMC simulation speed: %.2f photon/ms\n",
                cfg->nphoton, cfg->nphoton, (double)cfg->nphoton / kernel_runtime);
    fflush(cfg->flog);

    // ==================================================================
    // Save output
    // ==================================================================
    optixcfg.outputBuffer.download(optixcfg.outputHostBuffer, optixcfg.outputBufferSize);
    MMC_FPRINTF(cfg->flog, "transfer complete:        %d ms\n", GetTimeMillis() - tic0);
    fflush(cfg->flog);

    /* the launched and escaped energy sums of each pattern follow the output grids, see saveEnergy */
    float* energybuf = optixcfg.outputHostBuffer + optixcfg.outputBufferSize - 4 * cfg->srcnum;
    double energy[2] = {0.0, 0.0};

    if (cfg->energytot == NULL) {
        cfg->energytot = (double*)calloc(cfg->srcnum, sizeof(double));
        cfg->energyesc = (double*)calloc(cfg->srcnum, sizeof(double));
    }

    for (int j = 0; j < cfg->srcnum; j++) {
        double etot = (double)energybuf[j * 4] + energybuf[j * 4 + 2], eesc = (double)energybuf[j * 4 + 1] + energybuf[j * 4 + 3];
        cfg->energytot[j] += etot;
        cfg->energyesc[j] += eesc;
        energy[0] += etot;
        energy[1] += eesc;
    }

    // detected photons: the counter, then the records, as the OpenCL/CUDA hosts store them
    if (optixcfg.launchParams.issavedet) {
        size_t reclen = 2 * mesh->prop + 2 + optixcfg.launchParams.issaveexit * 6;
        unsigned int ndet = 0;
        CUDA_ASSERT(cudaMemcpy(&ndet, (void*)optixcfg.detBuffer.d_pointer(), sizeof(unsigned int), cudaMemcpyDeviceToHost));

        if (ndet > cfg->maxdetphoton) {
            MMC_FPRINTF(cfg->flog, S_RED "WARNING: the detected photon (%u) is more than what your have specified (%d), please use the -H option to specify a greater number\n" S_RESET, ndet, cfg->maxdetphoton);
            ndet = cfg->maxdetphoton;
        }

        float* hrec = (float*)malloc(sizeof(float) * (1 + reclen * ndet));
        CUDA_ASSERT(cudaMemcpy(hrec, (void*)optixcfg.detBuffer.d_pointer(), sizeof(float) * (1 + reclen * ndet), cudaMemcpyDeviceToHost));
        cfg->exportdetected = (float*)realloc(cfg->exportdetected, sizeof(float) * reclen * (cfg->detectedcount + ndet + 1));
        memcpy(cfg->exportdetected + reclen * cfg->detectedcount, hrec + 1, sizeof(float) * reclen * ndet);
        free(hrec);
        cfg->detectedcount += ndet;
        cfg->his.detected += ndet;
        MMC_FPRINTF(cfg->flog, "detected %u photons, total: %d\n", ndet, cfg->detectedcount);
    }

    MMC_FPRINTF(cfg->flog, "total simulated energy: %.2f\tabsorbed: %5.5f%%\n(loss due to initial specular reflection is excluded in the total)\n",
                energy[0], (energy[0] > 0.0) ? (energy[0] - energy[1]) / energy[0] * 100.f : 0.f);
    fflush(cfg->flog);

    int isrfforward_o = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE);
    size_t fieldlen_o = optixcfg.launchParams.crop0.w;

    for (size_t i = 0; i < fieldlen_o; i++) {
        // combine real double-buffer outputs
        #pragma omp atomic
        mesh->weight[i] += optixcfg.outputHostBuffer[i] +
                           optixcfg.outputHostBuffer[i + fieldlen_o];
    }

    /* RF forward: store imaginary part in exportadjoint */
    if (isrfforward_o) {
        if (cfg->exportadjoint == NULL) {
            cfg->exportadjoint = (float*)calloc(fieldlen_o, sizeof(float));
        }

        for (size_t i = 0; i < fieldlen_o; i++) {
            #pragma omp atomic
            cfg->exportadjoint[i] += optixcfg.outputHostBuffer[i + fieldlen_o * 2] +
                                     optixcfg.outputHostBuffer[i + fieldlen_o * 3];
        }
    }

    // ==================================================================
    // normalize output
    // ==================================================================
    if (cfg->isnormalized) {
        MMC_FPRINTF(cfg->flog, "normalizing raw data ...\t");
        fflush(cfg->flog);

        double sum_normalizer = 0.0;

        for (int j = 0; j < cfg->srcnum; j++) {
            double etot = cfg->energytot[j], eabs = cfg->energytot[j] - cfg->energyesc[j];
            double normalizer = mesh_normalize(mesh, cfg, eabs, etot, j);
            sum_normalizer += normalizer;
            MMCDEBUG(cfg, dlTime, (cfg->flog, "source %d\ttotal simulated energy: %f\tabsorbed: " S_BOLD "" S_BLUE "%5.5f%%" S_RESET "\tnormalizor=%g\n",
                                   j + 1, etot, (etot > 0.0) ? 100.f * eabs / etot : 0.f, normalizer));
        }

        cfg->his.normalizer = sum_normalizer / cfg->srcnum;

        /* multi-source mode: mesh_normalize scales the field of the first source slot, the
           other slots use the same factor, as the OpenCL/CUDA backends do */
        if (cfg->extrasrclen > cfg->srcnum && optixcfg.launchParams.srcid < 0) {
            size_t slotlen = (size_t)cfg->crop0.z * cfg->maxgate;

            for (size_t i = slotlen; i < fieldlen_o; i++) {
                mesh->weight[i] *= cfg->his.normalizer;

                if (isrfforward_o && cfg->exportadjoint) {
                    cfg->exportadjoint[i] *= cfg->his.normalizer;
                }
            }
        }
        MMC_FPRINTF(cfg->flog, "normalization complete:    %d ms\n",
                    GetTimeMillis() - tic0);
        fflush(cfg->flog);
    }

    // adjoint Jacobians from the fields of the source and detector slots, by the CUDA kernels
    mmc_cu_adjoint_grid(cfg, mesh, mesh->weight, fieldlen_o, isrfforward_o,
                        optixcfg.launchParams.savevoxtype ? (int*)optixcfg.voxtypeBuffer.d_pointer() : NULL, tic0);

#ifndef MCX_CONTAINER
    #pragma omp master
    {
        if (cfg->issave2pt && cfg->parentid == mpStandalone) {
            MMC_FPRINTF(cfg->flog, "saving data to file ...\t");
            mesh_saveweight(mesh, cfg, 0);
            MMC_FPRINTF(cfg->flog, "saving data complete : %d ms\n\n",
                        GetTimeMillis() - tic0);
            fflush(cfg->flog);
        }

        if (optixcfg.launchParams.issavedet && cfg->parentid == mpStandalone && cfg->exportdetected) {
            MMC_FPRINTF(cfg->flog, "saving detected photon data to file ...\t");
            cfg->his.totalphoton = cfg->nphoton;
            cfg->his.unitinmm = cfg->unitinmm;
            cfg->his.savedphoton = cfg->detectedcount;
            cfg->his.detected = cfg->detectedcount;
            cfg->his.colcount = 2 * cfg->his.maxmedia + (cfg->issaveexit > 0) * 6 + 2;
            cfg->his.seedbyte = 0;
            mcx_savedetphoton(cfg->exportdetected, NULL, cfg->detectedcount, 0, cfg);
            MMC_FPRINTF(cfg->flog, "saving detected photon data complete : %d ms\n\n", GetTimeMillis() - tic0);
            fflush(cfg->flog);
        }
    }
#endif

    // ==================================================================
    // Free memory
    // ==================================================================
    #pragma omp master
    {
        if (gpu) {
            free(gpu);
        }
    }

    clearOptixParams(&optixcfg);
}

/**
 * @brief reject the unsupported settings and set the launch parameters that stay the
 * same during a run, before the module is compiled
 */
void prepBoundParams(mcconfig* cfg, tetmesh* mesh, OptixParams* optixcfg) {
    if (cfg->method != rtBLBadouelGrid) {
        mcx_error(-1, "Optix MMC only supports dual grid mode", __FILE__, __LINE__);
    }

    /* reject the features that the OptiX kernel does not implement, rather than
       returning silently wrong outputs */
    if (cfg->seed == SEED_FROM_FILE || cfg->outputtype == otJacobian || cfg->outputtype == otWL ||
            cfg->outputtype == otWP || cfg->outputtype == otRFmus) {
        mcx_error(-1, "Optix MMC does not support photon replay or the replay outputs (-O J/L/P); use --compute cuda or opencl", __FILE__, __LINE__);
    }

    if (cfg->extrasrclen > 0 && cfg->srcnum > 1) {
        mcx_error(-1, "Optix MMC does not support multiple pattern sources sharing photons (SrcNum > 1); use --compute cuda or opencl", __FILE__, __LINE__);
    }

    if ((cfg->srctype == stPattern || cfg->srctype == stPattern3D) && cfg->srcpattern == NULL) {
        mcx_error(-1, "the pattern source has no pattern data", __FILE__, __LINE__);
    }

    if (cfg->issavedet && (cfg->issaveseed || cfg->ismomentum || cfg->isextdet || cfg->issaveexit > 1)) {
        mcx_error(-1, "Optix MMC saves detected photons without seeds (-q 1), momentum transfer, wide-field detectors or detector images (-x 2); use --compute cuda or opencl", __FILE__, __LINE__);
    }

    if (cfg->issaveref) {
        mcx_error(-1, "Optix MMC does not save the surface reflectance (-X 1); use --compute cuda or opencl", __FILE__, __LINE__);
    }

    if (cfg->debuglevel & dlTraj) {
        mcx_error(-1, "Optix MMC does not save photon trajectories (-D M); use --compute cuda or opencl", __FILE__, __LINE__);
    }

    if (cfg->mcmethod == mmMCML) {
        mcx_error(-1, "Optix MMC only supports the MCX (micro-Beer-Lambert) method, not -m 1", __FILE__, __LINE__);
    }

    if (cfg->isnodalmua || cfg->isnodalmusp) {
        mcx_error(-1, "Optix MMC does not support per-node optical properties (nodemua/nodemusp)", __FILE__, __LINE__);
    }

    int timeSteps = (int)((cfg->tend - cfg->tstart) / cfg->tstep + 0.5);

    if (timeSteps < 1) {
        mcx_error(-1, "There must be at least one time step.", __FILE__, __LINE__);
    }

    // set up optical properties
    if (mesh->prop + 1 > MAX_PROP) {
        mcx_error(-1, "Medium type count exceeds limit.", __FILE__, __LINE__);
    }

    /* the settings that stay the same during a run, which createModule binds as
       constants so that the compiler removes the code of the unused features */
    memset(&optixcfg->launchParams, 0, sizeof(optixcfg->launchParams));
    optixcfg->launchParams.srctype = cfg->srctype;
    optixcfg->launchParams.focus = cfg->srcdir.w;
    optixcfg->launchParams.srcnum = cfg->srcnum;
    optixcfg->launchParams.isreflect = cfg->isreflect;
    optixcfg->launchParams.outputtype = static_cast<int>(cfg->outputtype);
    optixcfg->launchParams.omega = cfg->omega;
    optixcfg->launchParams.srcid = (cfg->extrasrclen > 0 && cfg->srcid == 0) ? -1 : cfg->srcid;  /* < 0: a field per source slot, > 0: one source slot */
    optixcfg->launchParams.extrasrclen = cfg->extrasrclen;
    optixcfg->launchParams.issavedet = (cfg->issavedet && cfg->detnum > 0);
    optixcfg->launchParams.issaveexit = (cfg->issaveexit > 0);
    optixcfg->launchParams.detnum = cfg->detnum;
    optixcfg->launchParams.maxmedia = mesh->prop;
    optixcfg->launchParams.savevoxtype = (cfg->issave2pt && (cfg->outputtype == otAdjointMus || cfg->outputtype == otAdjointMusp
                                          || cfg->outputtype == otAdjointMuaMusp));

    if (cfg->e0 == -1) {
        optixcfg->launchParams.mediumid0 = INITIAL_MEDIUM_UNKNOWN;
    } else {
        optixcfg->launchParams.mediumid0 = MAX(mesh->type[cfg->e0 - 1], 0);  /* labels <= 0: the ambient medium */
    }
}

/**
 * @brief prepare launch parameters
 */
void prepLaunchParams(mcconfig* cfg, tetmesh* mesh, GPUInfo* gpu,
                      OptixParams* optixcfg) {

    for (int i = 0; i <= mesh->prop; ++i) {
        optixcfg->launchParams.medium[i].mua = mesh->med[i].mua;
        optixcfg->launchParams.medium[i].mus = mesh->med[i].mus;
        optixcfg->launchParams.medium[i].g = mesh->med[i].g;
        optixcfg->launchParams.medium[i].n = mesh->med[i].n;
    }

    // source setup
    optixcfg->launchParams.srcpos = make_float3(cfg->srcpos.x,
                                    cfg->srcpos.y,
                                    cfg->srcpos.z);
    optixcfg->launchParams.srcdir = make_float3(cfg->srcdir.x,
                                    cfg->srcdir.y,
                                    cfg->srcdir.z);
    optixcfg->launchParams.srcparam1 = make_float4(cfg->srcparam1.x,
                                       cfg->srcparam1.y,
                                       cfg->srcparam1.z,
                                       cfg->srcparam1.w);
    optixcfg->launchParams.srcparam2 = make_float4(cfg->srcparam2.x,
                                       cfg->srcparam2.y,
                                       cfg->srcparam2.z,
                                       cfg->srcparam2.w);

    // parameters of dual grid
    optixcfg->launchParams.nmin = make_float3(mesh->nmin.x,
                                  mesh->nmin.y,
                                  mesh->nmin.z);
    optixcfg->launchParams.nmax = make_float3(mesh->nmax.x - mesh->nmin.x,
                                  mesh->nmax.y - mesh->nmin.y,
                                  mesh->nmax.z - mesh->nmin.z);
    /* the output holds srcnum interleaved fields for a pattern source shared by srcnum patterns,
       or one field per source slot in the multi-source mode */
    int nslots = (cfg->extrasrclen > 0 && optixcfg->launchParams.srcid < 0) ? cfg->extrasrclen : 1;
    optixcfg->launchParams.crop0 = make_uint4(cfg->crop0.x,
                                   cfg->crop0.y,
                                   cfg->crop0.z,
                                   cfg->crop0.z * cfg->maxgate * cfg->srcnum * nslots);
    optixcfg->launchParams.dstep = 1.0f / cfg->steps.x;  /* the dual grid voxel size (--gridsize), not unitinmm */

    // time-gate settings
    optixcfg->launchParams.tstart = cfg->tstart;
    optixcfg->launchParams.tend = cfg->tend;
    optixcfg->launchParams.Rtstep = 1.0f / cfg->tstep;
    optixcfg->launchParams.maxgate = cfg->maxgate;

    optixcfg->launchParams.oneoverc0 = (float)(3.335640951981520e-12);

    // number of photons for each thread
    int totalthread = cfg->nthread;

    int gpuid, threadid = 0;
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    gpuid = cfg->deviceid[threadid] - 1;

    if (cfg->autopilot) {
        totalthread = gpu[gpuid].autothread;
    }

    optixcfg->launchWidth = totalthread;
    optixcfg->launchParams.threadphoton = cfg->nphoton / totalthread;
    optixcfg->launchParams.oddphoton =
        cfg->nphoton - optixcfg->launchParams.threadphoton * totalthread;

    // output buffer (single precision); 4x for RF (re1, re2, im1, im2)
    int isrfforward = (cfg->omega > 0.f && cfg->seed != SEED_FROM_FILE);
    optixcfg->outputBufferSize = optixcfg->launchParams.crop0.w * (isrfforward ? 4 : 2) + 4 * cfg->srcnum;  /* + the energy sums */
    optixcfg->outputHostBuffer = (float*)calloc(optixcfg->outputBufferSize, sizeof(float));
    optixcfg->outputBuffer.alloc_and_upload(optixcfg->outputHostBuffer,
                                            optixcfg->outputBufferSize);
    optixcfg->launchParams.outputbuffer = optixcfg->outputBuffer.d_pointer();

    // photon seed buffer
    if (cfg->seed > 0) {
        srand(cfg->seed);
    } else {
        srand(time(0));
    }

    uint4* hseed = (uint4*)malloc(sizeof(uint4) * totalthread);

    for (int i = 0; i < totalthread; ++i) {
        // separate statements: the evaluation order of function arguments is unspecified
        hseed[i].x = mmc_rand_uint32();
        hseed[i].y = mmc_rand_uint32();
        hseed[i].z = mmc_rand_uint32();
        hseed[i].w = mmc_rand_uint32();
    }

    optixcfg->seedBuffer.alloc_and_upload(hseed, totalthread);
    optixcfg->launchParams.seedbuffer = optixcfg->seedBuffer.d_pointer();

    // detected photons: detectors, a counter plus the records, and the per-thread scratch
    optixcfg->launchParams.maxdetphoton = cfg->maxdetphoton;

    if (optixcfg->launchParams.issavedet) {
        size_t reclen = 2 * mesh->prop + 2 + optixcfg->launchParams.issaveexit * 6;
        float4* hdet = (float4*)malloc(sizeof(float4) * cfg->detnum);

        for (uint i = 0; i < cfg->detnum; i++) {
            hdet[i] = make_float4(cfg->detpos[i].x, cfg->detpos[i].y, cfg->detpos[i].z, cfg->detpos[i].w);
        }

        optixcfg->detposBuffer.alloc_and_upload(hdet, cfg->detnum);
        free(hdet);

        float* hrec = (float*)calloc(1 + reclen * cfg->maxdetphoton, sizeof(float));
        optixcfg->detBuffer.alloc_and_upload(hrec, 1 + reclen * cfg->maxdetphoton);
        free(hrec);

        optixcfg->ppathBuffer.alloc(sizeof(float) * (2 * mesh->prop + 1) * totalthread);

        optixcfg->launchParams.detposbuffer = optixcfg->detposBuffer.d_pointer();
        optixcfg->launchParams.detbuffer = optixcfg->detBuffer.d_pointer();
        optixcfg->launchParams.ppathbuffer = optixcfg->ppathBuffer.d_pointer();
    }

    if (hseed) {
        free(hseed);
    }

    // pattern source: the intensities, srcnum per pixel, and for srcnum > 1 patterns sharing the
    // photons, the per-thread launch pixel and per-pattern energy sums
    optixcfg->launchParams.srcpattern = 0;
    optixcfg->launchParams.slotbuffer = 0;
    optixcfg->launchParams.energybuffer = 0;

    if ((cfg->srctype == stPattern || cfg->srctype == stPattern3D) && cfg->srcpattern) {
        size_t patlen = (cfg->srctype == stPattern) ? (size_t)(cfg->srcparam1.w * cfg->srcparam2.w)
                        : (size_t)(cfg->srcparam1.x * cfg->srcparam1.y * cfg->srcparam1.z);
        optixcfg->srcpatternBuffer.alloc_and_upload(cfg->srcpattern, patlen * cfg->srcnum);
        optixcfg->launchParams.srcpattern = optixcfg->srcpatternBuffer.d_pointer();
    }

    if (optixcfg->launchParams.savevoxtype) {
        optixcfg->voxtypeBuffer.alloc(sizeof(int) * cfg->crop0.z);
        CUDA_ASSERT(cudaMemset((void*)optixcfg->voxtypeBuffer.d_pointer(), 0, sizeof(int) * cfg->crop0.z));
        optixcfg->launchParams.voxtypebuffer = optixcfg->voxtypeBuffer.d_pointer();
    }

    if (cfg->srcnum > 1 || nslots > 1) {
        optixcfg->slotBuffer.alloc(sizeof(uint) * totalthread);
        optixcfg->launchParams.slotbuffer = optixcfg->slotBuffer.d_pointer();
    }

    if (cfg->srcnum > 1) {
        float* hen = (float*)calloc((size_t)totalthread * 4 * cfg->srcnum, sizeof(float));
        optixcfg->energyBuffer.alloc_and_upload(hen, (size_t)totalthread * 4 * cfg->srcnum);
        free(hen);
        optixcfg->launchParams.energybuffer = optixcfg->energyBuffer.d_pointer();
    }

    // upload extra source list if in multi-source mode
    if (cfg->extrasrclen > 0 && cfg->srcdata != NULL) {
        /* the kernel has no elements: replace the enclosing element of each source, found by
           mesh_init_srcdata_eid, by its medium, -1 to find the medium by ray tracing */
        ExtraSrc* hsrc = (ExtraSrc*)malloc(sizeof(ExtraSrc) * cfg->extrasrclen);
        memcpy(hsrc, cfg->srcdata, sizeof(ExtraSrc) * cfg->extrasrclen);

        for (int i = 0; i < cfg->extrasrclen; i++) {
            int eid = (int)hsrc[i].srcparam2.w;
            hsrc[i].srcparam2.w = (eid > 0 && eid <= (int)mesh->ne) ? (float)MAX(mesh->type[eid - 1], 0) : -1.f;
        }

        optixcfg->srcdataBuffer.alloc_and_upload(hsrc, cfg->extrasrclen);
        free(hsrc);
        optixcfg->launchParams.srcdatabuffer = optixcfg->srcdataBuffer.d_pointer();
    } else {
        optixcfg->launchParams.srcdatabuffer = 0;
    }

    // upload launch parameters to device
    optixcfg->launchParamsBuffer.alloc_and_upload(&optixcfg->launchParams, 1);
}

/**************************************************************************
 * helper functions for Optix pipeline creation
******************************************************************************/

/**
 * @brief initialize optix
 */
void initOptix() {
    cudaFree(0);
    OPTIX_CHECK(optixInit());
}

/**
 * @brief creates and configures a optix device context
 */
void createContext(mcconfig* cfg, OptixParams* optixcfg) {
    int gpuid, threadid = 0;

#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif

    gpuid = cfg->deviceid[threadid] - 1;

    if (gpuid < 0) {
        mcx_error(-1, "GPU ID must be non-zero", __FILE__, __LINE__);
    }

    CUDA_ASSERT(cudaSetDevice(gpuid));
    CUDA_ASSERT(cudaStreamCreate(&optixcfg->stream));

    cudaGetDeviceProperties(&optixcfg->deviceProps, gpuid);
    std::cout << "Running on device: " << optixcfg->deviceProps.name << std::endl;

    CUresult cuRes = cuCtxGetCurrent(&optixcfg->cudaContext);

    if (cuRes != CUDA_SUCCESS) {
        fprintf(stderr, "Error querying current context: error code %d\n", cuRes);
    }

    OptixDeviceContextOptions options = {};
    options.logCallbackFunction = [](unsigned int level, const char* tag, const char* message,
    void*) {
        std::cerr << "[" << std::setw( 2 ) << level
                  << "][" << std::setw( 12 ) << tag << "]: "
                  << message << "\n";
    };
#ifdef MMC_OPTIX_DEBUG  /* full OptiX debug pipeline, opt-in: also fails on some driver/SDK combinations */
    options.logCallbackLevel = 4;
    options.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;
#else
    options.logCallbackLevel = 0;
    options.validationMode = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_OFF;
#endif

    OPTIX_CHECK(optixDeviceContextCreate(optixcfg->cudaContext, &options,
                                         &optixcfg->optixContext));
}

/**
 * @brief creates the module that contains all programs
 */
void createModule(mcconfig* cfg, OptixParams* optixcfg, std::string ptxcode) {
    // moduleCompileOptions
    optixcfg->moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
#ifdef MMC_OPTIX_DEBUG
    optixcfg->moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_FULL;
    optixcfg->moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_0;
#else
    optixcfg->moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
    optixcfg->moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_3;
#endif

    // pipelineCompileOptions
    optixcfg->pipelineCompileOptions = {};
    optixcfg->pipelineCompileOptions.traversableGraphFlags =
        OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS;
    optixcfg->pipelineCompileOptions.usesMotionBlur     = false;
    optixcfg->pipelineCompileOptions.numPayloadValues   = 15; /* +1 for RF imaginary weight */
    optixcfg->pipelineCompileOptions.numAttributeValues = 2;  // for triangle
#ifdef MMC_OPTIX_DEBUG
    optixcfg->pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_DEBUG |
            OPTIX_EXCEPTION_FLAG_TRACE_DEPTH | OPTIX_EXCEPTION_FLAG_STACK_OVERFLOW;
#else
    optixcfg->pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
#endif
    optixcfg->pipelineCompileOptions.pipelineLaunchParamsVariableName = "gcfg";

    // pipelineLinkOptions
    optixcfg->pipelineLinkOptions.maxTraceDepth = 1;

    /* bind the launch parameters that stay the same during a run (set by prepBoundParams)
       as constants, so that the compiler removes the branches of the unused source types
       and features, as the OpenCL backend does with its JIT macros */
#define MMC_BIND(field) {offsetof(MMCParam, field), sizeof(optixcfg->launchParams.field), &optixcfg->launchParams.field, #field}
    OptixModuleCompileBoundValueEntry bound[] = {
        MMC_BIND(srctype), MMC_BIND(focus), MMC_BIND(srcnum), MMC_BIND(isreflect), MMC_BIND(outputtype),
        MMC_BIND(omega), MMC_BIND(srcid), MMC_BIND(extrasrclen), MMC_BIND(issavedet), MMC_BIND(issaveexit),
        MMC_BIND(detnum), MMC_BIND(maxmedia), MMC_BIND(mediumid0), MMC_BIND(savevoxtype)
    };
#undef MMC_BIND
    optixcfg->moduleCompileOptions.boundValues = bound;
    optixcfg->moduleCompileOptions.numBoundValues = sizeof(bound) / sizeof(bound[0]);

    char log[2048];
    size_t logsize = sizeof(log);
    OPTIX_CHECK(optixModuleCreateFromPTX(optixcfg->optixContext,
                                         &optixcfg->moduleCompileOptions,
                                         &optixcfg->pipelineCompileOptions,
                                         ptxcode.c_str(),
                                         ptxcode.size(),
                                         log,
                                         &logsize,
                                         &optixcfg->module
                                        ));

    if (logsize > 1) {
        std::cout << log << std::endl;
    }
}

/**
 * @brief set up ray generation programs
 */
void createRaygenPrograms(OptixParams* optixcfg) {
    optixcfg->raygenPGs.resize(1);

    OptixProgramGroupOptions pgOptions = {};
    OptixProgramGroupDesc pgDesc    = {};
    pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    pgDesc.raygen.module            = optixcfg->module;
    pgDesc.raygen.entryFunctionName = "__raygen__rg";

    char log[2048];
    size_t logsize = sizeof(log);
    OPTIX_CHECK(optixProgramGroupCreate(optixcfg->optixContext,
                                        &pgDesc,
                                        1,
                                        &pgOptions,
                                        log,
                                        &logsize,
                                        &optixcfg->raygenPGs[0]
                                       ));

    if (logsize > 1) {
        std::cout << log << std::endl;
    }
}

/**
 * @brief set up miss programs
 */
void createMissPrograms(OptixParams* optixcfg) {
    optixcfg->missPGs.resize(1);

    OptixProgramGroupOptions pgOptions = {};
    OptixProgramGroupDesc pgDesc    = {};
    pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_MISS;
    pgDesc.raygen.module            = optixcfg->module;
    pgDesc.raygen.entryFunctionName = "__miss__ms";

    char log[2048];
    size_t logsize = sizeof(log);
    OPTIX_CHECK(optixProgramGroupCreate(optixcfg->optixContext,
                                        &pgDesc,
                                        1,
                                        &pgOptions,
                                        log,
                                        &logsize,
                                        &optixcfg->missPGs[0]
                                       ));

    if (logsize > 1) {
        std::cout << log << std::endl;
    }
}

/**
 * @brief set up hitgroup programs
 */
void createHitgroupPrograms(OptixParams* optixcfg) {
    optixcfg->hitgroupPGs.resize(1);

    OptixProgramGroupOptions pgOptions = {};
    OptixProgramGroupDesc pgDesc    = {};
    pgDesc.kind                     = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    pgDesc.hitgroup.moduleCH            = optixcfg->module;
    pgDesc.hitgroup.entryFunctionNameCH = "__closesthit__ch";

    char log[2048];
    size_t logsize = sizeof(log);
    OPTIX_CHECK(optixProgramGroupCreate(optixcfg->optixContext,
                                        &pgDesc,
                                        1,
                                        &pgOptions,
                                        log,
                                        &logsize,
                                        &optixcfg->hitgroupPGs[0]
                                       ));

    if (logsize > 1) {
        std::cout << log << std::endl;
    }
}

/**
 * @brief set up acceleration structures
 */
OptixTraversableHandle buildAccel(tetmesh* mesh, OptixParams* optixcfg) {
    // ==================================================================
    // upload the model to the device
    // note: mesh->fnode needs to be float3
    // mesh->face needs to be uint3 (zero-indexed)
    // ==================================================================
    optixcfg->vertexBuffer.alloc_and_upload(mesh->fnode, mesh->nn);
    optixcfg->indexBuffer.alloc_and_upload(mesh->face, mesh->nface);

    OptixTraversableHandle asHandle {0};

    // ==================================================================
    // triangle inputs
    // ==================================================================
    OptixBuildInput triangleInput = {};
    triangleInput.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;

    // create local variables, because we need a *pointer* to the
    // device pointers
    CUdeviceptr d_vertices = optixcfg->vertexBuffer.d_pointer();
    CUdeviceptr d_indices  = optixcfg->indexBuffer.d_pointer();

    triangleInput.triangleArray.vertexFormat        = OPTIX_VERTEX_FORMAT_FLOAT3;
    triangleInput.triangleArray.vertexStrideInBytes = sizeof(float3);
    triangleInput.triangleArray.numVertices         = mesh->nn;
    triangleInput.triangleArray.vertexBuffers       = &d_vertices;

    triangleInput.triangleArray.indexFormat         = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    triangleInput.triangleArray.indexStrideInBytes  = sizeof(uint3);
    triangleInput.triangleArray.numIndexTriplets    = mesh->nface;
    triangleInput.triangleArray.indexBuffer         = d_indices;

    uint32_t triangleInputFlags[1] = { OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT };

    // in this example we have one SBT entry, and no per-primitive
    // materials:
    triangleInput.triangleArray.flags               = triangleInputFlags;
    triangleInput.triangleArray.numSbtRecords               = 1;
    triangleInput.triangleArray.sbtIndexOffsetBuffer        = 0;
    triangleInput.triangleArray.sbtIndexOffsetSizeInBytes   = 0;
    triangleInput.triangleArray.sbtIndexOffsetStrideInBytes = 0;

    // ==================================================================
    // BLAS setup
    // ==================================================================
    OptixAccelBuildOptions accelOptions = {};
    accelOptions.buildFlags             = OPTIX_BUILD_FLAG_NONE
                                          | OPTIX_BUILD_FLAG_ALLOW_COMPACTION
                                          ;
    accelOptions.motionOptions.numKeys  = 1;
    accelOptions.operation              = OPTIX_BUILD_OPERATION_BUILD;

    OptixAccelBufferSizes blasBufferSizes;
    OPTIX_CHECK(optixAccelComputeMemoryUsage
                (optixcfg->optixContext,
                 &accelOptions,
                 &triangleInput,
                 1,  // num_build_inputs
                 &blasBufferSizes
                ));

    // ==================================================================
    // prepare compaction
    // ==================================================================
    osc::CUDABuffer compactedSizeBuffer;
    compactedSizeBuffer.alloc(sizeof(uint64_t));

    OptixAccelEmitDesc emitDesc;
    emitDesc.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    emitDesc.result = compactedSizeBuffer.d_pointer();

    // ==================================================================
    // execute build (main stage)
    // ==================================================================
    osc::CUDABuffer tempBuffer;
    tempBuffer.alloc(blasBufferSizes.tempSizeInBytes);

    osc::CUDABuffer outputBuffer;
    outputBuffer.alloc(blasBufferSizes.outputSizeInBytes);

    OPTIX_CHECK(optixAccelBuild(optixcfg->optixContext,
                                optixcfg->stream,
                                &accelOptions,
                                &triangleInput,
                                1,
                                tempBuffer.d_pointer(),
                                tempBuffer.sizeInBytes,
                                outputBuffer.d_pointer(),
                                outputBuffer.sizeInBytes,
                                &asHandle,
                                &emitDesc, 1
                               ));
    CUDA_SYNC_CHECK();

    // ==================================================================
    // perform compaction
    // ==================================================================
    uint64_t compactedSize;
    compactedSizeBuffer.download(&compactedSize, 1);

    optixcfg->asBuffer.alloc(compactedSize);
    OPTIX_CHECK(optixAccelCompact(optixcfg->optixContext,
                                  optixcfg->stream,
                                  asHandle,
                                  optixcfg->asBuffer.d_pointer(),
                                  optixcfg->asBuffer.sizeInBytes,
                                  &asHandle));
    CUDA_SYNC_CHECK();

    // ==================================================================
    // clean up
    // ==================================================================
    outputBuffer.free(); // << the UNcompacted, temporary output buffer
    tempBuffer.free();
    compactedSizeBuffer.free();

    return asHandle;
}

/**
 * @brief assemble the pipeline of all programs
 */
void createPipeline(OptixParams* optixcfg) {
    std::vector<OptixProgramGroup> programGroups;

    for (auto pg : optixcfg->raygenPGs) {
        programGroups.push_back(pg);
    }

    for (auto pg : optixcfg->missPGs) {
        programGroups.push_back(pg);
    }

    for (auto pg : optixcfg->hitgroupPGs) {
        programGroups.push_back(pg);
    }

    char log[2048];
    size_t logsize = sizeof(log);
    OPTIX_CHECK(optixPipelineCreate(optixcfg->optixContext,
                                    &optixcfg->pipelineCompileOptions,
                                    &optixcfg->pipelineLinkOptions,
                                    programGroups.data(),
                                    (int)programGroups.size(),
                                    log,
                                    &logsize,
                                    &optixcfg->pipeline
                                   ));

    if (logsize > 1) {
        std::cout << log << std::endl;
    }
}

/**
 * @ set up the shader binding table
 */
void buildSBT(tetmesh* mesh, OptixParams* optixcfg) {
    // ==================================================================
    // build raygen records
    // ==================================================================
    std::vector<RaygenRecord> raygenRecords;

    for (size_t i = 0; i < optixcfg->raygenPGs.size(); i++) {
        RaygenRecord rec;
        OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->raygenPGs[i], &rec));
        rec.data = nullptr;
        raygenRecords.push_back(rec);
    }

    optixcfg->raygenRecordsBuffer.alloc_and_upload(raygenRecords);
    optixcfg->sbt.raygenRecord = optixcfg->raygenRecordsBuffer.d_pointer();

    // ==================================================================
    // build miss records
    // ==================================================================
    std::vector<MissRecord> missRecords;

    for (size_t i = 0; i < optixcfg->missPGs.size(); i++) {
        MissRecord rec;
        OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->missPGs[i], &rec));
        rec.data = nullptr; /* for now ... */
        missRecords.push_back(rec);
    }

    optixcfg->missRecordsBuffer.alloc_and_upload(missRecords);
    optixcfg->sbt.missRecordBase          = optixcfg->missRecordsBuffer.d_pointer();
    optixcfg->sbt.missRecordStrideInBytes = sizeof(MissRecord);
    optixcfg->sbt.missRecordCount         = (int)missRecords.size();

    // ==================================================================
    // build hitgroup records
    // ==================================================================
    std::vector<HitgroupRecord> hitgroupRecords;
    HitgroupRecord rec;
    // all meshes use the same code, so all same hit group
    OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->hitgroupPGs[0], &rec));

    // combine face normal + front + back into a float4 array
    float4* fnorm = (float4*)calloc(mesh->nface, sizeof(float4));

    for (size_t i = 0; i < mesh->nface; ++i) {
        uint media = (mesh->front[i] << 16) | (0xFFFF & mesh->back[i]);  /* the kernel reads 16 bits each */
        fnorm[i] = make_float4(mesh->fnorm[i].x, mesh->fnorm[i].y, mesh->fnorm[i].z,
                               *(float*)&media);
    }

    optixcfg->faceBuffer.alloc_and_upload(fnorm, mesh->nface);
    rec.data.fnorm = (float4*)optixcfg->faceBuffer.d_pointer();
    hitgroupRecords.push_back(rec);
    free(fnorm);

    optixcfg->hitgroupRecordsBuffer.alloc_and_upload(hitgroupRecords);
    optixcfg->sbt.hitgroupRecordBase          = optixcfg->hitgroupRecordsBuffer.d_pointer();
    optixcfg->sbt.hitgroupRecordStrideInBytes = sizeof(HitgroupRecord);
    optixcfg->sbt.hitgroupRecordCount         = (int)hitgroupRecords.size();
}

/**
 * @ Free allocated memory
 */
void clearOptixParams(OptixParams* optixcfg) {
    optixcfg->raygenRecordsBuffer.free();
    optixcfg->missRecordsBuffer.free();
    optixcfg->hitgroupRecordsBuffer.free();
    optixcfg->launchParamsBuffer.free();
    optixcfg->vertexBuffer.free();
    optixcfg->indexBuffer.free();
    optixcfg->faceBuffer.free();
    optixcfg->asBuffer.free();
    optixcfg->seedBuffer.free();
    optixcfg->outputBuffer.free();
    optixcfg->detposBuffer.free();
    optixcfg->detBuffer.free();
    optixcfg->ppathBuffer.free();
    optixcfg->srcdataBuffer.free();
    optixcfg->voxtypeBuffer.free();
    optixcfg->srcpatternBuffer.free();
    optixcfg->slotBuffer.free();
    optixcfg->energyBuffer.free();
    free(optixcfg->outputHostBuffer);

    optixPipelineDestroy(optixcfg->pipeline);
    optixProgramGroupDestroy(optixcfg->raygenPGs[0]);
    optixProgramGroupDestroy(optixcfg->missPGs[0]);
    optixProgramGroupDestroy(optixcfg->hitgroupPGs[0]);
    optixModuleDestroy(optixcfg->module);
    optixDeviceContextDestroy(optixcfg->optixContext);
}
