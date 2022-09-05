#include <cstdlib>
#include <iostream>
#include <time.h>
#include <cstring>
#include <optix_function_table_definition.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "mmc_optix_host.h"
#include "mmc_tictoc.h"
#include "incbin.h"

INCTXT(mmcShaderPtx, mmcShaderPtxSize, "mmc_optix_core.ptx");

/************************************************************************** In
this unit, we first launch a master thread and initialize the necessary data
structures.This include the command line options(cfg), tetrahedral mesh(mesh)
and the ray tracer precomputed data (tracer).
******************************************************************************/
#define CUDA_ASSERT(a)                                                         \
    mcx_cu_assess((a), __FILE__, __LINE__) ///< macro to report CUDA error
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
    } else {
        return 128;
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

void mmc_run_optix(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, 
    void (*progressfun)(float, void*), void* handle) {
    GPUInfo* gpuinfo = NULL;      /** gpuinfo: structure to store GPU information */
    unsigned int activedev = 0;   /** activedev: count of total active GPUs to be used */

    if (!(activedev = mcx_list_cu_gpu(cfg, &gpuinfo))) {
        mcx_error(-1, "No GPU device found\n", __FILE__, __LINE__);
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
        optix_run_simulation(cfg, mesh, tracer, gpuinfo, progressfun, handle);

#ifdef _OPENMP
    }
#endif
}

void optix_run_simulation(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, GPUInfo* gpu, 
    void (*progressfun)(float, void*), void* handle) {
    uint tic0 = StartTimer();
    // ==================================================================
    // prepare optix pipeline
    // ==================================================================
    OptixParams optixcfg;

    initOptix();

    createContext(cfg, &optixcfg);
    MMC_FPRINTF(cfg->flog, "optix init complete:  \t%d ms\n", GetTimeMillis() - tic0);
    fflush(cfg->flog);
      
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
    MMC_FPRINTF(cfg->flog, "lauching OptiX for time window [%.1fns %.1fns] ...\n"
                    , cfg->tstart * 1e9, cfg->tend * 1e9);
    fflush(cfg->flog);
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
    MMC_FPRINTF(cfg->flog, "kernel complete:  \t%d ms\nretrieving flux ... \t", 
        GetTimeMillis() - tic0);
    fflush(cfg->flog);

    // ==================================================================
    // Save output: combine two outputs into one
    // ==================================================================
    optixcfg.outputBuffer.download(optixcfg.outputHostBuffer, optixcfg.outputBufferSize);
    MMC_FPRINTF(cfg->flog, "transfer complete:        %d ms\n", GetTimeMillis() - tic0);
    fflush(cfg->flog);
    for (size_t i = 0; i < optixcfg.launchParams.crop0.w; i++) {
        #pragma omp atomic
        mesh->weight[i] += optixcfg.outputHostBuffer[i] +
            optixcfg.outputHostBuffer[i + optixcfg.launchParams.crop0.w];
    }

    // ==================================================================
    // Free memory
    // ==================================================================
    clearOptixParams(&optixcfg);
}

/**
 * @brief prepare launch parameters
 */
void prepLaunchParams(mcconfig* cfg, tetmesh* mesh, GPUInfo* gpu, 
    OptixParams *optixcfg) {
    if (cfg->method != rtBLBadouelGrid) {
        mcx_error(-1, "Optix MMC only supports dual grid mode", __FILE__, __LINE__);
    }

    int timeSteps = (int)((cfg->tend - cfg->tstart) / cfg->tstep + 0.5);
    if (timeSteps < 1) {
        mcx_error(-1, "There must be at least one time step.", __FILE__, __LINE__);
    }

    // set up optical properties
    if (mesh->prop + 1 > MAX_PROP) {
        mcx_error(-1, "Medium type count exceeds limit.", __FILE__, __LINE__);
    }
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

    // parameters of dual grid
    optixcfg->launchParams.nmin = make_float3(mesh->nmin.x,
                                              mesh->nmin.y,
                                              mesh->nmin.z);
    optixcfg->launchParams.crop0 = make_uint4(cfg->crop0.x,
                                               cfg->crop0.y,
                                               cfg->crop0.z,
                                               cfg->crop0.z * timeSteps);
    optixcfg->launchParams.dstep = 1.0f / cfg->unitinmm;

    // time-gate settings
    optixcfg->launchParams.tstart = cfg->tstart;
    optixcfg->launchParams.tend = cfg->tend;
    optixcfg->launchParams.Rtstep = 1.0f / cfg->tstep;
    optixcfg->launchParams.maxgate = cfg->maxgate;

    // init medium ID using element based
    optixcfg->launchParams.mediumid0 = mesh->type[cfg->e0-1];

    // number of photons for each thread
    int totalthread = cfg->nthread;

    int gpuid, threadid = 0;
#ifdef _OPENMP
    threadid = omp_get_thread_num();
#endif
    gpuid = cfg->deviceid[threadid] - 1;
    
    if (cfg->autopilot)
        totalthread = gpu[gpuid].autothread;

    optixcfg->launchWidth = totalthread;
    optixcfg->launchParams.threadphoton = cfg->nphoton / totalthread;
    optixcfg->launchParams.oddphoton = 
        cfg->nphoton - optixcfg->launchParams.threadphoton * totalthread;

    // output buffer (single precision)
    optixcfg->outputBufferSize = (optixcfg->launchParams.crop0.w << 1);
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
    uint4 *hseed = (uint4 *)malloc(sizeof(uint4) * totalthread);
    for (int i = 0; i < totalthread; ++i) {
        hseed[i] = make_uint4(rand(), rand(), rand(), rand());
    }
    optixcfg->seedBuffer.alloc_and_upload(hseed, totalthread);
    optixcfg->launchParams.seedbuffer = optixcfg->seedBuffer.d_pointer();

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
    if(cuRes != CUDA_SUCCESS) 
        fprintf(stderr, "Error querying current context: error code %d\n", cuRes);
      
    OPTIX_CHECK(optixDeviceContextCreate(optixcfg->cudaContext, 0, &optixcfg->optixContext));  
}

/**
 * @brief creates the module that contains all programs
 */
void createModule(mcconfig* cfg, OptixParams* optixcfg, std::string ptxcode) {
    // moduleCompileOptions
    optixcfg->moduleCompileOptions.maxRegisterCount = OPTIX_COMPILE_DEFAULT_MAX_REGISTER_COUNT;
#ifndef NDEBUG
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
    optixcfg->pipelineCompileOptions.numPayloadValues   = 14;
    optixcfg->pipelineCompileOptions.numAttributeValues = 2;  // for triangle
#ifndef NDEBUG
    optixcfg->pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_DEBUG |
        OPTIX_EXCEPTION_FLAG_TRACE_DEPTH | OPTIX_EXCEPTION_FLAG_STACK_OVERFLOW;
#else
    optixcfg->pipelineCompileOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
#endif
    optixcfg->pipelineCompileOptions.pipelineLaunchParamsVariableName = "gcfg";
      
    // pipelineLinkOptions
    optixcfg->pipelineLinkOptions.maxTraceDepth = 1;

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
    if (logsize > 1) std::cout << log << std::endl;
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
    if (logsize > 1) std::cout << log << std::endl;
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
    if (logsize > 1) std::cout << log << std::endl;
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
    if (logsize > 1) std::cout << log << std::endl;
}

/**
 * @brief set up acceleration structures
 */
OptixTraversableHandle buildAccel(tetmesh* mesh, OptixParams* optixcfg) {
    // ==================================================================
    // upload the model to the device
    // note: mesh->fnode needs to be float3 and 
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
    
    uint32_t triangleInputFlags[1] = { 0 }; // OPTIX_GEOMETRY_FLAG_NONE?
    
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
                                
                                &emitDesc,1
                                ));
    CUDA_SYNC_CHECK();
    
    // ==================================================================
    // perform compaction
    // ==================================================================
    uint64_t compactedSize;
    compactedSizeBuffer.download(&compactedSize,1);
    
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
    for (auto pg : optixcfg->raygenPGs)
      programGroups.push_back(pg);
    for (auto pg : optixcfg->missPGs)
      programGroups.push_back(pg);
    for (auto pg : optixcfg->hitgroupPGs)
      programGroups.push_back(pg);

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
    if (logsize > 1) std::cout << log << std::endl;
}

/**
 * @ set up the shader binding table
 */
void buildSBT(tetmesh* mesh, OptixParams* optixcfg) {
    // ==================================================================
    // build raygen records
    // ==================================================================
    std::vector<RaygenRecord> raygenRecords;
    for (size_t i = 0;i < optixcfg->raygenPGs.size();i++) {
      RaygenRecord rec;
      OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->raygenPGs[i],&rec));
      rec.data = nullptr;
      raygenRecords.push_back(rec);
    }
    optixcfg->raygenRecordsBuffer.alloc_and_upload(raygenRecords);
    optixcfg->sbt.raygenRecord = optixcfg->raygenRecordsBuffer.d_pointer();

    // ==================================================================
    // build miss records
    // ==================================================================
    std::vector<MissRecord> missRecords;
    for (size_t i = 0;i < optixcfg->missPGs.size();i++) {
      MissRecord rec;
      OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->missPGs[i],&rec));
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
    OPTIX_CHECK(optixSbtRecordPackHeader(optixcfg->hitgroupPGs[0],&rec));
    rec.data.node = (float3*)optixcfg->vertexBuffer.d_pointer();

    // combine face + front + back into a uint4 array
    uint4 *face = (uint4*)calloc(mesh->nface, sizeof(uint4));
    for (size_t i = 0; i < mesh->nface; ++i) {
        face[i] = make_uint4(mesh->face[i].x, mesh->face[i].y, mesh->face[i].z,
            (mesh->front[i] << 16) | (0xFF & mesh->back[i]));
    }
    optixcfg->faceBuffer.alloc_and_upload(face, mesh->nface);
    rec.data.face = (uint4*)optixcfg->faceBuffer.d_pointer();
    hitgroupRecords.push_back(rec);

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
    free(optixcfg->outputHostBuffer);
}