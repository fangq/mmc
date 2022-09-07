#include <cstring>
#include "mmc_cuda_query_gpu.h"

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