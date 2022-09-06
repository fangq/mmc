#ifndef _MMC_CUDA_QUERY_GPU_H
#define _MMC_CUDA_QUERY_GPU_H

#include "mmc_utils.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#define CUDA_ASSERT(a)                                                         \
    mcx_cu_assess((a), __FILE__, __LINE__) ///< macro to report CUDA error

void mcx_cu_assess(cudaError_t cuerr, const char* file, const int linenum);
int mcx_corecount(int v1, int v2);
int mcx_smxblock(int v1, int v2);
int mcx_list_cu_gpu(mcconfig* cfg, GPUInfo** info);

#endif