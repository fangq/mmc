#ifndef _MCEXTREME_GPU_LAUNCH_H
#define _MCEXTREME_GPU_LAUNCH_H

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#else
  #include <CL/cl.h>
#endif

#include "mcx_utils.h"

#ifdef  __cplusplus
extern "C" {
#endif

#define ABS(a)  ((a)<0?-(a):(a))

#define MCX_DEBUG_RNG       2                   /**< MCX debug flags */
#define MCX_DEBUG_MOVE      1
#define MCX_DEBUG_PROGRESS  2048

#define MIN(a,b)           ((a)<(b)?(a):(b))

#define OCL_ASSERT(x)  ocl_assess((x),__FILE__,__LINE__)

#define CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV           0x4000
#define CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV           0x4001
#define CL_DEVICE_REGISTERS_PER_BLOCK_NV                0x4002
#define CL_DEVICE_WARP_SIZE_NV                          0x4003
#define CL_DEVICE_GPU_OVERLAP_NV                        0x4004
#define CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV                0x4005
#define CL_DEVICE_INTEGRATED_MEMORY_NV                  0x4006

#define CL_DEVICE_BOARD_NAME_AMD                        0x4038
#define CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD             0x4040
#define CL_DEVICE_WAVEFRONT_WIDTH_AMD                   0x4043
#define CL_DEVICE_GFXIP_MAJOR_AMD                       0x404A
#define CL_DEVICE_GFXIP_MINOR_AMD                       0x404B

cl_platform_id mcx_list_gpu(mcconfig *cfg,unsigned int *activedev,cl_device_id *activedevlist,GPUInfo **info);
void ocl_assess(int cuerr,const char *file,const int linenum);

#ifdef  __cplusplus
}
#endif

#endif
