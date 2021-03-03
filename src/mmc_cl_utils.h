/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2021
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

#ifndef _MCEXTREME_GPU_LAUNCH_H
#define _MCEXTREME_GPU_LAUNCH_H

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#define CL_TARGET_OPENCL_VERSION 120

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

cl_platform_id mcx_list_cl_gpu(mcconfig *cfg,unsigned int *activedev,cl_device_id *activedevlist,GPUInfo **info);
void ocl_assess(int cuerr,const char *file,const int linenum);

#ifdef  __cplusplus
}
#endif

#endif
