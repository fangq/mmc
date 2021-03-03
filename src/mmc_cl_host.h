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

/***************************************************************************//**
\file    mcx_host_cl.hpp

\brief   Definitions of mmc high-level driver functions for OpenCL backend
*******************************************************************************/

#ifndef _MMC_HOSTCODE_CL_H
#define _MMC_HOSTCODE_CL_H

#include "mmc_cl_utils.h"
#include "simpmesh.h"
#include "tettracing.h"

#ifdef  __cplusplus
extern "C" {
#endif

#ifndef CL_MEM_LOCATION_HOST_NV
  #define CL_MEM_LOCATION_HOST_NV                     (1 << 0)
  typedef cl_bitfield         cl_mem_flags_NV;
#endif

#define RO_MEM             (CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR)
#define WO_MEM             (CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR)
#define RW_MEM             (CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR)
#define RW_PTR             (CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR)
#define NV_PIN             CL_MEM_LOCATION_HOST_NV

#define OCL_ASSERT(x)  ocl_assess((x),__FILE__,__LINE__)

#define RAND_SEED_WORD_LEN      4        //48 bit packed with 64bit length

typedef struct PRE_ALIGN(32) GPU_mcconfig{
  cl_float3 srcpos;
  cl_float3 srcdir;
  cl_float  tstart,tend;
  cl_uint   isreflect,issavedet,issaveexit,ismomentum,isatomic,isspecular;
  cl_float  Rtstep;
  cl_float  minenergy;
  cl_uint   maxdetphoton;
  cl_uint   maxmedia;
  cl_uint   detnum;
  cl_int    voidtime;
  cl_int    srctype;                    /**< type of the source */
  cl_float4 srcparam1;                  /**< source parameters set 1 */
  cl_float4 srcparam2;                  /**< source parameters set 2 */
  cl_uint   issaveref;     /**<1 save diffuse reflectance at the boundary voxels, 0 do not save*/
  cl_uint   maxgate;
  cl_uint   debuglevel;           /**< debug flags */
  cl_int    reclen;                 /**< record (4-byte per record) number per detected photon */
  cl_int    outputtype;
  cl_int    elemlen;
  cl_int    mcmethod;
  cl_int    method;
  cl_float  dstep;
  cl_float  focus;
  cl_int    nn, ne, nf;
  cl_float3 nmin;
  cl_float  nout;
  cl_uint   roulettesize;
  cl_int    srcnum;
  cl_int4   crop0;
  cl_int    srcelemlen;
  cl_float4 bary0;
  cl_int    e0;
  cl_int    isextdet;
  cl_int    framelen;
  cl_uint   nbuffer;
  cl_uint   maxpropdet;
  cl_uint   normbuf;
  cl_int    issaveseed;
  cl_uint   seed;
} MCXParam POST_ALIGN(32);

typedef struct POST_ALIGN(32) GPU_reporter{
  float  raytet;
} MCXReporter  POST_ALIGN(32);

void mmc_run_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, void (*progressfun)(float, void *),void *handle);

#ifdef  __cplusplus
}
#endif

#endif
