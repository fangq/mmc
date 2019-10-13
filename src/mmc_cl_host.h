/*******************************************************************************
**  Mesh-based Monte Carlo (MMC) -- this unit was ported from MCX
**
**  Monte Carlo eXtreme (MCX)  - GPU accelerated Monte Carlo 3D photon migration
**  Author: Qianqian Fang <q.fang at neu.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates," Biomed. Opt. Express, 1(1) 165-175 (2010)
**
**  (Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, 17(22) 20178-20190 (2009)
**
**  License: GPL v3, see LICENSE.txt for details
**
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

typedef struct GPU_mcconfig{
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
  cl_uint   buffermask;
  //cl_int    issaveseed;
} MCXParam __attribute__ ((aligned (32)));

typedef struct GPU_reporter{
  float  raytet;
} MCXReporter  __attribute__ ((aligned (32)));

void mmc_run_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);

#ifdef  __cplusplus
}
#endif

#endif
