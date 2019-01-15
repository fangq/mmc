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
\file    mcx_utils_cl.hpp

\brief   Definitions of mmc high-level driver functions for OpenCL backend
*******************************************************************************/

#ifndef _MMC_UTILITIES_CL_H
#define _MMC_UTILITIES_CL_H

#ifdef __APPLE__
  #include <OpenCL/opencl.h>
#else
  #include <CL/opencl.h>
#endif

#include "mcx_utils.h"

char *print_cl_errstring(cl_int err);
void ocl_assess(int cuerr,const char *file,const int linenum);
cl_platform_id mcx_list_gpu(mcconfig *cfg,unsigned int *activedev,cl_device_id *activedevlist);


#endif
