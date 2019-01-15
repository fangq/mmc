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

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define RO_MEM             (CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR)
#define WO_MEM             (CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR)
#define RW_MEM             (CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR)
#define RW_PTR             (CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR)

#define OCL_ASSERT(x)  ocl_assess((x),__FILE__,__LINE__)

typedef struct GPU_mcconfig{
        cl_uint nphoton;      /**<(total simulated photon number) we now use this to 
                             temporarily alias totalmove, as to specify photon
                             number is causing some troubles*/
        cl_int seed;         /**<random number generator seed*/
        cl_float4 srcpos;    /**<src position in mm*/
        cl_float4 srcdir;    /**<src normal direction*/
        cl_int srctype;      /**<src type: 0 - pencil beam, 1 - isotropic ... */
        cl_float4 srcparam1;       /**<source parameters set 1*/
        cl_float4 srcparam2;       /**<source parameters set 2*/
        cl_int voidtime;
        cl_float4 bary0;     /**<initial bary centric coordinates of the source*/
        cl_float tstart;     /**<start time in second*/
        cl_float tstep;      /**<time step in second*/
        cl_float tend;       /**<end time in second*/
        cl_int medianum;     /**<total types of media*/
        cl_int detnum;       /**<total detector numbers*/
        cl_int e0;       /**<total detector numbers*/
        cl_float detradius;  /**<detector radius*/
        cl_float sradius;    /**<source region radius: if set to non-zero, accumulation 
                            will not perform for dist<sradius; this can reduce
                            normalization error when using non-atomic write*/   
        cl_int maxgate;        /**<simultaneous recording gates*/
        cl_int replaydet;      /**<the detector id for which to replay the detected photons, start from 1*/
        cl_char isreflect;     /**<1 for reflecting photons at boundary,0 for exiting*/
        cl_char isnormalized;  /**<1 to normalize the fluence, 0 for raw fluence*/
        cl_char issavedet;     /**<1 to count all photons hits the detectors*/
        cl_char ismomentum;    /**<1 to save momentum transfer for detected photons, implies issavedet=1*/
        cl_char issaveexit;    /**<1 to save the exit position and vector of a detected photon, implies issavedet=1*/
        cl_char issave2pt;     /**<1 to save the 2-point distribution, 0 do not save*/
        cl_char isgpuinfo;     /**<1 to print gpu info when attach, 0 do not print*/
        cl_char isspecular;    /**<1 calculate the initial specular ref if outside the mesh, 0 do not calculate*/
        cl_char issaveseed;    /**<1 save the seed for a detected photon, 0 do not save*/
        cl_char method;        /**<0-Plucker 1-Havel, 2-Badouel, 3-branchless Badouel*/
        cl_char basisorder;    /**<0 to use piece-wise-constant basis for fluence, 1, linear*/
        cl_char outputtype;    /**<'X' output is flux, 'F' output is fluence, 'E' energy deposit*/
        cl_float roulettesize; /**<number of roulette for termination*/
        cl_float minenergy;    /**<minimum energy to propagate photon*/
        cl_float nout;         /**<refractive index for the domain outside the mesh*/
        cl_int isextdet;      /**<if 1, there is external wide-field detector (marked by -2 in the mesh)*/
        cl_uint debuglevel; /**<a flag to control the printing of the debug information*/
        cl_float unitinmm;     /**<define the length unit in mm*/
} gmcconfig __attribute__ ((aligned (16)));

typedef struct OCL_tetmesh{
        cl_int nn;      /**< number of nodes */
        cl_int ne;      /**< number of elements */
        cl_int prop;    /**< number of media */
        cl_int srcelemlen;	  /**< length of the elements that may contain the source*/
        cl_int detelemlen;	  /**< length of the elements that may contain the detector*/
} gtetmesh  __attribute__ ((aligned (16)));

typedef struct OCL_mmcdata{
        cl_mem gcfg, gmesh, gmedia;

	cl_mem *detpos;   /**<detector positions and radius, overwrite detradius*/
	cl_mem *srcpattern;	/**<source pattern*/
        cl_mem *node;/**< node coordinates */
        cl_mem *elem;  /**< element indices */
        cl_mem *elem2; /**< second order element, storing the last 6 elements of a 10-node tet */
        cl_mem *srcelem;  /**< candidate list of elements containing the source*/
        cl_mem *detelem;  /**< candidate list of elements containing a widefield detector*/
        cl_mem *type;  /**< element-based media index */
        cl_mem *facenb;/**< face neighbors, idx of the element sharing a face */
        cl_mem *weight;/**< volumetric fluence for all nodes at all time-gates */
        cl_mem *evol; /**< volume of an element */
        cl_mem *nvol; /**< veronio volume of a node */
	cl_mem *gvisitor;
} gmmcdata; 

int mmc_init_from_cmd(mcconfig *cfg, tetmesh *mesh, raytracer *tracer,int argc, char**argv);
int mmc_init_from_json(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, char *jcfg, char *jmesh);
int mmc_reset(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_cleanup(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_prep(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_run_cl(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);

#ifdef  __cplusplus
}
#endif

#endif
