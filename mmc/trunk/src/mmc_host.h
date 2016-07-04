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
\file    mcx_host.h

\brief   Definition of mmc high-level driver functions
*******************************************************************************/

#ifndef _MMC_HOSTCODE_H
#define _MMC_HOSTCODE_H

#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"

int mmc_init_from_cmd(mcconfig *cfg, tetmesh *mesh, raytracer *tracer,int argc, char**argv);
int mmc_init_from_json(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, char *jcfg, char *jmesh);
int mmc_reset(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_cleanup(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_prep(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);
int mmc_run_mp(mcconfig *cfg, tetmesh *mesh, raytracer *tracer);

#endif
