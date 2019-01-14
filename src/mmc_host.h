/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmc_host.h

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
