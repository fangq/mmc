/***************************************************************************//**
**  \mainpage Monte Carlo eXtreme - GPU accelerated Monte Carlo Photon Migration
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2009-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**  \li \c (\b Yu2018) Leiming Yu, Fanny Nina-Paravecino, David Kaeli, and Qianqian Fang,
**          "Scalable and massively parallel Monte Carlo photon transport
**           simulations for heterogeneous computing platforms," J. Biomed. Optics, 23(1), 010504, 2018.
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmcx_core.h

@brief   MMC GPU kernel header file
*******************************************************************************/

#ifndef _MMCX_HOSTCODE_H
#define _MMCX_HOSTCODE_H

#include "mcx_utils.h"
#include "simpmesh.h"

#define RAND_SEED_WORD_LEN      4        //48 bit packed with 64bit length

#ifdef  __cplusplus
extern "C" {
#endif

#define ABS(a)  ((a)<0?-(a):(a))

#define MCX_DEBUG_RNG       2                   /**< MCX debug flags */
#define MCX_DEBUG_MOVE      1
#define MCX_DEBUG_PROGRESS  2048

typedef unsigned char uchar;

void mmc_run_cu(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, void (*progressfun)(float, void *),void *handle);

#ifdef  __cplusplus
}
#endif

#endif

