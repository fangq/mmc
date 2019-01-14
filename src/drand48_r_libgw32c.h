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
\file    drand48_r_libgw32c.h

\brief   Windows 32 port of drand48_r random number generator from libgw2c
*******************************************************************************/

#ifndef _MMC_POSIX_RAND_LIBGW32C_H
#define _MMC_POSIX_RAND_LIBGW32C_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <ieee754.h>
#include <string.h>

//typedef unsigned long long uint64_t

typedef signed __int8 sint8;
typedef unsigned __int8 uint8;
typedef signed __int16 sint16;
typedef unsigned __int16 uint16;
typedef signed __int32 sint32;
typedef unsigned __int32 uint32;
typedef signed __int64 sint64;
typedef unsigned __int64 uint64;


/** 
\struct drand48_data drand48_r_libgw32c.h
\brief Data structure for communication with thread safe versions.

This type is to be regarded as opaque. It's only exported because users
have to allocate objects of this type.
*/
struct drand48_data
{
  unsigned short int __x[3]; /* Current state. */
  unsigned short int __old_x[3]; /* Old state. */
  unsigned short int __c; /* Additive const. in congruential formula. */
  unsigned short int __init; /* Flag for initializing. */
  uint64 __a; /* Factor in congruential formula. */
};

#ifdef  __cplusplus
extern "C" {
#endif


/* Global state for non-reentrant functions. */
int __drand48_iterate (unsigned short int xsubi[3], struct drand48_data *buffer );
int __erand48_r (unsigned short int xsubi[3],struct drand48_data *buffer, double *result);
int drand48_r (struct drand48_data *buffer, double *result);
int seed48_r (unsigned short int seed16v[3],struct drand48_data *buffer);
double erand48 ( unsigned short int xsubi[3]);


#ifdef  __cplusplus
}
#endif

#endif
