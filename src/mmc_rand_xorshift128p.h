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
\file    sfmt_rand.h

\brief   An interface to use the SFMT-19937 random number generator
*******************************************************************************/

#ifndef _MCEXTREME_SFMT_RAND_H
#define _MCEXTREME_SFMT_RAND_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef inlinefun
    #define inlinefun static inline
#endif

typedef unsigned long long RandType;
typedef unsigned int uint;

#if defined(_WIN32) || defined(__APPLE__)
    typedef unsigned long long ulong;
#endif

#define MCX_RNG_NAME       "xorshift128+ RNG"
#define RAND_BUF_LEN       2        //register arrays
#define RAND_SEED_WORD_LEN      4        //48 bit packed with 64bit length

inlinefun void rng_init(RandType t[RAND_BUF_LEN], RandType tnew[RAND_BUF_LEN], uint* n_seed, int idx);

inlinefun void rand_need_more(RandType t[RAND_BUF_LEN], RandType tbuf[RAND_BUF_LEN]);
// generate [0,1] random number for the next scattering length
inlinefun float rand_next_scatlen(RandType t[RAND_BUF_LEN]);
// generate [0,1] random number for the next arimuthal angle
inlinefun float rand_next_aangle(RandType t[RAND_BUF_LEN]);
// generate random number for the next zenith angle
inlinefun float rand_next_zangle(RandType t[RAND_BUF_LEN]);
inlinefun float rand_next_reflect(RandType t[RAND_BUF_LEN]);
inlinefun float rand_do_roulette(RandType t[RAND_BUF_LEN]);

#ifdef MMC_USE_SSE_MATH
    inlinefun void rand_next_aangle_sincos(RandType t[RAND_BUF_LEN], float* si, float* co);
    inlinefun float rand_next_scatlen_ps(RandType t[RAND_BUF_LEN]);
#endif

#endif
