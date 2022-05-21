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
\file    rng_common.h

\brief   random number generator (RNG) independent interface functions
*******************************************************************************/

#ifndef _MMC_RNG_COMMON_H
#define _MMC_RNG_COMMON_H

#define TWO_PI     (M_PI*2.0)
#define EPS        1e-6f
#define EPS2       0.1*EPS

//! generate [0,1] random number for the next scattering length
inlinefun float rand_next_scatlen(RandType t[RAND_BUF_LEN]) {
    return -logf(rand_uniform01(t) + EPS);
}
//! generate [0,1] random number for the next arimuthal angle
inlinefun float rand_next_aangle(RandType t[RAND_BUF_LEN]) {
    return rand_uniform01(t);
}
//! generate random number for the next zenith angle
inlinefun float rand_next_zangle(RandType t[RAND_BUF_LEN]) {
    return rand_uniform01(t);
}
//! generate random number for reflection test
inlinefun float rand_next_reflect(RandType t[RAND_BUF_LEN]) {
    return rand_uniform01(t);
}
//! generate random number for the next zenith angle
inlinefun float rand_do_roulette(RandType t[RAND_BUF_LEN]) {
    return rand_uniform01(t);
}
#ifdef MMC_USE_SSE_MATH      //! when SSE Math functions are used

#define MATH_BLOCK 8

//! generate [0,1] random number for the sin/cos of arimuthal angles
inlinefun void rand_next_aangle_sincos(RandType t[RAND_BUF_LEN], float* si, float* co) {
    static __thread V4SF sine[MATH_BLOCK], cosine[MATH_BLOCK];
    static __thread int pos = (MATH_BLOCK << 2);

    if (pos >= (MATH_BLOCK << 2)) {
        V4SF ran[MATH_BLOCK];
        int i, j;

        for (i = 0; i < MATH_BLOCK; i++)
            for (j = 0; j < 4; j++) {
                ran[i].f[j] = TWO_PI * rand_uniform01(t);
            }

        for (i = 0; i < MATH_BLOCK; i++) {
            sincos_ps(ran[i].v, &(sine[i].v), &(cosine[i].v));
        }

        pos = 0;
    }

    *si = sine[0].f[pos];
    *co = cosine[0].f[pos++];
}

//! generate [0,1] random number for the next scattering length
inlinefun float rand_next_scatlen_ps(RandType t[RAND_BUF_LEN]) {
    static __thread V4SF logval[MATH_BLOCK];
    static __thread int pos = (MATH_BLOCK << 2);
    float res;

    if (pos >= (MATH_BLOCK << 2)) {
        V4SF ran[MATH_BLOCK];
        int i, j;

        for (i = 0; i < MATH_BLOCK; i++)
            for (j = 0; j < 4; j++) {
                ran[i].f[j] = rand_uniform01(t);
            }

        for (i = 0; i < MATH_BLOCK; i++) {
            logval[i].v = log_ps(ran[i].v);
        }

        pos = 0;
    }

    res = ((logval[0].f[pos] != logval[0].f[pos]) ? LOG_RNG_MAX : (-logval[0].f[pos]));
    pos++;
    return res;
}
#endif

#endif
