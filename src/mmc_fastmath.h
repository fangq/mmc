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
\file    fastmath.h

\brief   A library for approximated math functions
*******************************************************************************/

#ifndef _MMC_FAST_MATH_FUNCTIONS_H
#define _MMC_FAST_MATH_FUNCTIONS_H

/**
 * @brief An approximated floating-point exponential function using 9th order Taylor's expansion
 *
 * GNU Glibc math function expf() is quite slow compared to Intel's MKL, we used an approximated
 * expf in the code when GCC is used to compile MMC.
 *
 * This function can only handle small x values
 */

static inline float fast_expf9(float x) {
    return (362880.f + x * (362880.f + x * (181440.f + x * (60480.f + x * (15120.f + x * (3024.f + x * (504.f + x * (72.f + x * (9.f + x))))))))) * 2.75573192e-6f;
}

/*
static inline float fast_log2(float x) {
        int i = (*(int *)&x);
        return (((i&0x7f800000)>>23)-0x7f)+(i&0x007fffff)/(float)0x800000;
}

static inline float fast_logf(float x) {
        return 0.693147180559945f*fast_log2(x);
}
*/

#endif
