/*****************************************************************//**
*            A library for approximated math functions               *
*                                                                    *
*  \author Qianqian Fang <q.fang at neu.edu>                         *
*                                                                    *
*********************************************************************/

/***************************************************************************//**
\file    fastmath.h

\brief   A library for approximated math functions
*******************************************************************************/

#ifndef _MMC_FAST_MATH_FUNCTIONS_H
#define _MMC_FAST_MATH_FUNCTIONS_H

static inline float fast_expf9(float x) {
        return (362880.f+x*(362880.f+x*(181440.f+x*(60480.f+x*(15120.f+x*(3024.f+x*(504.f+x*(72.f+x*(9.f+x)))))))))*2.75573192e-6f;
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
