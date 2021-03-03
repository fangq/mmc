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
\file    vector_types.h

\brief   Definitions of the basic short vector data structures
*******************************************************************************/

#ifndef _MMC_VECTOR_H
#define _MMC_VECTOR_H

#ifdef _MSC_VER
    #define PRE_ALIGN(x) __declspec(align(x))
    #define POST_ALIGN(x)
#else
    #define PRE_ALIGN(x)
    #define POST_ALIGN(x) __attribute__ ((aligned(x)))
#endif

/**
 \struct MMC_float4 vector_types.h
 \brief  floating-point quadraplet {x,y,z,w}

 the data structure is 16byte aligned to facilitate SSE operations
*/

typedef struct PRE_ALIGN(16) MMC_float4{
    float x,y,z,w;
} float4 POST_ALIGN(16);

/**
 \struct MMC_float3 vector_types.h
 \brief  floating-point triplet {x,y,z}

 if SSE is enabled, float3 is identical to float4
*/

#if defined(MMC_USE_SSE) || defined(USE_OPENCL)
 typedef struct MMC_float4 float3;
#else
 typedef struct MMC_float3{
    float x,y,z;
 } float3;
#endif

/**
 \struct MMC_int2 vector_types.h
 \brief  integer pair {ix,iy}
*/

typedef struct MMC_int2{
    int x,y;
} int2;

/**
 \struct MMC_int3 vector_types.h
 \brief  integer triplet {ix,iy,iz}
*/

typedef struct MMC_int3{
    int x,y,z;
} int3;

/**
 \struct MMC_int4 vector_types.h
 \brief  unsigned integer quadraplet {ix,iy,iz,iw}
*/
typedef struct PRE_ALIGN(16) MMC_int4{
    int x,y,z,w;
} int4 POST_ALIGN(16);

/**
 \struct MMC_uint3 vector_types.h
 \brief  unsigned integer triplet {ix,iy,iz}
*/
typedef struct MMC_uint3{
    unsigned int x,y,z;
} uint3;

/**
 \struct MMC_uint3 vector_types.h
 \brief  unsigned integer triplet {ix,iy,iz}
*/
typedef struct PRE_ALIGN(16) MMC_uint4{
    unsigned int x,y,z,w;
} uint4 POST_ALIGN(16);

/**
 \struct MMC_uint2 vector_types.h
 \brief  unsigned integer pair {ix,iy}
*/

typedef struct MMC_uint2{
    unsigned int x,y;
} uint2;


typedef unsigned int uint;
typedef unsigned char uchar;

#endif
