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

#ifndef _MCEXTREME_CONSTANT_H
#define _MCEXTREME_CONSTANT_H

#define DET_MASK           0xFFFF0000
#define MED_MASK           0x0000FFFF

#define MAX_PROP           4000

#define MCX_DEBUG_REC_LEN  6  /**<  number of floating points per position saved when -D M is used for trajectory */

#define MCX_SRC_PENCIL     0  /**<  default-Pencil beam src, no param */
#define MCX_SRC_ISOTROPIC  1  /**<  isotropic source, no param */
#define MCX_SRC_CONE       2  /**<  uniform cone, srcparam1.x=max zenith angle in rad */
#define MCX_SRC_GAUSSIAN   3  /**<  Gaussian beam, srcparam1.x=sigma */
#define MCX_SRC_PLANAR     4  /**<  quadrilateral src, vectors spanned by srcparam{1}.{x,y,z} */
#define MCX_SRC_PATTERN    5  /**<  same as above, load srcpattern as intensity */
#define MCX_SRC_FOURIER    6  /**<  same as above, srcparam1.w and 2.w defines the spatial freq in x/y */
#define MCX_SRC_ARCSINE    7  /**<  same as isotropic, but more photons near the pole dir */
#define MCX_SRC_DISK       8  /**<  uniform 2D disk along v */
#define MCX_SRC_FOURIERX   9  /**<  same as Fourier, except the v1/v2 and v are orthogonal */
#define MCX_SRC_FOURIERX2D 10 /**<  2D (sin(kx*x+phix)*sin(ky*y+phiy)+1)/2 */
#define MCX_SRC_ZGAUSSIAN  11 /**<  Gaussian zenith anglular distribution */
#define MCX_SRC_LINE       12 /**<  a non-directional line source */
#define MCX_SRC_SLIT       13 /**<  a collimated line source */
#define MCX_SRC_PENCILARRAY 14 /**<  a rectangular array of pencil beams */
#define MCX_SRC_PATTERN3D  15  /**<  a 3D pattern source, starting from srcpos, srcparam1.{x,y,z} define the x/y/z dimensions */

#define SAVE_DETID(a)         ((a)    & 0x1)   /**<  mask to save detector ID*/
#define SAVE_NSCAT(a)         ((a)>>1 & 0x1)   /**<  output partial scattering counts */
#define SAVE_PPATH(a)         ((a)>>2 & 0x1)   /**<  output partial path */
#define SAVE_MOM(a)           ((a)>>3 & 0x1)   /**<  output momentum transfer */
#define SAVE_PEXIT(a)         ((a)>>4 & 0x1)   /**<  save exit positions */
#define SAVE_VEXIT(a)         ((a)>>5 & 0x1)   /**<  save exit vector/directions */
#define SAVE_W0(a)            ((a)>>6 & 0x1)   /**<  save initial weight */

#define SET_SAVE_DETID(a)     ((a) | 0x1   )   /**<  mask to save detector ID*/
#define SET_SAVE_NSCAT(a)     ((a) | 0x1<<1)   /**<  output partial scattering counts */
#define SET_SAVE_PPATH(a)     ((a) | 0x1<<2)   /**<  output partial path */
#define SET_SAVE_MOM(a)       ((a) | 0x1<<3)   /**<  output momentum transfer */
#define SET_SAVE_PEXIT(a)     ((a) | 0x1<<4)   /**<  save exit positions */
#define SET_SAVE_VEXIT(a)     ((a) | 0x1<<5)   /**<  save exit vector/directions */
#define SET_SAVE_W0(a)        ((a) | 0x1<<6)   /**<  save initial weight */

#define UNSET_SAVE_DETID(a)     ((a) & ~(0x1)   )   /**<  mask to save detector ID*/
#define UNSET_SAVE_NSCAT(a)     ((a) & ~(0x1<<1))   /**<  output partial scattering counts */
#define UNSET_SAVE_PPATH(a)     ((a) & ~(0x1<<2))   /**<  output partial path */
#define UNSET_SAVE_MOM(a)       ((a) & ~(0x1<<3))   /**<  output momentum transfer */
#define UNSET_SAVE_PEXIT(a)     ((a) & ~(0x1<<4))   /**<  save exit positions */
#define UNSET_SAVE_VEXIT(a)     ((a) & ~(0x1<<5))   /**<  save exit vector/directions */
#define UNSET_SAVE_W0(a)        ((a) & ~(0x1<<6))   /**<  save initial weight */

#endif
