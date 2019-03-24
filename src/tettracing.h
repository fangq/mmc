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
\file    tettracing.h

\brief   Definition of the core ray-tracing functions
*******************************************************************************/

#ifndef _MMC_RAY_TRACING_H
#define _MMC_RAY_TRACING_H

#include "simpmesh.h"
#include "mcx_utils.h"

#define MAX_TRIAL          3          /**< number of fixes when a photon hits an edge/vertex */
#define FIX_PHOTON         1e-3f      /**< offset to the ray to avoid edge/vertex */

/***************************************************************************//**
\struct MMC_ray tettracing.h
\brief  Data structure associated with the current photon/ray

*******************************************************************************/   

typedef struct MMC_ray{
	float3 p0;                    /**< current photon position */
	float3 vec;                   /**< current photon direction vector */
	float3 pout;                  /**< the intersection position of the ray to the enclosing tet */
	float4 bary0;                 /**< the Barycentric coordinate of the intersection with the tet */
	int eid;                      /**< the index of the enclosing tet (starting from 1) */
	int faceid;                   /**< the index of the face at which ray intersects with tet */
	int isend;                    /**< if 1, the scattering event ends before reaching the intersection */
	int nexteid;                  /**< the index to the neighboring tet to be moved into */
	float weight;                 /**< photon current weight */
	float photontimer;            /**< the total time-of-fly of the photon */
	float slen0;                  /**< initial unitless scattering length = length*mus */
	float slen;                   /**< the remaining unitless scattering length = length*mus  */
	float Lmove;                  /**< last photon movement length */
	double Eabsorb;               /**< accummulated energy being absorbed */
	unsigned int photonid;        /**< index of the current photon */
	float *partialpath;           /**< pointer to the partial path data buffer */
        void  *photonseed;            /**< pointer to store the photon seed */
	float focus;                  /**< focal length of the source, if defined */
	unsigned int posidx;	      /**< launch position index of the photon for pattern source type */
	unsigned int oldidx;
	double oldweight;
} ray;

/***************************************************************************//**
\struct MMC_visitor tettracing.h
\brief  A structure that accumulates the statistics about the simulation

*******************************************************************************/  

typedef struct MMC_visitor{
	float raytet;                 /**< total number of ray-tet tests */
	float raytet0;                /**< total number of ray-tet tests outside of the object mesh */
	float rtstep;                 /**< reciprocal of the movement step - obsolete */
	int   detcount;               /**< number of total detected photons */
	int   bufpos;                 /**< the position of the detected photon buffer */
	int   reclen;                 /**< record (4-byte per record) number per detected photon */
	float *partialpath;           /**< pointer to the partial path data buffer */
	void  *photonseed;            /**< pointer to store the photon seed */
	double *launchweight;         /**< pointer to accumulated launched photon weight */
	double *absorbweight;         /**< pointer to accumulated absorbed photon weight */
	double *kahanc0;              /**< temp variable to enable Kahan summation to reduce round-off error */
	double *kahanc1;              /**< temp variable to enable Kahan summation to reduce round-off error */
} visitor;

#ifdef	__cplusplus
extern "C" {
#endif
void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout);
void fixphoton(float3 *p,float3 *nodes, int *ee);
void onephoton(size_t id,raytracer *tracer,tetmesh *mesh,mcconfig *cfg,RandType *ran,RandType *ran0, visitor *visit);
void launchphoton(mcconfig *cfg, ray *r, tetmesh *mesh, RandType *ran, RandType *ran0);
float reflectray(mcconfig *cfg,float3 *c0,raytracer *tracer,int *oldeid,int *eid,int faceid,RandType *ran);
void save_scatter_events(ray *r, tetmesh *mesh, mcconfig *cfg, visitor *visit);
void albedoweight(ray *r, tetmesh *mesh, mcconfig *cfg, visitor *visit);
void visitor_init(mcconfig *cfg, visitor* visit);
void visitor_clear(visitor* visit);

#ifdef MCX_CONTAINER
#ifdef __cplusplus
extern "C"
#endif
 int mcx_throw_exception(const int id, const char *msg, const char *filename, const int linenum);
#endif

#ifdef	__cplusplus
}
#endif

#endif
