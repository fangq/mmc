/***************************************************************************//**
**  \mainpage Monte Carlo eXtreme - GPU accelerated Monte Carlo Photon Migration \
**      -- OpenCL edition
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2009-2018
**
**  \section sref Reference:
**  \li \c (\b Yu2018) Leiming Yu, Fanny Nina-Paravecino, David Kaeli, and Qianqian Fang,
**         "Scalable and massively parallel Monte Carlo photon transport simulations 
**         for heterogeneous computing platforms," J. Biomed. Optics, 23(1), 010504 (2018)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

#ifdef __CUDA_ARCH__
  #define __constant const
  #define __private 
  #define __local 
  #define __global
  #define __kernel __global__
  #define get_global_id 

  inline __device__ __host__ float3 operator *(float3 a, float3 b){
      return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
  }
  inline __device__ __host__ float3 operator *(float f, float3 v){
      return make_float3(v.x*f, v.y*f, v.z*f);
  }
  inline __device__ __host__ float3 operator *(float3 v, float f){
      return make_float3(v.x*f, v.y*f, v.z*f);
  }
  inline __device__ __host__ void operator *=(float3 & b, float f){
      b.x *= f;
      b.y *= f;
      b.z *= f;
  }
  inline __device__ __host__ float3 operator +(float3 a, float3 b){
      return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
  }
  inline __device__ __host__ void operator +=(float3 & b, float3 a){
      b.x += a.x;
      b.y += a.y;
      b.z += a.z;
  }
  inline __device__ __host__ float3 operator -(float3 a, float3 b){
      return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
  }
  inline __device__ __host__ void operator -=(float3 & b, float3 a){
      b.x -= a.x;
      b.y -= a.y;
      b.z -= a.z;
  }
  inline __device__ __host__ float3 operator /(float3 v, float f){
      float inv = 1.0f / f;
      return v * inv;
  }
  inline __device__ __host__ void operator /=(float3 & b, float f){
      float inv = 1.0f / f;
      b.x *= inv;
      b.y *= inv;
      b.z *= inv;
  }

  inline __device__ __host__ float4 operator *(float4 a, float4 b){
      return make_float4(a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w);
  }
  inline __device__ __host__ float4 operator *(float f, float4 v){
      return make_float4(v.x*f, v.y*f, v.z*f, v.w*f);
  }
  inline __device__ __host__ float4 operator *(float4 v, float f){
      return make_float4(v.x*f, v.y*f, v.z*f, v.w*f);
  }
  inline __device__ __host__ void operator *=(float4 & b, float f){
      b.x *= f;
      b.y *= f;
      b.z *= f;
      b.w *= f;
  }
  inline __device__ __host__ float4 operator +(float4 a, float4 b){
      return make_float4(a.x+b.x, a.y+b.y, a.z+b.w, a.w+b.w);
  }
  inline __device__ __host__ void operator +=(float4 & b, float4 a){
      b.x += a.x;
      b.y += a.y;
      b.z += a.z;
      b.w += a.w;
  }
  inline __device__ __host__ float4 operator -(float4 a, float4 b){
      return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w);
  }
  inline __device__ __host__ void operator -=(float4 & b, float4 a){
      b.x -= a.x;
      b.y -= a.y;
      b.z -= a.z;
      b.w -= a.w;
  }
  inline __device__ __host__ float3 operator /(float3 a, float3 b){
      return make_float3(a.x/b.x, a.y/b.y, a.z/b.z);
  }
  inline __device__ __host__ float4 operator /(float4 v, float f){
      float inv = 1.0f / f;
      return v * inv;
  }
  inline __device__ __host__ void operator /=(float4 & b, float f){
      float inv = 1.0f / f;
      b.x *= inv;
      b.y *= inv;
      b.z *= inv;
      b.w *= inv;
  }
  inline __device__ __host__ float4 operator /(float4 a, float4 b){
      return make_float4(a.x/b.x, a.y/b.y, a.z/b.z, a.w/b.w);
  }

  inline __device__ __host__ float dot(float3 a, float3 b){
      return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  inline __device__ __host__ float dot(float4 a, float4 b){
      return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
  inline __device__ __host__ float clamp(float f, float a, float b){
      return max(a, min(f, b));
  }
  inline __device__ __host__ float3 clamp(float3 v, float a, float b){
      return make_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
  }

  inline __device__ __host__ float3 clamp(float3 v, float3 a, float3 b){
      return make_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
  }

  #define FL4(f) make_float4(f,f,f,f)
  #define FL3(f) make_float3(f,f,f)
  #define FLT_EPSILON   1.19209290E-07F
  #define atomicadd(a,b)  atomicAdd(a,b)
  #define atomic_add(a,b) atomicAdd(a,b)
  #define sincos(a,b)     sincosf(a,b)

#else
  #define FL4(f) (f)
  #define FL3(f) (f)
#endif

#ifdef MCX_SAVE_DETECTORS
  #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#endif

#ifdef USE_HALF
  #pragma OPENCL EXTENSION cl_khr_fp16 : enable
  #define FLOAT4VEC half4
  #define TOFLOAT4  convert_half4
#else
  #define FLOAT4VEC float4
  #define TOFLOAT4
#endif

#ifdef MCX_USE_NATIVE
  #define MCX_MATHFUN(fun)              native_##fun
  #define MCX_SINCOS(theta,osin,ocos)   {(osin)=native_sin(theta);(ocos)=native_cos(theta);} 
#else
  #define MCX_MATHFUN(fun)              fun
  #define MCX_SINCOS(theta,osin,ocos)   (ocos)=sincos((theta),&(osin))
#endif

#ifdef MCX_GPU_DEBUG
  #define MMC_FPRINTF(x)        printf x             // enable debugging in CPU mode
#else
  #define MMC_FPRINTF(x)        {}
#endif

#define R_PI               0.318309886183791f

#define ONE_PI             3.1415926535897932f     //pi
#define TWO_PI             6.28318530717959f       //2*pi

#define C0                 299792458000.f          //speed of light in mm/s
#define R_C0               3.335640951981520e-12f  //1/C0 in s/mm

#define EPS                FLT_EPSILON             //round-off limit
#define VERY_BIG           (1.f/FLT_EPSILON)       //a big number
#define JUST_ABOVE_ONE     1.0001f                 //test for boundary
#define SAME_VOXEL         -9999.f                 //scatter within a voxel
#define NO_LAUNCH          9999                    //when fail to launch, for debug
#define MAX_PROP           2000                     /*maximum property number*/
#define ID_UNDEFINED       0xFFFFFFFFU              /**< flag indicating the index is outside of the volume */

#define DET_MASK           0xFFFF0000
#define MED_MASK           0x0000FFFF
#define MAX_ACCUM          1000.f
#define R_MIN_MUS          1e9f
#define FIX_PHOTON         1e-3f      /**< offset to the ray to avoid edge/vertex */
#define MAX_TRIAL          3          /**< number of fixes when a photon hits an edge/vertex */

#define MCX_DEBUG_MOVE      1
#define MCX_DEBUG_PROGRESS  2048

#define MMC_UNDEFINED      (3.40282347e+38F)

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define F32N(a) ((a) & 0x80000000)          /**<  Macro to test if a floating point is negative */
#define F32P(a) ((a) ^ 0x80000000)          /**<  Macro to test if a floating point is positive */

typedef struct MMC_Ray{
	float3 p0;                    /**< current photon position */
	float3 vec;                   /**< current photon direction vector */
	float3 pout;                  /**< the intersection position of the ray to the enclosing tet */
	int eid;                      /**< the index of the enclosing tet (starting from 1) */
	int faceid;                   /**< the index of the face at which ray intersects with tet */
	int isend;                    /**< if 1, the scattering event ends before reaching the intersection */
	float weight;                 /**< photon current weight */
	float photontimer;            /**< the total time-of-fly of the photon */
	float slen;                   /**< the remaining unitless scattering length = length*mus  */
	float Lmove;                  /**< last photon movement length */
	uint oldidx;
	float oldweight;
	//int nexteid;                  /**< the index to the neighboring tet to be moved into */
	//float4 bary0;                 /**< the Barycentric coordinate of the intersection with the tet */
	//float slen0;                  /**< initial unitless scattering length = length*mus */
	//unsigned int photonid;        /**< index of the current photon */
	//unsigned int posidx;	      /**< launch position index of the photon for pattern source type */
} ray __attribute__ ((aligned (32)));


typedef struct MMC_Parameter {
  float3 srcpos;
  float3 srcdir;
  float  tstart,tend;
  uint   isreflect,issavedet,issaveexit,ismomentum,isatomic,isspecular;
  float  Rtstep;
  float  minenergy;
  uint   maxdetphoton;
  uint   maxmedia;
  uint   detnum;
  int    voidtime;
  int    srctype;                    /**< type of the source */
  float4 srcparam1;                  /**< source parameters set 1 */
  float4 srcparam2;                  /**< source parameters set 2 */
  uint   issaveref;     /**<1 save diffuse reflectance at the boundary voxels, 0 do not save*/
  uint   maxgate;
  uint   debuglevel;           /**< debug flags */
  int    reclen;                 /**< record (4-byte per record) number per detected photon */
  int    outputtype;
  int    elemlen;
  int    mcmethod;
  int    method;
  float  dstep;
  float  focus;
  int    nn, ne, nf;
  float3 nmin;
  float  nout;
  uint   roulettesize;
  int    srcnum;
  int4   crop0;
  int    srcelemlen;
  float4 bary0;
  int    e0;
  int    isextdet;
  int    framelen;
  uint   nbuffer;
  uint   buffermask;
  //int    issaveseed;
} MCXParam __attribute__ ((aligned (32)));

typedef struct MMC_Reporter{
  float  raytet;
} MCXReporter  __attribute__ ((aligned (32)));

typedef struct MMC_medium{
        float mua;                     /**<absorption coeff in 1/mm unit*/
        float mus;                     /**<scattering coeff in 1/mm unit*/
        float g;                       /**<anisotropy*/
        float n;                       /**<refractive index*/
} medium __attribute__ ((aligned (32)));

__constant int faceorder[]={1,3,2,0,-1};
__constant int ifaceorder[]={3,0,2,1};
//__constant int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
//__constant int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
__constant int out[4][3]={{0,3,1},{3,2,1},{0,2,3},{0,1,2}};
__constant int facemap[]={2,0,1,3};
__constant int ifacemap[]={1,2,0,3};

enum TDebugLevel {dlMove=1,dlTracing=2,dlBary=4,dlWeight=8,dlDist=16,dlTracingEnter=32,
                  dlTracingExit=64,dlEdge=128,dlAccum=256,dlTime=512,dlReflect=1024,
                  dlProgress=2048,dlExit=4096};

enum TRTMethod {rtPlucker, rtHavel, rtBadouel, rtBLBadouel, rtBLBadouelGrid};   
enum TMCMethod {mmMCX, mmMCML};

enum TSrcType {stPencil, stIsotropic, stCone, stGaussian, stPlanar,
               stPattern, stFourier, stArcSin, stDisk, stFourierX, 
               stFourier2D, stZGaussian, stLine, stSlit};
enum TOutputType {otFlux, otFluence, otEnergy, otJacobian, otWL, otWP};
enum TOutputFormat {ofASCII, ofBin, ofJSON, ofUBJSON};
enum TOutputDomain {odMesh, odGrid};

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define RAND_BUF_LEN       2        //register arrays
#define RAND_SEED_WORD_LEN      4        //48 bit packed with 64bit length
#define LOG_MT_MAX         22.1807097779182f
#define IEEE754_DOUBLE_BIAS     0x3FF0000000000000ul /* Added to exponent.  */

typedef ulong  RandType;

static float xorshift128p_nextf (__private RandType t[RAND_BUF_LEN]){
   union {
        ulong  i;
	float f[2];
	uint  u[2];
   } s1;
   const ulong s0 = t[1];
   s1.i = t[0];
   t[0] = s0;
   s1.i ^= s1.i << 23; // a
   t[1] = s1.i ^ s0 ^ (s1.i >> 18) ^ (s0 >> 5); // b, c
   s1.i = t[1] + s0;
   s1.u[0] = 0x3F800000U | (s1.u[0] >> 9);

   return s1.f[0] - 1.0f;
}

static void copystate(__local float *v1, __private float *v2, int len){
    for(int i=0;i<len;i++)
        v1[i]=v2[i];
}

static float rand_uniform01(__private RandType t[RAND_BUF_LEN]){
    return xorshift128p_nextf(t);
}

static void xorshift128p_seed (__global uint *seed,RandType t[RAND_BUF_LEN]){
    t[0] = (ulong)seed[0] << 32 | seed[1] ;
    t[1] = (ulong)seed[2] << 32 | seed[3];
}

static void gpu_rng_init(__private RandType t[RAND_BUF_LEN], __global uint *n_seed, int idx){
    xorshift128p_seed((n_seed+idx*RAND_SEED_WORD_LEN),t);
}

float rand_next_scatlen(__private RandType t[RAND_BUF_LEN]){
    return -MCX_MATHFUN(log)(rand_uniform01(t)+EPS);
}

#define rand_next_aangle(t)  rand_uniform01(t)
#define rand_next_zangle(t)  rand_uniform01(t)
#define rand_next_reflect(t) rand_uniform01(t)
#define rand_do_roulette(t)  rand_uniform01(t) 

#ifdef USE_ATOMIC

#ifndef __CUDA_ARCH__

// OpenCL float atomicadd hack:
// http://suhorukov.blogspot.co.uk/2011/12/opencl-11-atomic-operations-on-floating.html
// https://devtalk.nvidia.com/default/topic/458062/atomicadd-float-float-atomicmul-float-float-/

inline float atomicadd(volatile __global float* address, const float value){
    float old = value;
    while ((old = atomic_xchg(address, atomic_xchg(address, 0.0f)+old))!=0.0f);
    return old;
}

/*

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

inline double atomicadd(__global double *val, const double delta){
  union {
  double f;
  ulong  i;
  } old, new;

  do{
     old.f = *val;
     new.f = old.f + delta;
  } while (atom_cmpxchg((volatile __global ulong *)val, old.i, new.i) != old.i);
  return old.f;
}
*/
#endif

#endif

void clearpath(__local float *p, int len){
      uint i;
      for(i=0;i<len;i++)
      	   p[i]=0.f;
}

#ifdef MCX_SAVE_DETECTORS
uint finddetector(float3 *p0,__constant float4 *gdetpos,__constant MCXParam *gcfg){
      uint i;
      for(i=0;i<gcfg->detnum;i++){
      	if((gdetpos[i].x-p0[0].x)*(gdetpos[i].x-p0[0].x)+
	   (gdetpos[i].y-p0[0].y)*(gdetpos[i].y-p0[0].y)+
	   (gdetpos[i].z-p0[0].z)*(gdetpos[i].z-p0[0].z) < gdetpos[i].w){
	        return i+1;
	   }
      }
      return 0;
}

void savedetphoton(__global float *n_det,__global uint *detectedphoton,
                   __local float *ppath,float3 *p0,float3 *v,__constant float4 *gdetpos,
		   int extdetid, __constant MCXParam *gcfg){
      uint detid=(extdetid<0)? finddetector(p0,gdetpos,gcfg) : extdetid;
      if(detid){
	 uint baseaddr=atomic_inc(detectedphoton);
	 if(baseaddr<gcfg->maxdetphoton){
	    uint i;
	    baseaddr*=(gcfg->reclen+1);
	    n_det[baseaddr++]=detid;
	    for(i=0;i<(gcfg->maxmedia<<1);i++)
		n_det[baseaddr++]=ppath[i]; // save partial pathlength to the memory
	    for(i=0;i<gcfg->ismomentum*gcfg->maxmedia;i++)
		n_det[baseaddr++]=ppath[i+(gcfg->maxmedia<<1)]; // save partial pathlength to the memory
	    if(gcfg->issaveexit){
	        n_det[baseaddr++]=p0->x;
		n_det[baseaddr++]=p0->y;
		n_det[baseaddr++]=p0->z;
		n_det[baseaddr++]=v->x;
		n_det[baseaddr++]=v->y;
		n_det[baseaddr++]=v->z;
	    }
            n_det[baseaddr++]=ppath[gcfg->reclen-1]; // save partial pathlength to the memory
	 }
      }
}
#endif

/** 
 * \brief Branch-less Badouel-based SSE4 ray-tracer to advance photon by one step
 * 
 * this function uses Branch-less Badouel-based SSE4 ray-triangle intersection 
 * tests to advance photon by one step, see Fang2012. Both Badouel and 
 * Branch-less Badouel algorithms do not calculate the Barycentric coordinates
 * and can only store energy loss using 0-th order basis function. This function
 * is the fastest among the 4 ray-tracers.
 *
 * \param[in,out] r: the current ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

float branchless_badouel_raytet(ray *r, __constant MCXParam *gcfg,__constant int *elem,__global float *weight,
    int type, __constant int *facenb, __constant float4 *normal, __constant medium *med){

	float Lmin;
	float ww,totalloss=0.f;
	int tshift,faceidx=-1,eid;
	float4 T,S;
	union {
	    float f;
	    uint  i;
	} currweight;

	if(r->eid<=0) 
		return -1;

	eid=(r->eid-1)<<2;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;

	S = FL4(r->vec.x)*normal[eid]+FL4(r->vec.y)*normal[eid+1]+FL4(r->vec.z)*normal[eid+2];
	T = normal[eid+3] - (FL4(r->p0.x)*normal[eid]+FL4(r->p0.y)*normal[eid+1]+FL4(r->p0.z)*normal[eid+2]);
        T = -convert_float4_rte(isgreater(T,FL4(0.f))*2)*FL4(0.5f)*T;
	T = T/S;

        S = -convert_float4_rte(isgreater(S,FL4(0.f))*2)*FL4(0.5f);
        T =  S * T + (FL4(1.f)-S) * FL4(1e10f);

	eid=r->eid-1;

	Lmin=fmin(fmin(fmin(T.x,T.y),T.z),T.w);
	faceidx=((Lmin==1e10f)? 4: Lmin==T.x? 0: (Lmin==T.y? 1 : (Lmin==T.z ? 2 : 3)));
	r->faceid=faceorder[faceidx];

	if(r->faceid>=0 && Lmin>=0.f){
	    medium prop;

	    prop=med[type];
            currweight.f=r->weight;

	    r->Lmove=(prop.mus <= EPS) ? R_MIN_MUS : r->slen/prop.mus;
	    r->isend=(Lmin>r->Lmove);
	    r->Lmove=((r->isend) ? r->Lmove : Lmin);
	    r->pout=r->p0+FL3(Lmin)*r->vec;

	    if((int)((r->photontimer+r->Lmove*(prop.n*R_C0)-gcfg->tstart)*gcfg->Rtstep)>gcfg->maxgate-1){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(gcfg->tend-r->photontimer)/(prop.n*R_C0)-1e-4f;
	    }
            totalloss=MCX_MATHFUN(exp)(-prop.mua*r->Lmove);
            r->weight*=totalloss;

	    totalloss=1.f-totalloss;    /*remaining fraction*/
	    r->slen-=r->Lmove*prop.mus;
	    ww=currweight.f-r->weight;
            r->photontimer+=r->Lmove*(prop.n*R_C0);
/*
	    if(gcfg->outputtype==otWL || gcfg->outputtype==otWP)
		    tshift=MIN( ((int)(replaytime[r->photonid]*gcfg->Rtstep)), gcfg->maxgate-1 )*gcfg->framelen;
	    else
*/
                    tshift=MIN( ((int)((r->photontimer-gcfg->tstart)*gcfg->Rtstep)), gcfg->maxgate-1 )*gcfg->framelen;
            {
#ifndef MCX_SKIP_VOLUME
	       if(prop.mua>0.f){
	           if(gcfg->outputtype!=otEnergy && gcfg->outputtype!=otWP)
                      ww/=prop.mua;
	       }
  #ifndef USE_DMMC
               uint newidx=eid+tshift;
	       r->oldidx=(r->oldidx==ID_UNDEFINED)? newidx: r->oldidx;
	       if(newidx!=r->oldidx){
#ifndef DO_NOT_SAVE
    #ifdef USE_ATOMIC
                   atomicadd(weight+r->oldidx,r->oldweight);
    #else
                   weight[r->oldidx]+=r->oldweight;
    #endif
#endif
                   r->oldidx=newidx;
		   r->oldweight=ww;
               }else{
	           r->oldweight+=ww;
	       }
  #else
		   eid=(int)(r->Lmove*gcfg->dstep)+1;    // number of segments
		   eid=(eid<<1);
		   S.w=r->Lmove/eid;                     // segment length
	           T.w=MCX_MATHFUN(exp)(-prop.mua*S.w);  // segment loss
		   T.xyz =  r->vec * FL3(S.w);      // delta vector
		   S.xyz =  (r->p0 - gcfg->nmin) + (T.xyz * FL3(0.5f)); /*starting point*/
		   totalloss=(totalloss==0.f)? 0.f : (1.f-T.w)/totalloss;    // fraction of total loss per segment
		   S.w=ww;                               // S.w is now the current weight
                   for(faceidx=0; faceidx< eid; faceidx++){
		       int3 idx= convert_int3_rtn(S.xyz * FL3(gcfg->dstep));
		       idx = idx & (idx>=(int3)(0));
		       uint newidx=(idx.z*gcfg->crop0.y+idx.y*gcfg->crop0.x+idx.x)+tshift;
		       r->oldidx=(r->oldidx==ID_UNDEFINED)? newidx: r->oldidx;
		       if(newidx!=r->oldidx){
#ifndef DO_NOT_SAVE
    #ifdef USE_ATOMIC
		           atomicadd(weight+r->oldidx,r->oldweight);
    #else
                           weight[r->oldidx]+=r->oldweight;
    #endif
#endif
                           r->oldidx=newidx;
                           r->oldweight=S.w*totalloss;
		       }else
                           r->oldweight+=S.w*totalloss;
		       S.w*=T.w;
		       S.xyz += T.xyz;
                   }
  #endif
#endif
            }
	    r->p0=r->p0+FL3(r->Lmove)*r->vec;
	}
	return ((r->faceid==-2) ? 0.f : r->slen);
}



/**
 * @brief Calculate the reflection/transmission of a ray at an interface
 *
 * This function handles the reflection and transmission events at an interface
 * where the refractive indices mismatch.
 *
 * \param[in,out] cfg: simulation configuration structure
 * \param[in] c0: the current direction vector of the ray
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] oldeid: the index of the element the photon moves away from
 * \param[in] eid: the index of the element the photon about to move into
 * \param[in] faceid: index of the face through which the photon reflects/transmits
 * \param[in,out] ran: the random number generator states
 */

#ifdef MCX_DO_REFLECTION

float reflectray(__constant MCXParam *gcfg,float3 *c0, int *oldeid,int *eid,int faceid,__private RandType *ran, __constant int *type, __constant float4 *normal, __constant medium *med){
	/*to handle refractive index mismatch*/
        float3 pnorm={0.f,0.f,0.f};
	float Icos,Re,Im,Rtotal,tmp0,tmp1,tmp2,n1,n2;
	int offs=(*oldeid-1)<<2;

	faceid=ifaceorder[faceid];
	/*calculate the normal direction of the intersecting triangle*/
	pnorm.x=((__constant float*)&(normal[offs]))[faceid];
	pnorm.y=((__constant float*)&(normal[offs]))[faceid+4];
	pnorm.z=((__constant float*)&(normal[offs]))[faceid+8];

	/*pn pointing outward*/

	/*compute the cos of the incidence angle*/
        Icos=fabs(dot(*c0,pnorm));

	n1=((*oldeid!=*eid) ? med[type[*oldeid-1]].n : gcfg->nout);
	n2=((*eid>0) ? med[type[*eid-1]].n : gcfg->nout);

	tmp0=n1*n1;
	tmp1=n2*n2;
        tmp2=1.f-tmp0/tmp1*(1.f-Icos*Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

        if(tmp2>0.f){ /*if no total internal reflection*/
          Re=tmp0*Icos*Icos+tmp1*tmp2;      /*transmission angle*/
	  tmp2=MCX_MATHFUN(sqrt)(tmp2); /*to save one sqrt*/
          Im=2.f*n1*n2*Icos*tmp2;
          Rtotal=(Re-Im)/(Re+Im);     /*Rp*/
          Re=tmp1*Icos*Icos+tmp0*tmp2*tmp2;
          Rtotal=(Rtotal+(Re-Im)/(Re+Im))*0.5f; /*(Rp+Rs)/2*/
	  if(*oldeid==*eid) return Rtotal; /*initial specular reflection*/
	  if(rand_next_reflect(ran)<=Rtotal){ /*do reflection*/
              *c0+=(FL3(-2.f*Icos))*pnorm;
              //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(("R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal));
	      *eid=*oldeid; /*stay with the current element*/
	  }else if(gcfg->isspecular==2 && *eid==0){
              // if do transmission, but next neighbor is 0, terminate
          }else{                              /*do transmission*/
              *c0+=(FL3(-Icos))*pnorm;
	      *c0=(FL3(tmp2))*pnorm+FL3(n1/n2)*(*c0);
              //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(("Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal));
	  }
       }else{ /*total internal reflection*/
          *c0+=(FL3(-2.f*Icos))*pnorm;
	  *eid=*oldeid;
          //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(("V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f));
       }
       tmp0=MCX_MATHFUN(rsqrt)(dot(*c0,*c0));
       (*c0)*=FL3(tmp0);
       return 1.f;
}

#endif

/**
 * @brief Performing one scattering event of the photon
 *
 * This function updates the direction of the photon by performing a scattering calculation
 *
 * @param[in] g: anisotropy g
 * @param[out] dir: current ray direction vector
 * @param[out] ran: random number generator states
 * @param[out] cfg: the simulation configuration
 * @param[out] pmom: buffer to store momentum transfer data if needed
 */

float mc_next_scatter(float g, float3 *dir,__private RandType *ran, __constant MCXParam *gcfg, float *pmom){

    float nextslen;
    float sphi,cphi,tmp0,theta,stheta,ctheta,tmp1;
    float3 p;

    //random scattering length (normalized)
    nextslen=rand_next_scatlen(ran);

    tmp0=TWO_PI*rand_next_aangle(ran); //next arimuth angle
    MCX_SINCOS(tmp0,sphi,cphi);

    if(g>EPS){  //if g is too small, the distribution of theta is bad
	tmp0=(1.f-g*g)/(1.f-g+2.f*g*rand_next_zangle(ran));
	tmp0*=tmp0;
	tmp0=(1.f+g*g-tmp0)/(2.f*g);
	tmp0=clamp(tmp0,-1.f,1.f);

	theta=acos(tmp0);
	stheta=MCX_MATHFUN(sqrt)(1.f-tmp0*tmp0);
	//stheta=MCX_MATHFUN(sin)(theta);
	ctheta=tmp0;
    }else{
	theta=acos(2.f*rand_next_zangle(ran)-1.f);
    	MCX_SINCOS(theta,stheta,ctheta);
    }

    if( dir->z>-1.f+EPS && dir->z<1.f-EPS ) {
	tmp0=1.f-dir->z*dir->z;   //reuse tmp to minimize registers
	tmp1=MCX_MATHFUN(rsqrt)(tmp0);
	tmp1=stheta*tmp1;

	p.x=tmp1*(dir->x*dir->z*cphi - dir->y*sphi) + dir->x*ctheta;
	p.y=tmp1*(dir->y*dir->z*cphi + dir->x*sphi) + dir->y*ctheta;
	p.z=-tmp1*tmp0*cphi			    + dir->z*ctheta;
    }else{
	p.x=stheta*cphi;
	p.y=stheta*sphi;
	p.z=(dir->z>0.f)?ctheta:-ctheta;
    }
    if (gcfg->ismomentum)
        pmom[0]+=(1.f-ctheta);

    dir->x=p.x;
    dir->y=p.y;
    dir->z=p.z;
    return nextslen;
}

/** 
 * \brief Function to deal with ray-edge/ray-vertex intersections
 * 
 * when a photon is crossing a vertex or edge, (slightly) pull the
 * photon toward the center of the element and try again
 *
 * \param[in,out] p: current photon position
 * \param[in] nodes: pointer to the 4 nodes of the tet
 * \param[in] ee: indices of the 4 nodes ee=elem[eid]
 */

void fixphoton(float3 *p,__global float3 *nodes, __constant int *ee){
        float3 c0={0.f,0.f,0.f};
	int i;
        /*calculate element centroid*/
	for(i=0;i<4;i++)
		c0+=nodes[ee[i]-1];
	*p+=(c0*FL3(0.25f)-*p)*(FL3(FIX_PHOTON));
}



/**
 * @brief Launch a new photon
 *
 * This function launch a new photon using one of the dozen supported source forms.
 *
 * \param[in,out] cfg: simulation configuration structure
 * \param[in,out] r: the current ray
 * \param[in] mesh: the mesh data structure
 * \param[in,out] ran: the random number generator states
 */

void launchphoton(__constant MCXParam *gcfg, ray *r, __global float3 *node,__constant int *elem,__constant int *srcelem, __private RandType *ran){
        int canfocus=1;
        float3 origin=r->p0;

	r->slen=rand_next_scatlen(ran);
#if defined(MCX_SRC_PENCIL)
		if(r->eid>0)
		      return;
#elif defined(MCX_SRC_PLANAR) || defined(MCX_SRC_PATTERN) || defined(MCX_SRC_PATTERN3D) || defined(MCX_SRC_FOURIER) /*a rectangular grid over a plane*/
		  float rx=rand_uniform01(ran);
		  float ry=rand_uniform01(ran);
		  r->p0.x=gcfg->srcpos.x+rx*gcfg->srcparam1.x+ry*gcfg->srcparam2.x;
		  r->p0.y=gcfg->srcpos.y+rx*gcfg->srcparam1.y+ry*gcfg->srcparam2.y;
		  r->p0.z=gcfg->srcpos.z+rx*gcfg->srcparam1.z+ry*gcfg->srcparam2.z;
		  r->weight=1.f;
    #if defined(MCX_SRC_PATTERN) 
		  if(gcfg->srctype==stPattern){
		    int xsize=(int)gcfg->srcparam1.w;
		    int ysize=(int)gcfg->srcparam2.w;
		    //r->posidx=MIN((int)(ry*ysize),ysize-1)*xsize+MIN((int)(rx*xsize),xsize-1);
    #elif defined(MCX_SRC_FOURIER)  // need to prevent rx/ry=1 here
		    r->weight=(MCX_MATHFUN(cos)((floor(gcfg->srcparam1.w)*rx+floor(gcfg->srcparam2.w)*ry+gcfg->srcparam1.w-floor(gcfg->srcparam1.w))*TWO_PI)*(1.f-gcfg->srcparam2.w+floor(gcfg->srcparam2.w))+1.f)*0.5f;
    #endif
		origin.x+=(gcfg->srcparam1.x+gcfg->srcparam2.x)*0.5f;
		origin.y+=(gcfg->srcparam1.y+gcfg->srcparam2.y)*0.5f;
		origin.z+=(gcfg->srcparam1.z+gcfg->srcparam2.z)*0.5f;
#elif defined(MCX_SRC_FOURIERX) || defined(MCX_SRC_FOURIERX2D) // [v1x][v1y][v1z][|v2|]; [kx][ky][phi0][M], unit(v0) x unit(v1)=unit(v2)
		float rx=rand_uniform01(ran);
		float ry=rand_uniform01(ran);
		float4 v2=gcfg->srcparam1;
		v2.w*=MCX_MATHFUN(rsqrt)(gcfg->srcparam1.x*gcfg->srcparam1.x+gcfg->srcparam1.y*gcfg->srcparam1.y+gcfg->srcparam1.z*gcfg->srcparam1.z);
		v2.x=v2.w*(gcfg->srcdir.y*gcfg->srcparam1.z - gcfg->srcdir.z*gcfg->srcparam1.y);
		v2.y=v2.w*(gcfg->srcdir.z*gcfg->srcparam1.x - gcfg->srcdir.x*gcfg->srcparam1.z); 
		v2.z=v2.w*(gcfg->srcdir.x*gcfg->srcparam1.y - gcfg->srcdir.y*gcfg->srcparam1.x);
		r->p0.x=gcfg->srcpos.x+rx*gcfg->srcparam1.x+ry*v2.x;
		r->p0.y=gcfg->srcpos.y+rx*gcfg->srcparam1.y+ry*v2.y;
		r->p0.z=gcfg->srcpos.z+rx*gcfg->srcparam1.z+ry*v2.z;
    #if defined(MCX_SRC_FOURIERX2D)
			r->weight=(MCX_MATHFUN(sin)((gcfg->srcparam2.x*rx+gcfg->srcparam2.z)*TWO_PI)*MCX_MATHFUN(sin)((gcfg->srcparam2.y*ry+gcfg->srcparam2.w)*TWO_PI)+1.f)*0.5f; //between 0 and 1
    #else
			r->weight=(MCX_MATHFUN(cos)((gcfg->srcparam2.x*rx+gcfg->srcparam2.y*ry+gcfg->srcparam2.z)*TWO_PI)*(1.f-gcfg->srcparam2.w)+1.f)*0.5f; //between 0 and 1
    #endif
		origin.x+=(gcfg->srcparam1.x+v2.x)*0.5f;
		origin.y+=(gcfg->srcparam1.y+v2.y)*0.5f;
		origin.z+=(gcfg->srcparam1.z+v2.z)*0.5f;
#elif defined(MCX_SRC_DISK) || defined(MCX_SRC_GAUSSIAN) // uniform disk distribution or Gaussian-beam
		float sphi, cphi;
		float phi=TWO_PI*rand_uniform01(ran);
		sphi=MCX_MATHFUN(sin)(phi);	cphi=MCX_MATHFUN(cos)(phi);
		float r0;
    #if defined(MCX_SRC_DISK)
		    r0=MCX_MATHFUN(sqrt)(rand_uniform01(ran))*gcfg->srcparam1.x;
    #else
		if(fabs(gcfg->focus) < 1e-5f || fabs(gcfg->srcparam1.y) < 1e-5f)
		    r0=MCX_MATHFUN(sqrt)(-MCX_MATHFUN(log)((rand_uniform01(ran)))*gcfg->srcparam1.x;
		else{
		    float z0=gcfg->srcparam1.x*gcfg->srcparam1.x*M_PI/gcfg->srcparam1.y; //Rayleigh range
		    r0=MCX_MATHFUN(sqrt)(-MCX_MATHFUN(log)((rand_uniform01(ran))*(1.f+(gcfg->focus*gcfg->focus/(z0*z0))))*gcfg->srcparam1.x;
		}
    #endif

		if(gcfg->srcdir.z>-1.f+EPS && gcfg->srcdir.z<1.f-EPS){
		    float tmp0=1.f-gcfg->srcdir.z*gcfg->srcdir.z;
		    float tmp1=r0*MCX_MATHFUN(rsqrt)(tmp0);
		    r->p0.x=gcfg->srcpos.x+tmp1*(gcfg->srcdir.x*gcfg->srcdir.z*cphi - gcfg->srcdir.y*sphi);
		    r->p0.y=gcfg->srcpos.y+tmp1*(gcfg->srcdir.y*gcfg->srcdir.z*cphi + gcfg->srcdir.x*sphi);
		    r->p0.z=gcfg->srcpos.z-tmp1*tmp0*cphi;
		}else{
   		    r->p0.x+=r0*cphi;
		    r->p0.y+=r0*sphi;
		}
#elif defined(MCX_SRC_CONE) || defined(MCX_SRC_ISOTROPIC) || defined(MCX_SRC_ARCSINE) 
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=MCX_MATHFUN(sin)(ang);	cphi=MCX_MATHFUN(cos)(ang);
    #if defined(MCX_SRC_CONE) // a solid-angle section of a uniform sphere
		        do{
				ang=(gcfg->srcparam1.y>0) ? TWO_PI*rand_uniform01(ran) : acos(2.f*rand_uniform01(ran)-1.f); //sine distribution
		        }while(ang>gcfg->srcparam1.x);
    #else
			if(gcfg->srctype==stIsotropic) // uniform sphere
				ang=acos(2.f*rand_uniform01(ran)-1.f); //sine distribution
			else
				ang=M_PI*rand_uniform01(ran); //uniform distribution in zenith angle, arcsine
    #endif
		stheta=MCX_MATHFUN(sin)(ang);
		ctheta=MCX_MATHFUN(cos)(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
		canfocus=0;
                if(gcfg->srctype==stIsotropic)
                    if(r->eid>0)
                        return;
#elif defined(MCX_SRC_ZGAUSSIAN)
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=MCX_MATHFUN(sin)(ang);	cphi=MCX_MATHFUN(cos)(ang);
		ang=MCX_MATHFUN(sqrt)(-2.f*MCX_MATHFUN(log)((rand_uniform01(ran)))*(1.f-2.f*rand_uniform01(ran))*gcfg->srcparam1.x;
		stheta=MCX_MATHFUN(sin)(ang);
		ctheta=MCX_MATHFUN(cos)(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
		canfocus=0;
#elif defined(MCX_SRC_LINE) || defined(MCX_SRC_SLIT) 
	      float t=rand_uniform01(ran);
	      r->p0.x+=t*gcfg->srcparam1.x;
	      r->p0.y+=t*gcfg->srcparam1.y;
	      r->p0.z+=t*gcfg->srcparam1.z;

    #if defined(MCX_SRC_LINE)
	              float s,p;
		      t=1.f-2.f*rand_uniform01(ran);
		      s=1.f-2.f*rand_uniform01(ran);
		      p=MCX_MATHFUN(sqrt)(1.f-r->vec.x*r->vec.x-r->vec.y*r->vec.y)*(rand_uniform01(ran)>0.5f ? 1.f : -1.f);
		      float3 vv;
		      vv.x=r->vec.y*p-r->vec.z*s;
		      vv.y=r->vec.z*t-r->vec.x*p;
		      vv.z=r->vec.x*s-r->vec.y*t;
		      r->vec=vv;
		      //*((float3*)&(r->vec))=(float3)(r->vec.y*p-r->vec.z*s,r->vec.z*t-r->vec.x*p,r->vec.x*s-r->vec.y*t);
    #endif
              origin.x+=(gcfg->srcparam1.x)*0.5f;
              origin.y+=(gcfg->srcparam1.y)*0.5f;
              origin.z+=(gcfg->srcparam1.z)*0.5f;
              canfocus=(gcfg->srctype==stSlit);
#endif

        if(canfocus && gcfg->focus!=0.f){ // if beam focus is set, determine the incident angle
	        float Rn2;
	        origin.x+=gcfg->focus*r->vec.x;
		origin.y+=gcfg->focus*r->vec.y;
		origin.z+=gcfg->focus*r->vec.z;
		if(gcfg->focus<0.f){ // diverging beam
                     r->vec.x=r->p0.x-origin.x;
                     r->vec.y=r->p0.y-origin.y;
                     r->vec.z=r->p0.z-origin.z;
		}else{             // converging beam
                     r->vec.x=origin.x-r->p0.x;
                     r->vec.y=origin.y-r->p0.y;
                     r->vec.z=origin.z-r->p0.z;
		}
		Rn2=MCX_MATHFUN(rsqrt)(dot(r->vec,r->vec)); // normalize
		r->vec=r->vec*Rn2;
	}

        r->p0+=r->vec*EPS;

#if defined(MCX_SRC_PLANAR) || defined(MCX_SRC_PATTERN) || defined(MCX_SRC_PATTERN3D) || defined(MCX_SRC_FOURIER) || defined(MCX_SRC_FOURIERX) || defined(MCX_SRC_FOURIERX2D)	
	/*Caluclate intial element id and bary-centric coordinates for area sources - position changes everytime*/
	float3 vecS=FL3(0.f), vecAB, vecAC, vecN;
	int is,i,ea,eb,ec;
	float bary[4]={0.f};
	for(is=0;is<gcfg->srcelemlen;is++){
		int include = 1;
		constant int *elems=elem+(srcelem[is]-1)*gcfg->elemlen;
		for(i=0;i<4;i++){
			ea=elems[out[i][0]]-1;
			eb=elems[out[i][1]]-1;
			ec=elems[out[i][2]]-1;
			vecAB=node[eb]-node[ea];
			vecAC=node[ec]-node[ea];
			vecS=r->p0-node[ea];
			vecN=cross(vecAB,vecAC);
			bary[facemap[i]]=-dot(vecS,vecN);	
		}
		for(i=0;i<4;i++){
			if(bary[i]<-1e-4f){
				include = 0;
			}
		}
		if(include){
			r->eid=srcelem[is];
			float s=0.f;
			for(i=0;i<4;i++){s+=bary[i];}
			for(i=0;i<4;i++){
				if((bary[i]/s)<1e-4f)
					r->faceid=ifacemap[i]+1;
			}
			break;
		}
	}
#endif		
}


/**
 * @brief The core Monte Carlo function simulating a single photon (!!!Important!!!)
 *
 * This is the core Monte Carlo simulation function. It simulates the life-time
 * of a single photon packet, from launching to termination.
 *
 * \param[in] id: the linear index of the current photon, starting from 0.
 * \param[in] tracer: the ray-tracer aux data structure
 * \param[in] mesh: the mesh data structure
 * \param[in,out] ran: the random number generator states
 * \param[in,out] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

void onephoton(unsigned int id,__local float *ppath, __constant MCXParam *gcfg,__global float3 *node,__constant int *elem, __global float *weight,__global float *dref,
    __constant int *type, __constant int *facenb,  __constant int *srcelem, __constant float4 *normal, __constant medium *med,
    __global float *n_det, __global uint *detectedphoton, float *energytot, float *energyesc, __constant float4 *gdetpos, __private RandType *ran, int *raytet){

	int oldeid,fixcount=0;
	ray r={gcfg->srcpos,gcfg->srcdir,{MMC_UNDEFINED,0.f,0.f},gcfg->e0,0,0,1.f,0.f,0.f,0.f,ID_UNDEFINED,0.f};

	//r.photonid=id;

	/*initialize the photon parameters*/
        launchphoton(gcfg, &r, node, elem, srcelem, ran);
	*energytot+=r.weight;
#ifdef MCX_SAVE_DETECTORS
	if(gcfg->issavedet)
	    ppath[gcfg->reclen-1] = r.weight; /*last record in partialpath is the initial photon weight*/
#endif
	/*use Kahan summation to accumulate weight, otherwise, counter stops at 16777216*/
	/*http://stackoverflow.com/questions/2148149/how-to-sum-a-large-number-of-float-number*/

	while(1){  /*propagate a photon until exit*/
	    r.slen=branchless_badouel_raytet(&r, gcfg, elem, weight, type[r.eid-1], facenb, normal, med);
	    (*raytet)++;
	    if(r.pout.x==MMC_UNDEFINED){
	    	  if(r.faceid==-2) break; /*reaches the time limit*/
		  if(fixcount++<MAX_TRIAL){
			fixphoton(&r.p0,node,(__constant int *)(elem+(r.eid-1)*gcfg->elemlen));
			continue;
                  }
	    	  r.eid=ID_UNDEFINED;
        	  r.faceid=-1;
	    }
#ifdef MCX_SAVE_DETECTORS
	    if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
	            ppath[gcfg->maxmedia+type[r.eid-1]-1]+=r.Lmove;  /*second medianum block is the partial path*/
#endif
	    /*move a photon until the end of the current scattering path*/
	    while(r.faceid>=0 && !r.isend){
	            r.p0=r.pout;

		    oldeid=r.eid;
	    	    r.eid=((__constant int *)(facenb+(r.eid-1)*gcfg->elemlen))[r.faceid];
#ifdef MCX_DO_REFLECTION
		    if(gcfg->isreflect && (r.eid<=0 || (r.eid>0 && med[type[r.eid-1]].n != med[type[oldeid-1]].n ))){
			if(! (r.eid<=0 && med[type[oldeid-1]].n == gcfg->nout ))
			    reflectray(gcfg,&r.vec,&oldeid,&r.eid,r.faceid,ran,type,normal,med);
		    }
#endif
	    	    if(r.eid<=0) break;
		    /*when a photon enters the domain from the background*/
		    if(type[oldeid-1]==0 && type[r.eid-1]){
                        //if(gcfg->debuglevel&dlExit)
			    MMC_FPRINTF(("e %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			    	r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid));
                        if(!gcfg->voidtime)
                            r.photontimer=0.f;
                    }
		    /*when a photon exits the domain into the background*/
		    if(type[oldeid-1] && type[r.eid-1]==0){
                        //if(gcfg->debuglevel&dlExit)
		            MMC_FPRINTF(("x %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			        r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid));
			if(!gcfg->isextdet){
                            r.eid=0;
        		    break;
                        }
		    }
//		    if(r.eid==0 && med[type[oldeid-1]].n == gcfg->nout ) break;
	    	    if(r.pout.x!=MMC_UNDEFINED) // && (gcfg->debuglevel&dlMove))
	    		MMC_FPRINTF(("P %f %f %f %d %u %e\n",r.pout.x,r.pout.y,r.pout.z,r.eid,id,r.slen));

	    	    r.slen=branchless_badouel_raytet(&r,gcfg, elem, weight, type[r.eid-1], facenb, normal, med);
		    (*raytet)++;
#ifdef MCX_SAVE_DETECTORS
		    if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
			ppath[gcfg->maxmedia+type[r.eid-1]-1]+=r.Lmove;
#endif
		    if(r.faceid==-2) break;
		    fixcount=0;
		    while(r.pout.x==MMC_UNDEFINED && fixcount++<MAX_TRIAL){
		       fixphoton(&r.p0,node,(__constant int *)(elem+(r.eid-1)*gcfg->elemlen));
                       r.slen=branchless_badouel_raytet(&r,gcfg, elem, weight, type[r.eid-1], facenb, normal, med);
		       (*raytet)++;
#ifdef MCX_SAVE_DETECTORS
		       if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
	            		ppath[gcfg->maxmedia+type[r.eid-1]-1]+=r.Lmove;
#endif
		    }
        	    if(r.pout.x==MMC_UNDEFINED){
        		/*possibily hit an edge or miss*/
			r.eid=ID_UNDEFINED;
        		break;
        	    }
	    }
	    if(r.eid<=0 || r.pout.x==MMC_UNDEFINED) {
        	    //if(r.eid==0 && (gcfg->debuglevel&dlMove))
        		 MMC_FPRINTF(("B %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen));
		    if(r.eid!=ID_UNDEFINED){
                       //if(gcfg->debuglevel&dlExit)
        		 MMC_FPRINTF(("E %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			    r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid));
#ifdef MCX_SAVE_DETECTORS
                       if(gcfg->issavedet && gcfg->issaveexit){                                     /*when issaveexit is set to 1*/
                            copystate(ppath+(gcfg->reclen-7),(float *)&(r.p0),3);  /*columns 7-5 from the right store the exit positions*/
                            copystate(ppath+(gcfg->reclen-4),(float *)&(r.vec),3); /*columns 4-2 from the right store the exit dirs*/
                       }
#endif
#ifdef MCX_SAVE_DREF
                       if(gcfg->issaveref && r.eid<0 && dref){
		            int tshift=MIN( ((int)((r.photontimer-gcfg->tstart)*gcfg->Rtstep)), gcfg->maxgate-1 )*gcfg->nf;
                            dref[((-r.eid)-1) + tshift]+=r.weight;
		       }
#endif
		    }else if(r.faceid==-2 && (gcfg->debuglevel&dlMove)){
                         MMC_FPRINTF(("T %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen));
	    	    }else if(r.eid && r.faceid!=-2  && gcfg->debuglevel&dlEdge){
        		 MMC_FPRINTF(("X %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen));
		    }
#ifdef MCX_SAVE_DETECTORS
		    if(r.eid<0){
                       if(gcfg->isextdet && type[oldeid-1]==gcfg->maxmedia+1){
                          savedetphoton(n_det,detectedphoton,ppath,&(r.p0),&(r.vec),gdetpos,oldeid,gcfg);
                       }else{
                          savedetphoton(n_det,detectedphoton,ppath,&(r.p0),&(r.vec),gdetpos,-1,gcfg);
		       }
                       clearpath(ppath,gcfg->reclen);
		    }
#endif
	    	    break;  /*photon exits boundary*/
	    }
	    //if(gcfg->debuglevel&dlMove)
	        MMC_FPRINTF(("M %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen));
	    if(gcfg->minenergy>0.f && r.weight < gcfg->minenergy && (gcfg->tend-gcfg->tstart)*gcfg->Rtstep<=1.f){ /*Russian Roulette*/
		if(rand_do_roulette(ran)*gcfg->roulettesize<=1.f){
			r.weight*=gcfg->roulettesize;
                        //if(gcfg->debuglevel&dlWeight)
			    MMC_FPRINTF(("Russian Roulette bumps r.weight to %f\n",r.weight));
		}else
			break;
	    }
            float mom=0.f;
	    r.slen=mc_next_scatter(med[type[r.eid-1]].g,&r.vec,ran,gcfg,&mom);
#ifdef MCX_SAVE_DETECTORS
            if(gcfg->ismomentum && type[r.eid-1]>0)                     /*when ismomentum is set to 1*/
                  ppath[(gcfg->maxmedia<<1)+type[r.eid-1]-1]+=mom; /*the third medianum block stores the momentum transfer*/
            if(gcfg->issavedet)
	          ppath[type[r.eid-1]-1]+=1.f;                          /*the first medianum block stores the scattering event counts*/
#endif
	}
	*energyesc+=r.weight;
}

__kernel void mmc_main_loop(const int nphoton, const int ophoton, __constant MCXParam *gcfg,__local float *sharedmem,
    __global float3 *node,__constant int *elem,  __global float *weight, __global float *dref, __constant int *type, __constant int *facenb,  __constant int *srcelem, __constant float4 *normal, 
    __constant medium *med,  __constant float4 *gdetpos,__global float *n_det, __global uint *detectedphoton, 
    __global uint *n_seed, __global int *progress, __global float *energy, __global MCXReporter *reporter){
 
 	RandType t[RAND_BUF_LEN];
	int idx=get_global_id(0);
	gpu_rng_init(t,n_seed,idx);
        float  energyesc=0.f, energytot=0.f;
	int raytet=0;

	/*launch photons*/
	for(int i=0;i<nphoton+(idx<ophoton);i++){
	    onephoton(idx*nphoton+MIN(idx,ophoton)+i,sharedmem+get_local_id(0)*gcfg->reclen,gcfg,node,elem,
	        weight,dref,type,facenb,srcelem, normal,med,n_det,detectedphoton,&energytot,&energyesc,gdetpos,t,&raytet);
	}
	energy[idx<<1]=energyesc;
	energy[1+(idx<<1)]=energytot;

        if(gcfg->debuglevel & MCX_DEBUG_PROGRESS)
            atomic_inc(progress);
	atomicadd(&(reporter->raytet),raytet);
}
