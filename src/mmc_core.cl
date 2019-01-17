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
#define RAND_MAX           4294967295

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
#define OUTSIDE_VOLUME     0xFFFFFFFF              /**< flag indicating the index is outside of the volume */

#define DET_MASK           0xFFFF0000
#define MED_MASK           0x0000FFFF
#define NULL               0
#define MAX_ACCUM          1000.f

#define MCX_DEBUG_RNG       1                   /**< MCX debug flags */
#define MCX_DEBUG_MOVE      2
#define MCX_DEBUG_PROGRESS  4

#define MIN(a,b)           ((a)<(b)?(a):(b))

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
	float focus;                  /**< focal length of the source, if defined */
	unsigned int posidx;	      /**< launch position index of the photon for pattern source type */
} ray __attribute__ ((aligned (32)));

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
} visitor __attribute__ ((aligned (32)));

typedef struct KernelParams {
  float4 ps,c0;
  float  twin0,twin1,tmax;
  uint   save2pt,doreflect,dorefint,savedet;
  float  Rtstep;
  float  minenergy;
  uint   maxdetphoton;
  uint   maxmedia;
  uint   detnum;
  uint   idx1dorig;
  uint   mediaidorig;
  uint   blockphoton;
  uint   blockextra;
  int    voidtime;
  int    srctype;                    /**< type of the source */
  float4 srcparam1;                  /**< source parameters set 1 */
  float4 srcparam2;                  /**< source parameters set 2 */
  uint   maxvoidstep;
  uint   issaveexit;    /**<1 save the exit position and dir of a detected photon, 0 do not save*/
  uint   issaveref;     /**<1 save diffuse reflectance at the boundary voxels, 0 do not save*/
  uint   maxgate;
  uint   threadphoton;                  /**< how many photons to be simulated in a thread */
  uint   debuglevel;           /**< debug flags */
} MCXParam __attribute__ ((aligned (32)));

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define RAND_BUF_LEN       2        //register arrays
#define RAND_SEED_LEN      4        //48 bit packed with 64bit length
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

static void copystate(__private RandType t[RAND_BUF_LEN], __private RandType tnew[RAND_BUF_LEN]){
    tnew[0]=t[0];
    tnew[1]=t[1];
}

static float rand_uniform01(__private RandType t[RAND_BUF_LEN]){
    return xorshift128p_nextf(t);
}

static void xorshift128p_seed (__global uint *seed,RandType t[RAND_BUF_LEN]){
    t[0] = (ulong)seed[0] << 32 | seed[1] ;
    t[1] = (ulong)seed[2] << 32 | seed[3];
}

static void gpu_rng_init(__private RandType t[RAND_BUF_LEN], __global uint *n_seed, int idx){
    xorshift128p_seed((n_seed+idx*RAND_SEED_LEN),t);
}

float rand_next_scatlen(__private RandType t[RAND_BUF_LEN]){
    return -MCX_MATHFUN(log)(rand_uniform01(t)+EPS);
}

#define rand_next_aangle(t)  rand_uniform01(t)
#define rand_next_zangle(t)  rand_uniform01(t)
#define rand_next_reflect(t) rand_uniform01(t)
#define rand_do_roulette(t)  rand_uniform01(t) 

#ifdef USE_ATOMIC
// OpenCL float atomicadd hack:
// http://suhorukov.blogspot.co.uk/2011/12/opencl-11-atomic-operations-on-floating.html
// https://devtalk.nvidia.com/default/topic/458062/atomicadd-float-float-atomicmul-float-float-/

inline float atomicadd(volatile __global float* address, const float value){
    float old = value;
    while ((old = atomic_xchg(address, atomic_xchg(address, 0.0f)+old))!=0.0f);
    return old;
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

float branchless_badouel_raytet(ray *r,visitor *visit, __constant MCXParam *gcfg,__global float3 *node,__global int4 *elem,
    __global int *type, __global int4 *facenb, __global float3 *normal, __global float4 *med){

	float4 bary={1e10f,0.f,0.f,0.f};
	float Lp0=0.f,rc,currweight,dlen,ww,totalloss=0.f;
	int tshift,faceidx=-1,baseid,eid;
	float3 O,T,S;
	int3 P;

	r->p0.w=1.f;
	r->vec.w=0.f;
	eid=r->eid-1;
	baseid=eid<<2;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;

	T = normal[baseid+3] - ((float3)(r->p0.x)*normal[baseid]+(float3)(r->p0.y)*normal[baseid+1]+(float3)(r->p0.z)*normal[baseid+2]);
	S = ((float3)(r->vec.x)*normal[baseid]+(float3)(r->vec.y)*normal[baseid+1]+(float3)(r->vec.z)*normal[baseid+2]);
	T = T/S;

	O = isgreaterequal(S,(float3)(0.f));
	T = isless(S,(float3)(0.f))*(float3)(1e10f))+(O*T);

	bary.x=fmin(fmin(T.x,T.y),T.z);
	faceidx=(bary.x==T.x? 0: (bary.x==T.y? 1 : 2));
	r->faceid=faceorder[faceidx];

	if(r->faceid>=0 && bary.x>=0){
	    medium prop;
	    int *ee=(int *)(elem+eid*gcfg->elemlen);
	    prop=med[type[eid]];
	    rc=prop.n*R_C0;
            currweight=r->weight;

            r->nexteid=((int *)(facenb+eid*gcfg->elemlen))[r->faceid]; // if I use nexteid-1, the speed got slower, strange!

	    dlen=(prop.mus <= EPS) ? R_MIN_MUS : r->slen/prop.mus;
	    Lp0=bary.x;
	    r->isend=(Lp0>dlen);
	    r->Lmove=((r->isend) ? dlen : Lp0);

            O = r->vec;
	    S = r->p0;
	    T = (float3)(bary.x);
	    T = O*T;
	    T = T+S;
	    r->pout=T;

	    if((int)((r->photontimer+r->Lmove*rc-gcfg->tstart)*visit->rtstep)>=(int)((gcfg->tend-gcfg->tstart)*visit->rtstep)){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(gcfg->tend-r->photontimer)/(prop.n*R_C0)-1e-4f;
	    }
            if(gcfg->mcmethod==mmMCX){
	       totalloss=exp(-prop.mua*r->Lmove);
               r->weight*=totalloss;
            }
	    totalloss=1.f-totalloss;
	    if(gcfg->seed==SEED_FROM_FILE && gcfg->outputtype==otJacobian){
#ifdef __INTEL_COMPILER
		currweight=exp(-DELTA_MUA*r->Lmove);
#else
		currweight=fast_expf9(-DELTA_MUA*r->Lmove);
#endif
                currweight*=replayweight[r->photonid];
		currweight+=r->weight;
	    }else if(gcfg->seed==SEED_FROM_FILE && gcfg->outputtype==otWL){
		currweight=r->Lmove;
		currweight*=replayweight[r->photonid];
		currweight+=r->weight;
            }
	    r->slen-=r->Lmove*prop.mus;
	    if(gcfg->seed==SEED_FROM_FILE && gcfg->outputtype==otWP){
		if(r->slen0<EPS)
		    currweight=1;
		else
		    currweight=r->Lmove*prop.mus/r->slen0;
		currweight*=replayweight[r->photonid];
		currweight+=r->weight;
	    }
	    if(bary.x>=0.f){
	        int framelen=(gcfg->basisorder?gcfg->nn:gcfg->ne);
		if(gcfg->method==rtBLBadouelGrid)
		    framelen=gcfg->crop0.z;
		ww=currweight-r->weight;
        	r->photontimer+=r->Lmove*rc;

		if(gcfg->outputtype==otWL || gcfg->outputtype==otWP)
			tshift=MIN( ((int)(gcfg->replaytime[r->photonid]*visit->rtstep)), gcfg->maxgate-1 )*framelen;
		else
                	tshift=MIN( ((int)((r->photontimer-gcfg->tstart)*visit->rtstep)), gcfg->maxgate-1 )*framelen;

        	if(gcfg->debuglevel&dlAccum) MMC_FPRINTF(gcfg->flog,"A %f %f %f %e %d %e\n",
        	   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

		if(prop.mua>0.f){
		  r->Eabsorb+=ww;
		  if(gcfg->outputtype==otFlux || gcfg->outputtype==otJacobian)
                     ww/=prop.mua;
		}

	        T = (float3)(r->Lmove);
	        T = S+(O*T);
	        r->p0=T;

                if(gcfg->mcmethod==mmMCX){
		  if(!gcfg->basisorder){
		     if(gcfg->method==rtBLBadouel){
#ifdef USE_ATOMIC
                        if(gcfg->isatomic)
			    atomicadd(weight+eid+tshift,ww);
                        else
#endif
                            weight[eid+tshift]+=ww;
                     }else{
			    float dstep, segloss, w0;
			    int3 idx;
			    int i, seg=(int)(r->Lmove/gcfg->unitinmm)+1;
			    seg=(seg<<1);
			    dstep=r->Lmove/seg;
	                    segloss=exp(-prop.mua*dstep);
			    T =  O * (float3)(dstep)); /*step*/
			    O =  S - gcfg->nmin;
			    S =  O + (T * (float3)(0.5f)); /*starting point*/
			    dstep=1.f/gcfg->unitinmm;
			    totalloss=(totalloss==0.f)? 0.f : (1.f-segloss)/totalloss;
			    w0=ww;
                            for(i=0; i< seg; i++){
				idx= convert_int3_rtn(S * (float3)(dstep));
				weight[idx.z*gcfg->crop0.y+idx.y*gcfg->crop0.x+idx.x+tshift]+=w0*totalloss;
				w0*=segloss;
			        S += T;
                            }
			}
		  }else{
			int i;
                        ww*=1.f/3.f;
#ifdef USE_ATOMIC
                        if(gcfg->isatomic)
			    for(i=0;i<3;i++)
				atomicadd(weight+ee[out[faceidx][i]]-1+tshift, ww);
                        else
#endif
                            for(i=0;i<3;i++)
                                weight[ee[out[faceidx][i]]-1+tshift]+=ww;
		  }
		}
	    }
	}
	visit->raytet++;
        if(type[eid]==0)
                visit->raytet0++;

	r->p0.w=0.f;
	if(r->faceid==-2)
           return 0.f;

	return r->slen;
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

float reflectray(mcconfig *cfg,float3 *c0,raytracer *tracer,int *oldeid,int *eid,int faceid,RandType *ran){
	/*to handle refractive index mismatch*/
        float3 pnorm={0.f}, *pn=&pnorm;
	float Icos,Re,Im,Rtotal,tmp0,tmp1,tmp2,n1,n2;
	int offs=(*oldeid-1)<<2;

	faceid=ifaceorder[faceid];
	/*calculate the normal direction of the intersecting triangle*/
/*	if(gcfg->method==rtPlucker) { //Plucker ray-tracing
		pn=normal+(offs)+faceid;
	}else if(gcfg->method<rtBLBadouel){
		pn=moment+(offs+faceid)*3;
	}else 
*/
	  if(gcfg->method==rtBLBadouel || gcfg->method==rtBLBadouelGrid){
		pnorm.x=(&(normal[offs].x))[faceid];
		pnorm.y=(&(normal[offs].x))[faceid+4];
		pnorm.z=(&(normal[offs].x))[faceid+8];
	}
	/*pn pointing outward*/

	/*compute the cos of the incidence angle*/
        Icos=fabs(dot(c0,pn));

	n1=(*oldeid!=*eid) ? med[type[*oldeid-1]].n : gcfg->nout;
	n2=(*eid>0) ? med[type[*eid-1]].n : gcfg->nout;

	tmp0=n1*n1;
	tmp1=n2*n2;
        tmp2=1.f-tmp0/tmp1*(1.f-Icos*Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

        if(tmp2>0.f){ /*if no total internal reflection*/
          Re=tmp0*Icos*Icos+tmp1*tmp2;      /*transmission angle*/
	  tmp2=sqrt(tmp2); /*to save one sqrt*/
          Im=2.f*n1*n2*Icos*tmp2;
          Rtotal=(Re-Im)/(Re+Im);     /*Rp*/
          Re=tmp1*Icos*Icos+tmp0*tmp2*tmp2;
          Rtotal=(Rtotal+(Re-Im)/(Re+Im))*0.5f; /*(Rp+Rs)/2*/
	  if(*oldeid==*eid) return Rtotal; /*initial specular reflection*/
	  if(rand_next_reflect(ran)<=Rtotal){ /*do reflection*/
              c0+=((float3)(-2.f*Icos))*pn;
              //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(gcfg->flog,"R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal);
	      *eid=*oldeid; /*stay with the current element*/
	  }else if(gcfg->isspecular==2 && *eid==0){
              // if do transmission, but next neighbor is 0, terminate
          }else{                              /*do transmission*/
              c0+=((float3)(-Icos))*pn;
	      c0=((float3)(tmp2))*pn+(float3)(n1/n2)*c0;
              //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(gcfg->flog,"Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal);
	  }
       }else{ /*total internal reflection*/
          c0+=((float3)(-2.f*Icos))*pn;
	  *eid=*oldeid;
          //if(gcfg->debuglevel&dlReflect) MMC_FPRINTF(gcfg->flog,"V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f);
       }
       tmp0=rsqrt(dot(c0,c0));
       c0*=tmp0;
       return 1.f;
}

/**
 * @brief Performing one scattering event of the photon
 *
 * This function updates the direction of the photon by performing a scattering calculation
 *
 * @param[in] g: anisotropy g
 * @param[out] dir: current ray direction vector
 * @param[out] ran: random number generator states
 * @param[out] ran0: additional random number generator states
 * @param[out] cfg: the simulation configuration
 * @param[out] pmom: buffer to store momentum transfer data if needed
 */

float mc_next_scatter(float g, float3 *dir,RandType *ran, RandType *ran0, __constant MCXParam *gcfg, float *pmom){

    float nextslen;
    float sphi,cphi,tmp0,theta,stheta,ctheta,tmp1;
    float3 p;

    //random scattering length (normalized)
#ifdef MMC_USE_SSE_MATH
    nextslen=rand_next_scatlen_ps(ran);
#else
    nextslen=rand_next_scatlen(ran);
#endif

    //random arimuthal angle
#ifdef MMC_USE_SSE_MATH
    rand_next_aangle_sincos(ran,&sphi,&cphi);
#else
    tmp0=TWO_PI*rand_next_aangle(ran); //next arimuth angle
    MCX_SINCOS(tmp0,&sphi,&cphi);
#endif

    //Henyey-Greenstein Phase Function, "Handbook of Optical Biomedical Diagnostics",2002,Chap3,p234
    //see Boas2002

    if(g>EPS){  //if g is too small, the distribution of theta is bad
	tmp0=(1.f-g*g)/(1.f-g+2.f*g*rand_next_zangle(ran));
	tmp0*=tmp0;
	tmp0=(1.f+g*g-tmp0)/(2.f*g);

    	// when ran=1, CUDA will give me 1.000002 for tmp0 which produces nan later
    	if(tmp0> 1.f) tmp0=1.f;
        if(tmp0<-1.f) tmp0=-1.f;

	theta=acosf(tmp0);
	stheta=sqrt(1.f-tmp0*tmp0);
	//stheta=sin(theta);
	ctheta=tmp0;
    }else{
	theta=acosf(2.f*rand_next_zangle(ran)-1.f);
    	MCX_SINCOS(theta,&stheta,&ctheta);
    }

    if( dir->z>-1.f+EPS && dir->z<1.f-EPS ) {
	tmp0=1.f-dir->z*dir->z;   //reuse tmp to minimize registers
	tmp1=rsqrt(tmp0);
	tmp1=stheta*tmp1;

	p.x=tmp1*(dir->x*dir->z*cphi - dir->y*sphi) + dir->x*ctheta;
	p.y=tmp1*(dir->y*dir->z*cphi + dir->x*sphi) + dir->y*ctheta;
	p.z=-tmp1*tmp0*cphi			    + dir->z*ctheta;
    }else{
	p.x=stheta*cphi;
	p.y=stheta*sphi;
	p.z=(dir->z>0.f)?ctheta:-ctheta;
    }
    if (cfg->ismomentum)
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

void fixphoton(float3 *p,float3 *nodes, int *ee){
        float3 c0={0.f,0.f,0.f};
	int i;
        /*calculate element centroid*/
	for(i=0;i<4;i++)
		c0+=nodes[ee[i]-1];
	*p+=(c0*(float3)(0.25f)-*p)*((float3)(FIX_PHOTON))
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
 * \param[in,out] ran0: the additional random number generator states
 */

void launchphoton(MCXParam *gcfg, ray *r, __global float3 *node,__global int4 *elem, RandType *ran, RandType *ran0){
        int canfocus=1;
        float3 origin={r->p0.x,r->p0.y,r->p0.z};

	r->slen=rand_next_scatlen(ran);
	if(gcfg->srctype==stPencil){ // pencil beam, use the old workflow, except when eid is not given
		if(r->eid>0)
		      return;
	}else if(gcfg->srctype==stPlanar || gcfg->srctype==stPattern || gcfg->srctype==stFourier){
		  float rx=rand_uniform01(ran);
		  float ry=rand_uniform01(ran);
		  r->p0.x=gcfg->srcpos.x+rx*gcfg->srcparam1.x+ry*gcfg->srcparam2.x;
		  r->p0.y=gcfg->srcpos.y+rx*gcfg->srcparam1.y+ry*gcfg->srcparam2.y;
		  r->p0.z=gcfg->srcpos.z+rx*gcfg->srcparam1.z+ry*gcfg->srcparam2.z;
		  r->weight=1.f;
		  if(gcfg->srctype==stPattern){
		    int xsize=(int)gcfg->srcparam1.w;
		    int ysize=(int)gcfg->srcparam2.w;
		    r->posidx=MIN((int)(ry*ysize),ysize-1)*xsize+MIN((int)(rx*xsize),xsize-1);
		    if(gcfg->seed==SEED_FROM_FILE && (gcfg->outputtype==otWL || gcfg->outputtype==otWP)){ // replay mode currently doesn't support multiple source patterns
		    	r->weight=gcfg->srcpattern[MIN( (int)(ry*gcfg->srcparam2.w), (int)gcfg->srcparam2.w-1 )*(int)(gcfg->srcparam1.w)+MIN( (int)(rx*gcfg->srcparam1.w), (int)gcfg->srcparam1.w-1 )];
		    	replayweight[r->photonid] *= r->weight;
		    }
		}else if(gcfg->srctype==stFourier){
		    r->weight=(cos((floorf(gcfg->srcparam1.w)*rx+floorf(gcfg->srcparam2.w)*ry+gcfg->srcparam1.w-floorf(gcfg->srcparam1.w))*TWO_PI)*(1.f-gcfg->srcparam2.w+floorf(gcfg->srcparam2.w))+1.f)*0.5f;
		}
		origin.x+=(gcfg->srcparam1.x+gcfg->srcparam2.x)*0.5f;
		origin.y+=(gcfg->srcparam1.y+gcfg->srcparam2.y)*0.5f;
		origin.z+=(gcfg->srcparam1.z+gcfg->srcparam2.z)*0.5f;
	}else if(gcfg->srctype==stFourierX || gcfg->srctype==stFourier2D){
		float rx=rand_uniform01(ran);
		float ry=rand_uniform01(ran);
		float4 v2=gcfg->srcparam1;
		v2.w*=1.f/(sqrt(gcfg->srcparam1.x*gcfg->srcparam1.x+gcfg->srcparam1.y*gcfg->srcparam1.y+gcfg->srcparam1.z*gcfg->srcparam1.z));
		v2.x=v2.w*(gcfg->srcdir.y*gcfg->srcparam1.z - gcfg->srcdir.z*gcfg->srcparam1.y);
		v2.y=v2.w*(gcfg->srcdir.z*gcfg->srcparam1.x - gcfg->srcdir.x*gcfg->srcparam1.z); 
		v2.z=v2.w*(gcfg->srcdir.x*gcfg->srcparam1.y - gcfg->srcdir.y*gcfg->srcparam1.x);
		r->p0.x=gcfg->srcpos.x+rx*gcfg->srcparam1.x+ry*v2.x;
		r->p0.y=gcfg->srcpos.y+rx*gcfg->srcparam1.y+ry*v2.y;
		r->p0.z=gcfg->srcpos.z+rx*gcfg->srcparam1.z+ry*v2.z;
		if(gcfg->srctype==stFourier2D)
			r->weight=(sin((gcfg->srcparam2.x*rx+gcfg->srcparam2.z)*TWO_PI)*sin((gcfg->srcparam2.y*ry+gcfg->srcparam2.w)*TWO_PI)+1.f)*0.5f; //between 0 and 1
		else
			r->weight=(cos((gcfg->srcparam2.x*rx+gcfg->srcparam2.y*ry+gcfg->srcparam2.z)*TWO_PI)*(1.f-gcfg->srcparam2.w)+1.f)*0.5f; //between 0 and 1
		origin.x+=(gcfg->srcparam1.x+v2.x)*0.5f;
		origin.y+=(gcfg->srcparam1.y+v2.y)*0.5f;
		origin.z+=(gcfg->srcparam1.z+v2.z)*0.5f;
	}else if(gcfg->srctype==stDisk || gcfg->srctype==stGaussian){  // uniform disk and Gaussian beam
		float sphi, cphi;
		float phi=TWO_PI*rand_uniform01(ran);
		sphi=sin(phi);	cphi=cos(phi);
		float r0;
		if(gcfg->srctype==stDisk)
		    r0=sqrt(rand_uniform01(ran))*gcfg->srcparam1.x;
		else if(fabs(r->focus) < 1e-5f || fabs(gcfg->srcparam1.y) < 1e-5f)
		    r0=sqrt(-log(rand_uniform01(ran)))*gcfg->srcparam1.x;
		else{
		    float z0=gcfg->srcparam1.x*gcfg->srcparam1.x*M_PI/gcfg->srcparam1.y; //Rayleigh range
		    r0=sqrt(-log(rand_uniform01(ran))*(1.f+(r->focus*r->focus/(z0*z0))))*gcfg->srcparam1.x;
		}
		if(gcfg->srcdir.z>-1.f+EPS && gcfg->srcdir.z<1.f-EPS){
		    float tmp0=1.f-gcfg->srcdir.z*gcfg->srcdir.z;
		    float tmp1=r0/sqrt(tmp0);
		    r->p0.x=gcfg->srcpos.x+tmp1*(gcfg->srcdir.x*gcfg->srcdir.z*cphi - gcfg->srcdir.y*sphi);
		    r->p0.y=gcfg->srcpos.y+tmp1*(gcfg->srcdir.y*gcfg->srcdir.z*cphi + gcfg->srcdir.x*sphi);
		    r->p0.z=gcfg->srcpos.z-tmp1*tmp0*cphi;
		}else{
   		    r->p0.x+=r0*cphi;
		    r->p0.y+=r0*sphi;
		}
	}else if(gcfg->srctype==stCone || gcfg->srctype==stIsotropic || gcfg->srctype==stArcSin){
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=sin(ang);	cphi=cos(ang);
		if(gcfg->srctype==stCone){  // a solid-angle section of a uniform sphere
		        do{
				ang=(gcfg->srcparam1.y>0) ? TWO_PI*rand_uniform01(ran) : acosf(2.f*rand_uniform01(ran)-1.f); //sine distribution
		        }while(ang>gcfg->srcparam1.x);
		}else{
			if(gcfg->srctype==stIsotropic) // uniform sphere
				ang=acosf(2.f*rand_uniform01(ran)-1.f); //sine distribution
			else
				ang=M_PI*rand_uniform01(ran); //uniform distribution in zenith angle, arcsine
		}
		stheta=sin(ang);
		ctheta=cos(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
		canfocus=0;
                if(gcfg->srctype==stIsotropic)
                    if(r->eid>0)
                        return;
	}else if(gcfg->srctype==stZGaussian){
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=sin(ang);	cphi=cos(ang);
		ang=sqrt(-2.f*log(rand_uniform01(ran)))*(1.f-2.f*rand_uniform01(ran0))*gcfg->srcparam1.x;
		stheta=sin(ang);
		ctheta=cos(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
		canfocus=0;
	}else if(gcfg->srctype==stLine || gcfg->srctype==stSlit){
	      float t=rand_uniform01(ran);
	      r->p0.x+=t*gcfg->srcparam1.x;
	      r->p0.y+=t*gcfg->srcparam1.y;
	      r->p0.z+=t*gcfg->srcparam1.z;

              if(gcfg->srctype==stLine){
	              float s,p;
		      t=1.f-2.f*rand_uniform01(ran);
		      s=1.f-2.f*rand_uniform01(ran);
		      p=sqrt(1.f-r->vec.x*r->vec.x-r->vec.y*r->vec.y)*(rand_uniform01(ran)>0.5f ? 1.f : -1.f);
		      float3 vv;
		      vv.x=r->vec.y*p-r->vec.z*s;
		      vv.y=r->vec.z*t-r->vec.x*p;
		      vv.z=r->vec.x*s-r->vec.y*t;
		      r->vec=vv;
		      //*((float3*)&(r->vec))=(float3)(r->vec.y*p-r->vec.z*s,r->vec.z*t-r->vec.x*p,r->vec.x*s-r->vec.y*t);
	      }
              origin.x+=(gcfg->srcparam1.x)*0.5f;
              origin.y+=(gcfg->srcparam1.y)*0.5f;
              origin.z+=(gcfg->srcparam1.z)*0.5f;
              canfocus=(gcfg->srctype==stSlit);
        }

        if(canfocus && r->focus!=0.f){ // if beam focus is set, determine the incident angle
	        float Rn2;
	        origin.x+=r->focus*r->vec.x;
		origin.y+=r->focus*r->vec.y;
		origin.z+=r->focus*r->vec.z;
		if(r->focus<0.f){ // diverging beam
                     r->vec.x=r->p0.x-origin.x;
                     r->vec.y=r->p0.y-origin.y;
                     r->vec.z=r->p0.z-origin.z;
		}else{             // converging beam
                     r->vec.x=origin.x-r->p0.x;
                     r->vec.y=origin.y-r->p0.y;
                     r->vec.z=origin.z-r->p0.z;
		}
		Rn2=rsqrt(dot(&r->vec,&r->vec)); // normalize
		r->vec=r->vec*Rn2;
	}

        r->p0+=r->vec*EPS;

	/*Caluclate intial element id and bary-centric coordinates*/
	float3 vecS={0.f}, vecAB, vecAC, vecN;
	int is,i,ea,eb,ec;
	float bary[4]={0.f};
	for(is=0;is<srcelemlen;is++){
		int include = 1;
		int *elems=(int *)(elem+(srcelem[is]-1)*gcfg->elemlen);
		for(i=0;i<4;i++){
			ea=elems[out[i][0]]-1;
			eb=elems[out[i][1]]-1;
			ec=elems[out[i][2]]-1;
			vecAB=node[ea]-node[eb]; //repeated for all photons
			vecAC=node[ea]-node[ec]; //repeated for all photons
			vecS=node[ea]-r->p0;
			vecN=cross(vecAB,vecAC); //repeated for all photons, vecN can be precomputed
			bary[facemap[i]]=-dot(&vecS,&vecN);
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
			s=1.f/s;
			r->bary0.x=bary[0]*s;
			r->bary0.y=bary[1]*s;
			r->bary0.z=bary[2]*s;
			r->bary0.w=bary[3]*s;
			for(i=0;i<4;i++){
				if((bary[i]*s)<1e-4f)
					r->faceid=ifacemap[i]+1;
			}
			break;
		}
	}
	if(is==srcelemlen){
	    return;
	}
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
 * \param[in,out] ran0: the additional random number generator states
 * \param[in,out] cfg: simulation configuration structure
 * \param[out] visit: statistics counters of this thread
 */

void onephoton(unsigned int id, __constant MCXParam *gcfg,__global float3 *node,__global int4 *elem,
    __global int *type, __global int4 *facenb, __global float3 *normal, __global float4 *med,
                RandType *ran, RandType *ran0, visitor *visit){

	int oldeid,fixcount=0,exitdet=0;
        float mom;
	float kahany, kahant;
	ray r={gcfg->srcpos,{gcfg->srcdir.x,gcfg->srcdir.y,gcfg->srcdir.z},{MMC_UNDEFINED,0.f,0.f},gcfg->bary0,gcfg->e0,gcfg->dim.y-1,0,0,1.f,0.f,0.f,0.f,0.f,0.,0,NULL,NULL,gcfg->srcdir.w,0};

	__local float *ppath=sharedmem+get_local_id(0)*(visit->reclen-1);
	r.photonid=id;
/*
        if(gcfg->issavedet && gcfg->issaveseed){
                r.photonseed=(void*)calloc(1,(sizeof(RandType)*RAND_BUF_LEN));
		memcpy(r.photonseed,(void *)ran, (sizeof(RandType)*RAND_BUF_LEN));
        }
*/
	/*initialize the photon parameters*/
        launchphoton(gcfg, &r, node, elem, type, facenb, normal, med, ran, ran0);
	ppath[visit->reclen-2] = r.weight; /*last record in partialpath is the initial photon weight*/

	/*use Kahan summation to accumulate weight, otherwise, counter stops at 16777216*/
	/*http://stackoverflow.com/questions/2148149/how-to-sum-a-large-number-of-float-number*/
	int pidx;
	if(gcfg->srctype!=stPattern){
	    if(gcfg->seed==SEED_FROM_FILE && (gcfg->outputtype==otWL || gcfg->outputtype==otWP))
		visit->launchweight=replayweight[r.photonid];
	    else
		visit->launchweight=r.weight;
	}

	while(1){  /*propagate a photon until exit*/
	    r.slen=branchless_badouel_raytet(&r,visit, gcfg,node,elem, type, facenb, normal, med);
	    if(r.pout.x==MMC_UNDEFINED){
	    	  if(r.faceid==-2) break; /*reaches the time limit*/
		  if(fixcount++<MAX_TRIAL){
			fixphoton(&r.p0,node,(int *)(elem+(r.eid-1)*gcfg->elemlen));
			continue;
                  }
	    	  r.eid=-r.eid;
        	  r.faceid=-1;
	    }
	    if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
	            ppath[gcfg->prop-1+type[r.eid-1]]+=r.Lmove;  /*second medianum block is the partial path*/
	    /*move a photon until the end of the current scattering path*/
	    while(r.faceid>=0 && !r.isend){
	            r.p0=r.pout;

		    oldeid=r.eid;
	    	    r.eid=((int *)(facenb+(r.eid-1)*gcfg->elemlen))[r.faceid];

		    if(gcfg->isreflect && (r.eid==0 || med[type[r.eid-1]].n != med[type[oldeid-1]].n )){
			if(! (r.eid==0 && med[type[oldeid-1]].n == gcfg->nout ))
			    reflectray(cfg,&r.vec,tracer,&oldeid,&r.eid,r.faceid,ran);
		    }
	    	    if(r.eid==0) break;
		    /*when a photon enters the domain from the background*/
		    if(type[oldeid-1]==0 && type[r.eid-1]){
                        //if(gcfg->debuglevel&dlExit)
			    MMC_FPRINTF(gcfg->flog,"e %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			    	r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid);
                        if(!gcfg->voidtime)
                            r.photontimer=0.f;
                    }
		    /*when a photon exits the domain into the background*/
		    if(type[oldeid-1] && type[r.eid-1]==0){
                        //if(gcfg->debuglevel&dlExit)
		            MMC_FPRINTF(gcfg->flog,"x %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			        r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid);
			if(!gcfg->isextdet){
                            r.eid=0;
        		    break;
                        }
		    }
//		    if(r.eid==0 && med[type[oldeid-1]].n == gcfg->nout ) break;
	    	    //if(r.pout.x!=MMC_UNDEFINED && (gcfg->debuglevel&dlMove))
	    		MMC_FPRINTF(gcfg->flog,"P %f %f %f %d %u %e\n",r.pout.x,r.pout.y,r.pout.z,r.eid,id,r.slen);

	    	    r.slen=branchless_badouel_raytet(&r,visit,gcfg,node,elem, type, facenb, normal, med);
		    if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
			ppath[gcfg->prop-1+type[r.eid-1]]+=r.Lmove;
		    if(r.faceid==-2) break;
		    fixcount=0;
		    while(r.pout.x==MMC_UNDEFINED && fixcount++<MAX_TRIAL){
		       fixphoton(&r.p0,node,(int *)(elem+(r.eid-1)*gcfg->elemlen));
                       r.slen=branchless_badouel_raytet(&r,visit,gcfg,node,elem, type, facenb, normal, med);
		       if(gcfg->issavedet && r.Lmove>0.f && type[r.eid-1]>0)
	            		ppath[gcfg->prop-1+type[r.eid-1]]+=r.Lmove;
		    }
        	    if(r.pout.x==MMC_UNDEFINED){
        		/*possibily hit an edge or miss*/
			r.eid=-r.eid;
        		break;
        	    }
	    }
	    if(r.eid<=0 || r.pout.x==MMC_UNDEFINED) {
        	    //if(r.eid==0 && (gcfg->debuglevel&dlMove))
        		 MMC_FPRINTF(gcfg->flog,"B %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
		    if(r.eid==0){
                       //if(gcfg->debuglevel&dlExit)
        		 MMC_FPRINTF(gcfg->flog,"E %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			    r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid);
                       if(gcfg->issavedet && gcfg->issaveexit){                                     /*when issaveexit is set to 1*/
                            copystate(ppath+(visit->reclen-2-6),&(r.p0.x),3);  /*columns 7-5 from the right store the exit positions*/
                            copystate(ppath+(visit->reclen-2-3),&(r.vec.x),3); /*columns 4-2 from the right store the exit dirs*/
                       }
		    }else if(r.faceid==-2 && (gcfg->debuglevel&dlMove)){
                         MMC_FPRINTF(gcfg->flog,"T %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    	    }else if(r.eid && r.faceid!=-2  && gcfg->debuglevel&dlEdge){
        		 MMC_FPRINTF(gcfg->flog,"X %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
		    }
		    if(gcfg->issavedet && r.eid==0){
		       int i;
                       if(gcfg->detnum==0 && gcfg->isextdet && type[oldeid-1]==prop+1){
                          exitdet=oldeid;
                       }else
		         for(i=0;i<gcfg->detnum;i++){
        		    if((gcfg->detpos[i].x-r.p0.x)*(gcfg->detpos[i].x-r.p0.x)+
        		       (gcfg->detpos[i].y-r.p0.y)*(gcfg->detpos[i].y-r.p0.y)+
        		       (gcfg->detpos[i].z-r.p0.z)*(gcfg->detpos[i].z-r.p0.z) < gcfg->detpos[i].w*gcfg->detpos[i].w){
			          exitdet=i+1;
                		  break;
        		     }
		         }
		    }
	    	    break;  /*photon exits boundary*/
	    }
	    //if(gcfg->debuglevel&dlMove)
	        MMC_FPRINTF(gcfg->flog,"M %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    if(gcfg->minenergy>0.f && r.weight < gcfg->minenergy && (gcfg->tend-gcfg->tstart)*visit->rtstep<=1.f){ /*Russian Roulette*/
		if(rand_do_roulette(ran)*gcfg->roulettesize<=1.f){
			r.weight*=gcfg->roulettesize;
                        //if(gcfg->debuglevel&dlWeight)
			    MMC_FPRINTF(gcfg->flog,"Russian Roulette bumps r.weight to %f\n",r.weight);
		}else
			break;
	    }
            mom=0.f;
	    r.slen0=mc_next_scatter(med[type[r.eid-1]].g,&r.vec,ran,ran0,cfg,&mom);
	    r.slen=r.slen0;
            if(gcfg->ismomentum && type[r.eid-1]>0)                     /*when ismomentum is set to 1*/
                  ppath[(prop<<1)-1+type[r.eid-1]]+=mom; /*the third medianum block stores the momentum transfer*/
            ppath[type[r.eid-1]-1]++;                          /*the first medianum block stores the scattering event counts*/
	}
	if(gcfg->issavedet && exitdet>0){
		int offset=visit->bufpos*visit->reclen;
		if(visit->bufpos>=visit->detcount){
		    visit->detcount+=DET_PHOTON_BUF;
		}
		visit->partialpath[offset]=exitdet;
	        copystate(visit->partialpath+offset+1,ppath,(visit->reclen-1));
                if(gcfg->issaveseed)
	            copystate(visit->photonseed+visit->bufpos*(sizeof(RandType)*RAND_BUF_LEN),r.photonseed,(sizeof(RandType)*RAND_BUF_LEN));
		visit->bufpos++;
	}
/*
        if(r.photonseed)
		free(r.photonseed);
*/	
	if(gcfg->srctype!=stPattern){
	    kahany=r.Eabsorb-visit->kahanc1[0];
	    kahant=visit->absorbweight[0]+kahany;
	    visit->kahanc1[0]=(kahant-visit->absorbweight[0])-kahany;
	    visit->absorbweight[0]=kahant;
	}
/*else{
	    int psize = (int)gcfg->srcparam1.w * (int)gcfg->srcparam2.w;
#pragma omp critical
{
	    for(pidx=0;pidx<gcfg->srcnum;pidx++){
	    	kahany=r.Eabsorb*gcfg->srcpattern[pidx*psize+r.posidx]-visit->kahanc1[pidx];
	    	kahant=visit->absorbweight[pidx]+kahany;
	    	visit->kahanc1[pidx]=(kahant-visit->absorbweight[pidx])-kahany;
	    	visit->absorbweight[pidx]=kahant;
	    }
}
	}
*/
}

void mmc_main_loop(__constant MCXParam *gcfg,__global float3 *node,__global int4 *elem,
    __global int *type, __global int4 *facenb, __global float3 *normal, __global float4 *med){
        dt=GetTimeMillis();
        MMCDEBUG(cfg,dlTime,(gcfg->flog,"seed=%u\nsimulating ... ",gcfg->seed));

	visitor visit={0.f,0.f,1.f/gcfg->tstep,DET_PHOTON_BUF,0,0,0.f,0.f};
	visit.reclen=(2+((gcfg->ismomentum)>0))*mesh->prop+(gcfg->issaveexit>0)*6+2;
	rng_init(ran0,ran1,(unsigned int *)&(gcfg->seed),threadid);

	/*launch photons*/
	for(i=0;i<gcfg->nphoton;i++){
		visit.raytet=0.f;
		visit.raytet0=0.f;
		if(gcfg->seed==SEED_FROM_FILE)
		    Eabsorb+=onephoton(i,tracer,mesh,gcfg,((RandType *)gcfg->photonseed)+i*RAND_BUF_LEN,ran1,&visit);
		else
		    Eabsorb+=onephoton(i,tracer,mesh,gcfg,ran0,ran1,&visit);
		raytri+=visit.raytet;
		raytri0+=visit.raytet0;
	}
	return 0;
}