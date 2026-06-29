/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2025
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

#ifdef __NVCC__
#define __constant const
#define __private
#define __local
#define __global
#define __kernel __global__

/* CUDA: dispatch persistent per-photon state (RF complex weight, multi-source
 * launch-slot index, detphoton/partial-path bookkeeping) via template params
 * so ptxas can dead-code-eliminate the unused features. Each flag gates state
 * that lives in registers across the photon lifetime, so unlike the previous
 * NODAL_USE_MUA/MUSP templating (which only gated leaf-of-dataflow values)
 * these actually reduce peak register pressure.
 *
 *   IS_RF          - cfg.omega > 0 and not replay: track r.weight_im /
 *                    oldweight_im / complex deposit (saves ~8 regs when off).
 *   IS_MULTISRC    - cfg.srcnum > 1 (pattern source) or extrasrclen>0/srcid<=0
 *                    (adjoint multi-source): track r.posidx / per-slot ppath
 *                    (saves ~12 regs when off).
 *   SAVE_DETPHOTON - cfg.issavedet: maintain ppath partial-path-length array,
 *                    save exit position/dir for detected photons (saves ~6
 *                    regs when off).
 *
 * NODAL_USE_MUA/MUSP moved out of templates: on CUDA they're now runtime reads
 * from gcfg (their lifetime is too short to affect peak register pressure, so
 * templating gave no benefit); on OpenCL they remain build-time #defines
 * since the kernel is JIT-built per simulation anyway. */
#define MMC_TEMPLATE template <const int IS_RF, const int IS_MULTISRC, const int SAVE_DETPHOTON>
#define MMC_TARGS    <IS_RF, IS_MULTISRC, SAVE_DETPHOTON>
#define NODAL_USE_MUA   (GPU_PARAM(gcfg, isnodalmua))
#define NODAL_USE_MUSP  (GPU_PARAM(gcfg, isnodalmusp))

/* Cap mmc_main_loop's per-thread register footprint via launch bounds so the
 * kernel can fit MIN_BLOCKS resident blocks per SM at the default block size.
 * 64 threads/block x 16 blocks/SM = 1024 threads/SM, requiring <=64 regs/thread
 * on archs with 64K regs/SM (sm_52+). Benchmarked on TITAN V (Volta sm_70) with
 * dmmc-cube60: 122 regs (no annotation) -> 64 regs (lb=16) gives +45% throughput
 * in DMMC-grid mode and +15% in pure mesh ray-tracing mode despite ~460 B of
 * spills. Going to lb=20 (48 regs) collapses to 0.34x baseline.
 *
 * Build-time overrides:
 *   -DMMC_NO_LAUNCH_BOUNDS    fully disable the annotation (let ptxas pick regs)
 *   -DMMC_BLOCKSIZE=<N>       override the threads-per-block target (default 64)
 *   -DMMC_MIN_BLOCKS=<N>      override the min-blocks-per-SM target (default 16)
 */
#ifndef MMC_BLOCKSIZE
    #define MMC_BLOCKSIZE 64
#endif
#ifndef MMC_MIN_BLOCKS
    #define MMC_MIN_BLOCKS 16
#endif
#ifdef MMC_NO_LAUNCH_BOUNDS
    #define MMC_LAUNCH_BOUNDS
#else
    #define MMC_LAUNCH_BOUNDS __launch_bounds__(MMC_BLOCKSIZE, MMC_MIN_BLOCKS)
#endif

inline __device__ float3 cross(float3 a, float3 b) {
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
inline __device__ __host__ int get_global_id(int idx) {
    return (idx == 0) ? blockIdx.x * blockDim.x + threadIdx.x
           : ( (idx == 1) ? blockIdx.y * blockDim.y + threadIdx.y
               : blockIdx.z * blockDim.z + threadIdx.z);
}
inline __device__ __host__ int get_local_id(int idx) {
    return (idx == 0) ? threadIdx.x : ( (idx == 1) ? threadIdx.y : threadIdx.z );
}
inline __device__ __host__ int get_local_size(int idx) {
    return (idx == 0) ? blockDim.x : ( (idx == 1) ? blockDim.y : blockDim.z );
}
inline __device__ __host__ float3 operator *(float3 a, float3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}
inline __device__ __host__ float3 operator *(float f, float3 v) {
    return make_float3(v.x * f, v.y * f, v.z * f);
}
inline __device__ __host__ float3 operator *(float3 v, float f) {
    return make_float3(v.x * f, v.y * f, v.z * f);
}
inline __device__ __host__ void operator *=(float3& b, float f) {
    b.x *= f;
    b.y *= f;
    b.z *= f;
}
inline __device__ __host__ float3 operator +(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline __device__ __host__ void operator +=(float3& b, float3 a) {
    b.x += a.x;
    b.y += a.y;
    b.z += a.z;
}
inline __device__ __host__ float3 operator -(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline __device__ __host__ void operator -=(float3& b, float3 a) {
    b.x -= a.x;
    b.y -= a.y;
    b.z -= a.z;
}
inline __device__ __host__ float3 operator /(float3 v, float f) {
    float inv = 1.0f / f;
    return v * inv;
}
inline __device__ __host__ void operator /=(float3& b, float f) {
    float inv = 1.0f / f;
    b.x *= inv;
    b.y *= inv;
    b.z *= inv;
}
inline __device__ __host__ void operator *=(float3& b, float3 f) {
    b.x *= f.x;
    b.y *= f.y;
    b.z *= f.z;
}
inline __device__ __host__ float4 operator -(float4 a) {
    return make_float4(-a.x, -a.y, -a.z, -a.w);
}
inline __device__ __host__ float4 operator *(float4 a, float4 b) {
    return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}
inline __device__ __host__ float4 operator *(float f, float4 v) {
    return make_float4(v.x * f, v.y * f, v.z * f, v.w * f);
}
inline __device__ __host__ float4 operator *(float4 v, float f) {
    return make_float4(v.x * f, v.y * f, v.z * f, v.w * f);
}
inline __device__ __host__ void operator *=(float4& b, float f) {
    b.x *= f;
    b.y *= f;
    b.z *= f;
    b.w *= f;
}
inline __device__ __host__ float4 operator +(float4 a, float4 b) {
    return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
inline __device__ __host__ void operator +=(float4& b, float4 a) {
    b.x += a.x;
    b.y += a.y;
    b.z += a.z;
    b.w += a.w;
}
inline __device__ __host__ float4 operator -(float4 a, float4 b) {
    return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
inline __device__ __host__ void operator -=(float4& b, float4 a) {
    b.x -= a.x;
    b.y -= a.y;
    b.z -= a.z;
    b.w -= a.w;
}
inline __device__ __host__ float3 operator /(float3 a, float3 b) {
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}
inline __device__ __host__ float4 operator /(float4 v, float f) {
    float inv = 1.0f / f;
    return v * inv;
}
inline __device__ __host__ void operator /=(float4& b, float f) {
    float inv = 1.0f / f;
    b.x *= inv;
    b.y *= inv;
    b.z *= inv;
    b.w *= inv;
}
#ifdef __USE_FAST_MATH__
inline __device__ float4 operator /(float4 a, float4 b) {
    return make_float4(__fdividef(a.x, b.x), __fdividef(a.y, b.y), __fdividef(a.z, b.z), __fdividef(a.w, b.w));
}
#else
inline __device__ __host__ float4 operator /(float4 a, float4 b) {
    return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}
#endif

inline __device__ __host__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline __device__ __host__ float dot(float4 a, float4 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline __device__ __host__ float clamp(float f, float a, float b) {
    return max(a, min(f, b));
}
inline __device__ __host__ float3 clamp(float3 v, float a, float b) {
    return make_float3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b));
}

inline __device__ __host__ float3 clamp(float3 v, float3 a, float3 b) {
    return make_float3(clamp(v.x, a.x, b.x), clamp(v.y, a.y, b.y), clamp(v.z, a.z, b.z));
}
inline __device__ __host__ float4 dropbelow0(float4 a) {
    return make_float4((a.x > 0.f) ? a.x : 0.f, (a.y > 0.f) ? a.y : 0.f, (a.z > 0.f) ? a.z : 0.f, (a.w > 0.f) ? a.w : 0.f);
}
inline __device__ __host__ float4 ispositive(float4 a) {
    return make_float4((a.x > 0.f) ? 1.f : 0.f, (a.y > 0.f) ? 1.f : 0.f, (a.z > 0.f) ? 1.f : 0.f, (a.w > 0.f) ? 1.f : 0.f);
}
inline __device__ __host__ float4 convert_float4_rte(float4 v) {
    return make_float4(roundf(v.x), roundf(v.y), roundf(v.z), roundf(v.w));
}
#define FL4(f) make_float4(f,f,f,f)
#define FL3(f) make_float3(f,f,f)
#define FL4_3(f) make_float3(f.x,f.y,f.z)
#define FLT_EPSILON   1.19209290E-07F
#define atomicadd(a,b)  atomicAdd(a,b)
#define atomic_inc(x)   atomicAdd(x,1)

#ifdef MCX_USE_NATIVE
    #define MCX_MATHFUN(fun)              fun
#else
    #define MCX_MATHFUN(fun)              fun
#endif
#define MCX_SINCOS(theta,osin,ocos)   sincosf((theta),&(osin),&(ocos))

#else
typedef struct MMC_FLOAT3 {
    float x, y, z;
} FLOAT3;
#define FL4(f) (f)
#define FL3(f) (f)
#define FL4_3(f) (f.x,f.y,f.z)
#define __constant__  __constant
#define __device__
#ifndef NULL
    #define NULL 0
#endif

/* OpenCL: per-node mua/musp dispatch via JIT macros (-DMCX_NODAL_MUA /
 * -DMCX_NODAL_MUSP appended in mmc_cl_host.c). The kernel always declares
 * gnodemua/gnodemusp args but only dereferences them when the macro is set,
 * so non-recon callers can bind NULL.
 *
 * CUDA: NODAL_USE_MUA / NODAL_USE_MUSP read from gcfg at runtime (see the
 * #ifdef __NVCC__ block at the head of this file). Templating gave no
 * register benefit since the mua/musp lookup is a leaf in the dataflow. */
#ifdef MCX_NODAL_MUA
    #define NODAL_USE_MUA   1
#else
    #define NODAL_USE_MUA   0
#endif
#ifdef MCX_NODAL_MUSP
    #define NODAL_USE_MUSP  1
#else
    #define NODAL_USE_MUSP  0
#endif

/* OpenCL: persistent per-photon state dispatch via build-time JIT flags
 * appended by mmc_cl_host.c. IS_RF / IS_MULTISRC / SAVE_DETPHOTON are 0 or 1
 * compile-time constants here, identical in effect to CUDA template params. */
#ifndef IS_RF
    #define IS_RF           0
#endif
#ifndef IS_MULTISRC
    #define IS_MULTISRC     0
#endif
#ifndef SAVE_DETPHOTON
    #define SAVE_DETPHOTON  0
#endif

/* OpenCL has no function templates; the dispatch lives entirely in the
 * #defines above, so the template prefix and call-site arglist are empty. */
#define MMC_TEMPLATE
#define MMC_TARGS
#define MMC_LAUNCH_BOUNDS


#ifdef MCX_USE_NATIVE
    #define MCX_MATHFUN(fun)              native_##fun
    #define MCX_SINCOS(theta,osin,ocos)   {(osin)=native_sin(theta);(ocos)=native_cos(theta);}
#else
    #define MCX_MATHFUN(fun)              fun
    #define MCX_SINCOS(theta,osin,ocos)   (ocos)=sincos((theta),&(osin))
#endif

#ifdef MCX_SAVE_DETECTORS
    #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#endif

#endif

#ifdef USE_HALF
    #ifndef NVCC
        #pragma OPENCL EXTENSION cl_khr_fp16 : enable
    #endif
    #define FLOAT4VEC half4
    #define TOFLOAT4  convert_half4
#else
    #define FLOAT4VEC float4
    #define TOFLOAT4
#endif

#ifdef MCX_DEBUG
    #define GPUDEBUG(x)        printf x             // enable debugging in CPU mode
#else
    #define GPUDEBUG(x)        {}
#endif

#define R_PI               0.318309886183791f

#define ONE_PI             3.1415926535897932f     //pi

#define C0                 299792458000.f          //speed of light in mm/s
#define R_C0               3.335640951981520e-12f  //1/C0 in s/mm

#define JUST_ABOVE_ONE     1.0001f                 //test for boundary
#define JUST_BELOW_ONE     0.9998f                 /**< test for boundary */
#define SAME_VOXEL         -9999.f                 //scatter within a voxel
#define NO_LAUNCH          9999                    //when fail to launch, for debug

#ifndef __NVCC__
    #define ID_UNDEFINED       0x7FFFFFFFU              /**< flag indicating the index is outside of the volume */
    #define TWO_PI             6.28318530717959f       //2*pi
    #define EPS                FLT_EPSILON             //round-off limit
    #define VERY_BIG           (1.f/FLT_EPSILON)       //a big number
#endif

#define MAX_PROP           4000                     /*maximum property number*/
#define DET_MASK           0xFFFF0000
#define MED_MASK           0x0000FFFF
#define MAX_ACCUM          1000.f
#define R_MIN_MUS          1e9f
#define FIX_PHOTON         1e-3f      /**< offset to the ray to avoid edge/vertex */
#define MAX_TRIAL          3          /**< number of fixes when a photon hits an edge/vertex */

#define MCX_DEBUG_MOVE      1
#define MCX_DEBUG_PROGRESS  2048

#define MCX_DEBUG_REC_LEN   6  /**<  number of floating points per position saved when -D M is used for trajectory */

#define MMC_UNDEFINED      (3.40282347e+38F)
#define SEED_FROM_FILE      -999                         /**< special flag indicating to load seeds from history file */

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define F32N(a) ((a) & 0x80000000)          /**<  Macro to test if a floating point is negative */
#define F32P(a) ((a) ^ 0x80000000)          /**<  Macro to test if a floating point is positive */

typedef struct MMC_Ray {
    float3 p0;                    /**< current photon position */
    float3 vec;                   /**< current photon direction vector */
    float3 pout;                  /**< the intersection position of the ray to the enclosing tet */
    int eid;                      /**< the index of the enclosing tet (starting from 1) */
    int faceid;                   /**< the index of the face at which ray intersects with tet */
    int isend;                    /**< if 1, the scattering event ends before reaching the intersection */
    float weight;                 /**< photon current weight (real part; magnitude for roulette) */
    float weight_im;              /**< imaginary part of photon weight, persists across element crossings for RF */
    float photontimer;            /**< the total time-of-fly of the photon */
    float slen;                   /**< the remaining unitless scattering length = length*mus  */
    float Lmove;                  /**< last photon movement length */
    uint oldidx;
    float oldweight;
    float oldweight_im;           /**< accumulated imaginary fluence deposit in current voxel for RF forward mode */
    unsigned int posidx;          /**< launch position index of the photon for pattern source type; also used as source-slot index in multi-source adjoint mode */
    //int nexteid;                /**< the index to the neighboring tet to be moved into */
    //float4 bary0;               /**< the Barycentric coordinate of the intersection with the tet */
    float slen0;                  /**< initial unitless scattering length = length*mus */
    unsigned int photonid;        /**< index of the current photon */
} ray __attribute__ ((aligned (4)));


typedef struct MMC_Parameter {
    float4 srcparam1;             /**< source parameters set 1 */
    float4 srcparam2;             /**< source parameters set 2 */
    uint4  crop0;
    float4 bary0;
    float3 srcpos;
    float3 srcdir;
    float3 nmin;
    float  tstart, tend;
    uint   isreflect, issavedet, issaveexit, ismomentum, isatomic, isspecular;
    float  Rtstep;
    float  minenergy;
    uint   maxdetphoton;
    uint   maxmedia;
    uint   detnum;
    int    voidtime;
    int    srctype;               /**< type of the source */
    uint   issaveref;             /**<1 save diffuse reflectance at the boundary voxels, 0 do not save*/
    uint   maxgate;
    uint   debuglevel;            /**< debug flags */
    int    reclen;                /**< record (4-byte per record) number per detected photon, does not include srcnum buffer for photon-sharing */
    int    outputtype;
    int    elemlen;
    int    mcmethod;
    int    method;
    float  dstep;
    float  focus;
    int    nn, ne, nf;
    float  nout;
    float  roulettesize;
    int    srcnum;
    int    srcelemlen;
    int    e0;
    int    isextdet;
    uint   framelen;
    uint   maxpropdet;
    uint   normbuf;
    int    issaveseed;
    int    seed;
    uint   maxjumpdebug;          /**< max number of positions to be saved to save photon trajectory when -D M is used */
    float  omega;                  /**< RF modulation angular frequency (rad/s); 0 for CW */
    float  oneoverc0;              /**< 1/C0 = 3.335640951981520e-12 s/mm */
    int    srcid;                  /**< < 0 for multi-source mode (adjoint); >= 0 for single source */
    int    extrasrclen;            /**< number of extra sources packed into gmed[] after media */
    int    srcpropoffset;          /**< gmed[] index where extra sources start (= prop+1+isextdet) */
    uint   isnodalmua;             /**< 1: read mua per-element from gnodemua centroid (DOT recon); CUDA only */
    uint   isnodalmusp;            /**< 1: read musp per-element from gnodemusp centroid (RF DOT recon); CUDA only */
} MCXParam __attribute__ ((aligned (16)));

typedef struct MMC_Reporter {
    float  raytet;
    uint   jumpdebug;
} MCXReporter  __attribute__ ((aligned (4)));

/** Extra source entry for multi-source / adjoint-mode simulation */
#ifndef MCX_EXTRASRC_DEFINED
#define MCX_EXTRASRC_DEFINED
typedef struct MCX_ExtraSrc {
    float4 srcpos;      /**< position (x,y,z) and importance weight (w) */
    float4 srcdir;      /**< direction (x,y,z) and focal length (w) */
    float4 srcparam1;   /**< source parameters set 1: x=radius for disk source */
    float4 srcparam2;   /**< source parameters set 2 */
} ExtraSrc;
#endif

typedef struct MCX_medium {
    float mua;                    /**<absorption coeff in 1/mm unit*/
    float mus;                    /**<scattering coeff in 1/mm unit*/
    float g;                      /**<anisotropy*/
    float n;                      /**<refractive index*/
} Medium __attribute__ ((aligned (16)));

#ifdef __NVCC__
__constant__ MCXParam gcfg[1];
__constant__ Medium   gmed[MAX_PROP];
#define GPU_PARAM(a,b) (a->b)

#else

enum TBoundary {bcNoReflect, bcReflect, bcAbsorbExterior, bcMirror /*, bcCylic*/};

#ifndef USE_MACRO_CONST
    #define GPU_PARAM(a,b) (a->b)
#else
    #define GPU_PARAM(a,b) (a ## b)
#endif

#endif

__constant__ int faceorder[] = {1, 3, 2, 0, -1};
__constant__ int ifaceorder[] = {3, 0, 2, 1};
//__constant int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
//__constant int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
#if defined(MCX_SRC_PLANAR) || defined(MCX_SRC_PATTERN) || defined(MCX_SRC_PATTERN3D) || defined(MCX_SRC_FOURIER) || defined(MCX_SRC_FOURIERX) || defined(MCX_SRC_FOURIERX2D)
__constant__ int out[4][3] = {{0, 3, 1}, {3, 2, 1}, {0, 2, 3}, {0, 1, 2}};
__constant__ int facemap[] = {2, 0, 1, 3};
__constant__ int ifacemap[] = {1, 2, 0, 3};
#endif

#ifndef __NVCC__
enum TDebugLevel {dlMove = 1, dlTracing = 2, dlBary = 4, dlWeight = 8, dlDist = 16, dlTracingEnter = 32,
                  dlTracingExit = 64, dlEdge = 128, dlAccum = 256, dlTime = 512, dlReflect = 1024,
                  dlProgress = 2048, dlExit = 4096, dlTraj = 8192
                 };

enum TRTMethod {rtPlucker, rtHavel, rtBadouel, rtBLBadouel, rtBLBadouelGrid};
enum TMCMethod {mmMCX, mmMCML};

enum TSrcType {stPencil, stIsotropic, stCone, stGaussian, stPlanar,
               stPattern, stFourier, stArcSin, stDisk, stFourierX,
               stFourier2D, stZGaussian, stLine, stSlit
              };
enum TOutputType {otFlux, otFluence, otEnergy, otJacobian, otWL, otWP};
enum TOutputFormat {ofASCII, ofBin, ofJSON, ofUBJSON};
enum TOutputDomain {odMesh, odGrid};
typedef ulong  RandType;
#endif


#define RAND_BUF_LEN       2        //register arrays
#define RAND_SEED_WORD_LEN      4        //48 bit packed with 64bit length
#define LOG_MT_MAX         22.1807097779182f
#define IEEE754_DOUBLE_BIAS     0x3FF0000000000000ul /* Added to exponent.  */


__device__ static float xorshift128p_nextf (__private RandType t[RAND_BUF_LEN]) {
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
#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)
__device__ static void copystate(__local float* v1, __private float* v2, int len) {
    for (int i = 0; i < len; i++) {
        v1[i] = v2[i];
    }
}
#endif

__device__ static float rand_uniform01(__private RandType t[RAND_BUF_LEN]) {
    return xorshift128p_nextf(t);
}

__device__ static void xorshift128p_seed (__global uint* seed, RandType t[RAND_BUF_LEN]) {
    t[0] = (ulong)seed[0] << 32 | seed[1] ;
    t[1] = (ulong)seed[2] << 32 | seed[3];
}

__device__ static void gpu_rng_init(__private RandType t[RAND_BUF_LEN], __global uint* n_seed, int idx) {
    xorshift128p_seed((n_seed + idx * RAND_SEED_WORD_LEN), t);
}

__device__ float rand_next_scatlen(__private RandType t[RAND_BUF_LEN]) {
    return -MCX_MATHFUN(log)(rand_uniform01(t) + EPS);
}

#define rand_next_aangle(t)  rand_uniform01(t)
#define rand_next_zangle(t)  rand_uniform01(t)
#define rand_next_reflect(t) rand_uniform01(t)
#define rand_do_roulette(t)  rand_uniform01(t)

#ifdef USE_ATOMIC

#ifndef __NVCC__

#if defined(USE_NVIDIA_GPU) && !defined(USE_OPENCL_ATOMIC)
// float atomicadd on NVIDIA GPU via PTX
// https://stackoverflow.com/a/72049624/4271392

__device__ inline float atomicadd(volatile __global float* address, const float value) {
    float old;
    asm volatile(
        "atom.global.add.f32 %0, [%1], %2;"
        : "=f"(old)
        : "l"(address), "f"(value)
        : "memory"
    );
    return old;
}
#else
// OpenCL float atomicadd hack:
// http://suhorukov.blogspot.co.uk/2011/12/opencl-11-atomic-operations-on-floating.html
// https://devtalk.nvidia.com/default/topic/458062/atomicadd-float-float-atomicmul-float-float-/

__device__ inline float atomicadd(volatile __global float* address, const float value) {
    float old = value, orig;

    while ((old = atomic_xchg(address, (orig = atomic_xchg(address, 0.0f)) + old)) != 0.0f);

    return orig;
}
#endif

#endif

#endif

__device__ void clearpath(__local float* p, int len) {
    int i;

    for (i = 0; i < len; i++) {
        p[i] = 0.f;
    }
}

#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)
__device__ uint finddetector(float3* p0, __constant float4* gmed, __constant MCXParam* gcfg) {
    uint i;

    for (i = GPU_PARAM(gcfg, maxmedia) + 1 + GPU_PARAM(gcfg, isextdet); i < GPU_PARAM(gcfg, maxmedia) + 1 + GPU_PARAM(gcfg, isextdet) + GPU_PARAM(gcfg, detnum); i++) {
        if ((gmed[i].x - p0[0].x) * (gmed[i].x - p0[0].x) +
                (gmed[i].y - p0[0].y) * (gmed[i].y - p0[0].y) +
                (gmed[i].z - p0[0].z) * (gmed[i].z - p0[0].z) < gmed[i].w * gmed[i].w) {
            return i - GPU_PARAM(gcfg, maxmedia) - GPU_PARAM(gcfg, isextdet);
        }
    }

    return 0;
}

MMC_TEMPLATE
__device__ void savedetphoton(__global float* n_det, __global uint* detectedphoton,
                              __local float* ppath, ray* r, __constant Medium* gmed,
                              int extdetid, __constant MCXParam* gcfg, __global RandType* photonseed, RandType* initseed) {
    uint detid = (extdetid < 0) ? finddetector(&(r->p0), (__constant float4*)gmed, gcfg) : extdetid;

    if (detid) {
        uint baseaddr = atomic_inc(detectedphoton);

        if (baseaddr < GPU_PARAM(gcfg, maxdetphoton)) {
            uint i;
#if defined(MCX_SAVE_SEED) || defined(__NVCC__)
#ifdef __NVCC__

            if (GPU_PARAM(gcfg, issaveseed)) {
#endif

                for (i = 0; i < RAND_BUF_LEN; i++) {
                    photonseed[baseaddr * RAND_BUF_LEN + i] = initseed[i];
                }

#ifdef __NVCC__
            }

#endif

#endif
            baseaddr *= (GPU_PARAM(gcfg, reclen) + 1);

            /* Pack the launch slot index (1-based, matching mcx convention) into the
             * upper 16 bits of the detid column when running in multi-source mode.
             * Lower 16 bits keep the detector id; the upper bits stay zero for
             * single-source runs (extrasrclen==0) so existing readers are unaffected. */
            if (IS_MULTISRC && GPU_PARAM(gcfg, extrasrclen) > 0 && GPU_PARAM(gcfg, srcid) <= 0) {
                detid |= ((unsigned int)(r->posidx + 1u) << 16);
            }

            n_det[baseaddr++] = detid;

            for (i = 0; i < (GPU_PARAM(gcfg, maxmedia) << 1); i++) {
                n_det[baseaddr++] = ppath[i];    // save partial pathlength to the memory
            }

            for (i = 0; i < GPU_PARAM(gcfg, ismomentum)*GPU_PARAM(gcfg, maxmedia); i++) {
                n_det[baseaddr++] = ppath[i + (GPU_PARAM(gcfg, maxmedia) << 1)];    // save partial pathlength to the memory
            }

            if (GPU_PARAM(gcfg, issaveexit)) {
                n_det[baseaddr++] = r->p0.x;
                n_det[baseaddr++] = r->p0.y;
                n_det[baseaddr++] = r->p0.z;
                n_det[baseaddr++] = r->vec.x;
                n_det[baseaddr++] = r->vec.y;
                n_det[baseaddr++] = r->vec.z;
            }

            n_det[baseaddr++] = ppath[GPU_PARAM(gcfg, reclen) - 1]; // save partial pathlength to the memory
        }
    }
}
#endif

/**
 * @brief Saving photon trajectory data for debugging purposes
 * @param[in] p: the position/weight of the current photon packet
 * @param[in] id: the global index of the photon
 * @param[in] gdebugdata: pointer to the global-memory buffer to store the trajectory info
 */

__device__ void savedebugdata(ray* r, uint id, __global MCXReporter* reporter, __global float* gdebugdata, __constant MCXParam* gcfg) {
    uint pos = atomic_inc(&reporter->jumpdebug);

    if (pos < GPU_PARAM(gcfg, maxjumpdebug)) {
        pos *= MCX_DEBUG_REC_LEN;
        ((__global uint*)gdebugdata)[pos++] = id;
        gdebugdata[pos++] = r->p0.x;
        gdebugdata[pos++] = r->p0.y;
        gdebugdata[pos++] = r->p0.z;
        gdebugdata[pos++] = r->weight;
        ((__global uint*)gdebugdata)[pos++] = r->eid;
    }
}

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

MMC_TEMPLATE
__device__ float branchless_badouel_raytet(ray* r, __constant MCXParam* gcfg, __local float* ppath, __global int* elem, __global float* weight,
        int type, __global int* facenb, __global float4* normal, __constant Medium* gmed,
        __global float* gnodemua, __global float* gnodemusp,
        __global float* replayweight, __global float* replaytime) {

    float Lmin;
    float ww, totalloss = 0.f;
    int tshift, faceidx = -1, eid;
    float4 T, S;

    union {
        float f;
        uint  i;
    } currweight;

    if (r->eid <= 0) {
        return -1;
    }

    eid = (r->eid - 1) << 2;

    r->pout.x = MMC_UNDEFINED;
    r->faceid = -1;
    r->isend = 0;

    if (r->eid <= GPU_PARAM(gcfg, normbuf)) {
        eid += GPU_PARAM(gcfg, maxpropdet);
        S = ((r->vec.x) * ((__constant float4*)gmed)[eid]) + ((r->vec.y) * ((__constant float4*)gmed)[eid + 1]) + ((r->vec.z) * ((__constant float4*)gmed)[eid + 2]);
        T = ((__constant float4*)gmed)[eid + 3] - (((r->p0.x) * ((__constant float4*)gmed)[eid]) + ((r->p0.y) * ((__constant float4*)gmed)[eid + 1]) + ((r->p0.z) * ((__constant float4*)gmed)[eid + 2]));
    } else {
        S = ((r->vec.x) * normal[eid]) + ((r->vec.y) * normal[eid + 1]) + ((r->vec.z) * normal[eid + 2]);
        T = normal[eid + 3] - (((r->p0.x) * normal[eid]) + ((r->p0.y) * normal[eid + 1]) + ((r->p0.z) * normal[eid + 2]));
    }

#ifndef __NVCC__
    T = -convert_float4_rte(isgreater(T, FL4(0.f)) * 2) * FL4(0.5f) * T;
#endif
    T = T / S;

#ifndef __NVCC__
    S = -convert_float4_rte(isgreater(S, FL4(0.f)) * 2) * FL4(0.5f);
    T =  (S * T) + ((FL4(1.f) - S) * FL4(1e10f));
#else
    T = make_float4(S.x > 0.f ? T.x : 1e10f, S.y > 0.f ? T.y : 1e10f, S.z > 0.f ? T.z : 1e10f, S.w > 0.f ? T.w : 1e10f);
#endif
    eid = r->eid - 1;

    Lmin = fmin(fmin(fmin(T.x, T.y), T.z), T.w);
    faceidx = ((Lmin == 1e10f) ? 4 : Lmin == T.x ? 0 : (Lmin == T.y ? 1 : (Lmin == T.z ? 2 : 3)));
    r->faceid = faceorder[faceidx];

    if (r->faceid >= 0 && Lmin >= 0.f) {
        Medium prop;

        prop = gmed[type];

        /* Per-node mua/musp override (set by redbird-style DOT reconstruction).
         * Uses the element-centroid average of the four nodal values, which is
         * what the linear-FEM mass-matrix integration of a per-node-linear
         * property reduces to. NODAL_USE_MUA/MUSP are gcfg runtime reads on
         * CUDA, build-time #defines on OpenCL (kernel JIT-rebuilds per
         * simulation), so the unused branch is dead-code-eliminated on OpenCL
         * and a uniform-divergence branch on CUDA. */
        if (NODAL_USE_MUA) {
            __global int* eelocal = elem + (r->eid - 1) * GPU_PARAM(gcfg, elemlen);
            prop.mua = 0.25f * (gnodemua[eelocal[0] - 1] + gnodemua[eelocal[1] - 1]
                                + gnodemua[eelocal[2] - 1] + gnodemua[eelocal[3] - 1]);

            if (NODAL_USE_MUSP) {
                prop.mus = 0.25f * (gnodemusp[eelocal[0] - 1] + gnodemusp[eelocal[1] - 1]
                                    + gnodemusp[eelocal[2] - 1] + gnodemusp[eelocal[3] - 1]);
            }
        }

        currweight.f = r->weight;

        r->Lmove = (prop.mus <= EPS) ? R_MIN_MUS : r->slen / prop.mus;
        r->isend = (Lmin > r->Lmove);
        r->Lmove = ((r->isend) ? r->Lmove : Lmin);
        r->pout = r->p0 + FL3(Lmin) * r->vec;

        if ((int)((r->photontimer + r->Lmove * (prop.n * R_C0) - gcfg->tstart)*GPU_PARAM(gcfg, Rtstep)) > GPU_PARAM(gcfg, maxgate) - 1) { /*exit time window*/
            r->faceid = -2;
            r->pout.x = MMC_UNDEFINED;
            r->Lmove = (gcfg->tend - r->photontimer) / (prop.n * R_C0) - 1e-4f;
        }

        totalloss = MCX_MATHFUN(exp)(-prop.mua * r->Lmove);
        r->weight *= totalloss;

        totalloss = 1.f - totalloss; /*remaining fraction*/

        if (GPU_PARAM(gcfg, seed) == SEED_FROM_FILE) {
            if ((GPU_PARAM(gcfg, outputtype) == otWL) + (GPU_PARAM(gcfg, outputtype) == otJacobian)) {
                currweight.f = r->Lmove;
                currweight.f *= replayweight[r->photonid];
                currweight.f += r->weight;
            } else if (GPU_PARAM(gcfg, outputtype) == otWP) {
                if (r->slen0 < EPS) {
                    currweight.f = 1;
                } else {
                    currweight.f = r->Lmove * prop.mus / r->slen0;
                }

                currweight.f *= replayweight[r->photonid];
                currweight.f += r->weight;
            }
        }

        r->slen -= r->Lmove * prop.mus;
        ww = currweight.f - r->weight;
        r->photontimer += r->Lmove * (prop.n * R_C0);

        if ((GPU_PARAM(gcfg, outputtype) == otWL) + (GPU_PARAM(gcfg, outputtype) == otWP)) {
            tshift = MIN( ((int)(replaytime[r->photonid] * GPU_PARAM(gcfg, Rtstep))), GPU_PARAM(gcfg, maxgate) - 1 ) * GPU_PARAM(gcfg, framelen);
        } else {
            tshift = MIN( ((int)((r->photontimer - gcfg->tstart) * GPU_PARAM(gcfg, Rtstep))), GPU_PARAM(gcfg, maxgate) - 1 ) * GPU_PARAM(gcfg, framelen);
        }

        {
#ifndef MCX_SKIP_VOLUME

            if ((GPU_PARAM(gcfg, outputtype) != otEnergy) * (GPU_PARAM(gcfg, outputtype) != otWP) * (GPU_PARAM(gcfg, outputtype) != otWL)) {
                ww = (prop.mua < EPS) ? (currweight.f * r->Lmove) :
#ifdef __NVCC__
                     __fdividef(ww, prop.mua);
#else
                     (ww / prop.mua);
#endif
            }

#ifdef USE_BLBADOUEL
#ifdef __NVCC__

            if (GPU_PARAM(gcfg, method) == rtBLBadouel) {
#endif
                /* multi-source (adjoint) mode: offset newidx by source-slot × ne × maxgate */
                uint src_slot_offset = (IS_MULTISRC && GPU_PARAM(gcfg, srcid) < 0) ?
                                       (uint)(r->posidx) * GPU_PARAM(gcfg, ne) * GPU_PARAM(gcfg, maxgate) : 0u;
                uint newidx = eid + tshift + src_slot_offset;
                r->oldidx = (r->oldidx == ID_UNDEFINED) ? newidx : r->oldidx;

                /* RF forward (omega>0, no seed-replay): apply complex Beer-Lambert over
                 * this one Lmove step. The deposit fluence per step is
                 *     ∫₀^L w₀ e^{-(mua+i ω n/c0) s} ds
                 *         = (w₀ - w_new) / (mua + i ω n/c0)
                 * matching the per-segment formula already in the rtBLBadouelGrid block. */
                float bl_dep_im = 0.f;
                float bl_dep_re_rf = 0.f;

                if (IS_RF) {
                    float w0_re = currweight.f, w0_im = r->weight_im;
                    float att_re = totalloss < 1.f ? (1.f - totalloss) : 1.f; /* exp(-mua*Lmove) */
                    float phase  = GPU_PARAM(gcfg, omega) * prop.n * GPU_PARAM(gcfg, oneoverc0) * r->Lmove;
                    float cphi, sphi;
                    MCX_SINCOS(phase, sphi, cphi);
                    /* w_new = w₀ * exp(-mua*Lmove) * exp(-i*phase) */
                    float new_re = att_re * (w0_re * cphi + w0_im * sphi);
                    float new_im = att_re * (-w0_re * sphi + w0_im * cphi);
                    /* deposit = (w₀ - w_new) / (mua + i ω n/c0) */
                    float dw_re = w0_re - new_re;
                    float dw_im = w0_im - new_im;
                    float a_im  = GPU_PARAM(gcfg, omega) * prop.n * GPU_PARAM(gcfg, oneoverc0);
                    float a_mag2 = prop.mua * prop.mua + a_im * a_im;
                    bl_dep_re_rf = (a_mag2 > 0.f) ? (dw_re * prop.mua + dw_im * a_im) / a_mag2
                                   : (w0_re * r->Lmove);
                    bl_dep_im    = (a_mag2 > 0.f) ? (dw_im * prop.mua - dw_re * a_im) / a_mag2
                                   : (w0_im * r->Lmove);
                    /* Replace the real-only ww with the matched complex Re(deposit), and
                     * advance r->weight_im for the next step. r->weight (real) is already
                     * att*real(rotation), match it. */
                    ww = bl_dep_re_rf;
                    r->weight    = new_re;
                    r->weight_im = new_im;
                }

                if (newidx != r->oldidx) {
#ifndef DO_NOT_SAVE

                    if (r->oldweight > 0.f) {
                        if ((GPU_PARAM(gcfg, srctype) != stPattern) + (GPU_PARAM(gcfg, srcnum) == 1)) {
#ifdef USE_ATOMIC
                            float oldval = atomicadd(weight + r->oldidx, r->oldweight);

                            if (oldval > MAX_ACCUM) {
                                if (atomicadd(weight + r->oldidx, -oldval) < 0.0f) {
                                    atomicadd(weight + r->oldidx, oldval);
                                } else {
                                    atomicadd(weight + r->oldidx + gcfg->crop0.w, oldval);
                                }
                            }

                            /* RF imag part lives at +2*crop0.w (matches rtBLBadouelGrid). */
                            if (IS_RF) {
                                atomicadd(weight + r->oldidx + gcfg->crop0.w * 2, r->oldweight_im);
                            }

#else
                            weight[r->oldidx] += r->oldweight;

                            if (IS_RF) {
                                weight[r->oldidx + gcfg->crop0.w * 2] += r->oldweight_im;
                            }

#endif
                        } else if (GPU_PARAM(gcfg, srctype) == stPattern) {

                            for (int pidx = 0; pidx < GPU_PARAM(gcfg, srcnum); pidx++) {
#ifdef USE_ATOMIC
                                float oldval = atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx]);

                                if (oldval > MAX_ACCUM) {
                                    if (atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, -oldval) < 0.0f) {
                                        atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, oldval);
                                    } else {
                                        atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx + gcfg->crop0.w, oldval);
                                    }
                                }

#else
                                weight[r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx] += r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx];
#endif
                            }

                        }
                    }

#endif
                    r->oldidx = newidx;
                    r->oldweight = ww;
                    r->oldweight_im = bl_dep_im;
                } else {
                    r->oldweight    += ww;
                    r->oldweight_im += bl_dep_im;
                }

#ifndef DO_NOT_SAVE

                if (r->faceid == -2 || !r->isend) {

                    if ((GPU_PARAM(gcfg, srctype) != stPattern) + (GPU_PARAM(gcfg, srcnum) == 1)) {

#ifdef USE_ATOMIC
                        float oldval = atomicadd(weight + newidx, r->oldweight);

                        if (oldval > MAX_ACCUM) {
                            if (atomicadd(weight + newidx, -oldval) < 0.0f) {
                                atomicadd(weight + newidx, oldval);
                            } else {
                                atomicadd(weight + newidx + gcfg->crop0.w, oldval);
                            }
                        }

                        if (IS_RF) {
                            atomicadd(weight + newidx + gcfg->crop0.w * 2, r->oldweight_im);
                        }

#else
                        weight[newidx] += r->oldweight;

                        if (IS_RF) {
                            weight[newidx + gcfg->crop0.w * 2] += r->oldweight_im;
                        }

#endif
                    } else if (GPU_PARAM(gcfg, srctype) == stPattern) {

                        for (int pidx = 0; pidx < GPU_PARAM(gcfg, srcnum); pidx++) {
#ifdef USE_ATOMIC
                            float oldval = atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx]);

                            if (oldval > MAX_ACCUM) {
                                if (atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, -oldval) < 0.0f) {
                                    atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, oldval);
                                } else {
                                    atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx + gcfg->crop0.w, oldval);
                                }
                            }

#else
                            weight[newidx * GPU_PARAM(gcfg, srcnum) + pidx] += r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx];
#endif
                        }

                    }

                    r->oldweight    = 0.f;
                    r->oldweight_im = 0.f;
                }

#endif // for ifdef DO_NOT_SAVE
#ifdef __NVCC__
            }

#endif // for if (GPU_PARAM(gcfg, method) == rtBLBadouel
#endif // for ifdef USE_BLBADOUEL

#ifdef USE_DMMC
#ifdef __NVCC__

            if (GPU_PARAM(gcfg, method) == rtBLBadouelGrid) {
#endif
                eid = (int)(r->Lmove * GPU_PARAM(gcfg, dstep)) + 1; // number of segments
                eid = (eid << 1);
                S.w = r->Lmove / eid;                 // segment length
                T.w = MCX_MATHFUN(exp)(-prop.mua * S.w); // segment real decay
#ifndef __NVCC__
                T.xyz =  r->vec * FL3(S.w);      // delta vector
                S.xyz =  (r->p0 - gcfg->nmin) + (T.xyz * FL3(0.5f)); /*starting point*/
#else
                T =  make_float4(r->vec.x * S.w, r->vec.y * S.w, r->vec.z * S.w, T.w); // delta vector
                S =  make_float4((r->p0.x - gcfg->nmin.x) + T.x * 0.5f, (r->p0.y - gcfg->nmin.y) + T.y * 0.5f, (r->p0.z - gcfg->nmin.z) + T.z * 0.5f, S.w); /*starting point*/
#endif
                totalloss = (totalloss == 0.f) ? 0.f : (1.f - T.w) / totalloss; // fraction of total loss per segment
                S.w = ww;                             // S.w is now the current weight

                /* multi-source (adjoint) mode: offset newidx by source-slot × field-size-per-source */
                uint src_slot_offset = (IS_MULTISRC && GPU_PARAM(gcfg, srcid) < 0) ?
                                       (uint)(r->posidx) * gcfg->crop0.z * GPU_PARAM(gcfg, maxgate) : 0u;

                /* RF forward: complex weight per-segment state */
                float seg_w_re = currweight.f, seg_w_im = r->weight_im;
                float seg_decay_cos = 1.f, seg_decay_sin = 0.f;

                if (IS_RF) {
                    float phase = GPU_PARAM(gcfg, omega) * prop.n * GPU_PARAM(gcfg, oneoverc0) * (r->Lmove / (float)eid);
                    MCX_SINCOS(phase, seg_decay_sin, seg_decay_cos);
                }

                for (faceidx = 0; faceidx < eid; faceidx++) {
#ifndef __NVCC__
                    int3 idx = convert_int3_rtn(S.xyz * FL3((float)GPU_PARAM(gcfg, dstep)));
                    idx = idx & (idx >= (int3)(0));
#else
                    int3 idx = make_int3((S.x > 0.f) ? __float2int_rd(S.x * GPU_PARAM(gcfg, dstep)) : 0,
                                         (S.y > 0.f) ? __float2int_rd(S.y * GPU_PARAM(gcfg, dstep)) : 0,
                                         (S.z > 0.f) ? __float2int_rd(S.z * GPU_PARAM(gcfg, dstep)) : 0);
#endif
                    uint newidx = (idx.z * gcfg->crop0.y + idx.y * gcfg->crop0.x + idx.x) + tshift + src_slot_offset;
                    r->oldidx = (r->oldidx == ID_UNDEFINED) ? newidx : r->oldidx;

                    /* per-segment RF complex weight deposit: (w0-w)/(mua + i*omega*n/c0) */
                    float seg_deposit_re = S.w * totalloss;
                    float seg_deposit_im = 0.f;

                    if (IS_RF) {
                        float w0_re = seg_w_re, w0_im = seg_w_im;
                        /* complex Beer-Lambert: w *= exp(-mua*s) * exp(-i*omega*n/c0*s) */
                        float att = T.w;
                        float new_re = att * (w0_re * seg_decay_cos + w0_im * seg_decay_sin);
                        float new_im = att * (-w0_re * seg_decay_sin + w0_im * seg_decay_cos);
                        seg_w_re = new_re;
                        seg_w_im = new_im;
                        /* fluence deposit = (w0 - w) / (mua + i*omega*n/c0) */
                        float dw_re = w0_re - new_re;
                        float dw_im = w0_im - new_im;
                        float a_im = GPU_PARAM(gcfg, omega) * prop.n * GPU_PARAM(gcfg, oneoverc0);
                        float a_mag2 = prop.mua * prop.mua + a_im * a_im;
                        seg_deposit_re = (a_mag2 < EPS) ? (w0_re * S.w) : (dw_re * prop.mua + dw_im * a_im) / a_mag2;
                        seg_deposit_im = (a_mag2 < EPS) ? (w0_im * S.w) : (dw_im * prop.mua - dw_re * a_im) / a_mag2;
                    }

                    if (newidx != r->oldidx) {
#ifndef DO_NOT_SAVE

                        if ((GPU_PARAM(gcfg, srctype) != stPattern) + (GPU_PARAM(gcfg, srcnum) == 1) || GPU_PARAM(gcfg, srcid) < 0) {

#ifdef USE_ATOMIC
                            float oldval = atomicadd(weight + r->oldidx, r->oldweight);

                            if (oldval > MAX_ACCUM) {
                                if (atomicadd(weight + r->oldidx, -oldval) < 0.0f) {
                                    atomicadd(weight + r->oldidx, oldval);
                                } else {
                                    atomicadd(weight + r->oldidx + gcfg->crop0.w, oldval);
                                }
                            }

                            /* RF imaginary part into buffer at +2*crop0.w */
                            if (IS_RF) {
                                atomicadd(weight + r->oldidx + gcfg->crop0.w * 2, r->oldweight_im);
                            }

#else
                            weight[r->oldidx] += r->oldweight;

                            if (IS_RF) {
                                weight[r->oldidx + gcfg->crop0.w * 2] += r->oldweight_im;
                            }

#endif
                        } else if (GPU_PARAM(gcfg, srctype) == stPattern) {

                            for (int pidx = 0; pidx < GPU_PARAM(gcfg, srcnum); pidx++) {
#ifdef USE_ATOMIC
                                float oldval = atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx]);

                                if (oldval > MAX_ACCUM) {
                                    if (atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, -oldval) < 0.0f) {
                                        atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx, oldval);
                                    } else {
                                        atomicadd(weight + r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx + gcfg->crop0.w, oldval);
                                    }
                                }

#else
                                weight[r->oldidx * GPU_PARAM(gcfg, srcnum) + pidx] += r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx];
#endif
                            }

                        }

#endif // for ifdef DO_NOT_SAVE
                        r->oldidx = newidx;
                        r->oldweight = seg_deposit_re;
                        r->oldweight_im = seg_deposit_im;
                    } else {
                        r->oldweight += seg_deposit_re;
                        r->oldweight_im += seg_deposit_im;
                    }

#ifndef DO_NOT_SAVE

                    if (r->faceid == -2 || !r->isend) {

                        if ((GPU_PARAM(gcfg, srctype) != stPattern) + (GPU_PARAM(gcfg, srcnum) == 1) || GPU_PARAM(gcfg, srcid) < 0) {

#ifdef USE_ATOMIC
                            float oldval = atomicadd(weight + newidx, r->oldweight);

                            if (oldval > MAX_ACCUM) {
                                if (atomicadd(weight + newidx, -oldval) < 0.0f) {
                                    atomicadd(weight + newidx, oldval);
                                } else {
                                    atomicadd(weight + newidx + gcfg->crop0.w, oldval);
                                }
                            }

                            /* RF imaginary part */
                            if (IS_RF) {
                                atomicadd(weight + newidx + gcfg->crop0.w * 2, r->oldweight_im);
                            }

#else
                            weight[newidx] += r->oldweight;

                            if (IS_RF) {
                                weight[newidx + gcfg->crop0.w * 2] += r->oldweight_im;
                            }

#endif
                        } else if (GPU_PARAM(gcfg, srctype) == stPattern) {

                            for (int pidx = 0; pidx < GPU_PARAM(gcfg, srcnum); pidx++) {
#ifdef USE_ATOMIC
                                float oldval = atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx]);

                                if (oldval > MAX_ACCUM) {
                                    if (atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, -oldval) < 0.0f) {
                                        atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx, oldval);
                                    } else {
                                        atomicadd(weight + newidx * GPU_PARAM(gcfg, srcnum) + pidx + gcfg->crop0.w, oldval);
                                    }
                                }

#else
                                weight[newidx * GPU_PARAM(gcfg, srcnum) + pidx] += r->oldweight * ppath[GPU_PARAM(gcfg, reclen) + pidx];
#endif
                            }

                        }

                        r->oldweight = 0.f;
                        r->oldweight_im = 0.f;
                    }

#endif // for ifdef DO_NOT_SAVE

#ifndef __NVCC__
                    S.w *= T.w;
                    S.xyz += T.xyz;
#else
                    S = make_float4(S.x + T.x, S.y + T.y, S.z + T.z, S.w * T.w);
#endif
                }

                /* RF: persist complex photon weight across element crossings */
                if (IS_RF) {
                    r->weight_im = seg_w_im;
                    r->weight    = seg_w_re;  /* real part also phase-mixed; needed as entry for next element */
                }

#ifdef __NVCC__
            }

#endif
#endif // for USE_DMMC
#endif // for MCX_SKIP_VOLUME
        }

        r->p0 = r->p0 + FL3(r->Lmove) * r->vec;
    }

    return ((r->faceid == -2) ? 0.f : r->slen);
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

__device__ float reflectray(__constant MCXParam* gcfg, float3* c0, int* oldeid, int* eid, int faceid, __private RandType* ran, __global int* type, __global float4* normal, __constant Medium* gmed) {
    /*to handle refractive index mismatch*/
    float3 pnorm = {0.f, 0.f, 0.f};
    float Icos, Re, Im, Rtotal, tmp0, tmp1, tmp2, n1, n2;
    int offs = (*oldeid - 1) << 2;

    faceid = ifaceorder[faceid];
    /*calculate the normal direction of the intersecting triangle*/
    pnorm.x = ((__global float*) & (normal[offs]))[faceid];
    pnorm.y = ((__global float*) & (normal[offs]))[faceid + 4];
    pnorm.z = ((__global float*) & (normal[offs]))[faceid + 8];

    /*pn pointing outward*/

    /*compute the cos of the incidence angle*/
    Icos = fabs(dot(*c0, pnorm));

    n1 = ((*oldeid != *eid) ? gmed[type[*oldeid - 1]].n : GPU_PARAM(gcfg, nout));
    n2 = ((*eid > 0) ? gmed[type[*eid - 1]].n : GPU_PARAM(gcfg, nout));

    tmp0 = n1 * n1;
    tmp1 = n2 * n2;
    tmp2 = 1.f - tmp0 / tmp1 * (1.f - Icos * Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

    if (tmp2 > 0.f && !((*eid <= 0) * (GPU_PARAM(gcfg, isreflect) == bcMirror))) { /*if no total internal reflection*/
        Re = tmp0 * Icos * Icos + tmp1 * tmp2; /*transmission angle*/
        tmp2 = MCX_MATHFUN(sqrt)(tmp2); /*to save one sqrt*/
        Im = 2.f * n1 * n2 * Icos * tmp2;
        Rtotal = (Re - Im) / (Re + Im); /*Rp*/
        Re = tmp1 * Icos * Icos + tmp0 * tmp2 * tmp2;
        Rtotal = (Rtotal + (Re - Im) / (Re + Im)) * 0.5f; /*(Rp+Rs)/2*/

        if (*oldeid == *eid) {
            return Rtotal;    /*initial specular reflection*/
        }

        if (rand_next_reflect(ran) <= Rtotal) { /*do reflection*/
            *c0 += (FL3(-2.f * Icos)) * pnorm;
            //if(GPU_PARAM(gcfg,debuglevel)&dlReflect) GPUDEBUG(("R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal));
            *eid = *oldeid; /*stay with the current element*/
        } else if ((GPU_PARAM(gcfg, isspecular) == 2) * (*eid == 0)) {
            // if do transmission, but next neighbor is 0, terminate
        } else {                             /*do transmission*/
            *c0 += (FL3(-Icos)) * pnorm;
            *c0 = (FL3(tmp2)) * pnorm + FL3(n1 / n2) * (*c0);
            //if(GPU_PARAM(gcfg,debuglevel)&dlReflect) GPUDEBUG(("Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal));
        }
    } else { /*total internal reflection*/
        *c0 += (FL3(-2.f * Icos)) * pnorm;
        *eid = *oldeid;
        //if(GPU_PARAM(gcfg,debuglevel)&dlReflect) GPUDEBUG(("V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f));
    }

    tmp0 = MCX_MATHFUN(rsqrt)(dot(*c0, *c0));
    (*c0) *= FL3(tmp0);
    return 1.f;
}

#endif

__device__ void rotatevector(float3* dir, float stheta, float ctheta, float sphi, float cphi) {
    float3 p;
    float tmp0;

    if ( dir->z > -1.f + EPS && dir->z < 1.f - EPS ) {
        tmp0 = 1.f - dir->z * dir->z; //reuse tmp to minimize registers
        float tmp1 = MCX_MATHFUN(rsqrt)(tmp0);
        tmp1 = stheta * tmp1;

        p.x = tmp1 * (dir->x * dir->z * cphi - dir->y * sphi) + dir->x * ctheta;
        p.y = tmp1 * (dir->y * dir->z * cphi + dir->x * sphi) + dir->y * ctheta;
        p.z = -tmp1 * tmp0 * cphi            + dir->z * ctheta;
    } else {
        p.x = stheta * cphi;
        p.y = stheta * sphi;
        p.z = (dir->z > 0.f) ? ctheta : -ctheta;
    }

    tmp0 = MCX_MATHFUN(rsqrt)(p.x * p.x + p.y * p.y + p.z * p.z);

    dir->x = p.x * tmp0;
    dir->y = p.y * tmp0;
    dir->z = p.z * tmp0;
}

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

__device__ float mc_next_scatter(float g, float3* dir, __private RandType* ran, __constant MCXParam* gcfg, float* pmom) {

    float nextslen;
    float sphi, cphi, tmp0, theta, stheta, ctheta;

    //random scattering length (normalized)
    nextslen = rand_next_scatlen(ran);

    tmp0 = TWO_PI * rand_next_aangle(ran); //next arimuth angle
    MCX_SINCOS(tmp0, sphi, cphi);

    if (g > EPS) { //if g is too small, the distribution of theta is bad
        tmp0 = (1.f - g * g) / (1.f - g + 2.f * g * rand_next_zangle(ran));
        tmp0 *= tmp0;
        tmp0 = (1.f + g * g - tmp0) / (2.f * g);
        tmp0 = clamp(tmp0, -1.f, 1.f);

        theta = acos(tmp0);
        stheta = MCX_MATHFUN(sqrt)(1.f - tmp0 * tmp0);
        //stheta=MCX_MATHFUN(sin)(theta);
        ctheta = tmp0;
    } else {
        theta = acos(2.f * rand_next_zangle(ran) - 1.f);
        MCX_SINCOS(theta, stheta, ctheta);
    }

    rotatevector(dir, stheta, ctheta, sphi, cphi);

    if (GPU_PARAM(gcfg, ismomentum)) {
        pmom[0] += (1.f - ctheta);
    }

    return nextslen;
}

/**
 * \brief Function to deal with ray-edge/ray-vertex intersections
 *
 * when a photon is crossing a vertex or edge, (slightly) pull the
 * photon toward the center of the element and try again
 *
 * \param[in,out] p: current photon position
 * \param[in] node: pointer to the 4 nodes of the tet
 * \param[in] ee: indices of the 4 nodes ee=elem[eid]
 */

__device__ void fixphoton(float3* p, __global FLOAT3* node, __global int* ee) {
    float3 c0 = {0.f, 0.f, 0.f};
    int i;

    /*calculate element centroid*/
    for (i = 0; i < 4; i++) {
        c0.x += node[ee[i] - 1].x;
        c0.y += node[ee[i] - 1].y;
        c0.z += node[ee[i] - 1].z;
    }

    *p += (c0 * FL3(0.25f) - *p) * (FL3(FIX_PHOTON));
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

MMC_TEMPLATE
__device__ void launchnewphoton(__constant MCXParam* gcfg, ray* r, __global FLOAT3* node, __global int* elem, __global int* srcelem, __private RandType* ran, __global float* srcpattern, __constant Medium* gmed) {
    int canfocus = 1;
    float3 origin = r->p0;

    /* Multi-source / single-slot mode: pick a source slot from srcdata[] (packed in gmed[] after media).
     *   srcid <  0  : pick a slot uniformly at random (every photon may go to a different slot;
     *                 r->posidx routes its deposits to the per-slot field buffer offset).
     *   srcid >  0  : launch only from srcdata[srcid-1] (1-based selector; matches mcx parity).
     *                 Field buffer collapses to one slot, so r->posidx is forced to 0.
     */
    /* Use bitwise & / | (rather than logical && / ||) so the OpenCL JIT
     * doesn't warn when GPU_PARAM(gcfg, srcid) and ...extrasrclen are
     * constant-folded macros - same convention as the omega/seed/srctype
     * guards elsewhere in this file. */
    if (IS_MULTISRC && (GPU_PARAM(gcfg, extrasrclen) > 0)
            & (((GPU_PARAM(gcfg, srcid) < 0))
               | ((GPU_PARAM(gcfg, srcid) > 0) & (GPU_PARAM(gcfg, srcid) <= GPU_PARAM(gcfg, extrasrclen))))) {
        unsigned int slot;            /* slot index into srcdata[] (source geometry) */
        unsigned int outslot;         /* slot index into field buffer (output) */

        if (GPU_PARAM(gcfg, srcid) < 0) {
            /* uniformly select a source slot */
            slot = (unsigned int)(rand_uniform01(ran) * GPU_PARAM(gcfg, extrasrclen));

            if (slot >= (unsigned int)GPU_PARAM(gcfg, extrasrclen)) {
                slot = (unsigned int)GPU_PARAM(gcfg, extrasrclen) - 1u;
            }

            outslot = slot;
        } else {
            slot = (unsigned int)(GPU_PARAM(gcfg, srcid) - 1);    /* 1-based to 0-based */
            outslot = 0u;                                          /* collapse to one output slot */
        }

        /* Extra sources are packed into gmed[] starting at srcpropoffset;
         * each ExtraSrc occupies 4 Medium slots (both are 16-byte aligned structs of 4 floats) */
        __constant ExtraSrc* srcs = (__constant ExtraSrc*)(gmed + GPU_PARAM(gcfg, srcpropoffset));

        /* Use scalar field access to remain compatible with both CUDA and OpenCL */
        float src_pos_x = srcs[slot].srcpos.x;
        float src_pos_y = srcs[slot].srcpos.y;
        float src_pos_z = srcs[slot].srcpos.z;
        float src_dir_x = srcs[slot].srcdir.x;
        float src_dir_y = srcs[slot].srcdir.y;
        float src_dir_z = srcs[slot].srcdir.z;

        r->posidx = outslot;  /* record source slot for field-buffer indexing */
        r->p0.x = src_pos_x;
        r->p0.y = src_pos_y;
        r->p0.z = src_pos_z;
        origin   = r->p0;

        /* direction */
        r->vec.x = src_dir_x;
        r->vec.y = src_dir_y;
        r->vec.z = src_dir_z;

        /* disk source: apply uniform disk sampling for adjoint (detector) sources */
        {
            float radius = srcs[slot].srcparam1.x;

            if (radius > 0.f) {
                float phi_d = TWO_PI * rand_uniform01(ran);
                float r0    = MCX_MATHFUN(sqrt)(rand_uniform01(ran)) * radius;
                float sphi, cphi;
                MCX_SINCOS(phi_d, sphi, cphi);

                if (src_dir_z > -1.f + EPS && src_dir_z < 1.f - EPS) {
                    float tmp0 = 1.f - src_dir_z * src_dir_z;
                    float tmp1 = r0 * MCX_MATHFUN(rsqrt)(tmp0);
                    r->p0.x = r->p0.x + tmp1 * (src_dir_x * src_dir_z * cphi - src_dir_y * sphi);
                    r->p0.y = r->p0.y + tmp1 * (src_dir_y * src_dir_z * cphi + src_dir_x * sphi);
                    r->p0.z = r->p0.z - tmp1 * tmp0 * cphi;
                } else {
                    r->p0.x += r0 * cphi;
                    r->p0.y += r0 * sphi;
                }
            }
        }

        /* focal-point focusing (srcdir.w = focal length) */
        float focallen = srcs[slot].srcdir.w;

        if (focallen != 0.f) {
            canfocus = 0;
            origin = r->p0 + FL3(focallen) * r->vec;
        }

        r->weight = srcs[slot].srcpos.w; /* importance weight */
        /* Per-slot initial tet from srcdata[slot].srcparam2.w (pre-computed
         * host-side in mmclab.cpp after the adjoint-slot setup). For slot 0
         * this typically matches gcfg->e0; for slots 1+ in multi-source mode
         * it points to the tet containing each slot's own srcpos. Fallback to
         * gcfg->e0 if the host didn't fill it (e.g., legacy paths). */
        int slot_eid = (int)srcs[slot].srcparam2.w;
        r->eid = (slot_eid > 0) ? slot_eid : GPU_PARAM(gcfg, e0);
        r->slen = rand_next_scatlen(ran);
        return;
    }

    r->slen = rand_next_scatlen(ran);
#if defined(__NVCC__) || defined(MCX_SRC_PENCIL)
#ifdef __NVCC__

    if (GPU_PARAM(gcfg, srctype) == MCX_SRC_PENCIL) {
#endif

        if (r->eid > 0) {
            return;
        }

#endif
#if defined(__NVCC__) || defined(MCX_SRC_PLANAR) || defined(MCX_SRC_PATTERN) || defined(MCX_SRC_PATTERN3D) || defined(MCX_SRC_FOURIER) /*a rectangular grid over a plane*/
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_PLANAR || GPU_PARAM(gcfg, srctype) == MCX_SRC_PATTERN || GPU_PARAM(gcfg, srctype) == MCX_SRC_PATTERN3D || GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIER) {
#endif
        float rx = rand_uniform01(ran);
        float ry = rand_uniform01(ran);
        r->p0.x = gcfg->srcpos.x + rx * gcfg->srcparam1.x + ry * gcfg->srcparam2.x;
        r->p0.y = gcfg->srcpos.y + rx * gcfg->srcparam1.y + ry * gcfg->srcparam2.y;
        r->p0.z = gcfg->srcpos.z + rx * gcfg->srcparam1.z + ry * gcfg->srcparam2.z;
        r->weight = 1.f;
#if defined(__NVCC__) || defined(MCX_SRC_PATTERN)
#ifdef __NVCC__

        if (GPU_PARAM(gcfg, srctype) == MCX_SRC_PATTERN) {
#endif
            int xsize = (int)gcfg->srcparam1.w;
            int ysize = (int)gcfg->srcparam2.w;

            if (IS_MULTISRC) {
                r->posidx = MIN((int)(ry * JUST_BELOW_ONE * ysize), ysize - 1) * xsize + MIN((int)(rx * JUST_BELOW_ONE * xsize), xsize - 1);
                r->weight = (GPU_PARAM(gcfg, srcnum) > 1) ? 1.f : srcpattern[r->posidx];
            }

#endif
#if defined(__NVCC__) || defined(MCX_SRC_FOURIER)  // need to prevent rx/ry=1 here
#ifdef __NVCC__
        } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIER)
#endif
            r->weight = (MCX_MATHFUN(cos)((floor(gcfg->srcparam1.w) * rx + floor(gcfg->srcparam2.w) * ry + gcfg->srcparam1.w - floor(gcfg->srcparam1.w)) * TWO_PI) * (1.f - gcfg->srcparam2.w + floor(gcfg->srcparam2.w)) + 1.f) * 0.5f;

#endif
        origin.x += (gcfg->srcparam1.x + gcfg->srcparam2.x) * 0.5f;
        origin.y += (gcfg->srcparam1.y + gcfg->srcparam2.y) * 0.5f;
        origin.z += (gcfg->srcparam1.z + gcfg->srcparam2.z) * 0.5f;
#endif
#if defined(__NVCC__) || defined(MCX_SRC_FOURIERX) || defined(MCX_SRC_FOURIERX2D) // [v1x][v1y][v1z][|v2|]; [kx][ky][phi0][M], unit(v0) x unit(v1)=unit(v2)
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIERX || GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIERX2D) {
#endif
        float rx = rand_uniform01(ran);
        float ry = rand_uniform01(ran);
        float tmp = gcfg->srcparam1.w * MCX_MATHFUN(rsqrt)(gcfg->srcparam1.x * gcfg->srcparam1.x + gcfg->srcparam1.y * gcfg->srcparam1.y + gcfg->srcparam1.z * gcfg->srcparam1.z);
        r->p0.x = gcfg->srcpos.x + rx * gcfg->srcparam1.x + ry * tmp * (gcfg->srcdir.y * gcfg->srcparam1.z - gcfg->srcdir.z * gcfg->srcparam1.y);
        r->p0.y = gcfg->srcpos.y + rx * gcfg->srcparam1.y + ry * tmp * (gcfg->srcdir.z * gcfg->srcparam1.x - gcfg->srcdir.x * gcfg->srcparam1.z);
        r->p0.z = gcfg->srcpos.z + rx * gcfg->srcparam1.z + ry * tmp * (gcfg->srcdir.x * gcfg->srcparam1.y - gcfg->srcdir.y * gcfg->srcparam1.x);
#if defined(__NVCC__) || defined(MCX_SRC_FOURIERX2D)
#ifdef __NVCC__

        if (GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIERX2D)
#endif
            r->weight = (MCX_MATHFUN(sin)((gcfg->srcparam2.x * rx + gcfg->srcparam2.z) * TWO_PI) * MCX_MATHFUN(sin)((gcfg->srcparam2.y * ry + gcfg->srcparam2.w) * TWO_PI) + 1.f) * 0.5f; //between 0 and 1

#endif
#if defined(__NVCC__) || defined(MCX_SRC_FOURIERX)
#ifdef __NVCC__
        else
#endif
            r->weight = (MCX_MATHFUN(cos)((gcfg->srcparam2.x * rx + gcfg->srcparam2.y * ry + gcfg->srcparam2.z) * TWO_PI) * (1.f - gcfg->srcparam2.w) + 1.f) * 0.5f; //between 0 and 1

#endif
        origin.x += (gcfg->srcparam1.x + tmp * (gcfg->srcdir.y * gcfg->srcparam1.z - gcfg->srcdir.z * gcfg->srcparam1.y)) * 0.5f;
        origin.y += (gcfg->srcparam1.y + tmp * (gcfg->srcdir.z * gcfg->srcparam1.x - gcfg->srcdir.x * gcfg->srcparam1.z)) * 0.5f;
        origin.z += (gcfg->srcparam1.z + tmp * (gcfg->srcdir.x * gcfg->srcparam1.y - gcfg->srcdir.y * gcfg->srcparam1.x)) * 0.5f;
#endif
#if defined(__NVCC__) || defined(MCX_SRC_DISK) || defined(MCX_SRC_GAUSSIAN) // uniform disk distribution or Gaussian-beam
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_DISK || GPU_PARAM(gcfg, srctype) == MCX_SRC_GAUSSIAN) {
#endif
        float sphi, cphi;
        float phi = TWO_PI * rand_uniform01(ran);
        sphi = MCX_MATHFUN(sin)(phi);
        cphi = MCX_MATHFUN(cos)(phi);
        float r0;
#if defined(__NVCC__) || defined(MCX_SRC_DISK)
#ifdef __NVCC__

        if (GPU_PARAM(gcfg, srctype) == MCX_SRC_DISK) {
#endif
            r0 = MCX_MATHFUN(sqrt)(rand_uniform01(ran)) * gcfg->srcparam1.x;
#endif
#if defined(__NVCC__) || defined(MCX_SRC_GAUSSIAN)
#ifdef __NVCC__
        } else {
#endif

            if (fabs(GPU_PARAM(gcfg, focus)) < 1e-5f || fabs(gcfg->srcparam1.y) < 1e-5f) {
                r0 = MCX_MATHFUN(sqrt)(-MCX_MATHFUN(log)((rand_uniform01(ran)))) * gcfg->srcparam1.x;
            } else {
                float z0 = gcfg->srcparam1.x * gcfg->srcparam1.x * M_PI / gcfg->srcparam1.y; //Rayleigh range
                r0 = MCX_MATHFUN(sqrt)(-MCX_MATHFUN(log)((rand_uniform01(ran)) * (1.f + (GPU_PARAM(gcfg, focus) * GPU_PARAM(gcfg, focus) / (z0 * z0))))) * gcfg->srcparam1.x;
            }

#ifdef __NVCC__
        }

#endif
#endif

        if (gcfg->srcdir.z > -1.f + EPS && gcfg->srcdir.z < 1.f - EPS) {
            float tmp0 = 1.f - gcfg->srcdir.z * gcfg->srcdir.z;
            float tmp1 = r0 * MCX_MATHFUN(rsqrt)(tmp0);
            r->p0.x = gcfg->srcpos.x + tmp1 * (gcfg->srcdir.x * gcfg->srcdir.z * cphi - gcfg->srcdir.y * sphi);
            r->p0.y = gcfg->srcpos.y + tmp1 * (gcfg->srcdir.y * gcfg->srcdir.z * cphi + gcfg->srcdir.x * sphi);
            r->p0.z = gcfg->srcpos.z - tmp1 * tmp0 * cphi;
        } else {
            r->p0.x += r0 * cphi;
            r->p0.y += r0 * sphi;
        }

#endif
#if defined(__NVCC__) || defined(MCX_SRC_CONE) || defined(MCX_SRC_ISOTROPIC) || defined(MCX_SRC_ARCSINE)
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_CONE || GPU_PARAM(gcfg, srctype) == MCX_SRC_ISOTROPIC || GPU_PARAM(gcfg, srctype) == MCX_SRC_ARCSINE) {
#endif
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
        sphi = MCX_MATHFUN(sin)(ang);
        cphi = MCX_MATHFUN(cos)(ang);
#if defined(__NVCC__) || defined(MCX_SRC_CONE) // a solid-angle section of a uniform sphere
#ifdef __NVCC__

        if (GPU_PARAM(gcfg, srctype) == MCX_SRC_CONE) {
#endif

            do {
                ang = (gcfg->srcparam1.y > 0) ? TWO_PI * rand_uniform01(ran) : acos(2.f * rand_uniform01(ran) - 1.f); //sine distribution
            } while (ang > gcfg->srcparam1.x);

#endif
#if defined(__NVCC__) || defined(MCX_SRC_ISOTROPIC) || defined(MCX_SRC_ARCSINE)
#ifdef __NVCC__
        } else {
#endif

            if (GPU_PARAM(gcfg, srctype) == stIsotropic) { // uniform sphere
                ang = acos(2.f * rand_uniform01(ran) - 1.f);    //sine distribution
            } else {
                ang = M_PI * rand_uniform01(ran);    //uniform distribution in zenith angle, arcsine
            }

#ifdef __NVCC__
        }

#endif
#endif
        stheta = MCX_MATHFUN(sin)(ang);
        ctheta = MCX_MATHFUN(cos)(ang);
        r->vec.x = stheta * cphi;
        r->vec.y = stheta * sphi;
        r->vec.z = ctheta;
        canfocus = 0;

        if (GPU_PARAM(gcfg, srctype) == stIsotropic)
            if (r->eid > 0) {
                return;
            }

#endif
#if defined(__NVCC__) || defined(MCX_SRC_ZGAUSSIAN)
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_ZGAUSSIAN) {
#endif
        float ang, stheta, ctheta, sphi, cphi;
        ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
        sphi = MCX_MATHFUN(sin)(ang);
        cphi = MCX_MATHFUN(cos)(ang);
        ang = MCX_MATHFUN(sqrt)(-2.f * MCX_MATHFUN(log)((rand_uniform01(ran)))) * (1.f - 2.f * rand_uniform01(ran)) * gcfg->srcparam1.x;
        stheta = MCX_MATHFUN(sin)(ang);
        ctheta = MCX_MATHFUN(cos)(ang);
        r->vec.x = stheta * cphi;
        r->vec.y = stheta * sphi;
        r->vec.z = ctheta;
        canfocus = 0;
#endif
#if defined(__NVCC__) || defined(MCX_SRC_LINE) || defined(MCX_SRC_SLIT)
#ifdef __NVCC__
    } else if (GPU_PARAM(gcfg, srctype) == MCX_SRC_LINE || GPU_PARAM(gcfg, srctype) == MCX_SRC_SLIT) {
#endif
        float t = rand_uniform01(ran);
        r->p0.x += t * gcfg->srcparam1.x;
        r->p0.y += t * gcfg->srcparam1.y;
        r->p0.z += t * gcfg->srcparam1.z;

#if defined(__NVCC__) || defined(MCX_SRC_LINE)
#ifdef __NVCC__

        if (GPU_PARAM(gcfg, srctype) == MCX_SRC_LINE) {
#endif
            float s, p;
            t = 1.f - 2.f * rand_uniform01(ran);
            s = 1.f - 2.f * rand_uniform01(ran);
            p = MCX_MATHFUN(sqrt)(1.f - r->vec.x * r->vec.x - r->vec.y * r->vec.y) * (rand_uniform01(ran) > 0.5f ? 1.f : -1.f);
            float3 vv;
            vv.x = r->vec.y * p - r->vec.z * s;
            vv.y = r->vec.z * t - r->vec.x * p;
            vv.z = r->vec.x * s - r->vec.y * t;
            r->vec = vv;
            //*((float3*)&(r->vec))=(float3)(r->vec.y*p-r->vec.z*s,r->vec.z*t-r->vec.x*p,r->vec.x*s-r->vec.y*t);
#ifdef __NVCC__
        }

#endif
#endif
        origin.x += (gcfg->srcparam1.x) * 0.5f;
        origin.y += (gcfg->srcparam1.y) * 0.5f;
        origin.z += (gcfg->srcparam1.z) * 0.5f;
        canfocus = (GPU_PARAM(gcfg, srctype) == stSlit);
#ifdef __NVCC__
    }

#endif
#endif

    if (canfocus) {
        if (isnan(GPU_PARAM(gcfg, focus))) { // if beam focus is set, determine the incident angle
            float ang, stheta, ctheta, sphi, cphi;
            ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
            MCX_SINCOS(ang, sphi, cphi);
            ang = acos(2.f * rand_uniform01(ran) - 1.f); //sine distribution
            MCX_SINCOS(ang, stheta, ctheta);
            rotatevector(&(r->vec), stheta, ctheta, sphi, cphi);
        } else if ((GPU_PARAM(gcfg, focus) < 0.f) * (isinf(GPU_PARAM(gcfg, focus)))) { // lambertian (cosine distribution) if focal length is -inf
            float ang, stheta, ctheta, sphi, cphi;
            ang = TWO_PI * rand_uniform01(ran); //next arimuth angle
            MCX_SINCOS(ang, sphi, cphi);
            stheta = MCX_MATHFUN(sqrt)(rand_uniform01(ran));
            ctheta = MCX_MATHFUN(sqrt)(1.f - stheta * stheta);
            rotatevector(&(r->vec), stheta, ctheta, sphi, cphi);
        } else if (GPU_PARAM(gcfg, focus) != 0.f) { // if beam focus is set, determine the incident angle
            float Rn2;
            origin.x += GPU_PARAM(gcfg, focus) * r->vec.x;
            origin.y += GPU_PARAM(gcfg, focus) * r->vec.y;
            origin.z += GPU_PARAM(gcfg, focus) * r->vec.z;

            if (GPU_PARAM(gcfg, focus) < 0.f) { // diverging beam
                r->vec.x = r->p0.x - origin.x;
                r->vec.y = r->p0.y - origin.y;
                r->vec.z = r->p0.z - origin.z;
            } else {            // converging beam
                r->vec.x = origin.x - r->p0.x;
                r->vec.y = origin.y - r->p0.y;
                r->vec.z = origin.z - r->p0.z;
            }

            Rn2 = MCX_MATHFUN(rsqrt)(dot(r->vec, r->vec)); // normalize
            r->vec = r->vec * Rn2;
        }
    }

    r->p0 += r->vec * EPS;

#if defined(__NVCC__) || defined(MCX_SRC_PLANAR) || defined(MCX_SRC_PATTERN) || defined(MCX_SRC_PATTERN3D) || defined(MCX_SRC_FOURIER) || defined(MCX_SRC_FOURIERX) || defined(MCX_SRC_FOURIERX2D)
#ifdef __NVCC__

    if (GPU_PARAM(gcfg, srctype) == MCX_SRC_PLANAR || GPU_PARAM(gcfg, srctype) == MCX_SRC_PATTERN || GPU_PARAM(gcfg, srctype) == MCX_SRC_PATTERN3D || GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIER || GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIERX || GPU_PARAM(gcfg, srctype) == MCX_SRC_FOURIERX2D) {
#endif
        /*Caluclate intial element id and bary-centric coordinates for area sources - position changes everytime*/
        float3 vecS = FL3(0.f), vecAB, vecAC, vecN;
        int is, i, ea, eb, ec;
        float bary[4] = {0.f};

        for (is = 0; is < GPU_PARAM(gcfg, srcelemlen); is++) {
            int include = 1;
            __global int* elems = elem + (srcelem[is] - 1) * GPU_PARAM(gcfg, elemlen);

            for (i = 0; i < 4; i++) {
                ea = elems[out[i][0]] - 1;
                eb = elems[out[i][1]] - 1;
                ec = elems[out[i][2]] - 1;
                vecAB = FL4_3(node[eb]) - FL4_3(node[ea]);
                vecAC = FL4_3(node[ec]) - FL4_3(node[ea]);
                vecS = r->p0 - FL4_3(node[ea]);
                vecN = cross(vecAB, vecAC);
                bary[facemap[i]] = -dot(vecS, vecN);
            }

            for (i = 0; i < 4; i++) {
                if (bary[i] < -1e-4f) {
                    include = 0;
                }
            }

            if (include) {
                r->eid = srcelem[is];
                float s = 0.f;

                for (i = 0; i < 4; i++) {
                    s += bary[i];
                }

                for (i = 0; i < 4; i++) {
                    if ((bary[i] / s) < 1e-4f) {
                        r->faceid = ifacemap[i] + 1;
                    }
                }

                break;
            }
        }

#ifdef __NVCC__
    }

#endif
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

MMC_TEMPLATE
__device__ void onephoton(unsigned int id, __local float* ppath, __constant MCXParam* gcfg, __global FLOAT3* node, __global int* elem, __global float* weight, __global float* dref,
                          __global int* type, __global int* facenb,  __global int* srcelem, __global float4* normal, __constant Medium* gmed,
                          __global float* gnodemua, __global float* gnodemusp,
                          __global float* n_det, __global uint* detectedphoton, __local float* energytot, __local float* energyesc, __private RandType* ran, int* raytet, __global float* srcpattern,
                          __global float* replayweight, __global float* replaytime, __global RandType* photonseed, __global MCXReporter* reporter, __global float* gdebugdata) {

    int oldeid, fixcount = 0;
    ray r = {gcfg->srcpos, gcfg->srcdir, {MMC_UNDEFINED, 0.f, 0.f}, GPU_PARAM(gcfg, e0), 0, 0, 1.f, 0.f, 0.f, 0.f, 0.f, ID_UNDEFINED, 0.f, 0.f};
#if defined(MCX_SAVE_SEED) || defined(__NVCC__)
    RandType initseed[RAND_BUF_LEN] = {NULL};
#endif

    clearpath(ppath, (GPU_PARAM(gcfg, reclen) + (GPU_PARAM(gcfg, srcnum) > 1) * GPU_PARAM(gcfg, srcnum)));
    r.photonid = id;

#if defined(MCX_SAVE_SEED) || defined(__NVCC__)

#ifdef __NVCC__

    if (GPU_PARAM(gcfg, issaveseed)) {
#endif

        for (oldeid = 0; oldeid < RAND_BUF_LEN; oldeid++) {
            initseed[oldeid] = ran[oldeid];
        }

#ifdef __NVCC__
    }

#endif

#endif

    /*initialize the photon parameters*/
    launchnewphoton MMC_TARGS (gcfg, &r, node, elem, srcelem, ran, srcpattern, gmed);

#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)
#ifdef __NVCC__

    if (SAVE_DETPHOTON) {
#endif

        if ((GPU_PARAM(gcfg, srctype) != stPattern) + (GPU_PARAM(gcfg, srcnum) == 1)) {
            ppath[GPU_PARAM(gcfg, reclen) - 1] = r.weight; /*last record in partialpath is the initial photon weight*/
        } else if (GPU_PARAM(gcfg, srctype) == stPattern) {
            *((__local uint*)(ppath + GPU_PARAM(gcfg, reclen) - 1)) = r.posidx;
        }

#ifdef __NVCC__
    }

#endif

#endif

    if (!IS_MULTISRC || GPU_PARAM(gcfg, srcnum) == 1) {
        *energytot += r.weight;
    } else {
        for (oldeid = 0; oldeid < GPU_PARAM(gcfg, srcnum); oldeid++) {
            ppath[GPU_PARAM(gcfg, reclen) + oldeid] = srcpattern[r.posidx * GPU_PARAM(gcfg, srcnum) + oldeid];
            energytot[oldeid] += r.weight * ppath[GPU_PARAM(gcfg, reclen) + oldeid];
        }
    }

    if (GPU_PARAM(gcfg, debuglevel) & dlTraj) {
        savedebugdata(&r, id, reporter, gdebugdata, gcfg);
    }

    /*use Kahan summation to accumulate weight, otherwise, counter stops at 16777216*/
    /*http://stackoverflow.com/questions/2148149/how-to-sum-a-large-number-of-float-number*/

    while (1) { /*propagate a photon until exit*/
        r.slen = branchless_badouel_raytet MMC_TARGS (&r, gcfg, ppath, elem, weight, type[r.eid - 1], facenb, normal, gmed, gnodemua, gnodemusp, replayweight, replaytime);
        (*raytet)++;

        if (r.pout.x == MMC_UNDEFINED) {
            if (r.faceid == -2) {
                break;    /*reaches the time limit*/
            }

            if (fixcount++ < MAX_TRIAL) {
                fixphoton(&r.p0, node, (__global int*)(elem + (r.eid - 1)*GPU_PARAM(gcfg, elemlen)));
                continue;
            }

            r.eid = ID_UNDEFINED;
            r.faceid = -1;
        }

#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)

        if (SAVE_DETPHOTON && r.Lmove > 0.f && type[r.eid - 1] > 0) {
            ppath[GPU_PARAM(gcfg, maxmedia) + type[r.eid - 1] - 1] += r.Lmove;    /*second medianum block is the partial path*/
        }

#endif

        /*move a photon until the end of the current scattering path*/
        while (r.faceid >= 0 && !r.isend) {
            r.p0 = r.pout;

            oldeid = r.eid;
            r.eid = ((__global int*)(facenb + (r.eid - 1) * GPU_PARAM(gcfg, elemlen)))[r.faceid];
#ifdef MCX_DO_REFLECTION

            if (GPU_PARAM(gcfg, isreflect) && (r.eid <= 0 || (r.eid > 0 && gmed[type[r.eid - 1]].n != gmed[type[oldeid - 1]].n ))) {
                if (! (r.eid <= 0 && (((gmed[type[oldeid - 1]].n == GPU_PARAM(gcfg, nout)) * (GPU_PARAM(gcfg, isreflect) != (int)bcMirror)) + (GPU_PARAM(gcfg, isreflect) == (int)bcAbsorbExterior)) )) {
                    reflectray(gcfg, &r.vec, &oldeid, &r.eid, r.faceid, ran, type, normal, gmed);
                }
            }

#endif

            if (r.eid <= 0) {
                break;
            }

            /*when a photon enters the domain from the background*/
            if (type[oldeid - 1] == 0 && type[r.eid - 1]) {
                //if(GPU_PARAM(gcfg,debuglevel)&dlExit)
                GPUDEBUG(("e %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                          r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid));

                if (!GPU_PARAM(gcfg, voidtime)) {
                    r.photontimer = 0.f;
                }
            }

            /*when a photon exits the domain into the background*/
            if (type[oldeid - 1] && type[r.eid - 1] == 0) {
                //if(GPU_PARAM(gcfg,debuglevel)&dlExit)
                GPUDEBUG(("x %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                          r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid));

                if (!GPU_PARAM(gcfg, isextdet)) {
                    r.eid = 0;
                    break;
                }
            }

            //          if(r.eid==0 && gmed[type[oldeid-1]].n == GPU_PARAM(gcfg,nout) ) break;
            if (r.pout.x != MMC_UNDEFINED) { // && (GPU_PARAM(gcfg,debuglevel)&dlMove))
                GPUDEBUG(("P %f %f %f %d %u %e\n", r.pout.x, r.pout.y, r.pout.z, r.eid, id, r.slen));
            }

            r.slen = branchless_badouel_raytet MMC_TARGS (&r, gcfg, ppath, elem, weight, type[r.eid - 1], facenb, normal, gmed, gnodemua, gnodemusp, replayweight, replaytime);
            (*raytet)++;
#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)

            if (SAVE_DETPHOTON && r.Lmove > 0.f && type[r.eid - 1] > 0) {
                ppath[GPU_PARAM(gcfg, maxmedia) + type[r.eid - 1] - 1] += r.Lmove;
            }

#endif

            if (r.faceid == -2) {
                break;
            }

            fixcount = 0;

            while (r.pout.x == MMC_UNDEFINED && fixcount++ < MAX_TRIAL) {
                fixphoton(&r.p0, node, (__global int*)(elem + (r.eid - 1)*GPU_PARAM(gcfg, elemlen)));
                r.slen = branchless_badouel_raytet MMC_TARGS (&r, gcfg, ppath, elem, weight, type[r.eid - 1], facenb, normal, gmed, gnodemua, gnodemusp, replayweight, replaytime);
                (*raytet)++;
#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)

                if (SAVE_DETPHOTON && r.Lmove > 0.f && type[r.eid - 1] > 0) {
                    ppath[GPU_PARAM(gcfg, maxmedia) + type[r.eid - 1] - 1] += r.Lmove;
                }

#endif
            }

            if (r.pout.x == MMC_UNDEFINED) {
                /*possibily hit an edge or miss*/
                r.eid = ID_UNDEFINED;
                break;
            }
        }

        if (r.eid <= 0 || r.pout.x == MMC_UNDEFINED) {
            //if(r.eid==0 && (GPU_PARAM(gcfg,debuglevel)&dlMove))
            GPUDEBUG(("B %f %f %f %d %u %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen));

            if (r.eid != ID_UNDEFINED) {
                //if(GPU_PARAM(gcfg,debuglevel)&dlExit)
                GPUDEBUG(("E %f %f %f %f %f %f %f %d\n", r.p0.x, r.p0.y, r.p0.z,
                          r.vec.x, r.vec.y, r.vec.z, r.weight, r.eid));
#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)

                if (SAVE_DETPHOTON * GPU_PARAM(gcfg, issaveexit)) {                                 /*when issaveexit is set to 1*/
                    copystate(ppath + (GPU_PARAM(gcfg, reclen) - 7), (__private float*) & (r.p0), 3); /*columns 7-5 from the right store the exit positions*/
                    copystate(ppath + (GPU_PARAM(gcfg, reclen) - 4), (__private float*) & (r.vec), 3); /*columns 4-2 from the right store the exit dirs*/
                }

#endif
#ifdef MCX_SAVE_DREF

                if (GPU_PARAM(gcfg, issaveref) && r.eid < 0 && dref) {
                    int tshift = MIN( ((int)((r.photontimer - gcfg->tstart) * GPU_PARAM(gcfg, Rtstep))), GPU_PARAM(gcfg, maxgate) - 1 ) * GPU_PARAM(gcfg, nf);
                    dref[((-r.eid) - 1) + tshift] += r.weight;
                }

#endif
            } else if ((r.faceid == -2) * (GPU_PARAM(gcfg, debuglevel)&dlMove)) {
                GPUDEBUG(("T %f %f %f %d %u %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen));
            } else if ((r.eid && r.faceid != -2) * (GPU_PARAM(gcfg, debuglevel)&dlEdge)) {
                GPUDEBUG(("X %f %f %f %d %u %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen));
            }

#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)
#ifdef __NVCC__

            if (SAVE_DETPHOTON) {
#endif

                if (r.eid <= 0) {

#if defined(MCX_SAVE_SEED) || defined(__NVCC__)
                    savedetphoton MMC_TARGS (n_det, detectedphoton, ppath, &r, gmed, ((GPU_PARAM(gcfg, isextdet) && type[oldeid - 1] == GPU_PARAM(gcfg, maxmedia) + 1) ? oldeid : -1), gcfg, photonseed, initseed);
#else
                    savedetphoton MMC_TARGS (n_det, detectedphoton, ppath, &r, gmed, ((GPU_PARAM(gcfg, isextdet) && type[oldeid - 1] == GPU_PARAM(gcfg, maxmedia) + 1) ? oldeid : -1), gcfg, photonseed, NULL);
#endif
                }

#ifdef __NVCC__
            }

#endif

#endif
            break;  /*photon exits boundary*/
        }

        //if(GPU_PARAM(gcfg,debuglevel)&dlMove)
        GPUDEBUG(("M %f %f %f %d %u %e\n", r.p0.x, r.p0.y, r.p0.z, r.eid, id, r.slen));

        {
            /* For RF forward, use the magnitude of the complex weight for roulette threshold.
             * The real part can be negative when cos(phi_accumulated)<0, incorrectly triggering
             * roulette for photons that still carry significant energy.
             * MCX reference: p.w = sqrt(w_re^2+w_im^2) is always used for the roulette check. */
            float roulette_w = r.weight;

            if (IS_RF) {
                roulette_w = MCX_MATHFUN(sqrt)(r.weight * r.weight + r.weight_im * r.weight_im);
            }

            if ((GPU_PARAM(gcfg, minenergy) > 0.f) * (roulette_w < GPU_PARAM(gcfg, minenergy)) * ((gcfg->tend - gcfg->tstart)*GPU_PARAM(gcfg, Rtstep) <= 1.f)) { /*Russian Roulette*/
                if (rand_do_roulette(ran)*GPU_PARAM(gcfg, roulettesize) <= 1.f) {
                    r.weight *= GPU_PARAM(gcfg, roulettesize);

                    if (IS_RF) {
                        r.weight_im *= GPU_PARAM(gcfg, roulettesize);  /* scale imaginary weight too, consistent with MCX */
                    }

                    //if(GPU_PARAM(gcfg,debuglevel)&dlWeight)
                    GPUDEBUG(("Russian Roulette bumps r.weight to %f\n", r.weight));
                } else {
                    break;
                }
            }
        }

        float mom = 0.f;
        r.slen0 = mc_next_scatter(gmed[type[r.eid - 1]].g, &r.vec, ran, gcfg, &mom);
        r.slen = r.slen0;

        if (GPU_PARAM(gcfg, debuglevel) & dlTraj) {
            savedebugdata(&r, id, reporter, gdebugdata, gcfg);
        }

#if defined(MCX_SAVE_DETECTORS) || defined(__NVCC__)

        if (SAVE_DETPHOTON) {
            if (GPU_PARAM(gcfg, ismomentum) && type[r.eid - 1] > 0) {             /*when ismomentum is set to 1*/
                ppath[(GPU_PARAM(gcfg, maxmedia) << 1) + type[r.eid - 1] - 1] += mom;    /*the third medianum block stores the momentum transfer*/
            }

            if (SAVE_DETPHOTON) {
                ppath[type[r.eid - 1] - 1] += 1.f;    /*the first medianum block stores the scattering event counts*/
            }
        }

#endif
    }

    if (GPU_PARAM(gcfg, debuglevel) & dlTraj) {
        savedebugdata(&r, id, reporter, gdebugdata, gcfg);
    }

    /* RF forward: photon weight is complex; the escaped energy is the magnitude
     * |w|, which is phase-rotation invariant.  Tracking Re(w) only biases low as
     * the phase advances along the path.  MCX uses the magnitude (see p.w =
     * sqrt(w_re^2 + w_im^2) in its kernel), so MMC mirrors that here for omega>0.
     * Replay (SEED_FROM_FILE) is excluded because the imaginary track is not
     * driven in that path. */
    float esc_mag = IS_RF
                    ? sqrt(r.weight * r.weight + r.weight_im * r.weight_im)
                    : r.weight;

    if (!IS_MULTISRC || GPU_PARAM(gcfg, srcnum) == 1) {
        *energyesc += esc_mag;
    } else {
        for (oldeid = 0; oldeid < GPU_PARAM(gcfg, srcnum); oldeid++) {
            energyesc[oldeid] += esc_mag * ppath[GPU_PARAM(gcfg, reclen) + oldeid];
        }
    }
}

MMC_TEMPLATE
__kernel MMC_LAUNCH_BOUNDS void mmc_main_loop(const int nphoton, const int ophoton,
#ifndef __NVCC__
    __constant__ MCXParam* gcfg, __local float* sharedmem, __constant__ Medium* gmed,
#endif
        __global FLOAT3* node, __global int* elem,  __global float* weight, __global float* dref, __global int* type, __global int* facenb,  __global int* srcelem, __global float4* normal,
        __global float* gnodemua, __global float* gnodemusp,
        __global float* n_det, __global uint* detectedphoton,
        __global uint* n_seed, __global int* progress, __global float* energy, __global MCXReporter* reporter, __global float* srcpattern,
        __global float* replayweight, __global float* replaytime, __global RandType* replayseed, __global RandType* photonseed, __global float* gdebugdata) {

    RandType t[RAND_BUF_LEN];
    int idx = get_global_id(0);
    int raytet = 0;

#ifdef __NVCC__
    extern __shared__ float sharedmem[];
#endif

    if (GPU_PARAM(gcfg, seed) != SEED_FROM_FILE) {
        gpu_rng_init(t, n_seed, idx);
    }

    clearpath(sharedmem + get_local_id(0) * GPU_PARAM(gcfg, srcnum), GPU_PARAM(gcfg, srcnum));
    clearpath(sharedmem + (get_local_size(0) + get_local_id(0)) * GPU_PARAM(gcfg, srcnum), GPU_PARAM(gcfg, srcnum));

    /*launch photons*/
    for (int i = 0; i < nphoton + (idx < ophoton); i++) {
        if (GPU_PARAM(gcfg, seed) == SEED_FROM_FILE)
            for (int j = 0; j < RAND_BUF_LEN; j++) {
                t[j] = replayseed[(idx * nphoton + MIN(idx, ophoton) + i) * RAND_BUF_LEN + j];
            }

        onephoton MMC_TARGS (idx * nphoton + MIN(idx, ophoton) + i, sharedmem + get_local_size(0) * (GPU_PARAM(gcfg, srcnum) << 1) +
                             get_local_id(0) * (GPU_PARAM(gcfg, reclen) + (GPU_PARAM(gcfg, srcnum) > 1) * GPU_PARAM(gcfg, srcnum)), gcfg, node, elem,
                             weight, dref, type, facenb, srcelem, normal, gmed,
                             gnodemua, gnodemusp,
                             n_det, detectedphoton, sharedmem + get_local_id(0) * GPU_PARAM(gcfg, srcnum),
                             sharedmem + (get_local_size(0) + get_local_id(0)) * GPU_PARAM(gcfg, srcnum), t, &raytet,
                             srcpattern, replayweight, replaytime, photonseed, reporter, gdebugdata);
    }

    for (int i = 0; i < GPU_PARAM(gcfg, srcnum); i++) {
        energy[(idx << 1) * GPU_PARAM(gcfg, srcnum) + i] += sharedmem[(get_local_size(0) + get_local_id(0)) * GPU_PARAM(gcfg, srcnum) + i];
        energy[((idx << 1) + 1) * GPU_PARAM(gcfg, srcnum) + i] += sharedmem[get_local_id(0) * GPU_PARAM(gcfg, srcnum) + i];
    }

    if (GPU_PARAM(gcfg, debuglevel) & MCX_DEBUG_PROGRESS && progress) {
        atomic_inc(progress);
    }

    atomicadd(&(reporter->raytet), raytet);
}

/*============================================================================*/
/* Adjoint Jacobian helper kernels (grid / BLBadouelGrid mode)                */
/*============================================================================*/

/**
 * @brief CW time-gate summation for one source/detector slot at a given voxel
 *
 * Sums all time gates for a given source-slot index to obtain the CW fluence.
 * The field layout is: field[vox + (tgate + slot*maxgate)*dimxyz]
 */
#ifdef __NVCC__
__device__ static inline float mmc_cw_sum(const float* field, unsigned int vox, unsigned int slot,
        unsigned int maxgate, unsigned int dimxyz) {
#else
static inline float mmc_cw_sum(__global const float* field, unsigned int vox, unsigned int slot,
                               unsigned int maxgate, unsigned int dimxyz) {
#endif
    float sum = 0.f;

    for (unsigned int t = 0; t < maxgate; t++) {
        sum += field[vox + (size_t)(t + slot * maxgate) * dimxyz];
    }

    return sum;
}

/**
 * @brief 2nd-order finite-difference spatial gradient along one axis for adjoint Jacobian
 */
#ifdef __NVCC__
__device__ static inline float mmc_fd_grad(const float* field, unsigned int vox, unsigned int slot,
        unsigned int maxgate, unsigned int dimxyz,
        unsigned int i, unsigned int N, unsigned int stride) {
#else
static inline float mmc_fd_grad(__global const float* field, unsigned int vox, unsigned int slot,
                                unsigned int maxgate, unsigned int dimxyz,
                                unsigned int i, unsigned int N, unsigned int stride) {
#endif

    if (N <= 1) {
        return 0.f;
    }

    float f0 = mmc_cw_sum(field, vox, slot, maxgate, dimxyz);

    if (i == 0) {
        float fp1 = mmc_cw_sum(field, vox + stride, slot, maxgate, dimxyz);

        if (N == 2) {
            return fp1 - f0;
        }

        float fp2 = mmc_cw_sum(field, vox + 2 * stride, slot, maxgate, dimxyz);
        return (-3.f * f0 + 4.f * fp1 - fp2) * 0.5f;
    } else if (i == N - 1) {
        float fm1 = mmc_cw_sum(field, vox - stride, slot, maxgate, dimxyz);

        if (N == 2) {
            return f0 - fm1;
        }

        float fm2 = mmc_cw_sum(field, vox - 2 * stride, slot, maxgate, dimxyz);
        return (fm2 - 4.f * fm1 + 3.f * f0) * 0.5f;
    } else {
        float fp1 = mmc_cw_sum(field, vox + stride, slot, maxgate, dimxyz);
        float fm1 = mmc_cw_sum(field, vox - stride, slot, maxgate, dimxyz);
        return (fp1 - fm1) * 0.5f;
    }
}

/**
 * @brief Adjoint mua Jacobian kernel: J[vox,s,d] = phi_src[vox,s] * phi_det[vox,d]
 *
 * For RF mode (gfield_im != NULL):
 *   Re(J) = Re(phi_src)*Re(phi_det) - Im(phi_src)*Im(phi_det)
 *   Im(J) = Re(phi_src)*Im(phi_det) + Im(phi_src)*Re(phi_det)
 *
 * Output layout: gadjoint[vox + (s*Nd+d)*dimxyz]        (real)
 *                gadjoint[vox + (s*Nd+d)*dimxyz + Ns*Nd*dimxyz] (imag, RF only)
 */
#ifdef __NVCC__
__global__ void mmc_adjoint_kernel(float* gfield_re, float* gfield_im, float* gadjoint,
                                   unsigned int dimxyz, unsigned int maxgate,
                                   unsigned int Ns, unsigned int Nd) {
    unsigned int vox = blockIdx.x * blockDim.x + threadIdx.x;
#else
__kernel void mmc_adjoint_kernel(__global float* gfield_re, __global float* gfield_im,
                                 __global float* gadjoint,
                                 unsigned int dimxyz, unsigned int maxgate,
                                 unsigned int Ns, unsigned int Nd) {
    unsigned int vox = get_global_id(0);
#endif

    if (vox >= dimxyz) {
        return;
    }

    size_t adjointlen = (size_t)dimxyz * Ns * Nd;

    for (unsigned int s = 0; s < Ns; s++) {
        float cw_src_re = mmc_cw_sum(gfield_re, vox, s, maxgate, dimxyz);
        float cw_src_im = (gfield_im != 0) ? mmc_cw_sum(gfield_im, vox, s, maxgate, dimxyz) : 0.f;

        for (unsigned int d = 0; d < Nd; d++) {
            float cw_det_re = mmc_cw_sum(gfield_re, vox, Ns + d, maxgate, dimxyz);

            unsigned int out_idx = vox + (unsigned int)((size_t)(s * Nd + d) * dimxyz);
            gadjoint[out_idx] = cw_src_re * cw_det_re;

            if (gfield_im != 0) {
                float cw_det_im = mmc_cw_sum(gfield_im, vox, Ns + d, maxgate, dimxyz);
                gadjoint[out_idx] -= cw_src_im * cw_det_im;
                gadjoint[out_idx + (unsigned int)(adjointlen)] = cw_src_re * cw_det_im + cw_src_im * cw_det_re;
            }
        }
    }
}

/**
 * @brief Adjoint D-coefficient Jacobian kernel: J_D[vox,s,d] = nabla(phi_src) . nabla(phi_det)
 *
 * For RF mode:
 *   Re(J_D) = Re(grad_s).Re(grad_d) - Im(grad_s).Im(grad_d)
 *   Im(J_D) = Re(grad_s).Im(grad_d) + Im(grad_s).Re(grad_d)
 */
#ifdef __NVCC__
__global__ void mmc_adjoint_dcoeff_kernel(float* gfield_re, float* gfield_im, float* gadjoint,
        unsigned int dimxyz, unsigned int maxgate,
        unsigned int Ns, unsigned int Nd,
        unsigned int Nx, unsigned int Ny) {
    unsigned int vox = blockIdx.x * blockDim.x + threadIdx.x;
#else
__kernel void mmc_adjoint_dcoeff_kernel(__global float* gfield_re, __global float* gfield_im,
                                        __global float* gadjoint,
                                        unsigned int dimxyz, unsigned int maxgate,
                                        unsigned int Ns, unsigned int Nd,
                                        unsigned int Nx, unsigned int Ny) {
    unsigned int vox = get_global_id(0);
#endif

    if (vox >= dimxyz) {
        return;
    }

    unsigned int Nxy = Nx * Ny;
    unsigned int Nz  = dimxyz / Nxy;
    unsigned int ix  = vox % Nx;
    unsigned int iy  = (vox / Nx) % Ny;
    unsigned int iz  = vox / Nxy;

    size_t adjointlen = (size_t)dimxyz * Ns * Nd;

    for (unsigned int s = 0; s < Ns; s++) {
        float gsx_re = mmc_fd_grad(gfield_re, vox, s, maxgate, dimxyz, ix, Nx, 1);
        float gsy_re = mmc_fd_grad(gfield_re, vox, s, maxgate, dimxyz, iy, Ny, Nx);
        float gsz_re = mmc_fd_grad(gfield_re, vox, s, maxgate, dimxyz, iz, Nz, Nxy);
        float gsx_im = 0.f, gsy_im = 0.f, gsz_im = 0.f;

        if (gfield_im != 0) {
            gsx_im = mmc_fd_grad(gfield_im, vox, s, maxgate, dimxyz, ix, Nx, 1);
            gsy_im = mmc_fd_grad(gfield_im, vox, s, maxgate, dimxyz, iy, Ny, Nx);
            gsz_im = mmc_fd_grad(gfield_im, vox, s, maxgate, dimxyz, iz, Nz, Nxy);
        }

        for (unsigned int d = 0; d < Nd; d++) {
            float gdx_re = mmc_fd_grad(gfield_re, vox, Ns + d, maxgate, dimxyz, ix, Nx, 1);
            float gdy_re = mmc_fd_grad(gfield_re, vox, Ns + d, maxgate, dimxyz, iy, Ny, Nx);
            float gdz_re = mmc_fd_grad(gfield_re, vox, Ns + d, maxgate, dimxyz, iz, Nz, Nxy);

            unsigned int out_idx = vox + (unsigned int)((size_t)(s * Nd + d) * dimxyz);
            gadjoint[out_idx] = gsx_re * gdx_re + gsy_re * gdy_re + gsz_re * gdz_re;

            if (gfield_im != 0) {
                float gdx_im = mmc_fd_grad(gfield_im, vox, Ns + d, maxgate, dimxyz, ix, Nx, 1);
                float gdy_im = mmc_fd_grad(gfield_im, vox, Ns + d, maxgate, dimxyz, iy, Ny, Nx);
                float gdz_im = mmc_fd_grad(gfield_im, vox, Ns + d, maxgate, dimxyz, iz, Nz, Nxy);

                gadjoint[out_idx] -= gsx_im * gdx_im + gsy_im * gdy_im + gsz_im * gdz_im;
                gadjoint[out_idx + (unsigned int)(adjointlen)] =
                    gsx_re * gdx_im + gsy_re * gdy_im + gsz_re * gdz_im
                    + gsx_im * gdx_re + gsy_im * gdy_re + gsz_im * gdz_re;
            }
        }
    }
}

/**
 * @brief Mesh-mode adjoint Jacobian (full FEM form), one thread per element.
 *
 * Implements the rb_femjacobian formula on a tet mesh given nodal fluences
 * phi_s, phi_d. Per element t with 4 nodes ee[0..3] and volume Ve:
 *
 *   J_mua(t) = -0.1*Ve * [ Σ_i φ_s(ee[i])·φ_d(ee[i])
 *                          + 0.5 * Σ_{i<j} (φ_s(ee[i])·φ_d(ee[j]) + φ_s(ee[j])·φ_d(ee[i])) ]
 *
 *   J_D(t)   = -[ Σ_i deldotdel[t][diag_i] · φ_s(ee[i])·φ_d(ee[i])
 *                + Σ_{i<j} deldotdel[t][off_ij] · (φ_s(ee[i])·φ_d(ee[j]) + φ_s(ee[j])·φ_d(ee[i])) ]
 *
 * Nodal-fluence layout (matches mmc_adjoint_kernel's grid convention):
 *   field[node + (gate + slot*maxgate)*nn]
 *
 * Output:
 *   if isnodal == 0: g_jmua_e/g_jd_e[t + (s*Nd+d)*ne]
 *   if isnodal != 0: g_jmua_n/g_jd_n[n + (s*Nd+d)*nn] (atomic scatter, *0.25)
 *
 * RF: gfield_im non-NULL; real/imag Jacobian written to first/second adjointlen blocks.
 *
 * @param[in]  gfield_re   nodal real fluence [nn * maxgate * nsrcslots]
 * @param[in]  gfield_im   nodal imag fluence (RF only) or NULL
 * @param[in]  gelem       per-element node indices (1-based), shape [ne * elemlen]
 * @param[in]  gevol       per-element volume (float, length ne)
 * @param[in]  gdeldotdel  per-element ⟨∇φ_i·∇φ_j⟩*Ve, packed upper-triangle [ne*10]
 * @param[out] gjmua       J_mua output (NULL to skip)
 * @param[out] gjd         J_D   output (NULL to skip)
 * @param[in]  ne, nn      mesh sizes
 * @param[in]  maxgate     number of time gates
 * @param[in]  Ns, Nd      number of source and detector slots
 * @param[in]  elemlen     usually 4 for tets
 * @param[in]  isnodal     0=element-based output; nonzero=nodal output (atomic scatter)
 */
#ifdef __NVCC__
__global__ void mmc_adjoint_mesh_full_kernel(float* gfield_re, float* gfield_im,
        int* gelem, float* gevol, float* gdeldotdel,
        float* gjmua, float* gjd,
        unsigned int ne, unsigned int nn, unsigned int maxgate,
        unsigned int Ns, unsigned int Nd, unsigned int elemlen,
        int isnodal) {
    unsigned int t = blockIdx.x * blockDim.x + threadIdx.x;
#else
__kernel void mmc_adjoint_mesh_full_kernel(__global float* gfield_re, __global float* gfield_im,
        __global int* gelem, __global float* gevol, __global float* gdeldotdel,
        __global float* gjmua, __global float* gjd,
        unsigned int ne, unsigned int nn, unsigned int maxgate,
        unsigned int Ns, unsigned int Nd, unsigned int elemlen,
        int isnodal) {
    unsigned int t = get_global_id(0);
#endif

    if (t >= ne) {
        return;
    }

    /* upper-triangle pair indices into the packed 10-entry deldotdel row:
     * order is [00,01,02,03,11,12,13,22,23,33] */
    const int diag_idx[4] = {0, 4, 7, 9};
    const int off_idx[6]  = {1, 2, 3, 5, 6, 8};
    const int pair_a[6]   = {0, 0, 0, 1, 1, 2};
    const int pair_b[6]   = {1, 2, 3, 2, 3, 3};

    int ee[4];

    for (int k = 0; k < 4; k++) {
        ee[k] = gelem[t * elemlen + k] - 1;  /* 1-based -> 0-based */
    }

    float Ve = gevol[t];

    size_t adjointlen_e = (size_t)ne * Ns * Nd;
    size_t adjointlen_n = (size_t)nn * Ns * Nd;
    int isrf = (gfield_im != 0) ? 1 : 0;

    for (unsigned int s = 0; s < Ns; s++) {
        float phisr[4], phisi[4];

        for (int k = 0; k < 4; k++) {
            phisr[k] = mmc_cw_sum(gfield_re, ee[k], s, maxgate, nn);
            phisi[k] = isrf ? mmc_cw_sum(gfield_im, ee[k], s, maxgate, nn) : 0.f;
        }

        for (unsigned int d = 0; d < Nd; d++) {
            unsigned int slot = Ns + d;
            float phidr[4], phidi[4];

            for (int k = 0; k < 4; k++) {
                phidr[k] = mmc_cw_sum(gfield_re, ee[k], slot, maxgate, nn);
                phidi[k] = isrf ? mmc_cw_sum(gfield_im, ee[k], slot, maxgate, nn) : 0.f;
            }

            /* assemble Re/Im of (phi_s * phi_d) products at each node-pair (float arithmetic
             * is fine here — MC noise dominates well above single-precision roundoff) */
            float jmua_re = 0.f, jmua_im = 0.f;
            float jd_re   = 0.f, jd_im   = 0.f;

            /* diagonal terms: i == j */
            for (int i = 0; i < 4; i++) {
                float pre = phisr[i] * phidr[i] - phisi[i] * phidi[i];
                float pim = isrf ? (phisr[i] * phidi[i] + phisi[i] * phidr[i]) : 0.f;
                jmua_re += pre;
                jmua_im += pim;

                if (gjd) {
                    float w = gdeldotdel[(size_t)t * 10 + diag_idx[i]];
                    jd_re += w * pre;
                    jd_im += w * pim;
                }
            }

            /* off-diagonal pairs (i<j): use both orderings phi_s[i]*phi_d[j] + phi_s[j]*phi_d[i] */
            for (int p = 0; p < 6; p++) {
                int a = pair_a[p], b = pair_b[p];
                float pre = phisr[a] * phidr[b] + phisr[b] * phidr[a]
                            - phisi[a] * phidi[b] - phisi[b] * phidi[a];
                float pim = isrf ? (phisr[a] * phidi[b] + phisr[b] * phidi[a]
                                    + phisi[a] * phidr[b] + phisi[b] * phidr[a]) : 0.f;
                jmua_re += 0.5f * pre;
                jmua_im += 0.5f * pim;

                if (gjd) {
                    float w = gdeldotdel[(size_t)t * 10 + off_idx[p]];
                    jd_re += w * pre;
                    jd_im += w * pim;
                }
            }

            jmua_re *= -0.1f * Ve;
            jmua_im *= -0.1f * Ve;
            jd_re   *= -1.f;
            jd_im   *= -1.f;

            unsigned int sdpair = s * Nd + d;

            if (!isnodal) {
                /* element-based output: one write per (sd, t) */
                if (gjmua) {
                    size_t k = (size_t)t + (size_t)sdpair * ne;
                    gjmua[k] = jmua_re;

                    if (isrf) {
                        gjmua[k + adjointlen_e] = jmua_im;
                    }
                }

                if (gjd) {
                    size_t k = (size_t)t + (size_t)sdpair * ne;
                    gjd[k] = jd_re;

                    if (isrf) {
                        gjd[k + adjointlen_e] = jd_im;
                    }
                }
            } else {
                /* nodal output: scatter 0.25 * elem-value to each of 4 nodes (atomic) */
                float c = 0.25f;
                float jmua_re_f = jmua_re * c;
                float jmua_im_f = jmua_im * c;
                float jd_re_f   = jd_re   * c;
                float jd_im_f   = jd_im   * c;

                for (int k = 0; k < 4; k++) {
                    size_t base = (size_t)ee[k] + (size_t)sdpair * nn;

                    if (gjmua) {
                        atomicadd(gjmua + base, jmua_re_f);

                        if (isrf) {
                            atomicadd(gjmua + base + adjointlen_n, jmua_im_f);
                        }
                    }

                    if (gjd) {
                        atomicadd(gjd + base, jd_re_f);

                        if (isrf) {
                            atomicadd(gjd + base + adjointlen_n, jd_im_f);
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Mesh-mode adjoint Jacobian (nodal approximation), one thread per node.
 *
 * Implements the nodal-adjoint approximation from PhD thesis eq.
 * (3d3d:adjoint:nodal) / rbjacmuafast.m:
 *
 *   J_mua_n(n) = -nvol[n] * φ_s(n) * φ_d(n)
 *
 * Valid when the forward mesh is fine relative to the parameter mesh.
 * J_D is NOT defined in this approximation; use the full FEM kernel for J_D.
 *
 * @param[in]  gfield_re   nodal real fluence [nn * maxgate * nsrcslots]
 * @param[in]  gfield_im   nodal imag fluence (RF only) or NULL
 * @param[in]  gnvol       per-node Voronoi volume (float, length nn)
 * @param[out] gjmua       J_mua_n output, [nn * Ns * Nd] (CW) plus another such block for imag
 * @param[in]  nn          number of nodes
 * @param[in]  maxgate     number of time gates
 * @param[in]  Ns, Nd      number of source and detector slots
 */
#ifdef __NVCC__
__global__ void mmc_adjoint_mesh_nodal_kernel(float* gfield_re, float* gfield_im,
        float* gnvol, float* gjmua,
        unsigned int nn, unsigned int maxgate,
        unsigned int Ns, unsigned int Nd) {
    unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
#else
__kernel void mmc_adjoint_mesh_nodal_kernel(__global float* gfield_re, __global float* gfield_im,
        __global float* gnvol, __global float* gjmua,
        unsigned int nn, unsigned int maxgate,
        unsigned int Ns, unsigned int Nd) {
    unsigned int n = get_global_id(0);
#endif

    if (n >= nn) {
        return;
    }

    size_t adjointlen = (size_t)nn * Ns * Nd;
    float vol = gnvol[n];
    int isrf = (gfield_im != 0) ? 1 : 0;

    for (unsigned int s = 0; s < Ns; s++) {
        float cw_src_re = mmc_cw_sum(gfield_re, n, s, maxgate, nn);
        float cw_src_im = isrf ? mmc_cw_sum(gfield_im, n, s, maxgate, nn) : 0.f;

        for (unsigned int d = 0; d < Nd; d++) {
            unsigned int slot = Ns + d;
            float cw_det_re = mmc_cw_sum(gfield_re, n, slot, maxgate, nn);

            size_t out_idx = (size_t)n + (size_t)(s * Nd + d) * nn;

            if (isrf) {
                float cw_det_im = mmc_cw_sum(gfield_im, n, slot, maxgate, nn);
                gjmua[out_idx]              = -vol * (cw_src_re * cw_det_re - cw_src_im * cw_det_im);
                gjmua[out_idx + adjointlen] = -vol * (cw_src_re * cw_det_im + cw_src_im * cw_det_re);
            } else {
                gjmua[out_idx] = -vol * cw_src_re * cw_det_re;
            }
        }
    }
}
