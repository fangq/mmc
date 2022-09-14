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
\file    simpmesh.h

\brief   Definition of basic mesh data structures and inline vector operations
*******************************************************************************/

#ifndef _MMC_MESH_UNIT_H
#define _MMC_MESH_UNIT_H

#include <stdio.h>
#include <math.h>
#include "mmc_utils.h"

#ifdef MMC_USE_SSE
    #include <smmintrin.h>
#endif

#if !defined(USE_OPENCL) && !defined(__NVCC__)

    #ifdef MMC_SFMT
        #include "mmc_rand_sfmt.c"
    #elif defined MMC_XORSHIFT
        #include "mmc_rand_xorshift128p.c"
    #else
        #include "mmc_rand_posix.c"
    #endif

#elif defined(__NVCC__)
    #include "mmc_rand_xorshift128p.h"
#else
    #include "mmc_rand_xorshift128p.c"
#endif

#define MMC_UNDEFINED (3.40282347e+38F)
#define ID_UNDEFINED  0xFFFFFFFF

#define R_RAND_MAX (1.f/RAND_MAX)
#define TWO_PI     (M_PI*2.0)
#define EPS        1e-6f
#define LOG_MT_MAX 22.1807097779182f
#define R_MIN_MUS  1e9f
#define R_C0       3.335640951981520e-12f  //1/C0 in s/mm
#define DELTA_MUA  1e-4f
#define VERY_BIG   1e30f
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define MAX(a,b)  ((a)>(b)?(a):(b))

#define MESH_ERROR(a)  mesh_error((a),__FILE__,__LINE__)

/***************************************************************************//**
\struct MMC_mesh simpmesh.h
\brief  Basic FEM mesh data structrure

We define nodes, elements, optical property indices, and other related data
related to an FEM mesh.

*******************************************************************************/

typedef struct MMC_mesh {
    int nn;                /**< number of nodes */
    int ne;                /**< number of elements */
    int nf;                /**< number of surface triangles */
    int prop;              /**< number of media */
    int elemlen;           /**< number of nodes per element */
    MMCfloat3* node;       /**< node coordinates */
    int*  elem;            /**< element indices */
    int*  elem2;           /**< element indices */
    float* edgeroi;        /**< immc: vessel edge radii */
    float* faceroi;        /**< immc: face thicknesses */
    float* noderoi;        /**< immc: node radius */
    int*  srcelem;         /**< candidate list of elements containing the source*/
    int  srcelemlen;       /**< length of the elements that may contain the source*/
    int*  detelem;         /**< candidate list of elements containing a widefield detector*/
    int  detelemlen;       /**< length of the elements that may contain the detector*/
    int*  type;            /**< element-based media index */
    int*  facenb;          /**< face neighbors, idx of the element sharing a face */
    medium* med;           /**< optical property of different media */
    float* atte;           /**< precomputed attenuation for each media */
    double* weight;        /**< volumetric fluence for all nodes at all time-gates */
    double* dref;          /**< surface diffuse reflectance */
    float* evol;           /**< volume of an element */
    float* nvol;           /**< voronoi volume of a node */
    float4 nmin;           /**< lower-corner of the mesh bounding box */
    float4 nmax;           /**< upper-corner of the mesh bounding box */
    uint nface;            /**< number of triangular meshes */
    float3* fnode;         /**< triangular mesh nodes */
    uint3* face;           /**< triangular meshes */
    float3* fnorm;         /**< face normal: pointing from back to front */
    uint* front;           /**< front face medium */
    uint* back;            /**< back face medium */
} tetmesh;

/***************************************************************************//**
\struct MMC_raytracer simpmesh.h
\brief  Ray-tracer data structrure for pre-computed data

We define the precomputed data in a ray-tracer structure. For the case of
Plucker-based ray-tracing, this structure contains the displacement and
moment vectors for each edge in a tetrahedron.

*******************************************************************************/

typedef struct MMC_raytracer {
    tetmesh* mesh;          /**< link to the mesh structure */
    char method;            /**< 1 for Plucker-based ray-tracing, 0 for Havel */
    MMCfloat3* d;           /**< precomputed data: for Pluckers, this is displacement */
    MMCfloat3* m;           /**< precomputed data: for Pluckers, this is moment */
    MMCfloat3* n;           /**< precomputed data: for Pluckers, face norm */
} raytracer;
#ifdef  __cplusplus
extern "C" {
#endif
void mesh_init(tetmesh* mesh);
void mesh_init_from_cfg(tetmesh* mesh, mcconfig* cfg);
void mesh_loadnode(tetmesh* mesh, mcconfig* cfg);
void mesh_loadelem(tetmesh* mesh, mcconfig* cfg);
void mesh_loadfaceneighbor(tetmesh* mesh, mcconfig* cfg);
void mesh_loadmedia(tetmesh* mesh, mcconfig* cfg);
void mesh_loadelemvol(tetmesh* mesh, mcconfig* cfg);
void mesh_loadseedfile(tetmesh* mesh, mcconfig* cfg);

void mesh_clear(tetmesh* mesh);
float mesh_normalize(tetmesh* mesh, mcconfig* cfg, float Eabsorb, float Etotal, int pair);
void mesh_build(tetmesh* mesh);
void mesh_error(const char* msg, const char* file, const int linenum);
void mesh_filenames(const char* format, char* foutput, mcconfig* cfg);
void mesh_saveweight(tetmesh* mesh, mcconfig* cfg, int isref);
void mesh_savedetphoton(float* ppath, void* seeds, int count, int seedbyte, mcconfig* cfg);
void mesh_getdetimage(float* detmap, float* ppath, int count, mcconfig* cfg, tetmesh* mesh);
void mesh_savedetimage(float* detmap, mcconfig* cfg);
float mesh_getdetweight(int photonid, int colcount, float* ppath, mcconfig* cfg);
void mesh_srcdetelem(tetmesh* mesh, mcconfig* cfg);
void mesh_createdualmesh(tetmesh* mesh, mcconfig* cfg);
void mesh_loadroi(tetmesh* mesh, mcconfig* cfg);

void tracer_init(raytracer* tracer, tetmesh* mesh, char methodid);
void tracer_build(raytracer* tracer);
void tracer_prep(raytracer* tracer, mcconfig* cfg);
void tracer_clear(raytracer* tracer);

float mc_next_scatter(float g, MMCfloat3* dir, RandType* ran, RandType* ran0, mcconfig* cfg, float* pmom);
#ifdef MCX_CONTAINER
#ifdef __cplusplus
extern "C"
#endif
int mcx_throw_exception(const int id, const char* msg, const char* filename, const int linenum);
#endif

#ifdef  __cplusplus
}
#endif

static inline void vec_add(MMCfloat3* a, MMCfloat3* b, MMCfloat3* res) {
    res->x = a->x + b->x;
    res->y = a->y + b->y;
    res->z = a->z + b->z;
}

static inline void vec_diff(MMCfloat3* a, MMCfloat3* b, MMCfloat3* res) {
    res->x = b->x - a->x;
    res->y = b->y - a->y;
    res->z = b->z - a->z;
}

static inline void vec_mult(MMCfloat3* a, float sa, MMCfloat3* res) {
    res->x = sa * a->x;
    res->y = sa * a->y;
    res->z = sa * a->z;
}

static inline void vec_mult_add(MMCfloat3* a, MMCfloat3* b, float sa, float sb, MMCfloat3* res) {
    res->x = sb * b->x + sa * a->x;
    res->y = sb * b->y + sa * a->y;
    res->z = sb * b->z + sa * a->z;
}

static inline void vec_cross(MMCfloat3* a, MMCfloat3* b, MMCfloat3* res) {
    res->x = a->y * b->z - a->z * b->y;
    res->y = a->z * b->x - a->x * b->z;
    res->z = a->x * b->y - a->y * b->x;
}

static inline void mmc_sincosf(float x, float* sine, float* cosine) {
#if defined(__GNUC__) && defined(__linux__) && !defined(__clang__)
    __builtin_sincosf(x, sine, cosine);
#else
    *sine = sinf(x);
    *cosine = cosf(x);
#endif
}

//#ifndef MMC_USE_SSE
static inline float vec_dot(MMCfloat3* a, MMCfloat3* b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}/*
#else

#ifndef __SSE4_1__
static inline float vec_dot(MMCfloat3 *a,MMCfloat3 *b){
        float dot;
        __m128 na,nb,res;
        na=_mm_load_ps(&a->x);
        nb=_mm_load_ps(&b->x);
        res=_mm_mul_ps(na,nb);
        res=_mm_hadd_ps(res,res);
        res=_mm_hadd_ps(res,res);
        _mm_store_ss(&dot,res);
        return dot;
}
#else
static inline float vec_dot(MMCfloat3 *a,MMCfloat3 *b){
        float dot;
        __m128 na,nb,res;
        na=_mm_load_ps(&a->x);
        nb=_mm_load_ps(&b->x);
        res=_mm_dp_ps(na,nb,0x7f);
        _mm_store_ss(&dot,res);
        return dot;
}
#endif
#endif
*/

static inline float pinner(MMCfloat3* Pd, MMCfloat3* Pm, MMCfloat3* Ad, MMCfloat3* Am) {
    return vec_dot(Pd, Am) + vec_dot(Pm, Ad);
}


static inline float dist2(MMCfloat3* p0, MMCfloat3* p1) {
    return (p1->x - p0->x) * (p1->x - p0->x) + (p1->y - p0->y) * (p1->y - p0->y) + (p1->z - p0->z) * (p1->z - p0->z);
}

static inline float dist(MMCfloat3* p0, MMCfloat3* p1) {
    return sqrtf(dist2(p0, p1));
}

static inline float dist2d2(float* p0, float* p1) {
    return (p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1]);
}

static inline float mmc_rsqrtf(float a) {
#ifdef MMC_USE_SSE
    return _mm_cvtss_f32( _mm_rsqrt_ss( _mm_set_ss( a ) ) );
#else
    return 1.f / sqrtf(a);
#endif
}

#endif
