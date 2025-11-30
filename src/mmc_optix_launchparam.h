#ifndef _MMC_OPTIX_LAUNCHPARAM_H
#define _MMC_OPTIX_LAUNCHPARAM_H

#include <limits>

#define MAX_PROP           4000              /*maximum property number*/

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

constexpr unsigned int INITIAL_MEDIUM_UNKNOWN = std::numeric_limits<unsigned int>::max();
constexpr unsigned int AMBIENT_MEDIUM_ID = 0;
constexpr unsigned int DEAD_MEDIUM_ID = std::numeric_limits<unsigned int>::max() - 1;

/**
 * @brief struct for medium optical properties
 */
typedef struct __attribute__((aligned(16))) MCX_medium {
    float mua;                     /**<absorption coeff in 1/mm unit*/
    float mus;                     /**<scattering coeff in 1/mm unit*/
    float g;                       /**<anisotropy*/
    float n;                       /**<refractive index*/
} Medium;

/**
 * @brief struct for simulation configuration paramaters
 */
typedef struct __attribute__((aligned(16))) MMC_Parameter {
    OptixTraversableHandle gashandle;

    CUdeviceptr seedbuffer;             /**< rng seed for each thread */
    CUdeviceptr outputbuffer;

    int    srctype;
    float4 srcparam1;
    float4 srcparam2;

    float3 srcpos;
    float3 srcdir;
    float3 nmin;
    float3 nmax;
    uint4 crop0;
    float dstep;
    float tstart, tend;
    float Rtstep;
    int maxgate;
    unsigned int mediumid0;             /**< initial medium type */

    uint isreflect;
    int outputtype;

    int threadphoton;
    int oddphoton;

    Medium medium[MAX_PROP];
} MMCParam;

struct TriangleMeshSBTData {
    /**< x,y,z: face normal; w: front(upper 2-byte) and
     * back(lower 2-byte) face medium type */
    float4* fnorm;
};

#endif
