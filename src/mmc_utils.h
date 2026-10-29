/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2026
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
\file    mmc_utils.h

\brief   Definition of program options and problem domain configurations
*******************************************************************************/

#ifndef _MMC_UTILITIES_H
#define _MMC_UTILITIES_H

#include <stdio.h>
#include "mmc_vector_types.h"
#include "cjson/cJSON.h"

#ifdef _OPENMP                      ///< use multi-threading for running simulation on multiple GPUs
    #include <omp.h>
#endif

#ifdef _MSC_VER                     ///< aligns a variable; placed before its declaration, as MSVC requires
    #define MMC_ALIGN(x) __declspec(align(x))
#else
    #define MMC_ALIGN(x) __attribute__ ((aligned(x)))
#endif

#ifdef __NVCC__
    #define FLOAT3 float3           ///< FLOAT3 is the true float3; in OpenCL, float3 is actually float4; but CUDA supports true float3
#endif

#define MAX_FULL_PATH       2048                         /**< max characters in a full file name string */
#define MAX_PATH_LENGTH     1024                         /**< max characters in a full file name string */
#define MAX_SESSION_LENGTH  64                           /**< max session name length */
#define MAX_CHECKPOINT      16                           /**< max number of photon save points */
#define DET_PHOTON_BUF      100000                       /**< initialize number of detected photons */
#define SEED_FROM_FILE      -999                         /**< special flag indicating to load seeds from history file */
#define MIN(a,b)            ((a)<(b)?(a):(b))            /**< macro to get the min values of two numbers */
#define MMC_ERROR(id,msg)   mcx_error(id,msg,__FILE__,__LINE__)
#define MMC_INFO            -99999
#define MAX_DEVICE          256

#ifndef MCX_CONTAINER
    #define S_RED     "\x1b[31m"
    #define S_GREEN   "\x1b[32m"
    #define S_YELLOW  "\x1b[33m"
    #define S_BLUE    "\x1b[34m"
    #define S_MAGENTA "\x1b[35m"
    #define S_CYAN    "\x1b[36m"
    #define S_BOLD    "\x1b[1m"
    #define S_ITALIC  "\x1b[3m"
    #define S_RESET   "\x1b[0m"
#else
    #define S_RED
    #define S_GREEN
    #define S_YELLOW
    #define S_BLUE
    #define S_MAGENTA
    #define S_CYAN
    #define S_BOLD
    #define S_ITALIC
    #define S_RESET
#endif

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define MAX(a,b)  ((a)>(b)?(a):(b))
#define ABS(a)    ((a)<0?-(a):(a))                        //! Macro to calculate the absolute value

typedef double OutputType;

enum TDebugLevel {dlMove = 1, dlTracing = 2, dlBary = 4, dlWeight = 8, dlDist = 16, dlTracingEnter = 32,
                  dlTracingExit = 64, dlEdge = 128, dlAccum = 256, dlTime = 512, dlReflect = 1024,
                  dlProgress = 2048, dlExit = 4096, dlTraj = 8192
                 };

enum TRTMethod {rtPlucker, rtHavel, rtBadouel, rtBLBadouel, rtBLBadouelGrid};
enum TMCMethod {mmMCX, mmMCML};
enum TComputeBackend {cbSSE, cbOpenCL, cbCUDA, cbOptiX};

enum TSrcType {stPencil, stIsotropic, stCone, stGaussian, stPlanar,
               stPattern, stFourier, stArcSin, stDisk, stFourierX,
               stFourier2D, stZGaussian, stLine, stSlit,
               stPencilArray, stPattern3D, stHyperboloid, stRing
              };
enum TOutputType {otFlux, otFluence, otEnergy, otJacobian, otWL, otWP,
                  otRF, otRFmus,
                  otAdjoint, otAdjointDcoeff, otAdjointMus, otAdjointMusp, otAdjointMuaD, otAdjointMuaMusp
                 };
#define MCX_IS_ADJOINT_TYPE(t)      ((int)(t) >= (int)otAdjoint)
#define MCX_IS_DUAL_ADJOINT_TYPE(t) ((int)(t) >= (int)otAdjointMuaD)
enum TOutputFormat {ofASCII, ofBin, ofNifti, ofAnalyze, ofMC2, ofTX3, ofJNifti, ofBJNifti};
enum TOutputDomain {odMesh, odGrid};
enum TDeviceVendor {dvUnknown, dvNVIDIA, dvAMD, dvIntel, dvIntelGPU, dvAppleCPU, dvAppleGPU};
enum TMCXParent  {mpStandalone, mpMATLAB};
enum TBoundary {bcNoReflect, bcReflect, bcAbsorbExterior, bcMirror /*, bcCylic*/};
enum TRayHitType {htNone, htInOut, htOutIn, htNoHitIn, htNoHitOut};
enum TROIType {rtNone, rtEdge, rtNode, rtFace, rtInner};  /**< rtInner: a crossing inside the ROI union, no medium change */

enum TBJData {JDB_mixed, JDB_nulltype, JDB_noop, JDB_true, JDB_false,
              JDB_char, JDB_string, JDB_hp,
              JDB_int8, JDB_uint8, JDB_int16, JDB_uint16,
              JDB_int32, JDB_uint32, JDB_int64, JDB_uint64,
              JDB_half, JDB_single, JDB_double,
              JDB_array, JDB_object, JDB_numtypes
             };

/***************************************************************************//**
\struct MMC_medium mmc_utils.h

\brief  This structure defines the optical properties for each medium

For each optical medium, we use a 4-element structure to identify its optical
properties, i.e, the absorption coeff (mu_a)in 1/mm, the scattering coeff (mu_s)
in 1/mm, the refractive index (n) and anisotropy (g).
*******************************************************************************/

/**
 * \struct MMC_medium mmc_utils.h
 * \brief The structure to store optical properties
 * Four relevant optical properties are needed
 */

/** \brief Extra source entry for multi-source / adjoint-mode simulation */
#ifndef MCX_EXTRASRC_DEFINED
#define MCX_EXTRASRC_DEFINED
typedef struct MCX_ExtraSrc {
    float4 srcpos;      /**< position and launch-weight (w component = importance) */
    float4 srcdir;      /**< direction and focal length (w component) */
    float4 srcparam1;   /**< source parameters set 1: x=radius for disk source */
    float4 srcparam2;   /**< source parameters set 2 */
} ExtraSrc;
#endif

typedef struct MMC_medium {
    float mua;                     /**<absorption coeff in 1/mm unit*/
    float mus;                     /**<scattering coeff in 1/mm unit*/
    float g;                       /**<anisotropy*/
    float n;                       /**<refractive index*/
} medium;

/**
 * \struct MMC_history mmc_utils.h
 * \brief Filer header data structure to .mch files to store detected photon data
 * This header has a total of 256 bytes
 */

typedef struct MMC_history {
    char magic[4];                 /**< magic bits= 'M','C','X','H' */
    unsigned int  version;         /**< version of the mch file format */
    unsigned int  maxmedia;        /**< number of media in the simulation */
    unsigned int  detnum;          /**< number of detectors in the simulation */
    unsigned int  colcount;        /**< how many output files per detected photon */
    unsigned int  totalphoton;     /**< how many total photon simulated */
    unsigned int  detected;        /**< how many photons are detected (not necessarily all saved) */
    unsigned int  savedphoton;     /**< how many detected photons are saved in this file */
    float unitinmm;                /**< what is the voxel size of the simulation */
    unsigned int  seedbyte;        /**< how many bytes per RNG seed */
    float normalizer;              /**< what is the normalization factor */
    unsigned int  srcnum;          /**< number of sources in the simulation*/
    int respin;                    /**< if positive, repeat count so total photon=totalphoton*respin; if negative, total number is processed in respin subset */
    unsigned int  savedetflag;     /**<  */
    int reserved[2];               /**< reserved fields for future extension */
} history;


typedef struct MCXGPUInfo {
    char name[MAX_SESSION_LENGTH];
    int id;
    int devcount;
    int platformid;
    int major, minor;
    size_t globalmem, constmem, sharedmem;
    int regcount;
    int clock;
    int sm, core;
    size_t autoblock, autothread;
    int maxgate;
    int maxmpthread;               /**< maximum thread number per multi-processor */
    enum TDeviceVendor vendor;
} GPUInfo;

/***************************************************************************//**
\struct MMC_config mmc_utils.h

\brief  This structure defines the problem settings (domain, filenames, session)

*******************************************************************************/

typedef struct MMC_config {
    size_t nphoton;                /**<total simulated photon number, max: 2^64-1*/
    int nblocksize;                /**<thread block size*/
    int nthread;                   /**<num of total threads, multiple of 128*/
    int seed;                      /**<random number generator seed*/
    int e0;                        /**<initial element id*/
    float4 srcpos;                 /**<src position in mm*/
    float4 srcdir;                 /**<src normal direction*/
    int srctype;                   /**<src type: 0 - pencil beam, 1 - isotropic ... */
    float4 srcparam1;              /**<source parameters set 1*/
    float4 srcparam2;              /**<source parameters set 2*/
    float* srcpattern;             /**<source pattern*/
    int voidtime;                  /**<1 start counting photon time when moves inside 0 voxels; 0: count time only after enters non-zero voxel*/
    float4 bary0;                  /**<initial bary centric coordinates of the source*/
    float tstart;                  /**<start time in second*/
    float tstep;                   /**<time step in second*/
    float tend;                    /**<end time in second*/
    FLOAT3 steps;                  /**<voxel sizes along x/y/z in mm; FLOAT3 (not float3, which is 16 bytes in the SSE/OpenCL C units but 12 in the CUDA units) keeps one mcconfig layout across all units*/
    uint3 dim;                     /**<dim.x is the initial element number in MMC, dim.y is faceid*/
    uint4 crop0;                   /**<sub-volume for cache*/
    uint4 crop1;                   /**<the other end of the caching box*/
    int medianum;                  /**<total types of media*/
    int srcnum;                    /**<total number of sources, could be larger than 1 only with pattern illumination*/
    int detnum;                    /**<total detector numbers*/
    float detradius;               /**<detector radius*/
    float sradius;                 /**<source region radius: if set to non-zero, accumulation \
                                       will not perform for dist<sradius; this can reduce \
                                       normalization error when using non-atomic write*/
    medium* prop;                  /**<optical property mapping table*/
    float4* detpos;                /**<detector positions and radius, overwrite detradius*/
    float4* detdir;               /**<detector normal directions and focal lengths for adjoint mode (Nd x 4)*/
    unsigned int* sdmap;          /**<adjoint Jacobian channels: 0-based (source, detector) index pairs, 2*sdmaplen values; NULL for all Ns x Nd pairs*/
    unsigned int sdmaplen;        /**<number of channels (rows) in sdmap*/
    float4 detparam1;              /**<parameters set 1 for wide-field detector*/
    float4 detparam2;              /**<parameters set 2 for wide-feild detector*/
    float* detpattern;             /**<detector pattern*/
    float  minstep;                /**<accumulation step size*/
    int maxgate;                   /**<simultaneous recording gates*/
    int respin;                    /**<number of repeatitions*/
    int printnum;                  /**<number of printed threads (for debugging)*/
    int replaydet;                 /**<the detector id for which to replay the detected photons, start from 1 \
                                       0 for wide-field detection pattern*/
    unsigned char* vol;            /**<pointer to the volume*/
    char session[MAX_SESSION_LENGTH];/**<session id, a string*/
    char meshtag[MAX_SESSION_LENGTH];/**<a string to tag all input mesh files*/
    char isrowmajor;               /**<1 for C-styled array in vol, 0 for matlab-styled array*/
    char isreflect;                /**<1 for reflecting photons at boundary,0 for exiting*/
    char isref3;                   /**<1 considering maximum 3 ref. interfaces; 0 max 2 ref*/
    char isnormalized;             /**<1 to normalize the fluence, 0 for raw fluence*/
    char issavedet;                /**<1 to count all photons hits the detectors*/
    char ismomentum;               /**<1 to save momentum transfer for detected photons, implies issavedet=1*/
    char issaveexit;               /**<1 to save the exit position and vector of a detected photon, implies issavedet=1 \
                                       2 to save accumulated photon weight in frames of images*/
    char issave2pt;                /**<1 to save the 2-point distribution, 0 do not save*/
    char isgpuinfo;                /**<1 to print gpu info when attach, 0 do not print*/
    char isspecular;               /**<1 calculate the initial specular ref if outside the mesh, 0 do not calculate*/
    char issaveseed;               /**<1 save the seed for a detected photon, 0 do not save*/
    char issaveref;                /**<1 to save diffuse reflectance on surface, 0 no save*/
    char method;                   /**<0-Plucker 1-Havel, 2-Badouel, 3-branchless Badouel*/
    int implicit;                  /**<1 for edge- or node-based implicit MMC, 2 for face-based implicit MMC*/
    char basisorder;               /**<0 to use piece-wise-constant basis for fluence, 1, linear*/
    char outputtype;               /**<'X' output is flux, 'F' output is fluence, 'E' energy deposit*/
    char outputformat;             /**<'ascii' output is text, 'bin': binary, 'json': regular json, 'ubjson': universal binary json*/
    int  mcmethod;                 /**<0 use MCX-styled MC (micro-Beer-Lambert law), 1 use MCML-styled MC (Albedo-Weight)*/
    float roulettesize;            /**<number of roulette for termination*/
    float minenergy;               /**<minimum energy to propagate photon*/
    float omega;                  /**<modulation angular frequency (rad/s) for RF forward simulation; 0 for CW*/
    float nout;                    /**<refractive index for the domain outside the mesh*/
    int isextdet;                  /**<if 1, there is external wide-field detector (marked by -2 in the mesh)*/
    FILE* flog;                    /**<stream handle to print log information*/
    char rootpath[MAX_PATH_LENGTH];/**<a string to specify the root folder of the simulation*/
    unsigned int debuglevel;       /**<a flag to control the printing of the debug information*/
    int debugphoton;               /**<if negative, print debug info for all photons, otherwise, only print for the selected one*/
    float unitinmm;                /**<define the length unit in mm*/
    history his;                   /**<header info of the history file*/
    unsigned int checkpt[MAX_CHECKPOINT]; /**<a list of photon numbers at which a snapshot of the weights will be saved*/
    void* photonseed;              /**< pointer to the seeds of the replayed photon */
    float* replayweight;           /**< pointer to the detected photon weight array */
    float* replaytime;             /**< pointer to the detected photon time-of-fly array */
    char seedfile[MAX_PATH_LENGTH];/**<if the seed is specified as a file (mch), mcx will replay the photons*/
    char deviceid[MAX_DEVICE];
    float workload[MAX_DEVICE];
    char compileropt[MAX_PATH_LENGTH];
    char kernelfile[MAX_PATH_LENGTH];
    char* clsource;
    int parentid;
    int optlevel;
    unsigned int maxdetphoton;     /*anticipated maximum detected photons*/
    unsigned int maxjumpdebug;     /**<num of  photon scattering events to save when saving photon trajectory is enabled*/
    unsigned int debugdatalen;     /**<max number of photon trajectory position length*/
    double* exportfield;           /*memory buffer when returning the flux to external programs such as matlab*/
    unsigned char* exportseed;     /*memory buffer when returning the RNG seed to matlab*/
    float* exportdetected;         /*memory buffer when returning the partial length info to external programs such as matlab*/
    float* exportdebugdata;        /**<pointer to the buffer where the photon trajectory data are stored*/
    float* exportadjoint;         /**<float buffer for RF imaginary fluence (forward mode only; preserved through adjoint post-processing)*/
    float* exportjacob;           /**<float buffer for adjoint Jacobian output (separate from RF imaginary fluence)*/
    ExtraSrc* srcdata;            /**<multi-source list for adjoint mode; length = extrasrclen*/
    int extrasrclen;              /**<number of entries in srcdata (>0 for adjoint/multi-source mode)*/
    int srcid;                    /**<multi-source selector: 0=merge all (default), -1=run each slot separately,
                                       -2=append detectors as adjoint sources without forcing adjoint output,
                                       N>0=run only the N-th source*/
    char adjointmode;             /**<mesh-mode adjoint Jacobian formula: 0=full FEM (default, rb_femjacobian
                                       formula with -0.1*Ve diag + 0.5*off-diag), 1=nodal approximation
                                       (J_mua_n = -nvol*phi_s*phi_d, rbjacmuafast.m). Affects J_mua only;
                                       J_D always uses the full FEM form.*/
    float* nodemua;               /**<optional per-node mua array, length mesh.nn; when set together with
                                       isnodalmua=1 the kernel reads mua per element from the centroid of
                                       these node values instead of gmed[type]. Used by redbird/DOT recon.*/
    float* nodemusp;              /**<optional per-node musp array, length mesh.nn; required for RF/D-coeff
                                       reconstruction along with nodemua (isnodalmusp=1). g is fixed to the
                                       per-tissue gmed[type] value; n is also taken from gmed[type].*/
    char isnodalmua;              /**<1: enable per-node mua via cfg.nodemua (CW or RF DOT recon).*/
    char isnodalmusp;             /**<1: also enable per-node musp via cfg.nodemusp (RF DOT recon).
                                       Requires isnodalmua=1.*/
    double* energytot;             /**<total energy launched for each source, a buffer of length srcnum */
    double* energyesc;             /**<total energy escaped for each source, a buffer of length srcnum */
    unsigned int detectedcount;    /**<total number of detected photons*/
    unsigned int runtime;          /**<total simulation runtime in ms*/
    char autopilot;                /**<1 optimal setting for dedicated card, 2, for non dedicated card*/
    float normalizer;              /**<normalization factor*/
    unsigned int gpuid;            /**<positive integer denotes the 1st/2nd/... OpenCL or CUDA devices, 0xFFFFFFFF for CPU only*/
    int compute;                   /**<0: sse, 1: opencl or 2: cuda*/
    char isdumpjson;               /**<1 to save json */
    int  zipid;                    /**<data zip method "zlib","gzip","base64","lzip","lzma","lz4","lz4hc"*/
    unsigned int savedetflag;      /**<a flag to control the output fields of detected photon data*/
    uint mediabyte;                /**<not used*/
    char* shapedata;               /**<a pointer points to a string defining the JSON-formatted shape data*/
    char* extrajson;               /**<a pointer points to a string defining the extra JSON input data*/
    FLOAT3* node;                  /**<node x/y/z data loaded from json input*/
    unsigned int nodenum;          /**<total number of nodes*/
    int* elem;                     /**<tetrahedron node indices from json input*/
    unsigned int elemnum;          /**<total number of elem*/
    unsigned int elemlen;          /**<number of nodes per elem*/
    int* surf;                     /**<surface triangles from json input (Mesh.MeshSurf), 5 columns: 3 node
                                        indices (1-based), inner label, outer label; OptiX only*/
    unsigned int surfnum;          /**<total number of surface triangles*/
    enum TROIType roitype;         /**<mesh roi type: rtNone: unknown, rtEdge: edge-roi, rtNode: node-roi, rtFace: face-roi */
    float* roidata;                /**<mesh roi (edge/node/face) data*/
    size_t roidatalen;             /**<number of floats in roidata*/
    float* roiprim;                /**<edge/node iMMC primitives given by coordinates, rows {Ax,Ay,Az,Bx,By,Bz,r}: capsule A-B, sphere if B==A*/
    int nroiprim;                  /**<number of rows in roiprim*/
    char jsonfile[MAX_PATH_LENGTH];/**<if the seed is specified as a file (mch), mcx will replay the photons*/
} mcconfig;

#ifdef  __cplusplus
extern "C" {
#endif
void mcx_savedata(OutputType* dat, size_t len, mcconfig* cfg, int isref);
void mcx_savenii(OutputType* dat, size_t len, char* name, int type32bit, int outputformatid, mcconfig* cfg);
void mcx_error(const int id, const char* msg, const char* file, const int linenum);
void mcx_loadconfig(FILE* in, mcconfig* cfg);
void mcx_saveconfig(FILE* in, mcconfig* cfg);
void mcx_readconfig(char* fname, mcconfig* cfg);
void mcx_writeconfig(char* fname, mcconfig* cfg);
void mcx_initcfg(mcconfig* cfg);
void mcx_clearcfg(mcconfig* cfg);
void mcx_validatecfg(mcconfig* cfg);
void mmc_gpu_method(mcconfig* cfg);
void mcx_parsecmd(int argc, char* argv[], mcconfig* cfg);
void mcx_usage(char* exename, mcconfig* cfg);
void mcx_loadvolume(char* filename, mcconfig* cfg, int isbuf);
void mcx_normalize(float field[], float scale, int fieldlen);
int  mcx_readarg(int argc, char* argv[], int id, void* output, const char* type);
void mcx_printlog(mcconfig* cfg, char* str);
int  mcx_remap(char* opt);
int  mcx_lookupindex(char* key, const char* index);
int  mcx_keylookup(char* key, const char* table[]);
int  mcx_parsedebugopt(char* debugopt, const char* debugflag);
void mcx_progressbar(float percent);
int  mcx_loadjson(cJSON* root, mcconfig* cfg);
void mcx_version(mcconfig* cfg);
int  mcx_loadfromjson(char* jbuf, mcconfig* cfg);
void mcx_prep(mcconfig* cfg);
unsigned int mcx_jacchannels(mcconfig* cfg);
unsigned int* mcx_getsdmap(mcconfig* cfg);
void mcx_printheader(mcconfig* cfg);
void mcx_cleargpuinfo(GPUInfo** gpuinfo);
void mcx_convertcol2row(unsigned int** vol, uint3* dim);
void mcx_convertcol2row4d(unsigned int** vol, uint4* dim);
void mcx_savejdata(char* filename, mcconfig* cfg);
void mcx_savedetphoton(float* ppath, void* seeds, int count, int seedbyte, mcconfig* cfg);
int  mcx_jdataencode(void* vol,  int ndim, uint* dims, char* type, int byte, int zipid, void* obj, int isubj, int iscol, mcconfig* cfg);
int  mcx_jdatadecode(void** vol, int* ndim, uint* dims, int maxdim, char** type, cJSON* obj, mcconfig* cfg);
void mcx_convertrow2col(float* vol, uint3* dim);
void mcx_savejnii(void* vol, int ndim, uint* dims, float* voxelsize, char* name, int isfloat, int iscol, int elemsize, mcconfig* cfg);
void mcx_savebnii(void* vol, int ndim, uint* dims, float* voxelsize, char* name, int isfloat, int iscol, int elemsize, mcconfig* cfg);
void mcx_savefloatjnii(float* vol, int ndim, uint* dims, float* voxelsize, char* name, mcconfig* cfg);
void mcx_savefloatbnii(float* vol, int ndim, uint* dims, float* voxelsize, char* name, mcconfig* cfg);
void mcx_savejdet(float* ppath, void* seeds, uint count, int doappend, mcconfig* cfg);
void mcx_savejdb(float* ppath, void* seeds, uint count, int doappend, mcconfig* cfg);
void mcx_fflush(FILE* out);
void mmc_validate_config(mcconfig* cfg, float* detps, int dimdetps[2], int seedbyte);

#if defined(MCX_CONTAINER) && (defined(MATLAB_MEX_FILE) || defined(OCTAVE_API_VERSION_NUMBER))
#ifdef _OPENMP
#define MMC_FPRINTF(fp,...) {if(omp_get_thread_num()==0) {(fp==stderr) ? mexPrintf(__VA_ARGS__) : fprintf(fp,__VA_ARGS__);}}  /**< macro to print messages, calls mexPrint if inside MATLAB */
#else
#define MMC_FPRINTF(fp,...) {(fp==stderr) ? mexPrintf(__VA_ARGS__) : fprintf(fp,__VA_ARGS__);} /**< macro to print messages, calls mexPrint in MATLAB */
#endif
#elif defined(MCX_CONTAINER) && defined(PYBIND11_VERSION_MAJOR)
/* in the Python modules, the log of cfg.flog=2 (stderr, the default) goes to Python's sys.stdout, the
   role mexPrintf plays in MATLAB, so that it shows in Jupyter; only the thread that holds the GIL prints */
#ifdef _OPENMP
#define MMC_FPRINTF(fp,...) {if(omp_get_thread_num()==0) {(fp==stderr) ? mcx_python_printf(__VA_ARGS__) : fprintf(fp,__VA_ARGS__);}}
#else
#define MMC_FPRINTF(fp,...) {(fp==stderr) ? mcx_python_printf(__VA_ARGS__) : fprintf(fp,__VA_ARGS__);}
#endif
#else
#define MMC_FPRINTF(fp,...) fprintf(fp,__VA_ARGS__) /**< macro to print messages, calls fprintf in command line mode */
#endif

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_API_VERSION_NUMBER) || defined (HAVE_OCTAVE)
int mexPrintf(const char* format, ... );
#else
int mexPrintf(const char* format, ... );
#endif

#define MMCDEBUG(cfg,debugflag,outputstr)  {\
        if((cfg)->debuglevel & (debugflag)) {\
            MMC_FPRINTF outputstr ;\
            fflush((cfg)->flog);\
        }\
    }


#ifdef MCX_CONTAINER
#ifdef __cplusplus
extern "C"
#endif
int  mmc_throw_exception(const int id, const char* msg, const char* filename, const int linenum);
void mcx_matlab_flush(void);
void mcx_python_flush(void);
int  mcx_python_printf(const char* format, ...);
#endif

/**
 * An exception can not leave an OpenMP parallel region (std::terminate aborts MATLAB, Octave or
 * Python), so an error raised in one, by a photon or a GPU thread, is recorded and the exception
 * caught in that thread; the first one recorded is raised again once the region has ended.
 * The MATLAB/Python modules catch the exceptions of the photons of mmc_run_mp through
 * mmc_task_guard, which runs task(arg) in a try block and returns non-zero if it threw.
 */
extern int (*mmc_task_guard)(void (*task)(void*), void* arg);
void mmc_defer_error(const int id, const char* msg, const char* file, const int linenum);
int  mmc_error_pending(void);
void mmc_error_reset(void);
void mmc_raise_deferred(void);


#ifdef  __cplusplus
}
#endif

#include <stdlib.h>

#include <stdint.h>

/**
 * @brief The host random number generator that derives the per-thread seeds of the RNGs
 *
 * The xorshift128+ generator of the GPU kernels, rather than the C library rand(): rand() is
 * a different generator in glibc, macOS and MSVC, so that the same seed (-E) simulated other
 * photons on each OS, and its state is shared by all threads. The 128-bit state is filled
 * from the seed by splitmix64, so that small or consecutive seeds give unrelated states.
 */
typedef struct {
    uint64_t s[2];
} MMCHostRNG;

static inline void mmc_hostrng_seed(MMCHostRNG* rng, uint64_t seed) {
    int i;

    for (i = 0; i < 2; i++) {
        uint64_t z = (seed += ((uint64_t)0x9E3779B9U << 32 | 0x7F4A7C15U));
        z = (z ^ (z >> 30)) * ((uint64_t)0xBF58476DU << 32 | 0x1CE4E5B9U);
        z = (z ^ (z >> 27)) * ((uint64_t)0x94D049BBU << 32 | 0x133111EBU);
        rng->s[i] = z ^ (z >> 31);
    }
}

/**
 * @brief Return 32 random bits, the upper half of the next xorshift128+ output
 */
static inline unsigned int mmc_rand_uint32(MMCHostRNG* rng) {
    uint64_t s1 = rng->s[0];
    const uint64_t s0 = rng->s[1];

    rng->s[0] = s0;
    s1 ^= s1 << 23;
    rng->s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return (unsigned int)((rng->s[1] + s0) >> 32);
}

#endif
