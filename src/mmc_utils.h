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
    #define S_BOLD     "\x1b[1m"
    #define S_ITALIC   "\x1b[3m"
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


typedef double OutputType;

enum TDebugLevel {dlMove = 1, dlTracing = 2, dlBary = 4, dlWeight = 8, dlDist = 16, dlTracingEnter = 32,
                  dlTracingExit = 64, dlEdge = 128, dlAccum = 256, dlTime = 512, dlReflect = 1024,
                  dlProgress = 2048, dlExit = 4096
                 };

enum TRTMethod {rtPlucker, rtHavel, rtBadouel, rtBLBadouel, rtBLBadouelGrid};
enum TMCMethod {mmMCX, mmMCML};
enum TComputeBackend {cbSSE, cbOpenCL, cbCUDA, cbOptiX};

enum TSrcType {stPencil, stIsotropic, stCone, stGaussian, stPlanar,
               stPattern, stFourier, stArcSin, stDisk, stFourierX,
               stFourier2D, stZGaussian, stLine, stSlit
              };
enum TOutputType {otFlux, otFluence, otEnergy, otJacobian, otWL, otWP};
enum TOutputFormat {ofASCII, ofBin, ofNifti, ofAnalyze, ofMC2, ofTX3, ofJNifti, ofBJNifti};
enum TOutputDomain {odMesh, odGrid};
enum TDeviceVendor {dvUnknown, dvNVIDIA, dvAMD, dvIntel, dvIntelGPU};
enum TMCXParent  {mpStandalone, mpMATLAB};
enum TBoundary {bcNoReflect, bcReflect, bcAbsorbExterior, bcMirror /*, bcCylic*/};
enum TRayHitType {htNone, htInOut, htOutIn, htNoHitIn, htNoHitOut};
enum TROIType {rtNone, rtEdge, rtNode, rtFace};

enum TBJData {JDB_mixed, JDB_nulltype, JDB_noop, JDB_true, JDB_false,
              JDB_char, JDB_string, JDB_hp, JDB_int8, JDB_uint8, JDB_int16, JDB_int32,
              JDB_int64, JDB_single, JDB_double, JDB_array, JDB_object, JDB_numtypes,
              JDB_uint16 = 10, JDB_uint32, JDB_uint64
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
    int maxmpthread;  /**< maximum thread number per multi-processor */
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
    MMCfloat3 srcpos;              /**<src position in mm*/
    float4 srcdir;                 /**<src normal direction*/
    int srctype;                   /**<src type: 0 - pencil beam, 1 - isotropic ... */
    float4 srcparam1;          /**<source parameters set 1*/
    float4 srcparam2;          /**<source parameters set 2*/
    float* srcpattern;         /**<source pattern*/
    int voidtime;
    float4 bary0;                  /**<initial bary centric coordinates of the source*/
    float tstart;                  /**<start time in second*/
    float tstep;                   /**<time step in second*/
    float tend;                    /**<end time in second*/
    MMCfloat3 steps;               /**<voxel sizes along x/y/z in mm*/
    uint3 dim;                     /**<dim.x is the initial element number in MMC, dim.y is faceid*/
    uint4 crop0;                   /**<sub-volume for cache*/
    uint4 crop1;                   /**<the other end of the caching box*/
    int medianum;                  /**<total types of media*/
    int srcnum;            /**<total number of sources, could be larger than 1 only with pattern illumination*/
    int detnum;                    /**<total detector numbers*/
    float detradius;               /**<detector radius*/
    float sradius;                 /**<source region radius: if set to non-zero, accumulation \
                                           will not perform for dist<sradius; this can reduce \
                                           normalization error when using non-atomic write*/
    medium* prop;                  /**<optical property mapping table*/
    float4* detpos;                /**<detector positions and radius, overwrite detradius*/
    float4 detparam1;          /**<parameters set 1 for wide-field detector*/
    float4 detparam2;          /**<parameters set 2 for wide-feild detector*/
    float* detpattern;         /**<detector pattern*/
    float  minstep;                /**<accumulation step size*/
    int maxgate;                   /**<simultaneous recording gates*/
    int respin;                    /**<number of repeatitions*/
    int printnum;                  /**<number of printed threads (for debugging)*/
    int replaydet;                 /**<the detector id for which to replay the detected photons, start from 1
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
    char issaveexit;               /**<1 to save the exit position and vector of a detected photon, implies issavedet=1*/
    /**<2 to save accumulated photon weight in frames of images*/
    char issave2pt;                /**<1 to save the 2-point distribution, 0 do not save*/
    char isgpuinfo;                /**<1 to print gpu info when attach, 0 do not print*/
    char isspecular;               /**<1 calculate the initial specular ref if outside the mesh, 0 do not calculate*/
    char issaveseed;               /**<1 save the seed for a detected photon, 0 do not save*/
    char issaveref;                /**<1 to save diffuse reflectance on surface, 0 no save*/
    char isatomic;                 /**<1 use atomic operations for weight accumulation, 0 do not use*/
    char method;                   /**<0-Plucker 1-Havel, 2-Badouel, 3-branchless Badouel*/
    int implicit;              /**<1 for edge- or node-based implicit MMC, 2 for face-based implicit MMC*/
    char basisorder;               /**<0 to use piece-wise-constant basis for fluence, 1, linear*/
    char outputtype;               /**<'X' output is flux, 'F' output is fluence, 'E' energy deposit*/
    char outputformat;             /**<'ascii' output is text, 'bin': binary, 'json': regular json, 'ubjson': universal binary json*/
    int  mcmethod;                 /**<0 use MCX-styled MC (micro-Beer-Lambert law), 1 use MCML-styled MC (Albedo-Weight)*/
    float roulettesize;            /**<number of roulette for termination*/
    float minenergy;               /**<minimum energy to propagate photon*/
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
    char kernelfile[MAX_SESSION_LENGTH];
    char* clsource;
    int parentid;
    int optlevel;
    unsigned int maxdetphoton; /*anticipated maximum detected photons*/
    double* exportfield;     /*memory buffer when returning the flux to external programs such as matlab*/
    unsigned char* exportseed;     /*memory buffer when returning the RNG seed to matlab*/
    float* exportdetected;  /*memory buffer when returning the partial length info to external programs such as matlab*/
    double energytot, energyabs, energyesc;
    unsigned int detectedcount; /**<total number of detected photons*/
    unsigned int runtime;
    char autopilot;     /**<1 optimal setting for dedicated card, 2, for non dedicated card*/
    float normalizer;            /**<normalization factor*/
    unsigned int nbuffer;        /**<2^nbuffer is the number of buffers for accummulation*/
    unsigned int gpuid;
    int compute;
    char isdumpjson;             /**<1 to save json */
    int  zipid;                  /**<data zip method "zlib","gzip","base64","lzip","lzma","lz4","lz4hc"*/
    unsigned int savedetflag;    /**<a flag to control the output fields of detected photon data*/
    uint mediabyte;
    char* shapedata;    /**<a pointer points to a string defining the JSON-formatted shape data*/
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
void mcx_parsecmd(int argc, char* argv[], mcconfig* cfg);
void mcx_usage(char* exename, mcconfig* cfg);
void mcx_loadvolume(char* filename, mcconfig* cfg, int isbuf);
void mcx_normalize(float field[], float scale, int fieldlen);
int  mcx_readarg(int argc, char* argv[], int id, void* output, const char* type);
void mcx_printlog(mcconfig* cfg, char* str);
int  mcx_remap(char* opt);
int  mcx_lookupindex(char* key, const char* index);
int  mcx_keylookup(char* key, const char* table[]);
int  mcx_parsedebugopt(char* debugopt);
void mcx_progressbar(float percent, void* cfg);
int  mcx_loadjson(cJSON* root, mcconfig* cfg);
void mcx_version(mcconfig* cfg);
int  mcx_loadfromjson(char* jbuf, mcconfig* cfg);
void mcx_prep(mcconfig* cfg);
void mcx_printheader(mcconfig* cfg);
void mcx_cleargpuinfo(GPUInfo** gpuinfo);
void mcx_convertcol2row(unsigned int** vol, uint3* dim);
void mcx_convertcol2row4d(unsigned int** vol, uint4* dim);
void mcx_savejdata(char* filename, mcconfig* cfg);
int  mcx_jdataencode(void* vol,  int ndim, uint* dims, char* type, int byte, int zipid, void* obj, int isubj, mcconfig* cfg);
int  mcx_jdatadecode(void** vol, int* ndim, uint* dims, int maxdim, char** type, cJSON* obj, mcconfig* cfg);
void mcx_savejnii(OutputType* vol, int ndim, uint* dims, float* voxelsize, char* name, int isfloat, mcconfig* cfg);
void mcx_savebnii(OutputType* vol, int ndim, uint* dims, float* voxelsize, char* name, int isfloat, mcconfig* cfg);
void mcx_savejdet(float* ppath, void* seeds, uint count, int doappend, mcconfig* cfg);

#ifdef MCX_CONTAINER
#ifdef _OPENMP
#define MMC_FPRINTF(fp,...) {if(omp_get_thread_num()==0) mexPrintf(__VA_ARGS__);}
#else
#define MMC_FPRINTF(fp,...) mexPrintf(__VA_ARGS__)
#endif
#else
#define MMC_FPRINTF(fp,...) fprintf(fp,__VA_ARGS__)
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
#endif


#ifdef  __cplusplus
}
#endif

#endif
