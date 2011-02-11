/*******************************************************************************
**  Mesh-based Monte Carlo (MMC) -- this unit was ported from MCX
**
**  Monte Carlo eXtreme (MCX)  - GPU accelerated Monte Carlo 3D photon migration
**  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates," Biomed. Opt. Express, 1(1) 165-175 (2010)
**
**  (Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, 17(22) 20178-20190 (2009)
**
**  License: GPL v3, see LICENSE.txt for details
**
*******************************************************************************/

/***************************************************************************//**
\file    mcx_utils.h

\brief   Definition of program options and problem domain configurations
*******************************************************************************/

#ifndef _MMC_UTILITIES_H
#define _MMC_UTILITIES_H

#include <stdio.h>
#include <vector_types.h>

#define MAX_PROP            256
#define MAX_DETECTORS       256
#define MAX_PATH_LENGTH     1024
#define MAX_SESSION_LENGTH  256
#define MIN(a,b)            ((a)<(b)?(a):(b))

#define MMCDEBUG(cfg,debugflag,outputstr)  {\
				if((cfg)->debuglevel & (debugflag)) {\
					fprintf outputstr ;\
					fflush((cfg)->flog);\
				}\
                            }

enum TDebugLevel {dlMove=1,dlTracing=2,dlBary=4,dlWeight=8,dlDist=16,dlTracingEnter=32,
                  dlTracingExit=64,dlEdge=128,dlAccum=256,dlTime=512,dlReflect=1024,
                  dlProgress=2048,dlExit=4096};

enum TRTMethod {rtPlucker, rtHavel, rtBadouel, rtBLBadouel};
enum TSrcType {stPencil, stCone, stGaussian};


/***************************************************************************//**
\struct MMC_medium mcx_utils.h

\brief  This structure defines the optical properties for each medium

For each optical medium, we use a 4-element structure to identify its optical
properties, i.e, the absorption coeff (mu_a)in 1/mm, the scattering coeff (mu_s) 
in 1/mm, the refractive index (n) and anisotropy (g).
*******************************************************************************/  

typedef struct MMC_medium{
	float mua;        /**<absorption coeff in 1/mm unit*/
	float mus;        /**<scattering coeff in 1/mm unit*/
	float n;          /**<refractive index*/
	float g;          /**<anisotropy*/
} medium;


/***************************************************************************//**
\struct MMC_config mcx_utils.h
                                                                                                                                                                                    
\brief  This structure defines the problem settings (domain, filenames, session)
                                                                                                                                                                                    
*******************************************************************************/  

typedef struct MMC_config{
	int nphoton;      /**<(total simulated photon number) we now use this to 
	                     temporarily alias totalmove, as to specify photon
			     number is causing some troubles*/
        int nblocksize;   /**<thread block size*/
	int nthread;      /**<num of total threads, multiple of 128*/
	int seed;         /**<random number generator seed*/
	
	float3 srcpos;    /**<src position in mm*/
	float3 srcdir;    /**<src normal direction*/
	int srctype;	  /**<src type: 0 - pencil beam, 1 - cone beam */
	float4 srcparam;  /**<additional parameters for advanced sources */
	float4 bary0;     /**<initial bary centric coordinates of the source*/
	float tstart;     /**<start time in second*/
	float tstep;      /**<time step in second*/
	float tend;       /**<end time in second*/
	float3 steps;     /**<voxel sizes along x/y/z in mm*/

	uint3 dim;        /**<domain size*/
	uint3 crop0;      /**<sub-volume for cache*/
	uint3 crop1;      /**<the other end of the caching box*/
	int medianum;     /**<total types of media*/
	int detnum;       /**<total detector numbers*/
	float detradius;  /**<detector radius*/
        float sradius;    /**<source region radius: if set to non-zero, accumulation 
                            will not perform for dist<sradius; this can reduce
                            normalization error when using non-atomic write*/

	medium *prop;     /**<optical property mapping table*/
	float4 *detpos;   /**<detector positions and radius, overwrite detradius*/
	float  minstep;   /**<accumulation step size*/

	int maxgate;        /**<simultaneous recording gates*/
	int respin;         /**<number of repeatitions*/
	int printnum;       /**<number of printed threads (for debugging)*/

	unsigned char *vol; /**<pointer to the volume*/
	char session[MAX_SESSION_LENGTH]; /**<session id, a string*/
        char meshtag[MAX_PATH_LENGTH];   /**<a string to tag all input mesh files*/
	char isrowmajor;    /**<1 for C-styled array in vol, 0 for matlab-styled array*/
	char isreflect;     /**<1 for reflecting photons at boundary,0 for exiting*/
        char isref3;        /**<1 considering maximum 3 ref. interfaces; 0 max 2 ref*/
	char isnormalized;  /**<1 to normalize the fluence, 0 for raw fluence*/
	char issavedet;     /**<1 to count all photons hits the detectors*/
	char issave2pt;     /**<1 to save the 2-point distribution, 0 do not save*/
	char isgpuinfo;     /**<1 to print gpu info when attach, 0 do not print*/
	char method;        /**<0-Plucker 1-Havel, 2-Badouel, 3-branchless Badouel*/
	char basisorder;    /**<0 to use piece-wise-constant basis for fluence, 1, linear*/
	float roulettesize; /**<number of roulette for termination*/
        float minenergy;    /**<minimum energy to propagate photon*/
	float nout;         /**<refractive index for the domain outside the mesh*/
        FILE *flog;         /**<stream handle to print log information*/
        char rootpath[MAX_PATH_LENGTH]; /**<a string to specify the root folder of the simulation*/
        unsigned int debuglevel; /**<a flag to control the printing of the debug information*/
	float unitinmm;     /**<define the length unit in mm*/
} mcconfig;

#ifdef	__cplusplus
extern "C" {
#endif
void mcx_savedata(float *dat,int len,mcconfig *cfg);
void mcx_error(int id,char *msg);
void mcx_loadconfig(FILE *in, mcconfig *cfg);
void mcx_saveconfig(FILE *in, mcconfig *cfg);
void mcx_readconfig(char *fname, mcconfig *cfg);
void mcx_writeconfig(char *fname, mcconfig *cfg);
void mcx_initcfg(mcconfig *cfg);
void mcx_clearcfg(mcconfig *cfg);
void mcx_parsecmd(int argc, char* argv[], mcconfig *cfg);
void mcx_usage(char *exename);
void mcx_loadvolume(char *filename,mcconfig *cfg);
void mcx_normalize(float field[], float scale, int fieldlen);
int  mcx_readarg(int argc, char *argv[], int id, void *output,char *type);
void mcx_printlog(mcconfig *cfg, char *str);
int  mcx_remap(char *opt);
int  mcx_getmethodid(char *method);
int  mcx_getsrcid(char *srctype);
int  mcx_parsedebugopt(char *debugopt);
void mcx_progressbar(unsigned int n, unsigned int ntotal, mcconfig *cfg);

#ifdef	__cplusplus
}
#endif

#endif
