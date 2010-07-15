/*******************************************************************************
**  Mesh-based Monte Carlo (MMC) -- this unit was ported from MCX
**
**  Monte Carlo eXtreme (MCX)  - GPU accelerated Monte Carlo 3D photon migration
**  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates," Biomed. Opt. Express, (in press)
**
**  (Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, vol. 17, issue 22, pp. 20178-20190 (2009)
**
**  mcx_utils.h: configuration and command line option processing unit
**
**  License: GPL v3, see LICENSE.txt for details
**
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
                  dlProgress=2048};

typedef struct MCXMedium{
	float mua;
	float mus;
	float n;
	float g;
} Medium;  /*this order shall match prop.{xyzw} in mcx_main_loop*/

typedef struct MCXConfig{
	int nphoton;      /*(total simulated photon number) we now use this to 
	                     temporarily alias totalmove, as to specify photon
			     number is causing some troubles*/
	//int totalmove;   /* [depreciated] total move per photon*/
        int nblocksize;   /*thread block size*/
	int nthread;      /*num of total threads, multiple of 128*/
	int seed;         /*random number generator seed*/
	
	float3 srcpos;    /*src position in mm*/
	float3 srcdir;    /*src normal direction*/
	float tstart;     /*start time in second*/
	float tstep;      /*time step in second*/
	float tend;       /*end time in second*/
	float3 steps;     /*voxel sizes along x/y/z in mm*/
	
	uint3 dim;        /*domain size*/
	uint3 crop0;      /*sub-volume for cache*/
	uint3 crop1;      /*the other end of the caching box*/
	int medianum;     /*total types of media*/
	int detnum;       /*total detector numbers*/
	float detradius;  /*detector radius*/
        float sradius;    /*source region radius: if set to non-zero, accumulation 
                            will not perform for dist<sradius; this can reduce
                            normalization error when using non-atomic write*/

	Medium *prop;     /*optical property mapping table*/
	float4 *detpos;   /*detector positions and radius, overwrite detradius*/
	float  minstep;   /*accumulation step size*/

	int maxgate;        /*simultaneous recording gates*/
	int respin;         /*number of repeatitions*/
	int printnum;       /*number of printed threads (for debugging)*/

	unsigned char *vol; /*pointer to the volume*/
	char session[MAX_SESSION_LENGTH]; /*session id, a string*/
	char isrowmajor;    /*1 for C-styled array in vol, 0 for matlab-styled array*/
	char isreflect;     /*1 for reflecting photons at boundary,0 for exiting*/
        char isref3;        /*1 considering maximum 3 ref. interfaces; 0 max 2 ref*/
	char isnormalized;  /*1 to normalize the fluence, 0 for raw fluence*/
	char issavedet;     /*1 to count all photons hits the detectors*/
	char issave2pt;     /*1 to save the 2-point distribution, 0 do not save*/
	char isgpuinfo;     /*1 to print gpu info when attach, 0 do not print*/
	float roulettesize; /*number of roulette for termination*/
        float minenergy;    /*minimum energy to propagate photon*/
	float nout;         /*refractive index for the domain outside the mesh*/
        FILE *flog;         /*stream handle to print log information*/
        char rootpath[MAX_PATH_LENGTH];
        unsigned int debuglevel;
	float unitinmm;     /*define the length unit in mm*/
} Config;

#ifdef	__cplusplus
extern "C" {
#endif
void mcx_savedata(float *dat,int len,Config *cfg);
void mcx_error(int id,char *msg);
void mcx_loadconfig(FILE *in, Config *cfg);
void mcx_saveconfig(FILE *in, Config *cfg);
void mcx_readconfig(char *fname, Config *cfg);
void mcx_writeconfig(char *fname, Config *cfg);
void mcx_initcfg(Config *cfg);
void mcx_clearcfg(Config *cfg);
void mcx_parsecmd(int argc, char* argv[], Config *cfg);
void mcx_usage(char *exename);
void mcx_loadvolume(char *filename,Config *cfg);
void mcx_normalize(float field[], float scale, int fieldlen);
int  mcx_readarg(int argc, char *argv[], int id, void *output,char *type);
void mcx_printlog(Config *cfg, char *str);
int  mcx_remap(char *opt);
int  mcx_parsedebugopt(char *debugopt);
void mcx_progressbar(int n, int ntotal, Config *cfg);

#ifdef	__cplusplus
}
#endif

#endif
