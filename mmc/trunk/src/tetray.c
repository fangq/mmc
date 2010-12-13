/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, 17(22) 20178-20190 (2009).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    tetray.c

\brief   Main program of MMC
*******************************************************************************/

#include <stdlib.h>
#include <time.h>
#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"
#include "tictoc.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the 
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

int main(int argc, char**argv){
	Config cfg;
	tetmesh mesh;
	raytracer tracer;
	float rtstep;
	double Eabsorb=0.0;
	RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
        RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
	int i;
	unsigned int threadid=0,ncomplete=0,t0;

	t0=StartTimer();

        /** \subsection sinit Initialization */
        mcx_initcfg(&cfg);

        /** parse command line options to initialize the configurations */
        mcx_parsecmd(argc,argv,&cfg);
	
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"initizing ... "));

	mesh_init_from_cfg(&mesh,&cfg);
	tracer_init(&tracer,&mesh);
	
	rtstep=1.f/cfg.tstep;

	if(cfg.seed<0) cfg.seed=time(NULL);

        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\nsimulating ... ",GetTimeMillis()-t0));

	/***************************************************************************//**
	The master thread then spawn multiple work-threads depending on your
	OpenMP settings. By default, the total thread number (master + work) is 
	your total CPU core number. For example, if you have a dual-core CPU, 
	the total thread number is 2; if you have two quad-core CPUs, the total 
	thread number is 8. If you want to set the total thread number manually, 
	you need to set the OMP_NUM_THREADS environment variable. For example, 
	\c OMP_NUM_THREADS=3 sets the total thread number to 3.
	*******************************************************************************/

/** \subsection ssimu Parallel photon transport simulation */

#pragma omp parallel private(ran0,ran1,threadid)
{
#ifdef _OPENMP
	threadid=omp_get_thread_num();	
#endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg.seed),threadid);

	/*launch photons*/
#pragma omp for reduction(+:Eabsorb)
	for(i=0;i<cfg.nphoton;i++){
		Eabsorb+=onephoton(i,&tracer,&mesh,&cfg,rtstep,ran0,ran1);
		#pragma omp atomic
		   ncomplete++;

		if((cfg.debuglevel & dlProgress) && threadid==0)
			mcx_progressbar(ncomplete,cfg.nphoton,&cfg);
	}
}

        /** \subsection sreport Post simulation */

	if((cfg.debuglevel & dlProgress))
		mcx_progressbar(cfg.nphoton,cfg.nphoton,&cfg);
	MMCDEBUG(&cfg,dlProgress,(cfg.flog,"\n"));
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	tracer_clear(&tracer);

	if(cfg.isnormalized){
          fprintf(cfg.flog,"total simulated energy: %d\tabsorbed: %5.3f%%\tnormalizor=%f\n",
		cfg.nphoton,100.*Eabsorb/cfg.nphoton,mesh_normalize(&mesh,&cfg,Eabsorb,cfg.nphoton));
	}
	if(cfg.issave2pt){
		MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving data ..."));
		mesh_saveweight(&mesh,&cfg);
	}

        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

        /** \subsection sclean End the simulation */

	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
