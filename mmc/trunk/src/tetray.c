/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    tetray.c

\brief   << Main program of MMC >>
*******************************************************************************/

#include <stdlib.h>
#include <time.h>
#include <string.h>
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
	mcconfig cfg;
	tetmesh mesh;
	raytracer tracer;
	visitor master={0.f,0.f,0,0,NULL};
	double Eabsorb=0.0;
	RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
        RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
	unsigned int i;
	float raytri=0.f;
	unsigned int threadid=0,ncomplete=0,t0,dt;

	t0=StartTimer();

        /** \subsection sinit Initialization */
        mcx_initcfg(&cfg);

        /** parse command line options to initialize the configurations */
        mcx_parsecmd(argc,argv,&cfg);
	
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"initializing ... "));

	mesh_init_from_cfg(&mesh,&cfg);
	tracer_init(&tracer,&mesh,cfg.method);
	tracer_prep(&tracer,&cfg);

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
	visitor visit={0.f,1.f/cfg.tstep,DET_PHOTON_BUF,0,NULL};
	int buflen=(1+((cfg.ismomentum)>0))*mesh.prop+2;
	if(cfg.issavedet) 
	    visit.partialpath=(float*)calloc(visit.detcount*buflen,sizeof(float));
#ifdef _OPENMP
	threadid=omp_get_thread_num();	
#endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg.seed),threadid);

	/*launch photons*/
        #pragma omp for reduction(+:Eabsorb) reduction(+:raytri)
	for(i=0;i<cfg.nphoton;i++){
		visit.raytet=0.f;
		Eabsorb+=onephoton(i,&tracer,&mesh,&cfg,ran0,ran1,&visit);
		raytri+=visit.raytet;
		#pragma omp atomic
		   ncomplete++;

		if((cfg.debuglevel & dlProgress) && threadid==0)
			mcx_progressbar(ncomplete,&cfg);
	}
	if(cfg.issavedet){
	    #pragma omp atomic
		master.detcount+=visit.bufpos;
            #pragma omp barrier
	    if(threadid==0)
		master.partialpath=(float*)calloc(master.detcount*buflen,sizeof(float));
            #pragma omp barrier
            #pragma omp critical
            {
		memcpy(master.partialpath+master.bufpos*buflen,
		       visit.partialpath,visit.bufpos*buflen*sizeof(float));
		master.bufpos+=visit.bufpos;
            }
            #pragma omp barrier
	    free(visit.partialpath);
	}
}

        /** \subsection sreport Post simulation */

	if((cfg.debuglevel & dlProgress))
		mcx_progressbar(cfg.nphoton,&cfg);

	dt=GetTimeMillis()-t0;
	MMCDEBUG(&cfg,dlProgress,(cfg.flog,"\n"));
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",dt));
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"speed ...\t%.0f ray-tetrahedron tests\n",raytri));

	tracer_clear(&tracer);

	if(cfg.isnormalized){
          fprintf(cfg.flog,"total simulated energy: %d\tabsorbed: %5.5f%%\tnormalizor=%g\n",
		cfg.nphoton,100.f*Eabsorb/cfg.nphoton,mesh_normalize(&mesh,&cfg,Eabsorb,cfg.nphoton));
	}
	if(cfg.issave2pt){
		switch(cfg.outputtype){
		    case otFlux:   MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving flux ...")); break;
                    case otFluence:MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving fluence ...")); break;
                    case otEnergy: MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving energy deposit ...")); break;
		}
		mesh_saveweight(&mesh,&cfg);
	}
	if(cfg.issavedet){
		MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving detected photons ..."));
		mesh_savedetphoton(master.partialpath,master.bufpos,&cfg);
		free(master.partialpath);
	}
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

        /** \subsection sclean End the simulation */

	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
