/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Pluker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmc_host.c

\brief   << Driver program of MMC >>
*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "mmc_host.h"
#include "tictoc.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the 
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

int mmc_init_from_cmd(mcconfig *cfg, tetmesh *mesh, raytracer *tracer,int argc, char**argv){
        mcx_initcfg(cfg);
        mcx_parsecmd(argc,argv,cfg);
	MMCDEBUG(cfg,dlTime,(cfg->flog,"initializing from commands ... "));
	mesh_init_from_cfg(mesh,cfg);
        return 0;
}
int mmc_init_from_json(mcconfig *cfg, tetmesh *mesh, raytracer *tracer, char *jcfg, char *jmesh){
        mcx_initcfg(cfg);
	MMCDEBUG(cfg,dlTime,(cfg->flog,"initializing from JSON ... "));
	mcx_loadfromjson(jcfg,cfg);
	mesh_init_from_cfg(mesh,cfg); /*need to define mesh load from json*/
	tracer_init(tracer,mesh,cfg->method);
        return 0;
}
int mmc_reset(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	mmc_cleanup(cfg,mesh,tracer);
	mcx_initcfg(cfg);
	mesh_init(mesh);
        return 0;
}
int mmc_cleanup(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
	tracer_clear(tracer);
	mesh_clear(mesh);
        mcx_clearcfg(cfg);
        return 0;
}
int mmc_prep(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
        mcx_prep(cfg);
	tracer_init(tracer,mesh,cfg->method);
	tracer_prep(tracer,cfg);
        return 0;
}

int mmc_run_mp(mcconfig *cfg, tetmesh *mesh, raytracer *tracer){
        double Eabsorb=0.0;
	RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
        RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
        unsigned int i;
        float raytri=0.f,raytri0=0.f;
        unsigned int threadid=0,ncomplete=0,t0,dt, debuglevel;
	visitor master={0.f,0.f,0.f,0,0,0,NULL,NULL,0.f,0.f};

	t0=StartTimer();
	
#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
	cfg->issaveseed=0;
#endif
        dt=GetTimeMillis();
        MMCDEBUG(cfg,dlTime,(cfg->flog,"seed=%u\nsimulating ... ",cfg->seed));
        if(cfg->debugphoton>=0){
            debuglevel=cfg->debuglevel;
            cfg->debuglevel &= 0xFFFFEA00;
#ifdef _OPENMP
            omp_set_num_threads(1);
#endif
        }

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
	visitor visit={0.f,0.f,1.f/cfg->tstep,DET_PHOTON_BUF,0,0,NULL,NULL,0.f,0.f};
	visit.reclen=(2+((cfg->ismomentum)>0))*mesh->prop+(cfg->issaveexit>0)*6+2;
	if(cfg->issavedet){
	    if(cfg->issaveseed)
	        visit.photonseed=calloc(visit.detcount,(sizeof(RandType)*RAND_BUF_LEN));
	    visit.partialpath=(float*)calloc(visit.detcount*visit.reclen,sizeof(float));
	}
#ifdef _OPENMP
	threadid=omp_get_thread_num();
#endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg->seed),threadid);

	/*launch photons*/
        #pragma omp for reduction(+:Eabsorb) reduction(+:raytri,raytri0)
	for(i=0;i<cfg->nphoton;i++){
		visit.raytet=0.f;
		visit.raytet0=0.f;
                if(i==cfg->debugphoton)
                    cfg->debuglevel = debuglevel;

		if(cfg->seed==SEED_FROM_FILE)
		    Eabsorb+=onephoton(i,tracer,mesh,cfg,((RandType *)cfg->photonseed)+i*RAND_BUF_LEN,ran1,&visit);
		else
		    Eabsorb+=onephoton(i,tracer,mesh,cfg,ran0,ran1,&visit);
		raytri+=visit.raytet;
		raytri0+=visit.raytet0;

                if(i==cfg->debugphoton)
                    cfg->debuglevel &= 0xFFFFEA00;

		#pragma omp atomic
		   ncomplete++;

		if((cfg->debuglevel & dlProgress) && threadid==0)
			mcx_progressbar(ncomplete,cfg);
		if(cfg->issave2pt && cfg->checkpt[0])
			mesh_saveweightat(mesh,cfg,i+1);
	}

	#pragma omp atomic
		master.totalweight += visit.totalweight;

	if(cfg->issavedet){
	    #pragma omp atomic
		master.detcount+=visit.bufpos;
            #pragma omp barrier
	    if(threadid==0){
		master.partialpath=(float*)calloc(master.detcount*visit.reclen,sizeof(float));
	        if(cfg->issaveseed)
        	    master.photonseed=calloc(master.detcount,(sizeof(RandType)*RAND_BUF_LEN));
            }
            #pragma omp barrier
            #pragma omp critical
            {
		memcpy(master.partialpath+master.bufpos*visit.reclen,
		       visit.partialpath,visit.bufpos*visit.reclen*sizeof(float));
                if(cfg->issaveseed)
                    memcpy((unsigned char*)master.photonseed+master.bufpos*(sizeof(RandType)*RAND_BUF_LEN),
                            visit.photonseed,visit.bufpos*(sizeof(RandType)*RAND_BUF_LEN));
		master.bufpos+=visit.bufpos;
            }
            #pragma omp barrier
	    free(visit.partialpath);
            if(cfg->issaveseed && visit.photonseed)
                 free(visit.photonseed);
	}
}

        /** \subsection sreport Post simulation */

	if((cfg->debuglevel & dlProgress))
		mcx_progressbar(cfg->nphoton,cfg);

	dt=GetTimeMillis()-dt;
	MMCDEBUG(cfg,dlProgress,(cfg->flog,"\n"));
        MMCDEBUG(cfg,dlTime,(cfg->flog,"\tdone\t%d\n",dt));
        MMCDEBUG(cfg,dlTime,(cfg->flog,"speed ...\t%.2f photon/ms, %.0f ray-tetrahedron tests (%.0f were overhead)\n",(double)cfg->nphoton/dt,raytri,raytri0));
        if(cfg->issavedet)
           fprintf(cfg->flog,"detected %d photons\n",master.detcount);

	if(cfg->isnormalized){
          cfg->his.normalizer=mesh_normalize(mesh,cfg,Eabsorb,master.totalweight);
          fprintf(cfg->flog,"total simulated energy: %f\tabsorbed: %5.5f%%\tnormalizor=%g\n",
		master.totalweight,100.f*Eabsorb/master.totalweight,cfg->his.normalizer);
	}
	if(cfg->issave2pt){
		switch(cfg->outputtype){
		    case otFlux:   MMCDEBUG(cfg,dlTime,(cfg->flog,"saving flux ...")); break;
                    case otFluence:MMCDEBUG(cfg,dlTime,(cfg->flog,"saving fluence ...")); break;
                    case otEnergy: MMCDEBUG(cfg,dlTime,(cfg->flog,"saving energy deposit ...")); break;
		}
		mesh_saveweight(mesh,cfg);
	}
	if(cfg->issavedet){
		MMCDEBUG(cfg,dlTime,(cfg->flog,"saving detected photons ..."));
		mesh_savedetphoton(master.partialpath,master.photonseed,master.bufpos,(sizeof(RandType)*RAND_BUF_LEN),cfg);
		free(master.partialpath);
                if(cfg->issaveseed && master.photonseed)
                    free(master.photonseed);
	}
        MMCDEBUG(cfg,dlTime,(cfg->flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	return 0;
}

