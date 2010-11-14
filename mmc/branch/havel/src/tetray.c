/*******************************************************************************
**  Mesh-based Monte Carlo (MMC)
**
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
**  tetray.c: main program of MMC
**
**  License: GPL v3, see LICENSE.txt for details
**
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

#ifdef MMC_LOGISTIC
  #include "logistic_rand.h"
#else
  #include "posix_randr.h"
#endif

int main(int argc, char**argv){
	Config cfg;
	tetmesh mesh;
	tetplucker plucker;
	float rtstep,Eabsorb=0.f;
	RandType ran0[RAND_BUF_LEN],ran1[RAND_BUF_LEN];
	int i,threadid=0,ncomplete=0;
	unsigned int t0;

	t0=StartTimer();

        mcx_initcfg(&cfg);

        /* parse command line options to initialize the configurations */
        mcx_parsecmd(argc,argv,&cfg);
	
	MMCDEBUG(&cfg,dlTime,(cfg.flog,"initizing ... "));

	mesh_init(&mesh);
	mesh_loadnode(&mesh,&cfg);
	mesh_loadelem(&mesh,&cfg);
	mesh_loadfaceneighbor(&mesh,&cfg);
	mesh_loadmedia(&mesh,&cfg);
	mesh_loadelemvol(&mesh,&cfg);

	plucker_init(&plucker,&mesh,cfg.isplucker);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL||mesh.med==NULL)
		mesh_error("encountered error while loading mesh files");

	rtstep=1.f/cfg.tstep;

	if(cfg.seed<0) cfg.seed=time(NULL);

        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\nsimulating ... ",GetTimeMillis()-t0));

#pragma omp parallel private(ran0,ran1,threadid)
{
#ifdef _OPENMP
	threadid=omp_get_thread_num();	
#endif
	rng_init(ran0,ran1,(unsigned int *)&(cfg.seed),threadid);

	/*launch photons*/
#pragma omp for reduction(+:Eabsorb)
	for(i=0;i<cfg.nphoton;i++){
		Eabsorb+=onephoton(i,&plucker,&mesh,&cfg,rtstep,ran0,ran1);
		#pragma omp atomic
		   ncomplete++;

		if((cfg.debuglevel & dlProgress) && threadid==0)
			mcx_progressbar(ncomplete,cfg.nphoton,&cfg);
	}
}
	if((cfg.debuglevel & dlProgress))
		mcx_progressbar(cfg.nphoton,cfg.nphoton,&cfg);
	MMCDEBUG(&cfg,dlProgress,(cfg.flog,"\n"));
        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	plucker_clear(&plucker);

	if(cfg.isnormalized){
          fprintf(cfg.flog,"total simulated energy: %d\tabsorbed: %5.3f%%\tnormalizor=%f\n",
		cfg.nphoton,100.f*Eabsorb/cfg.nphoton,mesh_normalize(&mesh,&cfg,Eabsorb,cfg.nphoton));
	}
	if(cfg.issave2pt){
		MMCDEBUG(&cfg,dlTime,(cfg.flog,"saving data ..."));
		mesh_saveweight(&mesh,&cfg);
	}

        MMCDEBUG(&cfg,dlTime,(cfg.flog,"\tdone\t%d\n",GetTimeMillis()-t0));

	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
