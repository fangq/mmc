#include <stdlib.h>
#include <time.h>
#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"

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
	int i,threadid=0;

        mcx_initcfg(&cfg);

        /* parse command line options to initialize the configurations */
        mcx_parsecmd(argc,argv,&cfg);
	
	mesh_init(&mesh);
	mesh_loadnode(&mesh,&cfg);
	mesh_loadelem(&mesh,&cfg);
	mesh_loadfaceneighbor(&mesh,&cfg);
	mesh_loadmedia(&mesh,&cfg);
	mesh_loadelemvol(&mesh,&cfg);

	plucker_init(&plucker,&mesh);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL||mesh.med==NULL)
		mesh_error("not all files were loaded");

	rtstep=1.f/cfg.tstep;

	if(cfg.seed<0) cfg.seed=time(NULL);

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
	}
}
	if(cfg.isnormalized)
	  mesh_normalize(&mesh,&cfg,Eabsorb,cfg.nphoton);
	plucker_clear(&plucker);
	mesh_saveweight(&mesh,&cfg);
	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
