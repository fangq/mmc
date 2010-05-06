#include <stdlib.h>
#include <time.h>
#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"


int main(int argc, char**argv){
	Config cfg;
	tetmesh mesh;
	tetplucker plucker;
	float rtstep,Eescape=0.f;
	int i;
	
        mcx_initcfg(&cfg);

        // parse command line options to initialize the configurations
        mcx_parsecmd(argc,argv,&cfg);
	
	mesh_init(&mesh);
	mesh_loadnode(&mesh,&cfg);
	mesh_loadelem(&mesh,&cfg);
	mesh_loadfaceneighbor(&mesh,&cfg);
	mesh_loadmedia(&mesh,&cfg);
	mesh_loadnodevol(&mesh,&cfg);

	plucker_init(&plucker,&mesh);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL||mesh.med==NULL)
		mesh_error("not all files were loaded");

	if(cfg.seed<0)
		srand(time(NULL));
	else
		srand(cfg.seed);

	rtstep=1.f/cfg.tstep;
	/*launch photons*/

#pragma omp parallel for schedule(static) default(shared) reduction(+:Eescape)
	for(i=0;i<cfg.nphoton;i++){
		Eescape+=onephoton(&plucker,&mesh,&cfg,rtstep);
	}
	if(cfg.isnormalized)
	  mesh_normalize(&mesh,&cfg,cfg.nphoton-Eescape,cfg.nphoton);
	plucker_clear(&plucker);
	mesh_saveweight(&mesh,&cfg);
	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
