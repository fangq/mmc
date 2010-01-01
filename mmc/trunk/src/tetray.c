#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "simpmesh.h"
#include "tettracing.h"
#include "mcx_utils.h"


int main(int argc, char**argv){
	Config cfg;
	tetmesh mesh;
	tetplucker plucker;

	int faceid,i,eid,isend;
	int *enb;
	float dlen,leninit,weight;
	float3 pout,ptmp;
	float3 p0,p1,c0;
	
        mcx_initcfg(&cfg);

        // parse command line options to initialize the configurations
        mcx_parsecmd(argc,argv,&cfg);
	
	mesh_init(&mesh);
	mesh_loadnode(&mesh,&cfg);
	mesh_loadelem(&mesh,&cfg);
	mesh_loadfaceneighbor(&mesh,&cfg);
	mesh_loadmedia(&mesh,&cfg);

	plucker_init(&plucker,&mesh);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL||mesh.med==NULL)
		mesh_error("not all files were loaded");

	if(cfg.seed<0)
		srand(time(NULL));
	else
		srand(cfg.seed);

	/*launch photons*/
	for(i=0;i<cfg.nphoton;i++){
	    eid=cfg.dim.x;
	    weight=1.f;
	    /*initialize the photon position*/
	    leninit=rand01();
	    leninit=((leninit==0.f)?LOG_MT_MAX:(-log(leninit)));
	    memcpy(&p0,&cfg.srcpos,sizeof(p0));
	    memcpy(&c0,&cfg.srcdir,sizeof(p0));
	    vec_mult_add(&p0,&c0,1.0f,leninit,&p1);
	    dlen=dist2(&p0,&p1);

	    while(1){  /*propagate a photon until exit*/
		trackpos(&p0,&p1,&plucker,eid,&pout,&faceid,&weight,&isend);
		if(pout.x==QLIMIT){
			eid=0;
			faceid=-1;
		}
		/*move a photon until the end of the current scattering path*/
		while(faceid>=0 && !isend){
			memcpy((void *)&ptmp,(void *)&pout,sizeof(ptmp));

			enb=(int *)(&plucker.mesh->facenb[eid-1]);
			eid=enb[faceid];
			if(eid==0){
		    	    if(cfg.isrowmajor) fprintf(cfg.flog,"hit boundary, exit %d\n",i);
		    	    break;
			}
			if(pout.x!=QLIMIT&&cfg.isrowmajor){
				fprintf(cfg.flog,"ray passes at: %f %f %f %d\n",pout.x,pout.y,pout.z,eid);
			}
			trackpos(&ptmp,&p1,&plucker,eid,&pout,&faceid,&weight,&isend);
		}
		if(eid==0) {
			if(pout.x==QLIMIT)
                          fprintf(cfg.flog,"%d %d %f %f %f %f %f %f\n",i,eid,p0.x,p0.y,p0.z,p1.x,p1.y,p1.z);
			break;  /*photon exits boundary*/
		}
		memcpy((void *)&p0,(void *)&p1,sizeof(p0));
		if(cfg.isrowmajor) fprintf(cfg.flog,"ray exits at: %f %f %f %d %d\n",p0.x,p0.y,p0.z,eid,i);
		mc_next_scatter(mesh.med[mesh.type[eid]-1].g,mesh.med[mesh.type[eid]-1].musp,&p1,&c0);
	        dlen=dist2(&p0,&p1);
	    }
	}
	plucker_clear(&plucker);
	mesh_saveweight(&mesh,&cfg);
	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
