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
	float slen,leninit,weight;
	float3 pout;
	float3 p0,c0;
	
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

	/*launch photons*/
#pragma omp parallel for schedule(static, 30) default(shared) private(i,eid,weight,slen,leninit,p0,c0,pout,faceid,isend,enb)
	for(i=0;i<cfg.nphoton;i++){
	    eid=cfg.dim.x;
	    weight=1.f;
	    /*initialize the photon position*/
	    leninit=rand01();
	    leninit=((leninit==0.f)?LOG_MT_MAX:(-log(leninit)));
	    memcpy(&p0,&cfg.srcpos,sizeof(p0));
	    memcpy(&c0,&cfg.srcdir,sizeof(p0));
	    slen=leninit;

	    while(1){  /*propagate a photon until exit*/
		slen=trackpos(&p0,&c0,&plucker,eid,&pout,slen,&faceid,&weight,&isend,&cfg);
		if(pout.x==QLIMIT){
                      faceid=-1;
		}
		/*move a photon until the end of the current scattering path*/
		while(faceid>=0 && !isend){
			memcpy((void *)&p0,(void *)&pout,sizeof(p0));

			enb=(int *)(&plucker.mesh->facenb[eid-1]);
			eid=enb[faceid];
			if(eid==0){
		    	    if(cfg.debuglevel&dlMove) fprintf(cfg.flog,"hit boundary, exit %d\n",i);
		    	    break;
			}
			if(pout.x!=QLIMIT&&cfg.debuglevel&dlMove){
			    fprintf(cfg.flog,"passes at: %f %f %f %d\n",pout.x,pout.y,pout.z,eid);
			}
                        if(pout.x==QLIMIT){
                              /*possibily hit an edge or miss*/
                              break;
                        }
			slen=trackpos(&p0,&c0,&plucker,eid,&pout,slen,&faceid,&weight,&isend,&cfg);
		}
		if(eid==0) {
			if(pout.x==QLIMIT&&cfg.debuglevel&dlMove)
                             fprintf(cfg.flog,"hit edge or miss: %d %d %e\n",i,eid,dist(&p0,&pout));
			break;  /*photon exits boundary*/
		}
		if(cfg.debuglevel&dlMove) fprintf(cfg.flog,"moves to: %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,i,slen);
		slen=mc_next_scatter(mesh.med[mesh.type[eid]-1].g,mesh.med[mesh.type[eid]-1].mus,&c0,&cfg);
	    }
	}
	if(cfg.isnormalized)
	  mesh_normalize(&mesh);
	plucker_clear(&plucker);
	mesh_saveweight(&mesh,&cfg);
	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);

	return 0;
}
