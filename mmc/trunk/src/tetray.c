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

	int faceid,i,eid;
	int *enb;
	float leninit,weight;
	float3 pout,ptmp;
	float3 p0,p1; /*{32.0967f,26.4944f,12.6675f};*/
	float3 c0={-0.577350269189626f,-0.577350269189626f,0.577350269189626f};
	float3 psrc={41.5f,38.f,2.f};
	
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

	/*srand(time(NULL));*/
	srand(314159265);

	
	/*launch photons*/
	for(i=0;i<cfg->nphoton;i++){
	    eid=1;
	    weight=1.f;
	    /*initialize the photon position*/
	    leninit=rand01();
	    leninit=((leninit==0.f)?LOG_MT_MAX:(-log(leninit)));
	    memcpy(&p0,&psrc,sizeof(p0));
	    vec_mult_add(&psrc,&c0,1.0f,leninit,&p1);

	    while(1){  /*propagate a photon until exit*/
		/*possible overshoot when p0-p1 is shorter*/
		trackpos(&p0,&p1,&plucker,eid,&pout,&faceid,&weight); 

		/*move a photon until end of the current scattering len*/
		while(faceid>=0 && dist2(&p0,&p1)>dist2(&p0,&pout)){
			memcpy((void *)&ptmp,(void *)&pout,sizeof(ptmp));

			enb=(int *)(&plucker.mesh->facenb[eid-1]);
			eid=enb[faceid];
			if(eid==0) {
		    	    printf("hit boundary, exit %d\n",i);
		    	    break;
			}
			if(pout.x!=QLIMIT){
				printf("ray passes at: %f %f %f %d\n",pout.x,pout.y,pout.z,eid);
			}
			trackpos(&ptmp,&p1,&plucker,eid,&pout,&faceid,&weight);
		}
		if(eid==0) break;  /*photon exits boundary*/
		memcpy((void *)&p0,(void *)&p1,sizeof(p0));
		printf("ray exits at: %f %f %f %d\n",p0.x,p0.y,p0.z,eid);
		mc_next_scatter(mesh.med[mesh.type[eid]-1].g,mesh.med[mesh.type[eid]-1].musp,&p1);
	    }
	}
	plucker_clear(&plucker);
	mesh_clear(&mesh);
        mcx_clearcfg(&cfg);


	return 0;
}
