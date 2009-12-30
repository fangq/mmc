#include "simpmesh.h"
#include "tettracing.h"
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define PHOTON_COUNT 1

int main(int argc, char**argv){
	tetmesh mesh;
	tetplucker plucker;
	int faceid,i,eid;
	int *enb;
	float dlen,leninit,movelen,weight;
	float3 pout,ptmp;
	float3 p0,p1; /*{32.0967f,26.4944f,12.6675f};*/
	float3 c0={-0.577350269189626f,-0.577350269189626f,0.577350269189626f};
	float3 psrc={41.5f,38.f,2.f};
	
	if(argc<5) {
		mesh_error("format: tetray nodefile elemfile facenbfile propfile");
	}
	
	mesh_init(&mesh);
	mesh_loadnode(&mesh,argv[1]);
	mesh_loadelem(&mesh,argv[2]);
	mesh_loadfaceneighbor(&mesh,argv[3]);
	mesh_loadmedia(&mesh,argv[4]);

	plucker_init(&plucker,&mesh);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL||mesh.med==NULL)
		mesh_error("not all files were loaded");

	srand(time(NULL));
	
	/*launch photons*/
	for(i=0;i<PHOTON_COUNT;i++){
	    eid=1;
	    weight=1.f;
	    /*initialize the photon position*/
	    leninit=rand01();
	    leninit=((leninit==0.f)?LOG_MT_MAX:(-log(leninit)));
	    memcpy(&p0,&psrc,sizeof(p0));
	    vec_mult_add(&psrc,&c0,1.0f,leninit,&p1);

	    while(1){  /*propagate a photon until exit*/
		dlen=dist2(&p0,&p1);
		trackpos(&p0,&p1,&plucker,eid,&pout,&faceid,&weight);
		movelen=dist2(&p0,&pout);

		/*move a photon until end of the current scattering len*/
		while(faceid>=0&&dlen>movelen){
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


	return 0;
}
