#include "simpmesh.h"
#include "tettracing.h"
#include <string.h>
#include <stdlib.h>
/*#include <math.h>*/
#include <time.h>

int main(int argc, char**argv){
	tetmesh mesh;
	tetplucker plucker;
	int faceid,i,eid;
	int *enb;
	float dlen;
	float3 pout,ptmp;
	float3 p1={32.0967,26.4944,12.6675}; /*{41.5,37,2};*/
	float3 p0={41.5,38,2};
	
	if(argc<4) {
		mesh_error("format: tetray nodefile elemfile");
	}
	
	mesh_init(&mesh);
	mesh_loadnode(&mesh,argv[1]);
	mesh_loadelem(&mesh,argv[2]);
	mesh_loadfaceneighbor(&mesh,argv[3]);
	plucker_init(&plucker,&mesh);
	
	if(mesh.node==NULL||mesh.elem==NULL||mesh.facenb==NULL)
		mesh_error("not all files were loaded");

	srand(time(NULL));
	
	eid=1;
	for(i=0;i<10;i++){
	    dlen=dist2(&p0,&p1);
	    trackpos(&p0,&p1,&plucker,eid,&pout,&faceid);
	    while(faceid>=0&&dlen>dist2(&p0,&pout)){
		    memcpy((void *)&ptmp,(void *)&pout,sizeof(ptmp));
		    enb=(int *)(&plucker.mesh->facenb[eid-1]);
		    eid=enb[faceid];
		    if(eid==0) {
		    	printf("hit boundary, exit");
		    	break;
		    }
		    if(pout.x!=QLIMIT){
			    printf("ray exits at: %f %f %f %d\n",pout.x,pout.y,pout.z,eid);
		    }
		    trackpos(&ptmp,&p1,&plucker,eid,&pout,&faceid);
	    }
	    if(eid==0) break;
	    memcpy((void *)&p0,(void *)&p1,sizeof(p0));
	    p1.x+=2.*rand()/RAND_MAX;
	    p1.y+=2.*rand()/RAND_MAX;
	    p1.z+=2.*rand()/RAND_MAX;
	}
	plucker_clear(&plucker);
	mesh_clear(&mesh);


	return 0;
}
