#include "tettracing.h"

void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
	pout->x=w->x*p1->x+w->y*p2->x+w->z*p3->x;
	pout->y=w->x*p1->y+w->y*p2->y+w->z*p3->y;
	pout->z=w->x*p1->z+w->y*p2->z+w->z*p3->z;
}
float getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
	float3 ptmp;
	float Rv;
	Rv=1.f/(w1+w2+w3);
	ptmp.x=w1*Rv;
	ptmp.y=w2*Rv;
	ptmp.z=w3*Rv;
	interppos(&ptmp,p1,p2,p3,pout);
	return Rv;
}

void trackpos(float3 *p0,float3 *p1,tetplucker *plucker,int eid /*start from 1*/, 
              float3 *pout, int *faceid, float *weight, int *isend){
	float3 pvec, pcrx;
	float3 pin;
	int *ee;
	int i;
	float w[6],Rv;
	int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
	int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
	int faceorder[]={1,3,2,0};
	
	if(plucker->mesh==NULL || plucker->d==NULL||eid<=0||eid>plucker->mesh->ne) 
		return;

	*faceid=-1;
	*isend=0;
	vec_diff(p0,p1,&pvec);
	vec_cross(p0,p1,&pcrx);
	ee=(int *)(plucker->mesh->elem+eid-1);
	for(i=0;i<6;i++){
		w[i]=pinner(&pvec,&pcrx,plucker->d+(eid-1)*6+i,plucker->m+(eid-1)*6+i);
	}
	pin.x=QLIMIT;
	pout->x=QLIMIT;
	for(i=0;i<4;i++){
	    if(w[fc[i][0]]||w[fc[i][1]]||w[fc[i][2]]){
	        if(i>=2) w[fc[i][1]]=-w[fc[i][1]]; // can not go back
		if(pin.x==QLIMIT&&w[fc[i][0]]<=0 && w[fc[i][1]]<=0 && w[fc[i][2]]>=0){
			// f_enter
			getinterp(-w[fc[i][0]],-w[fc[i][1]],w[fc[i][2]],
				plucker->mesh->node+ee[nc[i][0]]-1,
				plucker->mesh->node+ee[nc[i][1]]-1,
				plucker->mesh->node+ee[nc[i][2]]-1,&pin);
		}else if(pout->x==QLIMIT&&w[fc[i][0]]>=0 && w[fc[i][1]]>=0 && w[fc[i][2]]<=0){
			// f_leave
			Rv=getinterp(w[fc[i][0]],w[fc[i][1]],-w[fc[i][2]],
				plucker->mesh->node+ee[nc[i][0]]-1,
				plucker->mesh->node+ee[nc[i][1]]-1,
				plucker->mesh->node+ee[nc[i][2]]-1,pout);
			*faceid=faceorder[i];
			*isend=(dist2(p0,pout)>dist2(p0,p1));

			if(*isend)
			    *weight*=exp(-plucker->mesh->med[plucker->mesh->type[eid-1]-1].mua*dist(p0,p1));
			else
			    *weight*=exp(-plucker->mesh->med[plucker->mesh->type[eid-1]-1].mua*dist(p0,pout));

			Rv*=(*weight);
			plucker->mesh->weight[ee[0]-1]+=*weight*0.25f;
			plucker->mesh->weight[ee[1]-1]+=*weight*0.25f;
			plucker->mesh->weight[ee[2]-1]+=*weight*0.25f;
			plucker->mesh->weight[ee[3]-1]+=*weight*0.25f;

/*			plucker->mesh->weight[ee[nc[i][0]]-1]+=Rv*w[fc[i][0]];
			plucker->mesh->weight[ee[nc[i][1]]-1]+=Rv*w[fc[i][1]];
			plucker->mesh->weight[ee[nc[i][2]]-1]+=-Rv*w[fc[i][2]];*/
		}
	    }
	}
}
