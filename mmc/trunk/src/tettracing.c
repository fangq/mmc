#include <string.h>
#include <stdlib.h>
#include "tettracing.h"

void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
	pout->x=w->x*p1->x+w->y*p2->x+w->z*p3->x;
	pout->y=w->x*p1->y+w->y*p2->y+w->z*p3->y;
	pout->z=w->x*p1->z+w->y*p2->z+w->z*p3->z;
}
void getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
        pout->x=w1*p1->x+w2*p2->x+w3*p3->x;
        pout->y=w1*p1->y+w2*p2->y+w3*p3->y;
        pout->z=w1*p1->z+w2*p2->z+w3*p3->z;
}
/*
  when a photon is crossing a vertex or edge, pull the photon
  to the center of the element (slightly) and try again
*/
void fixphoton(float3 *p,float3 *nodes, int *ee){
        float3 c0={0.f,0.f,0.f};
	int i;
        /*calculate element centroid*/
	for(i=0;i<4;i++)
		vec_add(&c0,nodes+ee[i]-1,&c0);
	p->x+=(c0.x*0.25f-p->x)*FIX_PHOTON;
        p->y+=(c0.y*0.25f-p->y)*FIX_PHOTON;
        p->z+=(c0.z*0.25f-p->z)*FIX_PHOTON;
}
/*
  p0 and p1 ony determine the direction, slen determines the length
  so, how long is the vector p0->p1 does not matter. infact, the longer
  the less round off error when computing the Plucker coordinates.
*/
float trackpos(float3 *p0,float3 *pvec, tetplucker *plucker,int eid /*start from 1*/, 
              float3 *pout, float slen, int *faceid, float *weight, 
	      int *isend,float *photontimer, float rtstep, Config *cfg){
	float3 pcrx={0.f},p1={0.f};
	float3 pin={0.f};
	medium *prop;
	int *ee;
	int i;
	float w[6],Rv,ww,oldweight=*weight,ratio=0.f,newdlen=0.f,dlen=0.f; /*dlen is the physical distance*/
	int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
	int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
	int faceorder[]={1,3,2,0};
        float bary[2][4]={{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}};

	if(plucker->mesh==NULL || plucker->d==NULL||eid<=0||eid>plucker->mesh->ne) 
		return -1;

	*faceid=-1;
	*isend=0;
	vec_add(p0,pvec,&p1);
	vec_cross(p0,&p1,&pcrx);
	ee=(int *)(plucker->mesh->elem+eid-1);
	prop=plucker->mesh->med+(plucker->mesh->type[eid-1]-1);

#pragma unroll(6)
	for(i=0;i<6;i++){
		w[i]=pinner(pvec,&pcrx,plucker->d+(eid-1)*6+i,plucker->m+(eid-1)*6+i);
		/*if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"%f ",w[i]);*/
	}
	if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"%d \n",eid);

	pin.x=QLIMIT;
	pout->x=QLIMIT;
	for(i=0;i<4;i++){
	    if(w[fc[i][0]]||w[fc[i][1]]||w[fc[i][2]]){
		if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"testing face [%d]\n",i);

	        if(i>=2) w[fc[i][1]]=-w[fc[i][1]]; // can not go back
		if(pin.x==QLIMIT&&w[fc[i][0]]>=0.f && w[fc[i][1]]>=0.f && w[fc[i][2]]<=0.f){
			// f_enter
                        if(cfg->debuglevel&dlTracingEnter) fprintf(cfg->flog,"ray enters face %d[%d] of %d\n",i,faceorder[i],eid);

                        Rv=1.f/(-w[fc[i][0]]-w[fc[i][1]]+w[fc[i][2]]);
                        bary[0][nc[i][0]]=-w[fc[i][0]]*Rv;
                        bary[0][nc[i][1]]=-w[fc[i][1]]*Rv;
                        bary[0][nc[i][2]]=w[fc[i][2]]*Rv;

			getinterp(bary[0][nc[i][0]],bary[0][nc[i][1]],bary[0][nc[i][2]],
				plucker->mesh->node+ee[nc[i][0]]-1,
				plucker->mesh->node+ee[nc[i][1]]-1,
				plucker->mesh->node+ee[nc[i][2]]-1,&pin);

                        if(cfg->debuglevel&dlTracingEnter) fprintf(cfg->flog,"entrance point %f %f %f\n",pin.x,pin.y,pin.z);

		}else if(pout->x==QLIMIT&&w[fc[i][0]]<=0.f && w[fc[i][1]]<=0.f && w[fc[i][2]]>=0.f){
			// f_leave
                        if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"ray exits face %d[%d] of %d\n",i,faceorder[i],eid);

                        Rv=1.f/(w[fc[i][0]]+w[fc[i][1]]-w[fc[i][2]]);
                        bary[1][nc[i][0]]=w[fc[i][0]]*Rv;
                        bary[1][nc[i][1]]=w[fc[i][1]]*Rv;
                        bary[1][nc[i][2]]=-w[fc[i][2]]*Rv;

			getinterp(bary[1][nc[i][0]],bary[1][nc[i][1]],bary[1][nc[i][2]],
				plucker->mesh->node+ee[nc[i][0]]-1,
				plucker->mesh->node+ee[nc[i][1]]-1,
				plucker->mesh->node+ee[nc[i][2]]-1,pout);

                        if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"exit point %f %f %f\n",pout->x,pout->y,pout->z);

			newdlen=dist(p0,pout);
			/*ratio=newdlen;*/
			dlen=slen/prop->mus;
			*faceid=faceorder[i];
			*isend=(newdlen>dlen);
			newdlen=((*isend) ? dlen : newdlen);
		}
	    }
	}
        if(pin.x!=QLIMIT && pout->x!=QLIMIT){
		*photontimer+=newdlen*prop->n*R_C0;
		if(*photontimer>=cfg->tend){ /*exit time window*/
		   *faceid=-2;
	           pout->x=QLIMIT;
		   return 0.f;
		}
		*weight*=exp(-prop->mua*newdlen);
		slen-=newdlen*prop->mus;
                p0->x+=newdlen*pvec->x;
                p0->y+=newdlen*pvec->y;
                p0->z+=newdlen*pvec->z;
                if(cfg->debuglevel&dlWeight) fprintf(cfg->flog,"update weight to %f and path end %d \n",*weight,*isend);

		int tshift=(int)((*photontimer-cfg->tstart)*rtstep)*plucker->mesh->nn;
                ww=(oldweight-(*weight))*0.5f;
		/*ratio/=dist(&pin,pout);*/
                if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"barycentric [%f %f %f %f] [%f %f %f %f]\n",
                      bary[0][0],bary[0][1],bary[0][2],bary[0][3],bary[1][0],bary[1][1],bary[1][2],bary[1][3]);

                if(cfg->debuglevel&dlDist) fprintf(cfg->flog,"distances pin-p0: %f p0-pout: %f pin-pout: %f/%f p0-p1: %f\n",
                      dist(&pin,p0),dist(p0,pout),dist(&pin,pout),dist(&pin,p0)+dist(p0,pout)-dist(&pin,pout),dlen);
#pragma unroll(4)
                for(i=0;i<4;i++)
		     plucker->mesh->weight[ee[i]-1+tshift]+=ww*(bary[0][i]+bary[1][i]);
/*		     plucker->mesh->weight[ee[i]-1+tshift]+=ww*(ratio*bary[0][i]+(1.f-ratio)*bary[1][i]);*/
        }
	return slen;
}

float onephoton(int id,tetplucker *plucker,tetmesh *mesh,Config *cfg,float rtstep,RandType *ran, RandType *ran0){
	int faceid,eid,isend,fixcount=0;
	int *enb;
	float slen,weight,photontimer;
	float3 pout;
	float3 p0,c0;

	/*initialize the photon parameters*/
	eid=cfg->dim.x;
	weight=1.f;
	photontimer=0.f;
	slen=rand_next_scatlen(ran);
	memcpy(&p0,&(cfg->srcpos),sizeof(p0));
	memcpy(&c0,&(cfg->srcdir),sizeof(p0));

	while(1){  /*propagate a photon until exit*/
	    slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,rtstep,cfg);
	    if(pout.x==QLIMIT){
	    	  if(faceid==-2) break; /*reaches time gate*/
		  if(fixcount++<MAX_TRIAL){
			fixphoton(&p0,mesh->node,(int *)(mesh->elem+eid-1));
			continue;
                  }
	    	  eid=-eid;
        	  faceid=-1;
	    }
	    /*move a photon until the end of the current scattering path*/
	    while(faceid>=0 && !isend){
	    	    memcpy((void *)&p0,(void *)&pout,sizeof(p0));

	    	    enb=(int *)(mesh->facenb+eid-1);
	    	    eid=enb[faceid];
	    	    if(eid==0){
	    		if(cfg->debuglevel&dlMove) fprintf(cfg->flog,"hit boundary, exit %d\n",id);
	    		break;
	    	    }
	    	    if(pout.x!=QLIMIT && (cfg->debuglevel&dlMove)){
	    		fprintf(cfg->flog,"pass at: %f %f %f %d\n",pout.x,pout.y,pout.z,eid);
	    	    }
	    	    slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,rtstep,cfg);
		    if(faceid==-2) break;
		    fixcount=0;
		    while(pout.x==QLIMIT && fixcount++<MAX_TRIAL){
			fixphoton(&p0,mesh->node,(int *)(mesh->elem+eid-1));
                        slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,rtstep,cfg);
		    }
        	    if(pout.x==QLIMIT){
        		/*possibily hit an edge or miss*/
			eid=-eid;
        		break;
        	    }
	    }
	    if(eid<=0 || pout.x==QLIMIT) {
        	    if(eid==0 && (cfg->debuglevel&dlMove))
        		 fprintf(cfg->flog,"hit boundary: %d %d %f %f %f\n",id,eid,p0.x,p0.y,p0.z);
		    else if(faceid==-2 && (cfg->debuglevel&dlMove))
                         fprintf(cfg->flog,"time window ends: %d %d %f %f %f\n",id,eid,p0.x,p0.y,p0.z);
	    	    else if(eid && faceid!=-2  && cfg->debuglevel&dlEdge)
        		 fprintf(cfg->flog,"hit edge or vertex: %d %d %f %f %f\n",id,eid,p0.x,p0.y,p0.z);
	    	    break;  /*photon exits boundary*/
	    }
	    if(cfg->debuglevel&dlMove) fprintf(cfg->flog,"move to: %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,id,slen);
	    slen=mc_next_scatter(mesh->med[mesh->type[eid]-1].g,mesh->med[mesh->type[eid]-1].mus,&c0,ran,ran0,cfg);
	}
	return weight;
}
