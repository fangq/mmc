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
	      int *isend,float *photontimer, float *Eabsorb, float rtstep, Config *cfg){
	float3 pcrx={0.f},p1={0.f};
	float3 pin={0.f};
	medium *prop;
	int *ee;
	int i,tshift;
	float w[6],Rv,ww,currweight,ratio=0.f,dlen=0.f,rc; /*dlen is the physical distance*/
	const int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
	const int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
	const int faceorder[]={1,3,2,0};
        float bary[2][4]={{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}};
	float Lp0=0.f,Lio=0.f,Lmove=0.f,atte;

	if(plucker->mesh==NULL || plucker->d==NULL||eid<=0||eid>plucker->mesh->ne) 
		return -1;

	*faceid=-1;
	*isend=0;
	vec_add(p0,pvec,&p1);
	vec_cross(p0,&p1,&pcrx);
	ee=(int *)(plucker->mesh->elem+eid-1);
	prop=plucker->mesh->med+(plucker->mesh->type[eid-1]-1);
	atte=plucker->mesh->atte[plucker->mesh->type[eid-1]-1];
	rc=prop->n*R_C0;
        currweight=*weight;

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
                        if(pout->x!=QLIMIT) break;

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

			Lp0=dist(p0,pout);
			dlen=slen/prop->mus;
			*faceid=faceorder[i];
			*isend=(Lp0>dlen);
			Lmove=((*isend) ? dlen : Lp0);
			if(pin.x!=QLIMIT) break;
		}
	    }
	}
        if(pin.x!=QLIMIT && pout->x!=QLIMIT){
		if(*photontimer+Lmove*rc>=cfg->tend){ /*exit time window*/
		   *faceid=-2;
	           pout->x=QLIMIT;
		   Lmove=(cfg->tend-*photontimer)/(prop->n*R_C0)-1e-4f;
		}
		*weight*=exp(-prop->mua*Lmove);
		slen-=Lmove*prop->mus;
                p0->x+=Lmove*pvec->x;
                p0->y+=Lmove*pvec->y;
                p0->z+=Lmove*pvec->z;
                if(cfg->debuglevel&dlWeight) fprintf(cfg->flog,"update weight to %f and path end %d \n",*weight,*isend);

		Lio=dist(&pin,pout);

                if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"Y [%f %f %f %f] [%f %f %f %f]\n",
                      bary[0][0],bary[0][1],bary[0][2],bary[0][3],bary[1][0],bary[1][1],bary[1][2],bary[1][3]);

                if(cfg->debuglevel&dlDist) fprintf(cfg->flog,"D %f p0-pout: %f pin-pout: %f/%f p0-p1: %f\n",
                      dist(&pin,p0),Lp0,Lio,dist(&pin,p0)+dist(p0,pout)-dist(&pin,pout),dlen);

#ifndef SIMPLE_TRACING
		if(Lio>EPS){
		    Lio=1.f/Lio;
		    for(dlen=0;dlen<Lmove;dlen+=cfg->minstep){
			ratio=(Lp0-dlen)*Lio;
			ww=currweight;
			currweight*=atte;
			ww-=(currweight>*weight)?currweight:*weight;
			*photontimer+=((currweight>*weight)?cfg->minstep:Lmove-dlen)*rc;
			*Eabsorb+=ww;
			ww/=prop->mua; /*accumulate fluence*volume, not energy deposit*/

			tshift=(int)((*photontimer-cfg->tstart)*rtstep)*plucker->mesh->nn;

			/*ww will have the volume effect. volume of each nodes will be divided in the end*/

	                if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
			   p0->x-(Lmove-dlen)*pvec->x,p0->y-(Lmove-dlen)*pvec->y,p0->z-(Lmove-dlen)*pvec->z,ww,eid,dlen);

			for(i=0;i<4;i++)
				plucker->mesh->weight[ee[i]-1+tshift]+=ww*(ratio*bary[0][i]+(1.f-ratio)*bary[1][i]);
		    }
		}
#else
                if(Lio>EPS){
	                ratio=(Lp0-Lmove)/Lio;
        	        ww=currweight-*weight;
			*Eabsorb+=ww;
			ww/=prop->mua;
#pragma unroll(4)
                	for(i=0;i<4;i++)
                     		plucker->mesh->weight[ee[i]-1+tshift]+=ww*(ratio*bary[0][i]+(1.f-ratio)*bary[1][i]);
		}
#endif
		if(*faceid==-2)
		   return 0.f;
        }
	return slen;
}

float onephoton(int id,tetplucker *plucker,tetmesh *mesh,Config *cfg,float rtstep,RandType *ran, RandType *ran0){
	int faceid,eid,isend,fixcount=0;
	int *enb;
	float slen,weight,photontimer,Eabsorb;
	float3 pout;
	float3 p0,c0;

	/*initialize the photon parameters*/
	eid=cfg->dim.x;
	weight=1.f;
	photontimer=0.f;
	Eabsorb=0.f;
	slen=rand_next_scatlen(ran);
	memcpy(&p0,&(cfg->srcpos),sizeof(p0));
	memcpy(&c0,&(cfg->srcdir),sizeof(p0));

	while(1){  /*propagate a photon until exit*/
	    slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,&Eabsorb,rtstep,cfg);
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
	    	    if(eid==0) break;
	    	    if(pout.x!=QLIMIT && (cfg->debuglevel&dlMove))
	    		fprintf(cfg->flog,"P %f %f %f %d %d %f\n",pout.x,pout.y,pout.z,eid,id,slen);

	    	    slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,&Eabsorb,rtstep,cfg);
		    if(faceid==-2) break;
		    fixcount=0;
		    while(pout.x==QLIMIT && fixcount++<MAX_TRIAL){
			fixphoton(&p0,mesh->node,(int *)(mesh->elem+eid-1));
                        slen=trackpos(&p0,&c0,plucker,eid,&pout,slen,&faceid,&weight,&isend,&photontimer,&Eabsorb,rtstep,cfg);
		    }
        	    if(pout.x==QLIMIT){
        		/*possibily hit an edge or miss*/
			eid=-eid;
        		break;
        	    }
	    }
	    if(eid<=0 || pout.x==QLIMIT) {
        	    if(eid==0 && (cfg->debuglevel&dlMove))
        		 fprintf(cfg->flog,"B %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,id,slen);
		    else if(faceid==-2 && (cfg->debuglevel&dlMove))
                         fprintf(cfg->flog,"T %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,id,slen);
	    	    else if(eid && faceid!=-2  && cfg->debuglevel&dlEdge)
        		 fprintf(cfg->flog,"X %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,id,slen);
	    	    break;  /*photon exits boundary*/
	    }
	    if(cfg->debuglevel&dlMove) fprintf(cfg->flog,"M %f %f %f %d %d %f\n",p0.x,p0.y,p0.z,eid,id,slen);
	    slen=mc_next_scatter(mesh->med[mesh->type[eid]-1].g,&c0,ran,ran0,cfg);
	}
	return Eabsorb;
}
