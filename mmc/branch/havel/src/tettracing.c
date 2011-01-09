/*******************************************************************************
**  Mesh-based Monte Carlo (MMC)
**
**  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates," Biomed. Opt. Express, 1(1) 165-175 (2010)
**
**  (Fang2009) Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon 
**          Migration in 3D Turbid Media Accelerated by Graphics Processing 
**          Units," Optics Express, 17(22) 20178-20190 (2009)
**
**  License: GPL v3, see LICENSE.txt for details
**
*******************************************************************************/

/***************************************************************************//**
\file    tettracing.c

\brief   Core unit for Plucker-coordinate-based ray-tracing
*******************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "tettracing.h"

#ifdef MMC_USE_SSE
#include <smmintrin.h>
__m128 int_coef;
#endif

#define F32N(a) ((a) & 0x80000000)
#define F32P(a) ((a) ^ 0x80000000)

const int fc[4][3]={{0,4,2},{3,5,4},{2,5,1},{1,3,0}};
const int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}};
const int faceorder[]={1,3,2,0,-1};
const int ifaceorder[]={3,0,2,1};
const int fm[]={2,0,1,3};
const char maskmap[16]={4,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3};

inline void interppos(float3 *w,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
	pout->x=w->x*p1->x+w->y*p2->x+w->z*p3->x;
	pout->y=w->x*p1->y+w->y*p2->y+w->z*p3->y;
	pout->z=w->x*p1->z+w->y*p2->z+w->z*p3->z;
}
inline void getinterp(float w1,float w2,float w3,float3 *p1,float3 *p2,float3 *p3,float3 *pout){
        pout->x=w1*p1->x+w2*p2->x+w3*p3->x;
        pout->y=w1*p1->y+w2*p2->y+w3*p3->y;
        pout->z=w1*p1->z+w2*p2->z+w3*p3->z;
}

/*
  when a photon is crossing a vertex or edge, (slightly) pull the
  photon toward the center of the element and try again
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

float plucker_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){

	float3 pcrx={0.f},p1={0.f};
	float3 pin={0.f};
	medium *prop;
	int *ee;
	int i,tshift,eid;
	float w[6],Rv,ww,currweight,ratio=0.f,dlen=0.f,rc; /*dlen is the physical distance*/
	unsigned int *wi=(unsigned int*)w;
        float bary[2][4]={{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}};
	float Lp0=0.f,Lio=0.f,atte;

	if(tracer->mesh==NULL || tracer->d==NULL||r->eid<=0||r->eid>tracer->mesh->ne) 
		return -1;

	eid=r->eid-1;
	r->faceid=-1;
	r->isend=0;
	vec_add(&(r->p0),&(r->vec),&p1);
	vec_cross(&(r->p0),&p1,&pcrx);
	ee=(int *)(tracer->mesh->elem+eid);
	prop=tracer->mesh->med+(tracer->mesh->type[eid]-1);
	atte=tracer->mesh->atte[tracer->mesh->type[eid]-1];
	rc=prop->n*R_C0;
        currweight=r->weight;

	for(i=0;i<6;i++)
		w[i]=pinner(&(r->vec),&pcrx,tracer->d+(eid)*6+i,tracer->m+(eid)*6+i);

	if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"%d \n",eid);

	pin.x=MMC_UNDEFINED;
	r->pout.x=MMC_UNDEFINED;
	for(i=0;i<4;i++){
		if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"testing face [%d]\n",i);

	        if(i>=2) w[fc[i][1]]=-w[fc[i][1]]; // can not go back
		if(F32P(wi[fc[i][0]]) & F32P(wi[fc[i][1]]) & F32N(wi[fc[i][2]])){
			// f_enter
                        if(cfg->debuglevel&dlTracingEnter) fprintf(cfg->flog,"ray enters face %d[%d] of %d\n",i,faceorder[i],eid);

                        Rv=1.f/(w[fc[i][0]]+w[fc[i][1]]-w[fc[i][2]]);
                        bary[0][nc[i][0]]=w[fc[i][0]]*Rv;
                        bary[0][nc[i][1]]=w[fc[i][1]]*Rv;
                        bary[0][nc[i][2]]=-w[fc[i][2]]*Rv;

			getinterp(bary[0][nc[i][0]],bary[0][nc[i][1]],bary[0][nc[i][2]],
				tracer->mesh->node+ee[nc[i][0]]-1,
				tracer->mesh->node+ee[nc[i][1]]-1,
				tracer->mesh->node+ee[nc[i][2]]-1,&pin);

                        if(cfg->debuglevel&dlTracingEnter) fprintf(cfg->flog,"entrance point %f %f %f\n",pin.x,pin.y,pin.z);
                        if(r->pout.x!=MMC_UNDEFINED) break;

		}else if(F32N(wi[fc[i][0]]) & F32N(wi[fc[i][1]]) & F32P(wi[fc[i][2]])){
			// f_leave
                        if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"ray exits face %d[%d] of %d\n",i,faceorder[i],eid);

                        Rv=1.f/(-w[fc[i][0]]-w[fc[i][1]]+w[fc[i][2]]);
                        bary[1][nc[i][0]]=-w[fc[i][0]]*Rv;
                        bary[1][nc[i][1]]=-w[fc[i][1]]*Rv;
                        bary[1][nc[i][2]]=w[fc[i][2]]*Rv;

			getinterp(bary[1][nc[i][0]],bary[1][nc[i][1]],bary[1][nc[i][2]],
				tracer->mesh->node+ee[nc[i][0]]-1,
				tracer->mesh->node+ee[nc[i][1]]-1,
				tracer->mesh->node+ee[nc[i][2]]-1,&(r->pout));

                        if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"exit point %f %f %f\n",r->pout.x,r->pout.y,r->pout.z);

			Lp0=dist(&(r->p0),&(r->pout));
			dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
			r->faceid=faceorder[i];
#ifdef MMC_USE_SSE
		{
			int *enb=(int *)(tracer->mesh->facenb+eid);
			int nexteid=(enb[r->faceid]-1)*6;
	                if(nexteid){
        	                _mm_prefetch((char *)&((tracer->m+nexteid)->x),_MM_HINT_T0);
                	        _mm_prefetch((char *)&((tracer->m+(nexteid)*6+4)->x),_MM_HINT_T0);
                        	_mm_prefetch((char *)&((tracer->d+(nexteid)*6)->x),_MM_HINT_T0);
                                _mm_prefetch((char *)&((tracer->d+(nexteid)*6+4)->x),_MM_HINT_T0);
	                }
		}
#endif
			r->isend=(Lp0>dlen);
			r->Lmove=((r->isend) ? dlen : Lp0);
			if(pin.x!=MMC_UNDEFINED) break;
		}
	}
	(visit->raytet)++;
        if(pin.x!=MMC_UNDEFINED && r->pout.x!=MMC_UNDEFINED){
		if(r->photontimer+r->Lmove*rc>=cfg->tend){ /*exit time window*/
		   r->faceid=-2;
	           r->pout.x=MMC_UNDEFINED;
		   r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
		}
#ifdef __INTEL_COMPILER
		r->weight*=expf(-prop->mua*r->Lmove);
#else
		r->weight*=fast_expf9(-prop->mua*r->Lmove);
#endif
		r->slen-=r->Lmove*prop->mus;
                r->p0.x+=r->Lmove*r->vec.x;
                r->p0.y+=r->Lmove*r->vec.y;
                r->p0.z+=r->Lmove*r->vec.z;
                if(cfg->debuglevel&dlWeight) fprintf(cfg->flog,"update weight to %f and path end %d \n",r->weight,r->isend);

		if(!cfg->basisorder){
                        ww=currweight-r->weight;
                        r->Eabsorb+=ww;
                        r->photontimer+=r->Lmove*rc;
                        tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->ne;
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
			Lio=dist(&pin,&(r->pout));

	                if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"Y [%f %f %f %f] [%f %f %f %f]\n",
        	              bary[0][0],bary[0][1],bary[0][2],bary[0][3],bary[1][0],bary[1][1],bary[1][2],bary[1][3]);

                	if(cfg->debuglevel&dlDist) fprintf(cfg->flog,"D %f r->p0-r->pout: %f pin-r->pout: %f/%f r->p0-p1: %f\n",
	                      dist(&pin,&(r->p0)),Lp0,Lio,dist(&pin,&(r->p0))+dist(&(r->p0),&(r->pout))-dist(&pin,&(r->pout)),dlen);

#ifdef MAXSTEP_ADDITION
		if(Lio>EPS){
		    Lio=1.f/Lio;
		    for(dlen=0;dlen<r->Lmove;dlen+=cfg->minstep){
			ratio=(Lp0-dlen)*Lio;
			ww=currweight;
			currweight*=atte;
			ww-=(currweight>r->weight)?currweight:r->weight;
			r->photontimer+=((currweight>r->weight)?cfg->minstep:r->Lmove-dlen)*rc;
			r->Eabsorb+=ww;
			ww/=prop->mua; /*accumulate fluence*volume, not energy deposit*/

			tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->nn;

			/*ww has the volume effect. volume of each nodes will be divided in the end*/

	                if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
			   r->p0.x-(r->Lmove-dlen)*r->vec.x,r->p0.y-(r->Lmove-dlen)*r->vec.y,r->p0.z-(r->Lmove-dlen)*r->vec.z,ww,eid,dlen);

			for(i=0;i<4;i++)
				tracer->mesh->weight[ee[i]-1+tshift]+=ww*(ratio*bary[0][i]+(1.f-ratio)*bary[1][i]);
		    }
		}
#else
                if(Lio>EPS){
	                ratio=(Lp0-r->Lmove*0.5f)/Lio; /*add the weight to the middle-point*/
        	        ww=currweight-r->weight;
			r->Eabsorb+=ww;
			ww/=prop->mua;
                        r->photontimer+=r->Lmove*rc;
                        tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->nn;

                        if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
                           r->p0.x-(r->Lmove*0.5f)*r->vec.x,r->p0.y-(r->Lmove*0.5f)*r->vec.y,r->p0.z-(r->Lmove*0.5f)*r->vec.z,ww,eid,dlen);

//#pragma unroll(4)
                	for(i=0;i<4;i++)
                     		tracer->mesh->weight[ee[i]-1+tshift]+=ww*(ratio*bary[0][i]+(1.f-ratio)*bary[1][i]);
		}
#endif
		}
		if(r->faceid==-2)
		   return 0.f;
        }
	return r->slen;
}

#ifdef MMC_USE_SSE

inline __m128 rcp_nr(const __m128 a){
    const __m128 r = _mm_rcp_ps(a);
    return _mm_sub_ps(_mm_add_ps(r, r),_mm_mul_ps(_mm_mul_ps(r, a), r));
}

inline int havel_sse4(float3 *vecN, float3 *bary, const __m128 o,const __m128 d){
    const __m128 n = _mm_load_ps(&vecN->x);
    const __m128 det = _mm_dp_ps(n, d, 0x7f);
    if(_mm_movemask_ps(det) & 1)
        return 0;
    const __m128 dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
    const __m128 oldt = _mm_load_ss(&bary->x);
                
    if((_mm_movemask_ps(_mm_xor_ps(dett, _mm_sub_ss(_mm_mul_ss(oldt, det), dett)))&1) == 0){
        const __m128 detp = _mm_add_ps(_mm_mul_ps(o, det),_mm_mul_ps(dett, d));
        const __m128 detu = _mm_dp_ps(detp,_mm_load_ps(&((vecN+1)->x)), 0xf1);
                    
        if((_mm_movemask_ps(_mm_xor_ps(detu,_mm_sub_ss(det, detu)))&1) == 0){
            const __m128 detv = _mm_dp_ps(detp, _mm_load_ps(&((vecN+2)->x)), 0xf1);
  
            if((_mm_movemask_ps(_mm_xor_ps(detv,_mm_sub_ss(det, _mm_add_ss(detu, detv))))&1) == 0){
                const __m128 inv_det = rcp_nr(det);
                _mm_store_ss(&bary->x, _mm_mul_ss(dett, inv_det));
                _mm_store_ss(&bary->y, _mm_mul_ss(detu, inv_det));
                _mm_store_ss(&bary->z, _mm_mul_ss(detv, inv_det));
                //_mm_store_ps(&pout->x, _mm_mul_ps(detp,_mm_shuffle_ps(inv_det, inv_det, 0)));
                return 1;
            }
        }
    }                   
    return 0;
}

float havel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){

	float3 bary={1e10f,0.f,0.f,0.f};
	float barypout[4] __attribute__ ((aligned(16)));;
	float atte,rc,currweight,dlen,ww,Lp0;
	int i,j,k,tshift,*enb=NULL,*nextenb=NULL,eid;
	__m128 O,T,S;

	r->p0.w=1.f;
	r->vec.w=0.f;
	eid=r->eid-1;

	O = _mm_load_ps(&(r->p0.x));
	T = _mm_load_ps(&(r->vec.x));

	medium *prop;
	int *ee=(int *)(tracer->mesh->elem+eid);;
	prop=tracer->mesh->med+(tracer->mesh->type[eid]-1);
	atte=tracer->mesh->atte[tracer->mesh->type[eid]-1];
	rc=prop->n*R_C0;
        currweight=r->weight;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;
	for(i=0;i<4;i++)
	   if(havel_sse4(tracer->m+eid*12+i*3,&bary,O,T)){

		r->faceid=faceorder[i];
	   	dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
		Lp0=bary.x;
		r->isend=(Lp0>dlen);
		r->Lmove=((r->isend) ? dlen : Lp0);

                if(!r->isend){
		    enb=(int *)(tracer->mesh->facenb+eid);
                    r->nexteid=enb[r->faceid];
		    if(r->nexteid){
		        nextenb=(int *)(tracer->m+(r->nexteid-1)*12);
        		_mm_prefetch((char *)(nextenb),_MM_HINT_T0);
		        _mm_prefetch((char *)(nextenb+4*16),_MM_HINT_T0);
        		_mm_prefetch((char *)(nextenb+8*16),_MM_HINT_T0);
			nextenb=(int *)(tracer->mesh->elem+r->nexteid-1);
		    }
		}
	        S = _mm_set1_ps(bary.x);
	        S = _mm_mul_ps(S, T);
	        S = _mm_add_ps(S, O);
	        _mm_store_ps(&(r->pout.x),S);

		if(r->photontimer+r->Lmove*rc>=cfg->tend){ /*exit time window*/
		   r->faceid=-2;
	           r->pout.x=MMC_UNDEFINED;
		   r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
		}

		r->weight*=fast_expf9(-prop->mua*r->Lmove);
		r->slen-=r->Lmove*prop->mus;
	        if(bary.x==0.f){
			break;
		}
		ww=currweight-r->weight;
		r->Eabsorb+=ww;
		ww/=prop->mua;
                r->photontimer+=r->Lmove*rc;
		tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->nn-1;

                if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
                   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

	        T = _mm_mul_ps(T,_mm_set1_ps(r->Lmove));
	        S = _mm_add_ps(O, T);
	        _mm_store_ps(&(r->p0.x),S);

		barypout[nc[i][0]]=bary.y;
		barypout[nc[i][1]]=bary.z;
		barypout[nc[i][2]]=1.f-bary.y-bary.z;
		barypout[fm[i]]=0.f;

		//if(cfg->debuglevel&dlBary) 
		//    fprintf(cfg->flog,"barypout=[%f %f %f %f]\n",barypout[0],barypout[1],barypout[2],barypout[3]);

		T=_mm_load_ps(barypout);        /* bary centric at pout */
		O=_mm_load_ps(&(r->bary0.x)); /* bary centric at p0 */

		dlen=r->Lmove/bary.x;           /* normalized moving vector */

		//if(cfg->debuglevel&dlBary) 
		 //   fprintf(cfg->flog,"moved Lmove=%f from %f\n",r->Lmove, bary.x);

		if(r->isend)                    /* S is the bary centric for the moved photon */
		    S=_mm_add_ps(_mm_mul_ps(T,_mm_set1_ps(dlen)),_mm_mul_ps(O,_mm_set1_ps(1.f-dlen)));
		else
		    S=T;
		//if(cfg->debuglevel&dlBary) 
		//    fprintf(cfg->flog,"old bary0=[%f %f %f %f]\n",r->bary0.x,r->bary0.y,r->bary0.z,r->bary0.w);

		_mm_store_ps(barypout,S);
		if(nextenb && enb){
		    float *barynext=&(r->bary0.x);
		    _mm_store_ps(barynext,_mm_set1_ps(0.f));
		    for(j=0;j<3;j++)
		      for(k=0;k<4;k++){
		    	if(ee[nc[i][j]]==nextenb[k]){
		    	    barynext[k]=barypout[nc[i][j]];
			    break;
		    	}
		      }
		    //if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"[%d %d %d %d],[%d %d %d %d] - ",
		    //       ee[0],ee[1],ee[2],ee[3],nextenb[0],nextenb[1],nextenb[2],nextenb[3]);
		    //if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"[%f %f %f %f],[%f %f %f %f]\n",barypout[0],barypout[1],barypout[2],barypout[3],
		    //       barynext[0],barynext[1],barynext[2],barynext[3]);
		}else
		    _mm_store_ps(&(r->bary0.x),S);

		if(cfg->debuglevel&dlBary)
		    fprintf(cfg->flog,"new bary0=[%f %f %f %f]\n",r->bary0.x,r->bary0.y,r->bary0.z,r->bary0.w);
		T=_mm_mul_ps(_mm_add_ps(O,S),_mm_set1_ps(ww*0.5f));
		_mm_store_ps(barypout,T);

		for(j=0;j<4;j++)
		    tracer->mesh->weight[ee[j]+tshift]+=barypout[j];
		break;
	   }
	visit->raytet++;
	r->p0.w=0.f;
	if(r->faceid==-2)
           return 0.f;

	return r->slen;
}

float badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){

	float3 bary={1e10f,0.f,0.f,0.f};
	float Lp0=0.f,atte,rc,currweight,dlen,ww,t[4]={1e10f,1e10f,1e10f,1e10f};
	int i,tshift,faceidx=-1,eid;

	r->p0.w=1.f;
	r->vec.w=0.f;
	eid=r->eid-1;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;

	const __m128 o = _mm_load_ps(&(r->p0.x));
	const __m128 d = _mm_load_ps(&(r->vec.x));

	for(i=0;i<4;i++){
	    __m128 det,dett,n;
	    __m128 inv_det;
	    n=_mm_load_ps(&(tracer->m[3*((eid<<2)+i)].x));
	    det = _mm_dp_ps(n, d, 0x7f);
	    if(_mm_movemask_ps(det) & 1){
	  	 continue;
	    }
	    dett = _mm_dp_ps(_mm_mul_ps(int_coef, n), o, 0xff);
	    inv_det = rcp_nr(det);
	    _mm_store_ss(t+i, _mm_mul_ss(dett, inv_det));
	    if(t[i]<bary.x){
	       bary.x=t[i];
	       faceidx=i;
	       r->faceid=faceorder[i];
	    }
	}
/*	for(i=0;i<4;i++){
	    float det,dett;
	    float3 *n=tracer->n+eid*12+i*3;
	    det = vec_dot(n,r->vec);
	    if(det<0.f){
	  	 continue;
	    }
	    dett =n->w-vec_dot(n,r->p0);
	    t[i]=dett/det;
	    if(t[i]<bary.x){
	       bary.x=t[i];
	       faceidx=i;
	       r->faceid=faceorder[i];
	    }
	}*/
	if(r->faceid>=0){
	    medium *prop;
	    int *ee=(int *)(tracer->mesh->elem+eid);;
	    prop=tracer->mesh->med+(tracer->mesh->type[eid]-1);
	    atte=tracer->mesh->atte[tracer->mesh->type[eid]-1];
	    rc=prop->n*R_C0;
            currweight=r->weight;

	    dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
	    Lp0=bary.x;
	    r->isend=(Lp0>dlen);
	    r->Lmove=((r->isend) ? dlen : Lp0);

	    r->pout.x=r->p0.x+bary.x*r->vec.x;
	    r->pout.y=r->p0.y+bary.x*r->vec.y;
	    r->pout.z=r->p0.z+bary.x*r->vec.z;
	    if(r->photontimer+r->Lmove*rc>=cfg->tend){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
	    }

	    r->weight*=fast_expf9(-prop->mua*r->Lmove);
	    r->slen-=r->Lmove*prop->mus;
	    if(bary.x>=0.f){
		ww=currweight-r->weight;
		r->Eabsorb+=ww;
        	r->photontimer+=r->Lmove*rc;
		tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*
	             (cfg->basisorder?tracer->mesh->nn:tracer->mesh->ne);
        	if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
        	   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

        	r->p0.x+=r->Lmove*r->vec.x;
       		r->p0.y+=r->Lmove*r->vec.y;
        	r->p0.z+=r->Lmove*r->vec.z;
		if(!cfg->basisorder){
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
			ww/=prop->mua;
			tracer->mesh->weight[ee[nc[faceidx][0]]-1+tshift]+=ww*(1.f/3.f);
			tracer->mesh->weight[ee[nc[faceidx][1]]-1+tshift]+=ww*(1.f/3.f);
			tracer->mesh->weight[ee[nc[faceidx][2]]-1+tshift]+=ww*(1.f/3.f);
		}
	    }
	}
	visit->raytet++;
	r->p0.w=0.f;
	if(r->faceid==-2)
           return 0.f;

	return r->slen;
}

float branchless_badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){

	float3 bary={1e10f,0.f,0.f,0.f};
	float Lp0=0.f,atte,rc,currweight,dlen,ww;
	int tshift,faceidx=-1,baseid,eid;
	__m128 O,T,S;

	r->p0.w=1.f;
	r->vec.w=0.f;
	eid=r->eid-1;
	baseid=eid<<2;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;

	const __m128 Nx=_mm_load_ps(&(tracer->n[baseid].x));
	const __m128 Ny=_mm_load_ps(&(tracer->n[baseid+1].x));
	const __m128 Nz=_mm_load_ps(&(tracer->n[baseid+2].x));
	const __m128 dd=_mm_load_ps(&(tracer->n[baseid+3].x));

	O = _mm_set1_ps(r->p0.x);
	T = _mm_mul_ps(Nx,O);
	O = _mm_set1_ps(r->p0.y);
	T = _mm_add_ps(T, _mm_mul_ps(Ny,O));
	O = _mm_set1_ps(r->p0.z);
	T = _mm_add_ps(T, _mm_mul_ps(Nz,O));
	T = _mm_sub_ps(dd, T);

	O = _mm_set1_ps(r->vec.x);
	S = _mm_mul_ps(Nx,O);
	O = _mm_set1_ps(r->vec.y);
	S = _mm_add_ps(S, _mm_mul_ps(Ny,O));
	O = _mm_set1_ps(r->vec.z);
	S = _mm_add_ps(S, _mm_mul_ps(Nz,O));

	//T = _mm_mul_ps(T, rcp_nr(S));
	T = _mm_div_ps(T, S);

	O = _mm_cmpge_ps(S,_mm_set1_ps(0.f));
	T = _mm_add_ps(_mm_andnot_ps(O,_mm_set1_ps(1e10f)),_mm_and_ps(O,T));
	S = _mm_movehl_ps(T, T);
	O = _mm_min_ps(T, S);
	S = _mm_shuffle_ps(O, O, _MM_SHUFFLE(1,1,1,1));
	O = _mm_min_ss(O, S);
	
	_mm_store_ss(&(bary.x),O);
	faceidx=maskmap[_mm_movemask_ps(_mm_cmpeq_ps(T,_mm_set1_ps(bary.x)))];
	r->faceid=faceorder[faceidx];

	if(r->faceid>=0){
	    medium *prop;
	    int *enb, nexteid, *ee=(int *)(tracer->mesh->elem+eid);;
	    prop=tracer->mesh->med+(tracer->mesh->type[eid]-1);
	    atte=tracer->mesh->atte[tracer->mesh->type[eid]-1];
	    rc=prop->n*R_C0;
            currweight=r->weight;

            enb=(int *)(tracer->mesh->facenb+eid);
            nexteid=enb[r->faceid]; // if I use nexteid-1, the speed got slower, strange!
	    if(nexteid){
            	    _mm_prefetch((char *)&(tracer->n[nexteid<<2].x),_MM_HINT_T0);
	    }

	    dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
	    Lp0=bary.x;
	    r->isend=(Lp0>dlen);
	    r->Lmove=((r->isend) ? dlen : Lp0);

            O = _mm_load_ps(&(r->vec.x));
	    S = _mm_load_ps(&(r->p0.x));
	    T = _mm_set1_ps(bary.x);
	    T = _mm_mul_ps(O, T);
	    T = _mm_add_ps(T, S);
	    _mm_store_ps(&(r->pout.x),T);

	    if(r->photontimer+r->Lmove*rc>=cfg->tend){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
	    }

	    r->weight*=fast_expf9(-prop->mua*r->Lmove);
	    r->slen-=r->Lmove*prop->mus;
	    if(bary.x>=0.f){
		ww=currweight-r->weight;
		r->Eabsorb+=ww;
        	r->photontimer+=r->Lmove*rc;
		tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*
	             (cfg->basisorder?tracer->mesh->nn:tracer->mesh->ne);
        	if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %f\n",
        	   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

	        T = _mm_set1_ps(r->Lmove);
	        T = _mm_add_ps(S, _mm_mul_ps(O, T));
	        _mm_store_ps(&(r->p0.x),T);

		if(!cfg->basisorder){
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
			ww/=prop->mua;
			tracer->mesh->weight[ee[nc[faceidx][0]]-1+tshift]+=ww*(1.f/3.f);
			tracer->mesh->weight[ee[nc[faceidx][1]]-1+tshift]+=ww*(1.f/3.f);
			tracer->mesh->weight[ee[nc[faceidx][2]]-1+tshift]+=ww*(1.f/3.f);
		}
	    }
	}
	visit->raytet++;
	r->p0.w=0.f;
	if(r->faceid==-2)
           return 0.f;

	return r->slen;
}
#else

float havel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){
	mcx_error("wrong option, please recompile with SSE4 enabled");
}
float badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){
	mcx_error("wrong option, please recompile with SSE4 enabled");
}
float branchless_badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){
	mcx_error("wrong option, please recompile with SSE4 enabled");
}
#endif

float onephoton(int id,raytracer *tracer,tetmesh *mesh,mcconfig *cfg,
                RandType *ran, RandType *ran0, visitor *visit){

	int oldeid,fixcount=0;
	int *enb;
	const float int_coef_arr[4] = { -1.f, -1.f, -1.f, 1.f };
	float (*tracercore)(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit);

	ray r={cfg->srcpos,cfg->srcdir,{MMC_UNDEFINED,0.f,0.f},cfg->bary0,cfg->dim.x,-1,0,0,1.f,0.f,0.f,0.f};
	
	tracercore=havel_raytet;

	if(cfg->isplucker==1)
	    tracercore=plucker_raytet;
	else if(cfg->isplucker==2)
	    tracercore=badouel_raytet;
	else if(cfg->isplucker==3)
	    tracercore=branchless_badouel_raytet;

	/*initialize the photon parameters*/
	r.slen=rand_next_scatlen(ran);
	memcpy(&r.p0,&(cfg->srcpos),sizeof(r.p0));
	memcpy(&r.vec,&(cfg->srcdir),sizeof(r.p0));
#ifdef MMC_USE_SSE
	int_coef = _mm_load_ps(int_coef_arr);
#endif
	while(1){  /*propagate a photon until exit*/
	    r.slen=(*tracercore)(&r,tracer,cfg,visit);
	    if(r.pout.x==MMC_UNDEFINED){
	    	  if(r.faceid==-2) break; /*reaches the time limit*/
		  if(fixcount++<MAX_TRIAL){
			fixphoton(&r.p0,mesh->node,(int *)(mesh->elem+r.eid-1));
			continue;
                  }
	    	  r.eid=-r.eid;
        	  r.faceid=-1;
	    }
	    /*move a photon until the end of the current scattering path*/
	    while(r.faceid>=0 && !r.isend){
	    	    memcpy((void *)&r.p0,(void *)&r.pout,sizeof(r.p0));

	    	    enb=(int *)(mesh->facenb+r.eid-1);
		    oldeid=r.eid;
	    	    r.eid=enb[r.faceid];

		    if(cfg->isreflect && (r.eid==0 || mesh->med[mesh->type[r.eid-1]-1].n != mesh->med[mesh->type[oldeid-1]-1].n )){
			if(! (r.eid==0 && mesh->med[mesh->type[oldeid-1]-1].n == cfg->nout ))
			    r.weight*=reflectray(cfg,&r.vec,tracer,&oldeid,&r.eid,r.faceid,ran);
		    }

	    	    if(!cfg->isreflect && r.eid==0) break;
		    if(r.eid==0 && mesh->med[mesh->type[oldeid-1]-1].n == cfg->nout ) break;
	    	    if(r.pout.x!=MMC_UNDEFINED && (cfg->debuglevel&dlMove))
	    		fprintf(cfg->flog,"P %f %f %f %d %d %f\n",r.pout.x,r.pout.y,r.pout.z,r.eid,id,r.slen);

	    	    r.slen=(*tracercore)(&r,tracer,cfg,visit);
		    if(r.faceid==-2) break;
		    fixcount=0;
		    while(r.pout.x==MMC_UNDEFINED && fixcount++<MAX_TRIAL){
		       fixphoton(&r.p0,mesh->node,(int *)(mesh->elem+r.eid-1));
                       r.slen=(*tracercore)(&r,tracer,cfg,visit);
		    }
        	    if(r.pout.x==MMC_UNDEFINED){
        		/*possibily hit an edge or miss*/
			r.eid=-r.eid;
        		break;
        	    }
	    }
	    if(r.eid<=0 || r.pout.x==MMC_UNDEFINED) {
        	    if(r.eid==0 && (cfg->debuglevel&dlMove))
        		 fprintf(cfg->flog,"B %f %f %f %d %d %f\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
		    else if(r.faceid==-2 && (cfg->debuglevel&dlMove))
                         fprintf(cfg->flog,"T %f %f %f %d %d %f\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    	    else if(r.eid && r.faceid!=-2  && cfg->debuglevel&dlEdge)
        		 fprintf(cfg->flog,"X %f %f %f %d %d %f\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    	    break;  /*photon exits boundary*/
	    }
	    if(cfg->debuglevel&dlMove) fprintf(cfg->flog,"M %f %f %f %d %d %f\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    if(cfg->minenergy>0.f && r.weight < cfg->minenergy && (cfg->tend-cfg->tstart)*visit->rtstep<=1.f){ /*Russian Roulette*/
		if(rand_do_roulette(ran)*cfg->roulettesize<=1.f)
			r.weight*=cfg->roulettesize;
                        if(cfg->debuglevel&dlWeight) fprintf(cfg->flog,"Russian Roulette bumps r.weight to %f\n",r.weight);
		else
			break;
	    }
	    r.slen=mc_next_scatter(mesh->med[mesh->type[r.eid-1]-1].g,&r.vec,ran,ran0,cfg);
	}
	return r.Eabsorb;
}
inline float mmc_rsqrtf(float a){
#ifdef MMC_USE_SSE
        return _mm_cvtss_f32( _mm_rsqrt_ss( _mm_set_ss( a ) ) );
#else
	return 1.f/sqrtf(a);
#endif
}
float reflectray(mcconfig *cfg,float3 *c0,raytracer *tracer,int *oldeid,int *eid,int faceid,RandType *ran){
	/*to handle refractive index mismatch*/
        float3 pnorm;
	float Icos,Re,Im,Rtotal,tmp0,tmp1,tmp2,n1,n2;

	/*calculate the normal direction of the intersecting triangle*/
        vec_cross(tracer->d+(*oldeid-1)*6+fc[ifaceorder[faceid]][0],
                  tracer->d+(*oldeid-1)*6+fc[ifaceorder[faceid]][1],&pnorm);
	tmp0=mmc_rsqrtf(vec_dot(&pnorm,&pnorm));
	tmp0*=(ifaceorder[faceid]>=2) ? -1.f : 1.f;
	pnorm.x*=tmp0;
        pnorm.y*=tmp0;
        pnorm.z*=tmp0;
	
	/*compute the cos of the incidence angle*/
        Icos=fabs(vec_dot(c0,&pnorm));

	n1=tracer->mesh->med[tracer->mesh->type[*oldeid-1]-1].n;
	n2=(*eid>0) ? tracer->mesh->med[tracer->mesh->type[*eid-1]-1].n : cfg->nout;

	tmp0=n1*n1;
	tmp1=n2*n2;
        tmp2=1.f-tmp0/tmp1*(1.f-Icos*Icos); /*1-[n1/n2*sin(si)]^2 = cos(ti)^2*/

        if(tmp2>0.f){ /*if no total internal reflection*/
          Re=tmp0*Icos*Icos+tmp1*tmp2;      /*transmission angle*/
	  tmp2=sqrt(tmp2); /*to save one sqrt*/
          Im=2.f*n1*n2*Icos*tmp2;
          Rtotal=(Re-Im)/(Re+Im);     /*Rp*/
          Re=tmp1*Icos*Icos+tmp0*tmp2*tmp2;
          Rtotal=(Rtotal+(Re-Im)/(Re+Im))*0.5f; /*(Rp+Rs)/2*/
	  if(*eid==0 || rand_next_reflect(ran)<=Rtotal){ /*do reflection*/
              vec_mult_add(&pnorm,c0,2.f*Icos,1.f,c0);
              if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal);
	      return (*eid==0 ? *eid=*oldeid, Rtotal : 1.f); /*if reflect at boundary, use ref coeff, otherwise, use probability*/
          }else{                              /*do transmission*/
              vec_mult_add(&pnorm,c0, Icos,1.f,c0);
              vec_mult_add(&pnorm,c0,-tmp2,n1/n2,c0);
              if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal);
              return 1.f;
	  }
       }else{ /*total internal reflection*/
          vec_mult_add(&pnorm,c0,2.f*Icos,1.f,c0);
	  *eid=*oldeid;
          if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f);
          return 1.f;
       }
}
