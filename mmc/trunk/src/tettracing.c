/*******************************************************************************
**  Mesh-based Monte Carlo (MMC)
**
**  Author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
**
**  Reference:
**  (Fang2010) Qianqian Fang, "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Pluker Coordinates," Biomed. Opt. Express, 1(1) 165-175 (2010)
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
const int nc[4][3]={{3,0,1},{3,1,2},{2,0,3},{1,0,2}}; /*node orders for each face, in clock-wise orders*/
extern const int out[4][3]; /*defined in simpmesh.c, node orders for each face, in counter-clock-wise orders*/
extern const int facemap[4];
extern const int ifacemap[4];

const int faceorder[]={1,3,2,0,-1};
const int ifaceorder[]={3,0,2,1};
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
	medium *prop;
	int *ee;
	int i,tshift,eid,faceidx=-1;
	float w[6],Rv,ww,currweight,dlen=0.f,rc; /*dlen is the physical distance*/
	unsigned int *wi=(unsigned int*)w;
        float baryout[4]={0.f,0.f,0.f,0.f},*baryp0=&(r->bary0.x);
	float Lp0=0.f,ratio;

	if(tracer->mesh==NULL || tracer->d==NULL||r->eid<=0||r->eid>tracer->mesh->ne) 
		return -1;

	eid=r->eid-1;
	r->faceid=-1;
	r->isend=0;
	r->Lmove=0.f;
	vec_add(&(r->p0),&(r->vec),&p1);
	vec_cross(&(r->p0),&p1,&pcrx);
	ee=(int *)(tracer->mesh->elem+eid);
	prop=tracer->mesh->med+(tracer->mesh->type[eid]);
	rc=prop->n*R_C0;
        currweight=r->weight;

#ifdef MMC_USE_SSE
	{
	    __m128 O=_mm_load_ps(&(r->vec.x));
	    __m128 M,D,T=_mm_load_ps(&(pcrx.x));
	    for(i=0;i<6;i++){
	        D=_mm_load_ps(&(tracer->d[(eid)*6+i].x));
	        M=_mm_load_ps(&(tracer->m[(eid)*6+i].x));
		w[i]=_mm_cvtss_f32(_mm_add_ss(_mm_dp_ps(O,M,0x7F),_mm_dp_ps(T,D,0x7F)));
	    }
	}
#else
	for(i=0;i<6;i++)
		w[i]=pinner(&(r->vec),&pcrx,tracer->d+(eid)*6+i,tracer->m+(eid)*6+i);
#endif
	//exit(1);
	if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"%d \n",eid);

	r->pout.x=MMC_UNDEFINED;
	for(i=0;i<4;i++){
		//if(cfg->debuglevel&dlTracing) fprintf(cfg->flog,"testing face [%d]\n",i);

	        if(i>=2) w[fc[i][1]]=-w[fc[i][1]]; // can not go back
		if(F32N(wi[fc[i][0]]) & F32N(wi[fc[i][1]]) & F32P(wi[fc[i][2]])){
			// f_leave
                        //if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"ray exits face %d[%d] of %d\n",i,faceorder[i],eid);

                        Rv=1.f/(-w[fc[i][0]]-w[fc[i][1]]+w[fc[i][2]]);
                        baryout[nc[i][0]]=-w[fc[i][0]]*Rv;
                        baryout[nc[i][1]]=-w[fc[i][1]]*Rv;
                        baryout[nc[i][2]]=w[fc[i][2]]*Rv;

			getinterp(baryout[nc[i][0]],baryout[nc[i][1]],baryout[nc[i][2]],
				tracer->mesh->node+ee[nc[i][0]]-1,
				tracer->mesh->node+ee[nc[i][1]]-1,
				tracer->mesh->node+ee[nc[i][2]]-1,&(r->pout));

                        //if(cfg->debuglevel&dlTracingExit) fprintf(cfg->flog,"exit point %f %f %f\n",r->pout.x,r->pout.y,r->pout.z);

			Lp0=dist(&(r->p0),&(r->pout));
			dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
			faceidx=i;
			r->faceid=faceorder[i];

		       {
			int *enb=(int *)(tracer->mesh->facenb+eid);
			r->nexteid=enb[r->faceid];
#ifdef MMC_USE_SSE
			int nexteid=(r->nexteid-1)*6;
	                if(r->nexteid){
        	                _mm_prefetch((char *)&((tracer->m+nexteid)->x),_MM_HINT_T0);
                	        _mm_prefetch((char *)&((tracer->m+(nexteid)*6+4)->x),_MM_HINT_T0);
                        	_mm_prefetch((char *)&((tracer->d+(nexteid)*6)->x),_MM_HINT_T0);
                                _mm_prefetch((char *)&((tracer->d+(nexteid)*6+4)->x),_MM_HINT_T0);
	                }
#endif
		       }
			r->isend=(Lp0>dlen);
			r->Lmove=((r->isend) ? dlen : Lp0);
			break;
		}
	}
	visit->raytet++;
        if(r->pout.x!=MMC_UNDEFINED){
                if((int)((r->photontimer+r->Lmove*rc-cfg->tstart)*visit->rtstep)>=(int)((cfg->tend-cfg->tstart)*visit->rtstep)){ /*exit time window*/
		   r->faceid=-2;
	           r->pout.x=MMC_UNDEFINED;
		   r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
		}
		if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otJacobian){
#ifdef __INTEL_COMPILER
		    currweight=expf(-DELTA_MUA*r->Lmove);
#else
		    currweight=fast_expf9(-DELTA_MUA*r->Lmove);
#endif
                    currweight*=cfg->replayweight[r->photoid];
                    r->weight=0.f;
		}else if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otTaylor){
		    currweight=r->Lmove;
		    currweight*=cfg->replayweight[r->photoid];
		    r->weight=0.f;
                }else{
#ifdef __INTEL_COMPILER
		    r->weight*=expf(-prop->mua*r->Lmove);
#else
		    r->weight*=fast_expf9(-prop->mua*r->Lmove);
#endif
		}
		r->slen-=r->Lmove*prop->mus;
                r->p0.x+=r->Lmove*r->vec.x;
                r->p0.y+=r->Lmove*r->vec.y;
                r->p0.z+=r->Lmove*r->vec.z;
                //if(cfg->debuglevel&dlWeight) fprintf(cfg->flog,"update weight to %f and path end %d \n",r->weight,r->isend);

		if(!cfg->basisorder){
                        ww=currweight-r->weight;
                        r->Eabsorb+=ww;
                        r->photontimer+=r->Lmove*rc;
                        tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->ne;
#pragma omp atomic
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
	                if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"Y [%f %f %f %f]\n",
        	              baryout[0],baryout[1],baryout[2],baryout[3]);

                	if(Lp0>EPS){
                        	r->photontimer+=r->Lmove*rc;
        	        	ww=currweight-r->weight;
				if(prop->mua>0.f){
				  ratio=r->Lmove/Lp0;
				  r->Eabsorb+=ww;
				  if(cfg->outputtype!=otEnergy) ww/=prop->mua;
                        	  tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->nn;

                        	  if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %e\n",
                        	     r->p0.x-(r->Lmove*0.5f)*r->vec.x,r->p0.y-(r->Lmove*0.5f)*r->vec.y,r->p0.z-(r->Lmove*0.5f)*r->vec.z,ww,eid,dlen);

				  ww*=0.5f;
				  if(r->isend)
                  		    for(i=0;i<4;i++)
                     			baryout[i]=(1.f-ratio)*baryp0[i]+ratio*baryout[i];
                		  for(i=0;i<4;i++)
#pragma omp atomic
                     			tracer->mesh->weight[ee[i]-1+tshift]+=ww*(baryp0[i]+baryout[i]);
				}
				if(r->isend){
					memcpy(baryp0,baryout,sizeof(float4));
				}else{
		        		if(r->nexteid && faceidx>=0){
		        		    int j,k,*nextenb=(int *)(tracer->mesh->elem+r->nexteid-1);
					    memset(baryp0,0,sizeof(float4));
					    for(j=0;j<3;j++)
					      for(k=0;k<4;k++){
		    				if(ee[nc[faceidx][j]]==nextenb[k]){
		    				    baryp0[k]=baryout[nc[faceidx][j]];
						    break;
		    				}
					      }
					}
				}
			}
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
	float barypout[4] __attribute__ ((aligned(16)));
	float rc,currweight,dlen,ww,Lp0;
	int i,j,k,tshift,*enb=NULL,*nextenb=NULL,eid;
	__m128 O,T,S;

	r->p0.w=1.f;
	r->vec.w=0.f;
	eid=r->eid-1;

	O = _mm_load_ps(&(r->p0.x));
	T = _mm_load_ps(&(r->vec.x));

	medium *prop;
	int *ee=(int *)(tracer->mesh->elem+eid);
	prop=tracer->mesh->med+(tracer->mesh->type[eid]);
	rc=prop->n*R_C0;
        currweight=r->weight;

	r->pout.x=MMC_UNDEFINED;
	r->faceid=-1;
	r->isend=0;
	r->Lmove=0.f;
	
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

		if((int)((r->photontimer+r->Lmove*rc-cfg->tstart)*visit->rtstep)>=(int)((cfg->tend-cfg->tstart)*visit->rtstep)){ /*exit time window*/
		   r->faceid=-2;
	           r->pout.x=MMC_UNDEFINED;
		   r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
		}
		if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otJacobian){
#ifdef __INTEL_COMPILER
		    currweight=expf(-DELTA_MUA*r->Lmove);
#else
		    currweight=fast_expf9(-DELTA_MUA*r->Lmove);
#endif
                    currweight*=cfg->replayweight[r->photoid];
                    r->weight=0.f;
                }else if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otTaylor){
		    currweight=r->Lmove;
		    currweight*=cfg->replayweight[r->photoid];
		    r->weight=0.f;
                }else{
#ifdef __INTEL_COMPILER
		    r->weight*=expf(-prop->mua*r->Lmove);
#else
		    r->weight*=fast_expf9(-prop->mua*r->Lmove);
#endif
		}
		r->slen-=r->Lmove*prop->mus;
	        if(bary.x==0.f){
			break;
		}
		ww=currweight-r->weight;
                r->photontimer+=r->Lmove*rc;

		if(prop->mua>0.f){
		  r->Eabsorb+=ww;
		  if(cfg->outputtype!=otEnergy) ww/=prop->mua;
		}
                tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*tracer->mesh->nn-1;

                if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %e\n",
                   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

	        T = _mm_mul_ps(T,_mm_set1_ps(r->Lmove));
	        S = _mm_add_ps(O, T);
	        _mm_store_ps(&(r->p0.x),S);

		barypout[out[i][0]]=1.f-bary.y-bary.z;
		barypout[out[i][1]]=bary.y;
		barypout[out[i][2]]=bary.z;
		barypout[facemap[i]]=0.f;
	        if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"Y [%f %f %f %f]\n",
        	      barypout[0],barypout[1],barypout[2],barypout[3]);

		T=_mm_load_ps(barypout);        /* bary centric at pout */
		O=_mm_load_ps(&(r->bary0.x));   /* bary centric at p0 */

		dlen=r->Lmove/bary.x;           /* normalized moving length */

                if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %e\n",
                    r->p0.x-(r->Lmove*0.5f)*r->vec.x,r->p0.y-(r->Lmove*0.5f)*r->vec.y,r->p0.z-(r->Lmove*0.5f)*r->vec.z,ww,eid,dlen);

		if(r->isend)                    /* S is the bary centric for the photon after move */
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
		    	if(ee[out[i][j]]==nextenb[k]){
		    	    barynext[k]=barypout[out[i][j]];
			    break;
		    	}
		      }
		    //if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"[%d %d %d %d],[%d %d %d %d] - ",
		    //       ee[0],ee[1],ee[2],ee[3],nextenb[0],nextenb[1],nextenb[2],nextenb[3]);
		    //if(cfg->debuglevel&dlBary) fprintf(cfg->flog,"[%f %f %f %f],[%f %f %f %f]\n",barypout[0],barypout[1],barypout[2],barypout[3],
		    //       barynext[0],barynext[1],barynext[2],barynext[3]);
		}else
		    _mm_store_ps(&(r->bary0.x),S);

		//if(cfg->debuglevel&dlBary)
		//    fprintf(cfg->flog,"new bary0=[%f %f %f %f]\n",r->bary0.x,r->bary0.y,r->bary0.z,r->bary0.w);
		if(prop->mua>0.f){
		  T=_mm_mul_ps(_mm_add_ps(O,S),_mm_set1_ps(ww*0.5f));
		  _mm_store_ps(barypout,T);

		  for(j=0;j<4;j++)
#pragma omp atomic
		    tracer->mesh->weight[ee[j]+tshift]+=barypout[j];

		}
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
	float Lp0=0.f,rc,currweight,dlen,ww,t[4]={1e10f,1e10f,1e10f,1e10f};
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
	    prop=tracer->mesh->med+(tracer->mesh->type[eid]);
	    rc=prop->n*R_C0;
            currweight=r->weight;

	    dlen=(prop->mus <= EPS) ? R_MIN_MUS : r->slen/prop->mus;
	    Lp0=bary.x;
	    r->isend=(Lp0>dlen);
	    r->Lmove=((r->isend) ? dlen : Lp0);

	    r->pout.x=r->p0.x+bary.x*r->vec.x;
	    r->pout.y=r->p0.y+bary.x*r->vec.y;
	    r->pout.z=r->p0.z+bary.x*r->vec.z;
	    if((int)((r->photontimer+r->Lmove*rc-cfg->tstart)*visit->rtstep)>=(int)((cfg->tend-cfg->tstart)*visit->rtstep)){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
	    }
	    if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otJacobian){
#ifdef __INTEL_COMPILER
		currweight=expf(-DELTA_MUA*r->Lmove);
#else
		currweight=fast_expf9(-DELTA_MUA*r->Lmove);
#endif
                currweight*=cfg->replayweight[r->photoid];
                r->weight=0.f;
            }else if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otTaylor){
		currweight=r->Lmove;
		currweight*=cfg->replayweight[r->photoid];
		r->weight=0.f;
            }else{
#ifdef __INTEL_COMPILER
		r->weight*=expf(-prop->mua*r->Lmove);
#else
		r->weight*=fast_expf9(-prop->mua*r->Lmove);
#endif
	    }
	    r->slen-=r->Lmove*prop->mus;
	    if(bary.x>=0.f){
		ww=currweight-r->weight;
		r->Eabsorb+=ww;
        	r->photontimer+=r->Lmove*rc;
		tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*
	             (cfg->basisorder?tracer->mesh->nn:tracer->mesh->ne);
        	if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %e\n",
        	   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

        	r->p0.x+=r->Lmove*r->vec.x;
       		r->p0.y+=r->Lmove*r->vec.y;
        	r->p0.z+=r->Lmove*r->vec.z;
		if(!cfg->basisorder){
#pragma omp atomic
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
			if(cfg->outputtype!=otEnergy) ww/=prop->mua;
                        ww*=1.f/3.f;
			for(i=0;i<3;i++)
#pragma omp atomic
				tracer->mesh->weight[ee[out[faceidx][i]]-1+tshift]+=ww;
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
	float Lp0=0.f,rc,currweight,dlen,ww;
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
	    int *enb, *ee=(int *)(tracer->mesh->elem+eid);;
	    prop=tracer->mesh->med+(tracer->mesh->type[eid]);
	    rc=prop->n*R_C0;
            currweight=r->weight;

            enb=(int *)(tracer->mesh->facenb+eid);
            r->nexteid=enb[r->faceid]; // if I use nexteid-1, the speed got slower, strange!
	    if(r->nexteid){
            	    _mm_prefetch((char *)&(tracer->n[(r->nexteid-1)<<2].x),_MM_HINT_T0);
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

	    if((int)((r->photontimer+r->Lmove*rc-cfg->tstart)*visit->rtstep)>=(int)((cfg->tend-cfg->tstart)*visit->rtstep)){ /*exit time window*/
	       r->faceid=-2;
	       r->pout.x=MMC_UNDEFINED;
	       r->Lmove=(cfg->tend-r->photontimer)/(prop->n*R_C0)-1e-4f;
	    }
	    if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otJacobian){
#ifdef __INTEL_COMPILER
		currweight=expf(-DELTA_MUA*r->Lmove);
#else
		currweight=fast_expf9(-DELTA_MUA*r->Lmove);
#endif
                currweight*=cfg->replayweight[r->photoid];
                r->weight=0.f;
            }else if(cfg->seed==SEED_FROM_FILE && cfg->outputtype==otTaylor){
		currweight=r->Lmove;
		currweight*=cfg->replayweight[r->photoid];
		r->weight=0.f;
            }else{
#ifdef __INTEL_COMPILER
		r->weight*=expf(-prop->mua*r->Lmove);
#else
		r->weight*=fast_expf9(-prop->mua*r->Lmove);
#endif
	    }
	    r->slen-=r->Lmove*prop->mus;
	    if(bary.x>=0.f){
		ww=currweight-r->weight;
		r->Eabsorb+=ww;
        	r->photontimer+=r->Lmove*rc;
		tshift=(int)((r->photontimer-cfg->tstart)*visit->rtstep)*
	             (cfg->basisorder?tracer->mesh->nn:tracer->mesh->ne);
        	if(cfg->debuglevel&dlAccum) fprintf(cfg->flog,"A %f %f %f %e %d %e\n",
        	   r->p0.x,r->p0.y,r->p0.z,bary.x,eid+1,dlen);

	        T = _mm_set1_ps(r->Lmove);
	        T = _mm_add_ps(S, _mm_mul_ps(O, T));
	        _mm_store_ps(&(r->p0.x),T);
		if(!cfg->basisorder){
#pragma omp atomic
			tracer->mesh->weight[eid+tshift]+=ww;
		}else{
			int i;
			if(cfg->outputtype!=otEnergy) ww/=prop->mua;
                        ww*=1.f/3.f;
			for(i=0;i<3;i++)
#pragma omp atomic
				tracer->mesh->weight[ee[out[faceidx][i]]-1+tshift]+=ww;
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
	MMC_ERROR(-6,"wrong option, please recompile with SSE4 enabled");
	return MMC_UNDEFINED;
}
float badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){
	MMC_ERROR(-6,"wrong option, please recompile with SSE4 enabled");
	return MMC_UNDEFINED;
}
float branchless_badouel_raytet(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit){
	MMC_ERROR(-6,"wrong option, please recompile with SSE4 enabled");
	return MMC_UNDEFINED;
}
#endif

float onephoton(unsigned int id,raytracer *tracer,tetmesh *mesh,mcconfig *cfg,
                RandType *ran, RandType *ran0, visitor *visit){

	int oldeid,fixcount=0,exitdet=0;
	int *enb;
        float mom;
	ray r={cfg->srcpos,cfg->srcdir,{MMC_UNDEFINED,0.f,0.f},cfg->bary0,cfg->dim.x,cfg->dim.y-1,0,0,1.f,0.f,0.f,0.f,0.,0,NULL,NULL};

	float (*engines[4])(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit)=
	       {plucker_raytet,havel_raytet,badouel_raytet,branchless_badouel_raytet};
	float (*tracercore)(ray *r, raytracer *tracer, mcconfig *cfg, visitor *visit);

	r.partialpath=(float*)calloc(visit->reclen-1,sizeof(float));
	r.photoid=id;
	r.partialpath[visit->reclen-2] = r.weight;
	visit->accumu_weight += r.weight;

        if(cfg->issavedet && cfg->issaveseed){
                r.photonseed=(void*)calloc(1,sizeof(RandType));
		memcpy(r.photonseed,(void *)ran, sizeof(RandType));
        }
	tracercore=engines[0];
	if(cfg->method>=0 && cfg->method<4)
	    tracercore=engines[(int)(cfg->method)];
	else
	    MMC_ERROR(-6,"specified ray-tracing algorithm is not defined");

	/*initialize the photon parameters*/
        launchphoton(cfg, &r, mesh, ran, ran0);

/*
	if(cfg->srctype==stIsotropic){
                mom=0.f;
		r.slen=mc_next_scatter(0,&r.vec,ran,ran0,cfg,&mom);
                if(cfg->ismomentum)
                   r.partialpath[mesh->prop+mesh->type[r.eid-1]]+=mom;
	}else{
		r.slen=rand_next_scatlen(ran);
	}

	if(cfg->isspecular && r.faceid>=0 && mesh->med[mesh->type[r.eid-1]].n != cfg->nout && cfg->srctype==stPencil){
	    float Rspecular=reflectray(cfg,&r.vec,tracer,&r.eid,&r.eid,faceorder[r.faceid],ran);
	    if(Rspecular<1.f)
	       r.weight*=(1.f-Rspecular);
	    else
	       return 0.f;
	}
*/
#ifdef MMC_USE_SSE
	const float int_coef_arr[4] = { -1.f, -1.f, -1.f, 1.f };
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
	    if(cfg->issavedet && r.Lmove>0.f && mesh->type[r.eid-1]>0)
	            r.partialpath[mesh->type[r.eid-1]]+=r.Lmove;
	    /*move a photon until the end of the current scattering path*/
	    while(r.faceid>=0 && !r.isend){
	    	    memcpy((void *)&r.p0,(void *)&r.pout,sizeof(r.p0));

	    	    enb=(int *)(mesh->facenb+r.eid-1);
		    oldeid=r.eid;
	    	    r.eid=enb[r.faceid];

		    if(cfg->isreflect && (r.eid==0 || mesh->med[mesh->type[r.eid-1]].n != mesh->med[mesh->type[oldeid-1]].n )){
			if(! (r.eid==0 && mesh->med[mesh->type[oldeid-1]].n == cfg->nout ))
			    reflectray(cfg,&r.vec,tracer,&oldeid,&r.eid,r.faceid,ran);
		    }
	    	    if(r.eid==0) break;
//		    if(r.eid==0 && mesh->med[mesh->type[oldeid-1]].n == cfg->nout ) break;
	    	    if(r.pout.x!=MMC_UNDEFINED && (cfg->debuglevel&dlMove))
	    		fprintf(cfg->flog,"P %f %f %f %d %u %e\n",r.pout.x,r.pout.y,r.pout.z,r.eid,id,r.slen);

	    	    r.slen=(*tracercore)(&r,tracer,cfg,visit);
		    if(cfg->issavedet && r.Lmove>0.f && mesh->type[r.eid-1]>0)
			r.partialpath[mesh->type[r.eid-1]]+=r.Lmove;
		    if(r.faceid==-2) break;
		    fixcount=0;
		    while(r.pout.x==MMC_UNDEFINED && fixcount++<MAX_TRIAL){
		       fixphoton(&r.p0,mesh->node,(int *)(mesh->elem+r.eid-1));
                       r.slen=(*tracercore)(&r,tracer,cfg,visit);
		       if(cfg->issavedet && r.Lmove>0.f && mesh->type[r.eid-1]>0)
	            		r.partialpath[mesh->type[r.eid-1]]+=r.Lmove;
		    }
        	    if(r.pout.x==MMC_UNDEFINED){
        		/*possibily hit an edge or miss*/
			r.eid=-r.eid;
        		break;
        	    }
	    }
	    if(r.eid<=0 || r.pout.x==MMC_UNDEFINED) {
        	    if(r.eid==0 && (cfg->debuglevel&dlMove))
        		 fprintf(cfg->flog,"B %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
		    else if(r.eid==0){
                       if(cfg->debuglevel&dlExit)
        		 fprintf(cfg->flog,"E %f %f %f %f %f %f %f %d\n",r.p0.x,r.p0.y,r.p0.z,
			    r.vec.x,r.vec.y,r.vec.z,r.weight,r.eid);
                       if(cfg->issavedet && cfg->issaveexit){
                            memcpy(r.partialpath+(visit->reclen-1-6),&(r.p0.x),sizeof(float)*3);
                            memcpy(r.partialpath+(visit->reclen-1-3),&(r.vec.x),sizeof(float)*3);
                       }
		    }else if(r.faceid==-2 && (cfg->debuglevel&dlMove))
                         fprintf(cfg->flog,"T %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    	    else if(r.eid && r.faceid!=-2  && cfg->debuglevel&dlEdge)
        		 fprintf(cfg->flog,"X %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
		    if(cfg->issavedet && r.eid==0){
		       int i;
		       for(i=0;i<cfg->detnum;i++){
        		  if((cfg->detpos[i].x-r.p0.x)*(cfg->detpos[i].x-r.p0.x)+
        		     (cfg->detpos[i].y-r.p0.y)*(cfg->detpos[i].y-r.p0.y)+
        		     (cfg->detpos[i].z-r.p0.z)*(cfg->detpos[i].z-r.p0.z) < cfg->detpos[i].w*cfg->detpos[i].w){
			          exitdet=i+1;
                		  break;
        		     }
		       }
		    }
	    	    break;  /*photon exits boundary*/
	    }
	    if(cfg->debuglevel&dlMove) fprintf(cfg->flog,"M %f %f %f %d %u %e\n",r.p0.x,r.p0.y,r.p0.z,r.eid,id,r.slen);
	    if(cfg->minenergy>0.f && r.weight < cfg->minenergy && (cfg->tend-cfg->tstart)*visit->rtstep<=1.f){ /*Russian Roulette*/
		if(rand_do_roulette(ran)*cfg->roulettesize<=1.f)
			r.weight*=cfg->roulettesize;
                        if(cfg->debuglevel&dlWeight)
			    fprintf(cfg->flog,"Russian Roulette bumps r.weight to %f\n",r.weight);
		else
			break;
	    }
            mom=0.f;
	    r.slen=mc_next_scatter(mesh->med[mesh->type[r.eid-1]].g,&r.vec,ran,ran0,cfg,&mom);
            if(cfg->ismomentum && mesh->type[r.eid-1]>0)
                  r.partialpath[mesh->prop+mesh->type[r.eid-1]]+=mom;
            r.partialpath[0]++;
	}
	if(cfg->issavedet && exitdet>0){
		int offset=visit->bufpos*visit->reclen;
		if(visit->bufpos>=visit->detcount){
		    visit->detcount+=DET_PHOTON_BUF;
		    visit->partialpath=(float *)realloc(visit->partialpath,
				visit->detcount*visit->reclen*sizeof(float));
		    if(cfg->issaveseed)
	                    visit->photonseed=realloc(visit->photonseed,visit->detcount*sizeof(RandType));
		}
		visit->partialpath[offset]=exitdet;
	        memcpy(visit->partialpath+offset+1,r.partialpath,(visit->reclen-1)*sizeof(float));
                if(cfg->issaveseed)
	            memcpy(visit->photonseed+visit->bufpos*sizeof(RandType),r.photonseed,sizeof(RandType));
		visit->bufpos++;
	}
	free(r.partialpath);
        if(r.photonseed)
		free(r.photonseed);
	return r.Eabsorb;
}

float reflectray(mcconfig *cfg,float3 *c0,raytracer *tracer,int *oldeid,int *eid,int faceid,RandType *ran){
	/*to handle refractive index mismatch*/
        float3 pnorm={0.f}, *pn=&pnorm;
	float Icos,Re,Im,Rtotal,tmp0,tmp1,tmp2,n1,n2;
	int offs=(*oldeid-1)<<2;

	faceid=ifaceorder[faceid];
	/*calculate the normal direction of the intersecting triangle*/
	if(cfg->method==rtPlucker) { //Plucker ray-tracing
		pn=tracer->n+(offs)+faceid;
	}else if(cfg->method<rtBLBadouel){
		pn=tracer->m+(offs+faceid)*3;
	}else if(cfg->method==rtBLBadouel){
		pnorm.x=(&(tracer->n[offs].x))[faceid];
		pnorm.y=(&(tracer->n[offs].x))[faceid+4];
		pnorm.z=(&(tracer->n[offs].x))[faceid+8];
	}
	/*pn pointing outward*/

	/*compute the cos of the incidence angle*/
        Icos=fabs(vec_dot(c0,pn));

	n1=(*oldeid!=*eid) ? tracer->mesh->med[tracer->mesh->type[*oldeid-1]].n : cfg->nout;
	n2=(*eid>0) ? tracer->mesh->med[tracer->mesh->type[*eid-1]].n : cfg->nout;

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
	  if(*oldeid==*eid) return Rtotal; /*initial specular reflection*/
	  if(rand_next_reflect(ran)<=Rtotal){ /*do reflection*/
              vec_mult_add(pn,c0,-2.f*Icos,1.f,c0);
              //if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"R %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,Rtotal);
	      *eid=*oldeid; /*stay with the current element*/
	      return 1.f;
	  }else if(cfg->isspecular==2 && *eid==0){   /*if do transmission, but next neighbor is 0, terminate*/
	      return 1.f;
          }else{                              /*do transmission*/
              vec_mult_add(pn,c0,-Icos,1.f,c0);
              vec_mult_add(pn,c0,tmp2,n1/n2,c0);
              //if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"Z %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f-Rtotal);
              return 1.f;
	  }
       }else{ /*total internal reflection*/
          vec_mult_add(pn,c0,-2.f*Icos,1.f,c0);
	  *eid=*oldeid;
          //if(cfg->debuglevel&dlReflect) fprintf(cfg->flog,"V %f %f %f %d %d %f\n",c0->x,c0->y,c0->z,*eid,*oldeid,1.f);
          return 1.f;
       }
}

void launchphoton(mcconfig *cfg, ray *r, tetmesh *mesh, RandType *ran, RandType *ran0){
	/*a rectangular grid over a plane*/
	if(cfg->srctype==stPencil){ // pencil beam, use the old workflow
		r->slen=rand_next_scatlen(ran);
		return;
	}else if(cfg->srctype==stPlanar || cfg->srctype==stPattern || cfg->srctype==stFourier){
		float rx=rand_uniform01(ran);
		float ry=rand_uniform01(ran);
		r->p0.x=cfg->srcpos.x+rx*cfg->srcparam1.x+ry*cfg->srcparam2.x;
		r->p0.y=cfg->srcpos.y+rx*cfg->srcparam1.y+ry*cfg->srcparam2.y;
		r->p0.z=cfg->srcpos.z+rx*cfg->srcparam1.z+ry*cfg->srcparam2.z;
		r->weight=1.f;
		if(cfg->srctype==stPattern){
			r->weight=cfg->srcpattern[(int)(ry*(1.f-EPS)*cfg->srcparam2.w)*(int)(cfg->srcparam1.w)+(int)(rx*(1.f-EPS)*cfg->srcparam1.w)];
		}else if(cfg->srctype==stFourier){
			r->weight=(cosf((floorf(cfg->srcparam1.w)*rx+floorf(cfg->srcparam2.w)*ry+cfg->srcparam1.w-floorf(cfg->srcparam1.w))*TWO_PI)*(1.f-cfg->srcparam2.w+floorf(cfg->srcparam2.w))+1.f)*0.5f;
		}
	}else if(cfg->srctype==stFourierX || cfg->srctype==stFourier2D){
		float rx=rand_uniform01(ran);
		float ry=rand_uniform01(ran);
		float4 v2=cfg->srcparam1;
		v2.w*=1.f/(sqrt(cfg->srcparam1.x*cfg->srcparam1.x+cfg->srcparam1.y*cfg->srcparam1.y+cfg->srcparam1.z*cfg->srcparam1.z));
		v2.x=v2.w*(cfg->srcdir.y*cfg->srcparam1.z - cfg->srcdir.z*cfg->srcparam1.y);
		v2.y=v2.w*(cfg->srcdir.z*cfg->srcparam1.x - cfg->srcdir.x*cfg->srcparam1.z); 
		v2.z=v2.w*(cfg->srcdir.x*cfg->srcparam1.y - cfg->srcdir.y*cfg->srcparam1.x);
		r->p0.x=cfg->srcpos.x+rx*cfg->srcparam1.x+ry*v2.x;
		r->p0.y=cfg->srcpos.y+rx*cfg->srcparam1.y+ry*v2.y;
		r->p0.z=cfg->srcpos.z+rx*cfg->srcparam1.z+ry*v2.z;
		if(cfg->srctype==stFourier2D)
			r->weight=(sinf((cfg->srcparam2.x*rx+cfg->srcparam2.z)*TWO_PI)*sinf((cfg->srcparam2.y*ry+cfg->srcparam2.w)*TWO_PI)+1.f)*0.5f; //between 0 and 1
		else
			r->weight=(cosf((cfg->srcparam2.x*rx+cfg->srcparam2.y*ry+cfg->srcparam2.z)*TWO_PI)*(1.f-cfg->srcparam2.w)+1.f)*0.5f; //between 0 and 1
	}else if(cfg->srctype==stDisk){  // uniform disk distribution
		float sphi, cphi;
		float phi=TWO_PI*rand_uniform01(ran);
		sphi=sinf(phi);	cphi=cosf(phi);
		float r0=sqrt(rand_uniform01(ran))*cfg->srcparam1.x;
		if(cfg->srcdir.z>-1.f+EPS && cfg->srcdir.z<1.f-EPS){
		    float tmp0=1.f-cfg->srcdir.z*cfg->srcdir.z;
		    float tmp1=r0/sqrt(tmp0);
		    r->p0.x=cfg->srcpos.x+tmp1*(cfg->srcdir.x*cfg->srcdir.z*cphi - cfg->srcdir.y*sphi);
		    r->p0.y=cfg->srcpos.y+tmp1*(cfg->srcdir.y*cfg->srcdir.z*cphi + cfg->srcdir.x*sphi);
		    r->p0.z=cfg->srcpos.z-tmp1*tmp0*cphi;
		}else{
   		    r->p0.x+=r0*cphi;
		    r->p0.y+=r0*sphi;
		}
	}else if(cfg->srctype==stCone || cfg->srctype==stIsotropic || cfg->srctype==stArcSin){
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=sinf(ang);	cphi=cosf(ang);
		if(cfg->srctype==stCone){  // a solid-angle section of a uniform sphere
		        do{
				ang=(cfg->srcparam1.y>0) ? TWO_PI*rand_uniform01(ran) : acosf(2.f*rand_uniform01(ran)-1.f); //sine distribution
		        }while(ang>cfg->srcparam1.x);
		}else{
			if(cfg->srctype==stIsotropic) // uniform sphere
				ang=acosf(2.f*rand_uniform01(ran)-1.f); //sine distribution
			else
				ang=M_PI*rand_uniform01(ran); //uniform distribution in zenith angle, arcsine
		}
		stheta=sinf(ang);
		ctheta=cosf(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
	}else if(cfg->srctype==stGaussian){
		float ang,stheta,ctheta,sphi,cphi;
		ang=TWO_PI*rand_uniform01(ran); //next arimuth angle
		sphi=sinf(ang);	cphi=cosf(ang);
		ang=sqrt(-2.f*log(rand_uniform01(ran)))*(1.f-2.f*rand_uniform01(ran0))*cfg->srcparam1.x;
		stheta=sinf(ang);
		ctheta=cosf(ang);
		r->vec.x=stheta*cphi;
		r->vec.y=stheta*sphi;
		r->vec.z=ctheta;
	}

	/*Caluclate intial element id and bary-centric coordinates*/
	float3 vecS={0.f}, *nodes=mesh->node, vecAB, vecAC, vecN;
	int is,i,ea,eb,ec;
	float bary[4]={0.f};
	for(is=0;is<mesh->srcelemlen;is++){
		int include = 1;
		int *elems=(int *)(mesh->elem+(mesh->srcelem[is])-1);
		for(i=0;i<4;i++){
			ea=elems[out[i][0]]-1;
			eb=elems[out[i][1]]-1;
			ec=elems[out[i][2]]-1;
			vec_diff(&nodes[ea],&nodes[eb],&vecAB); //repeated for all photons
			vec_diff(&nodes[ea],&nodes[ec],&vecAC); //repeated for all photons
			vec_diff(&nodes[ea],&(r->p0),&vecS);
			vec_cross(&vecAB,&vecAC,&vecN); //repeated for all photons, vecN can be precomputed
			bary[facemap[i]]=-vec_dot(&vecS,&vecN);
		}
		for(i=0;i<4;i++){
			if(bary[i]<-1e-4f){
				include = 0;
			}
		}
		if(include){
			r->eid=mesh->srcelem[is];
			float s=0.f;
			for(i=0;i<4;i++){s+=bary[i];}
			r->bary0.x=bary[0]/s;
			r->bary0.y=bary[1]/s;
			r->bary0.z=bary[2]/s;
			r->bary0.w=bary[3]/s;
			for(i=0;i<4;i++){
				if((bary[i]/s)<1e-4f)
					r->faceid=ifacemap[i]+1;
			}
			break;
		}
	}
	if(is==mesh->srcelemlen){
//		fprintf(stdout,"initial element does not enclose the source!");
//		fprintf(stdout,"source position [%e %e %e] \n",r->p0.x,r->p0.y,r->p0.z);
//		fprintf(stdout,"bary centric volume [%e %e %e %e] \n",bary[0],bary[1],bary[2],bary[3]);
		mesh_error("initial element does not enclose the source!");
	}
}


