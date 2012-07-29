/*******************************************************************************
**  Mesh-based Monte Carlo (MMC) -- this unit was ported from MCX
**
**  Monte Carlo eXtreme (MCX)  - GPU accelerated Monte Carlo 3D photon migration
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
\file    mcx_utils.c

\brief   mcconfiguration and command line option processing unit
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/ioctl.h>
#include "mcx_utils.h"

#define FIND_JSON_KEY(id,idfull,parent,fallback,val) \
                    ((tmp=cJSON_GetObjectItem(parent,id))==0 ? \
                                ((tmp=cJSON_GetObjectItem(root,idfull))==0 ? fallback : tmp->val) \
                     : tmp->val)

#define FIND_JSON_OBJ(id,idfull,parent) \
                    ((tmp=cJSON_GetObjectItem(parent,id))==0 ? \
                                ((tmp=cJSON_GetObjectItem(root,idfull))==0 ? NULL : tmp) \
                     : tmp)

#define MMC_ASSERT(id)   mcx_assert(id,__FILE__,__LINE__)


const char shortopt[]={'h','E','f','n','t','T','s','a','g','b','B','D',
                 'd','r','S','e','U','R','l','L','I','o','u','C','M',
                 'i','V','O','m','F','\0'};
const char *fullopt[]={"--help","--seed","--input","--photon",
                 "--thread","--blocksize","--session","--array",
                 "--gategroup","--reflect","--reflect3","--debug","--savedet",
                 "--repeat","--save2pt","--minenergy",
                 "--normalize","--skipradius","--log","--listgpu",
                 "--printgpu","--root","--unitinmm","--continuity",
                 "--method","--interactive","--specular","--outputtype",
                 "--momentum","--outputformat",""};

const char debugflag[]={'M','C','B','W','D','I','O','X','A','T','R','P','E','\0'};
const char raytracing[]={'p','h','b','s','\0'};
const char outputtype[]={'x','f','e','\0'};
const char *outputformat[]={"ascii","bin","json","ubjson",""};
const char *srctypeid[]={"pencil","isotropic","cone","gaussian",""};

void mcx_initcfg(mcconfig *cfg){
     cfg->medianum=0;
     cfg->detnum=0;
     cfg->dim.x=0;
     cfg->dim.y=0;
     cfg->dim.z=0;
     cfg->nblocksize=128;
     cfg->nphoton=0;
     cfg->nthread=0;
     cfg->seed=0;
     cfg->isrowmajor=1; /* default is C array*/
     cfg->maxgate=1;
     cfg->isreflect=0;
     cfg->isref3=0;
     cfg->isnormalized=1;
     cfg->issavedet=0;
     cfg->respin=1;
     cfg->issave2pt=1;
     cfg->isgpuinfo=0;
     cfg->basisorder=1;
#ifndef MMC_USE_SSE
     cfg->method=0;
#else
     cfg->method=1;
#endif
     cfg->prop=NULL;
     cfg->detpos=NULL;
     cfg->vol=NULL;
     cfg->session[0]='\0';
     cfg->meshtag[0]='\0';
     cfg->minenergy=1e-6f;
     cfg->flog=stdout;
     cfg->sradius=0.f;
     cfg->rootpath[0]='\0';
     cfg->debuglevel=0;
     cfg->minstep=1.f;
     cfg->roulettesize=10.f;
     cfg->nout=1.f;
     cfg->unitinmm=1.f;
     cfg->srctype=0;
     cfg->isspecular=0;
     cfg->outputtype=otFlux;
     cfg->outputformat=ofASCII;
     cfg->ismomentum=0;
     cfg->his.version=1;
     cfg->his.unitinmm=1.f;
     memcpy(cfg->his.magic,"MCXH",4);

     memset(&(cfg->bary0),0,sizeof(float4));
     memset(&(cfg->srcparam),0,sizeof(float4));
     memset(cfg->checkpt,0,sizeof(unsigned int)*MAX_CHECKPOINT);
}

void mcx_clearcfg(mcconfig *cfg){
     if(cfg->medianum)
     	free(cfg->prop);
     if(cfg->detnum)
     	free(cfg->detpos);
     if(cfg->dim.x && cfg->dim.y && cfg->dim.z)
        free(cfg->vol);
     if(cfg->flog)
        fclose(cfg->flog);
     mcx_initcfg(cfg);
}

void mcx_savedata(float *dat,int len,mcconfig *cfg){
     FILE *fp;
     char name[MAX_PATH_LENGTH];
     sprintf(name,"%s.mc2",cfg->session);
     fp=fopen(name,"wb");
     fwrite(dat,sizeof(float),len,fp);
     fclose(fp);
}

void mcx_printlog(mcconfig *cfg, char *str){
     if(cfg->flog>0){ /*stdout is 1*/
         fprintf(cfg->flog,"%s\n",str);
     }
}

void mcx_normalize(float field[], float scale, int fieldlen){
     int i;
     for(i=0;i<fieldlen;i++){
         field[i]*=scale;
     }
}

void mcx_error(const int id,const char *msg,const char *file,const int linenum){
     fprintf(stdout,"\nMMC ERROR(%d):%s in unit %s:%d\n",id,msg,file,linenum);
#ifdef MCX_CONTAINER
     mmc_throw_exception(id,msg,file,linenum);
#else
     exit(id);
#endif
}

void mcx_assert(const int ret,const char *file,const int linenum){
     if(!ret) mcx_error(ret,"assert error",file,linenum);
}

void mcx_readconfig(char *fname, mcconfig *cfg){
     if(fname[0]==0){
     	mcx_loadconfig(stdin,cfg);
        if(cfg->session[0]=='\0'){
		strcpy(cfg->session,"default");
	}
     }else{
        FILE *fp=fopen(fname,"rt");
        if(fp==NULL) MMC_ERROR(-2,"can not load the specified config file");
        if(strstr(fname,".json")!=NULL){
            char *jbuf;
            int len;
            cJSON *jroot;

            fclose(fp);
            fp=fopen(fname,"rb");
            fseek (fp, 0, SEEK_END);
            len=ftell(fp)+1;
            jbuf=(char *)malloc(len);
            rewind(fp);
            if(fread(jbuf,len-1,1,fp)!=1)
                MMC_ERROR(-2,"reading input file is terminated");
            jbuf[len-1]='\0';
            jroot = cJSON_Parse(jbuf);
            if(jroot){
                mcx_loadjson(jroot,cfg);
                cJSON_Delete(jroot);
            }else{
                char *ptrold, *ptr=(char*)cJSON_GetErrorPtr();
                if(ptr) ptrold=strstr(jbuf,ptr);
                fclose(fp);
                if(ptr && ptrold){
                   char *offs=(ptrold-jbuf>=50) ? ptrold-50 : jbuf;
                   while(offs<ptrold){
                      fprintf(stderr,"%c",*offs);
                      offs++;
                   }
                   fprintf(stderr,"<error>%.50s\n",ptrold);
                }
                free(jbuf);
                MMC_ERROR(-9,"invalid JSON input file");
            }
            free(jbuf);
        }else{
	    mcx_loadconfig(fp,cfg); 
        }
        fclose(fp);
        if(cfg->session[0]=='\0'){
		strncpy(cfg->session,fname,MAX_SESSION_LENGTH);
	}
     }
}

int mcx_loadjson(cJSON *root, mcconfig *cfg){
     int i;
     cJSON *Mesh, *Optode, *Forward, *Session, *tmp, *subitem;
     char comment[MAX_PATH_LENGTH];
     Mesh    = cJSON_GetObjectItem(root,"Mesh");
     Optode  = cJSON_GetObjectItem(root,"Optode");
     Session = cJSON_GetObjectItem(root,"Session");
     Forward = cJSON_GetObjectItem(root,"Forward");

     if(Mesh){
        strncpy(cfg->meshtag, FIND_JSON_KEY("MeshID","Mesh.MeshID",Mesh,(MMC_ERROR(-1,"You must specify mesh files"),""),valuestring), MAX_PATH_LENGTH);
        cfg->dim.x=FIND_JSON_KEY("InitElem","Mesh.InitElem",Mesh,(MMC_ERROR(-1,"InitElem must be given"),0.0),valueint);
        if(cfg->rootpath[0]){
#ifdef WIN32
           sprintf(comment,"%s\\%s",cfg->rootpath,cfg->meshtag);
#else
           sprintf(comment,"%s/%s",cfg->rootpath,cfg->meshtag);
#endif
           strncpy(cfg->meshtag,comment,MAX_PATH_LENGTH);
        }
        cfg->unitinmm=FIND_JSON_KEY("LengthUnit","Mesh.LengthUnit",Mesh,1.0,valuedouble);
     }
     if(Optode){
        cJSON *dets, *src=FIND_JSON_OBJ("Source","Optode.Source",Optode);
        if(src){
           subitem=FIND_JSON_OBJ("Pos","Optode.Source.Pos",src);
           if(subitem){
              cfg->srcpos.x=subitem->child->valuedouble;
              cfg->srcpos.y=subitem->child->next->valuedouble;
              cfg->srcpos.z=subitem->child->next->next->valuedouble;
           }
           subitem=FIND_JSON_OBJ("Dir","Optode.Source.Dir",src);
           if(subitem){
              cfg->srcdir.x=subitem->child->valuedouble;
              cfg->srcdir.y=subitem->child->next->valuedouble;
              cfg->srcdir.z=subitem->child->next->next->valuedouble;
           }
           subitem=FIND_JSON_OBJ("Type","Optode.Source.Type",src);
           if(subitem){
              cfg->srctype=mcx_keylookup(subitem->valuestring,srctypeid);
           }
           subitem=FIND_JSON_OBJ("Param","Optode.Source.Param",src);
           if(subitem && cJSON_GetArraySize(subitem)==4){
              cfg->srcparam.x=subitem->child->valuedouble;
              cfg->srcparam.y=subitem->child->next->valuedouble;
              cfg->srcparam.z=subitem->child->next->next->valuedouble;
              cfg->srcparam.w=subitem->child->next->next->next->valuedouble;
           }
        }
        dets=FIND_JSON_OBJ("Detector","Optode.Detector",Optode);
        if(dets){
           cJSON *det=dets->child;
           if(det){
             cfg->detnum=cJSON_GetArraySize(dets);
             cfg->detpos=(float4*)malloc(sizeof(float4)*cfg->detnum);
             for(i=0;i<cfg->detnum;i++){
               cJSON *pos=dets, *rad=NULL;
               rad=FIND_JSON_OBJ("R","Optode.Detector.R",det);
               if(cJSON_GetArraySize(det)==2){
                   pos=FIND_JSON_OBJ("Pos","Optode.Detector.Pos",det);
               }
               if(pos){
	           cfg->detpos[i].x=pos->child->valuedouble;
                   cfg->detpos[i].y=pos->child->next->valuedouble;
	           cfg->detpos[i].z=pos->child->next->next->valuedouble;
               }
               if(rad){
                   cfg->detpos[i].w=rad->valuedouble;
               }
               det=det->next;
               if(det==NULL) break;
             }
           }
        }
     }
     if(Session){
        char val[1];
        cJSON *ck;
        if(cfg->seed==0)      cfg->seed=FIND_JSON_KEY("RNGSeed","Session.RNGSeed",Session,-1,valueint);
        if(cfg->nphoton==0)   cfg->nphoton=FIND_JSON_KEY("Photons","Session.Photons",Session,0,valueint);
        if(cfg->session[0]=='\0') strncpy(cfg->session, FIND_JSON_KEY("ID","Session.ID",Session,"default",valuestring), MAX_SESSION_LENGTH);

        if(!cfg->isreflect)   cfg->isreflect=FIND_JSON_KEY("DoMismatch","Session.DoMismatch",Session,cfg->isreflect,valueint);
        if(cfg->issave2pt)    cfg->issave2pt=FIND_JSON_KEY("DoSaveVolume","Session.DoSaveVolume",Session,cfg->issave2pt,valueint);
        if(cfg->isnormalized) cfg->isnormalized=FIND_JSON_KEY("DoNormalize","Session.DoNormalize",Session,cfg->isnormalized,valueint);
        if(!cfg->issavedet)   cfg->issavedet=FIND_JSON_KEY("DoPartialPath","Session.DoPartialPath",Session,cfg->issavedet,valueint);
        if(!cfg->isspecular)  cfg->isspecular=FIND_JSON_KEY("DoSpecular","Session.DoSpecular",Session,cfg->isspecular,valueint);
        if(!cfg->ismomentum)  cfg->ismomentum=FIND_JSON_KEY("DoDCS","Session.DoDCS",Session,cfg->ismomentum,valueint);
        if(cfg->basisorder)   cfg->basisorder=FIND_JSON_KEY("BasisOrder","Session.BasisOrder",Session,cfg->basisorder,valueint);
        if(!cfg->outputformat)  cfg->outputformat=mcx_keylookup(FIND_JSON_KEY("OutputFormat","Session.OutputFormat",Session,"ascii",valuestring),outputformat);
        if(cfg->outputformat<0)
		MMC_ERROR(-2,"the specified output format is not recognized");

        if(cfg->debuglevel==0)
           cfg->debuglevel=mcx_parsedebugopt(FIND_JSON_KEY("DebugFlag","Session.DebugFlag",Session,"",valuestring));
        strncpy(val,FIND_JSON_KEY("RayTracer","Session.RayTracer",Session,raytracing+cfg->method,valuestring),1);
        if(mcx_lookupindex(val, raytracing)){
		MMC_ERROR(-2,"the specified ray-tracing method is not recognized");
	}
	cfg->method=val[0];
        strncpy(val,FIND_JSON_KEY("OutputType","Session.OutputType",Session,outputtype+cfg->outputtype,valuestring),1);
        if(mcx_lookupindex(val, outputtype)){
                MMC_ERROR(-2,"the specified output data type is not recognized");
        }
	cfg->outputtype=val[0];
        ck=FIND_JSON_OBJ("Checkpoints","Session.Checkpoints",Session);
        if(ck){
            int num=MIN(cJSON_GetArraySize(ck),MAX_CHECKPOINT);
            ck=ck->child;
            for(i=0;i<num;i++){
               cfg->checkpt[i]=ck->valueint;
               ck=ck->next;
               if(ck==NULL) break;
            }
        }
     }
     if(Forward){
        cfg->tstart=FIND_JSON_KEY("T0","Forward.T0",Forward,0.0,valuedouble);
        cfg->tend  =FIND_JSON_KEY("T1","Forward.T1",Forward,0.0,valuedouble);
        cfg->tstep =FIND_JSON_KEY("Dt","Forward.Dt",Forward,0.0,valuedouble);
        cfg->nout=FIND_JSON_KEY("N0","Forward.N0",Forward,cfg->nout,valuedouble);

        cfg->maxgate=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);
     }
     if(cfg->meshtag[0]=='\0')
         MMC_ERROR(-1,"You must specify mesh files");
     if(cfg->dim.x==0)
	 MMC_ERROR(-1,"InitElem must be given");
     return 0;
}

void mcx_writeconfig(char *fname, mcconfig *cfg){
     if(fname[0]==0)
     	mcx_saveconfig(stdout,cfg);
     else{
     	FILE *fp=fopen(fname,"wt");
	if(fp==NULL) MMC_ERROR(-2,"can not write to the specified config file");
	mcx_saveconfig(fp,cfg);     
	fclose(fp);
     }
}

void mcx_loadconfig(FILE *in, mcconfig *cfg){
     int i,gates,srctype,itmp;
     char comment[MAX_PATH_LENGTH],*comm, strtypestr[MAX_SESSION_LENGTH]={'\0'};
     
     if(in==stdin)
     	fprintf(stdout,"Please specify the total number of photons: [1000000]\n\t");
     MMC_ASSERT(fscanf(in,"%d", &(i) )==1); 
     if(cfg->nphoton==0) cfg->nphoton=i;
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(in==stdin)
     	fprintf(stdout,"%d\nPlease specify the random number generator seed: [123456789]\n\t",cfg->nphoton);
     if(cfg->seed==0)
        MMC_ASSERT(fscanf(in,"%d", &(cfg->seed) )==1);
     else
        MMC_ASSERT(fscanf(in,"%d", &itmp )==1);
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(in==stdin)
     	fprintf(stdout,"%d\nPlease specify the position of the source: [10 10 5]\n\t",cfg->seed);
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->srcpos.x),&(cfg->srcpos.y),&(cfg->srcpos.z) )==3);
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(in==stdin)
     	fprintf(stdout,"%f %f %f\nPlease specify the normal direction of the source fiber: [0 0 1]\n\t",
                                   cfg->srcpos.x,cfg->srcpos.y,cfg->srcpos.z);
     //cfg->srcpos.x--;cfg->srcpos.y--;cfg->srcpos.z--; /*convert to C index, grid center*/
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->srcdir.x),&(cfg->srcdir.y),&(cfg->srcdir.z) )==3);
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(in==stdin)
     	fprintf(stdout,"%f %f %f\nPlease specify the time gates in seconds (start end and step) [0.0 1e-9 1e-10]\n\t",
                                   cfg->srcdir.x,cfg->srcdir.y,cfg->srcdir.z);
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->tstart),&(cfg->tend),&(cfg->tstep) )==3);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	fprintf(stdout,"%f %f %f\nPlease specify the path to the volume binary file:\n\t",
                                   cfg->tstart,cfg->tend,cfg->tstep);
     if(cfg->tstart>cfg->tend || cfg->tstep==0.f){
         MMC_ERROR(-9,"incorrect time gate settings");
     }
     gates=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);
     /*if(cfg->maxgate>gates)*/
	 cfg->maxgate=gates;

     MMC_ASSERT(fscanf(in,"%s", cfg->meshtag)==1);
     if(cfg->rootpath[0]){
#ifdef WIN32
         sprintf(comment,"%s\\%s",cfg->rootpath,cfg->meshtag);
#else
         sprintf(comment,"%s/%s",cfg->rootpath,cfg->meshtag);
#endif
         strncpy(cfg->meshtag,comment,MAX_PATH_LENGTH);
     }
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	fprintf(stdout,"%s\nPlease specify the index to the tetrahedral element enclosing the source [start from 1]:\n\t",cfg->meshtag);
     MMC_ASSERT(fscanf(in,"%d", &(cfg->dim.x))==1);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     cfg->minstep=MIN(cfg->steps.x,cfg->steps.y);
     cfg->minstep=MIN(cfg->minstep,cfg->steps.z);

     if(in==stdin)
     	fprintf(stdout,"Please specify the total number of detectors and fiber diameter (in mm):\n\t");
     MMC_ASSERT(fscanf(in,"%d %f", &(cfg->detnum), &(cfg->detradius))==2);
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(in==stdin)
     	fprintf(stdout,"%d %f\n",cfg->detnum,cfg->detradius);
     cfg->detpos=(float4*)malloc(sizeof(float4)*cfg->detnum);
     if(cfg->issavedet)
        cfg->issavedet=(cfg->detpos>0);
     for(i=0;i<cfg->detnum;i++){
        if(in==stdin)
		fprintf(stdout,"Please define detector #%d: x,y,z (in mm): [5 5 5 1]\n\t",i);
     	MMC_ASSERT(fscanf(in, "%f %f %f", &(cfg->detpos[i].x),&(cfg->detpos[i].y),&(cfg->detpos[i].z))==3);
        //cfg->detpos[i].x--;cfg->detpos[i].y--;cfg->detpos[i].z--;  /*convert to C index*/
        comm=fgets(comment,MAX_PATH_LENGTH,in);
        if(in==stdin)
		fprintf(stdout,"%f %f %f\n",cfg->detpos[i].x,cfg->detpos[i].y,cfg->detpos[i].z);
     }
     if(in==stdin)
     	fprintf(stdout,"Please specify the source type[pencil|cone|gaussian]:\n\t");
     if(fscanf(in,"%s", strtypestr)==1 && strtypestr[0]){
        srctype=mcx_keylookup(strtypestr,srctypeid);
	if(srctype==-1)
	   MMC_ERROR(-6,"the specified source type is not supported");
        if(srctype>0){
           comm=fgets(comment,MAX_PATH_LENGTH,in);
	   cfg->srctype=srctype;
	   if(in==stdin)
     	      fprintf(stdout,"Please specify the source parameters (4 floating-points):\n\t");
           MMC_ASSERT(fscanf(in, "%f %f %f %f", &(cfg->srcparam.x),
	          &(cfg->srcparam.y),&(cfg->srcparam.z),&(cfg->srcparam.w))==4);
	}else
	   return;
     }else
        return;
}

void mcx_saveconfig(FILE *out, mcconfig *cfg){
     int i;

     fprintf(out,"%d\n", (cfg->nphoton) ); 
     fprintf(out,"%d\n", (cfg->seed) );
     fprintf(out,"%f %f %f\n", (cfg->srcpos.x),(cfg->srcpos.y),(cfg->srcpos.z) );
     fprintf(out,"%f %f %f\n", (cfg->srcdir.x),(cfg->srcdir.y),(cfg->srcdir.z) );
     fprintf(out,"%f %f %f\n", (cfg->tstart),(cfg->tend),(cfg->tstep) );
     fprintf(out,"%f %d %d %d\n", (cfg->steps.x),(cfg->dim.x),(cfg->crop0.x),(cfg->crop1.x));
     fprintf(out,"%f %d %d %d\n", (cfg->steps.y),(cfg->dim.y),(cfg->crop0.y),(cfg->crop1.y));
     fprintf(out,"%f %d %d %d\n", (cfg->steps.z),(cfg->dim.z),(cfg->crop0.z),(cfg->crop1.z));
     fprintf(out,"%d", (cfg->medianum));
     for(i=0;i<cfg->medianum;i++){
     	fprintf(out, "%f %f %f %f\n", (cfg->prop[i].mus),(cfg->prop[i].g),(cfg->prop[i].mua),(cfg->prop[i].n));
     }
     fprintf(out,"%d", (cfg->detnum));
     for(i=0;i<cfg->detnum;i++){
     	fprintf(out, "%f %f %f %f\n", (cfg->detpos[i].x),(cfg->detpos[i].y),(cfg->detpos[i].z),(cfg->detpos[i].w));
     }
}

void mcx_loadvolume(char *filename,mcconfig *cfg){
     int datalen,res;
     FILE *fp=fopen(filename,"rb");
     if(fp==NULL){
     	     MMC_ERROR(-5,"the specified binary volume file does not exist");
     }
     if(cfg->vol){
     	     free(cfg->vol);
     	     cfg->vol=NULL;
     }
     datalen=cfg->dim.x*cfg->dim.y*cfg->dim.z;
     cfg->vol=(unsigned char*)malloc(sizeof(unsigned char)*datalen);
     res=fread(cfg->vol,sizeof(unsigned char),datalen,fp);
     fclose(fp);
     if(res!=datalen){
     	 MMC_ERROR(-6,"file size does not match specified dimensions");
     }
}

int mcx_parsedebugopt(char *debugopt){
    char *c=debugopt,*p;
    int debuglevel=0;

    while(*c){
       p=strchr(debugflag, ((*c<='z' && *c>='a') ? *c-'a'+'A' : *c) );
       if(p!=NULL)
	  debuglevel |= (1 << (p-debugflag));
       c++;
    }
    return debuglevel;
}

void mcx_progressbar(unsigned int n, mcconfig *cfg){
    unsigned int percentage, j,colwidth=79;
    static unsigned int oldmarker=0xFFFFFFFF;

#ifdef TIOCGWINSZ 
    struct winsize ttys;
    ioctl(0, TIOCGWINSZ, &ttys);
    colwidth=ttys.ws_col;
#endif
    
    percentage=(float)n*(colwidth-18)/cfg->nphoton;

    if(percentage != oldmarker){
        oldmarker=percentage;
    	for(j=0;j<colwidth;j++)     fprintf(stdout,"\b");
    	fprintf(stdout,"Progress: [");
    	for(j=0;j<percentage;j++)      fprintf(stdout,"=");
    	fprintf(stdout,(percentage<colwidth-18) ? ">" : "=");
    	for(j=percentage;j<colwidth-18;j++) fprintf(stdout," ");
    	fprintf(stdout,"] %3d%%",percentage*100/(colwidth-18));
	    fflush(stdout);
    }
}

int mcx_readarg(int argc, char *argv[], int id, void *output,char *type){
     /*
         when a binary option is given without a following number (0~1), 
         we assume it is 1
     */
     if(strcmp(type,"bool")==0 && (id>=argc-1||(argv[id+1][0]<'0'||argv[id+1][0]>'9'))){
	*((char*)output)=1;
	return id;
     }
     if(id<argc-1){
         if(strcmp(type,"bool")==0)
             *((char*)output)=atoi(argv[id+1]);
         else if(strcmp(type,"char")==0)
             *((char*)output)=argv[id+1][0];
	 else if(strcmp(type,"int")==0)
             *((int*)output)=atoi(argv[id+1]);
	 else if(strcmp(type,"float")==0)
             *((float*)output)=atof(argv[id+1]);
	 else if(strcmp(type,"string")==0)
	     strcpy((char *)output,argv[id+1]);
     }else{
     	 MMC_ERROR(-1,"incomplete input");
     }
     return id+1;
}
int mcx_remap(char *opt){
    int i=0;
    while(shortopt[i]!='\0'){
	if(strcmp(opt,fullopt[i])==0){
		opt[1]=shortopt[i];
		opt[2]='\0';
		return 0;
	}
	i++;
    }
    return 1;
}
int mcx_lookupindex(char *key, const char *index){
    int i=0;
    while(index[i]!='\0'){
        if(tolower(*key)==index[i]){
                *key=i;
                return 0;
        }
        i++;
    }
    return 1;
}
int mcx_keylookup(char *key, const char *table[]){
    int i=0;
    while(key[i]){
        key[i]=tolower(key[i]);
	i++;
    }
    i=0;
    while(table[i]!='\0'){
	if(strcmp(key,table[i])==0){
		return i;
	}
	i++;
    }
    return -1;
}
void mcx_parsecmd(int argc, char* argv[], mcconfig *cfg){
     int i=1,isinteractive=1,issavelog=0;
     char filename[MAX_PATH_LENGTH]={0};
     char logfile[MAX_PATH_LENGTH]={0};
     float np=0.f;

     if(argc<=1){
     	mcx_usage(argv[0]);
     	exit(0);
     }
     while(i<argc){
     	    if(argv[i][0]=='-'){
		if(argv[i][1]=='-'){
			if(mcx_remap(argv[i])){
				MMC_ERROR(-2,"unsupported verbose option");
			}
		}
	        switch(argv[i][1]){
		     case 'h':
		                mcx_usage(argv[0]);
				exit(0);
		     case 'i':
				if(filename[0]){
					MMC_ERROR(-2,"you can not specify both interactive mode and config file");
				}
		     		isinteractive=1;
				break;
		     case 'f': 
		     		isinteractive=0;
		     	        i=mcx_readarg(argc,argv,i,filename,"string");
				break;
		     case 'n':
		     	        i=mcx_readarg(argc,argv,i,&np,"float");
				cfg->nphoton=(int)np;
		     	        break;
		     case 't':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->nthread),"int");
		     	        break;
                     case 'T':
                               	i=mcx_readarg(argc,argv,i,&(cfg->nblocksize),"int");
                               	break;
		     case 's':
		     	        i=mcx_readarg(argc,argv,i,cfg->session,"string");
		     	        break;
		     case 'a':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isrowmajor),"bool");
		     	        break;
		     case 'g':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->maxgate),"int");
		     	        break;
		     case 'b':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isreflect),"bool");
		     	        break;
                     case 'B':
                                i=mcx_readarg(argc,argv,i,&(cfg->isref3),"bool");
                               	break;
		     case 'd':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->issavedet),"bool");
		     	        break;
		     case 'm':
		                i=mcx_readarg(argc,argv,i,&(cfg->ismomentum),"bool");
				if (cfg->ismomentum) cfg->issavedet=1;
				break;			       
		     case 'C':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->basisorder),"bool");
		     	        break;
		     case 'V':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isspecular),"bool");
		     	        break;
		     case 'r':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->respin),"int");
		     	        break;
		     case 'S':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->issave2pt),"bool");
		     	        break;
                     case 'e':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->minenergy),"float");
                                break;
		     case 'U':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isnormalized),"bool");
		     	        break;
		     case 'E':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->seed),"int");
		     	        break;
                     case 'F':
                                if(i>=argc)
                                        MMC_ERROR(-1,"incomplete input");
                                if((cfg->outputformat=mcx_keylookup(argv[++i], outputformat))<0)
                                        MMC_ERROR(-2,"the specified output data type is not recognized");
                                break;
                     case 'O':
                                i=mcx_readarg(argc,argv,i,&(cfg->outputtype),"char");
				if(mcx_lookupindex(&(cfg->outputtype), outputtype)){
                                        MMC_ERROR(-2,"the specified output data type is not recognized");
                                }
                                break;
                     case 'M':
                                i=mcx_readarg(argc,argv,i,&(cfg->method),"char");
				if(mcx_lookupindex(&(cfg->method), raytracing)){
					MMC_ERROR(-2,"the specified ray-tracing method is not recognized");
				}
                                break;
                     case 'R':
                                i=mcx_readarg(argc,argv,i,&(cfg->sradius),"float");
                                break;
                     case 'u':
                                i=mcx_readarg(argc,argv,i,&(cfg->unitinmm),"float");
                                break;
                     case 'l':
                                issavelog=1;
                                break;
		     case 'L':
                                cfg->isgpuinfo=2;
		                break;
		     case 'I':
                                cfg->isgpuinfo=1;
		                break;
		     case 'o':
		     	        i=mcx_readarg(argc,argv,i,cfg->rootpath,"string");
		     	        break;
                     case 'D':
				if(i+1<argc && isalpha(argv[i+1][0]) )
					cfg->debuglevel=mcx_parsedebugopt(argv[++i]);
				else
	                                i=mcx_readarg(argc,argv,i,&(cfg->debuglevel),"int");
                                break;
                     default:
				MMC_ERROR(-1,"unsupported command line option");
		}
	    }
	    i++;
     }
     if(issavelog && cfg->session){
          sprintf(logfile,"%s.log",cfg->session);
          cfg->flog=fopen(logfile,"wt");
          if(cfg->flog==NULL){
		cfg->flog=stdout;
		fprintf(cfg->flog,"unable to save to log file, will print from stdout\n");
          }
     }
     if(cfg->isgpuinfo!=2){ /*print gpu info only*/
       if(isinteractive){
          mcx_readconfig("",cfg);
       }else{
     	  mcx_readconfig(filename,cfg);
       }
     }
}

void mcx_usage(char *exename){
     printf("\
###############################################################################\n\
#                         Mesh-based Monte Carlo (MMC)                        #\n\
#     Copyright (c) 2010,2011 Qianqian Fang <fangq at nmr.mgh.harvard.edu>    #\n\
#                                                                             #\n\
#    Martinos Center for Biomedical Imaging, Massachusetts General Hospital   #\n\
###############################################################################\n\
$MMC $Rev::e72b60d $ Last Commit $Date::2011-12-20 13:58:34$ by $Author:: fangq$\n\
###############################################################################\n\
\n\
usage: %s <param1> <param2> ...\n\
where possible parameters include (the first item in [] is the default value)\n\
 -i 	       (--interactive) interactive mode\n\
 -s sessionid  (--session)     a string used to tag all output file names\n\
 -f config     (--input)       read config from a file\n\
 -n [0.|float] (--photon)      total photon number, max allowed value is 2^32-1\n\
 -b [0|1]      (--reflect)     1 do reflection at int&ext boundaries, 0 no ref.\n\
 -e [0.|float] (--minenergy)   minimum energy level to trigger Russian roulette\n\
 -U [1|0]      (--normalize)   1 to normalize the fluence to unitary,0 save raw\n\
 -d [0|1]      (--savedet)     1 to save photon info at detectors,0 not to save\n\
 -m [0|1]      (--momentum)    1 to save photon momentum transfer,0 not to save\n\
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save\n\
 -C [1|0]      (--basisorder)  1 piece-wise-linear basis for fluence,0 constant\n\
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh\n\
 -O [X|XFE]    (--outputtype)  X - output flux, F - fluence, E - energy deposit\n\
 -F format     (--outputformat)'ascii', 'bin' (in 'double'), 'json' or 'ubjson'\n\
 -u [1.|float] (--unitinmm)    define the length unit in mm for the mesh\n\
 -h            (--help)        print this message\n\
 -l            (--log)         print messages to a log file instead\n\
 -E [0|int]    (--seed)        set random-number-generator seed\n\
 -M [%c|PHBS]   (--method)      choose ray-tracing algorithm (only use 1 letter)\n\
                               P - Plucker-coordinate ray-tracing algorithm\n\
			       H - Havel's SSE4 ray-tracing algorithm\n\
			       B - partial Badouel's method (used by TIM-OS)\n\
			       S - branch-less Badouel's method with SSE\n\
 -D [0|int]    (--debug)       print debug information (you can use an integer\n\
  or                           or a string by combining the following flags)\n\
 -D [''|MCBWDIOXATRPE]         1 M  photon movement info\n\
                               2 C  print ray-polygon testing details\n\
                               4 B  print Bary-centric coordinates\n\
                               8 W  print photon weight changes\n\
                              16 D  print distances\n\
                              32 I  entering a triangle\n\
                              64 O  exiting a triangle\n\
                             128 X  hitting an edge\n\
                             256 A  accumulating weights to the mesh\n\
                             512 T  timing information\n\
                            1024 R  debugging reflection\n\
                            2048 P  show progress bar\n\
                            4096 E  exit photon info\n\
      combine multiple items by using a string, or add selected numbers together\n\
example:\n\
       %s -n 1000000 -f input.inp -s test -b 0 -D TP\n",exename,
#ifdef MMC_USE_SSE
'H',
#else
'P',
#endif
exename);
}
