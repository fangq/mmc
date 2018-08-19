/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mcx_utils.c

\brief   MC simulation settings and command line option processing unit
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/ioctl.h>
#include "mcx_utils.h"

/**
 * Macro to load JSON keys
 */

#define FIND_JSON_KEY(id,idfull,parent,fallback,val) \
                    ((tmp=cJSON_GetObjectItem(parent,id))==0 ? \
                                ((tmp=cJSON_GetObjectItem(root,idfull))==0 ? fallback : tmp->val) \
                     : tmp->val)

/**
 * Macro to load JSON object
 */

#define FIND_JSON_OBJ(id,idfull,parent) \
                    ((tmp=cJSON_GetObjectItem(parent,id))==0 ? \
                                ((tmp=cJSON_GetObjectItem(root,idfull))==0 ? NULL : tmp) \
                     : tmp)

/**
 * Macro to include unit name and line number in the error message
 */
#define MMC_ASSERT(id)   mcx_assert(id,__FILE__,__LINE__)


/**
 * Short command line options
 * If a short command line option is '-' that means it only has long/verbose option.
 * Array terminates with '\0'.
 */

const char shortopt[]={'h','E','f','n','t','T','s','a','g','b','D',
                 'd','r','S','e','U','R','l','L','I','o','u','C','M',
                 'i','V','O','-','F','q','x','P','k','v','m','-','-','-','\0'};
		 
/**
 * Long command line options
 * The length of this array must match the length of shortopt[], terminates with ""
 */

const char *fullopt[]={"--help","--seed","--input","--photon",
                 "--thread","--blocksize","--session","--array",
                 "--gategroup","--reflect","--debug","--savedet",
                 "--repeat","--save2pt","--minenergy",
                 "--normalize","--skipradius","--log","--listgpu",
                 "--printgpu","--root","--unitinmm","--basisorder",
                 "--method","--interactive","--specular","--outputtype",
                 "--momentum","--outputformat","--saveseed","--saveexit",
                 "--replaydet","--voidtime","--version","--mc","--atomic",
                 "--debugphoton",""};

/**
 * Debug flags
 * R: debug random number generator
 * M: record photon movement and trajectory
 * P: show progress bar
 */

const char debugflag[]={'M','C','B','W','D','I','O','X','A','T','R','P','E','\0'};

/**
 * Selecting mesh-based ray-tracing algorithm:
 * p: Plucker-based ray-tracer, see Fang2010
 * h: Havel-based SSE4 ray-tracer, see Fang2012
 * b: Badouel ray-tracing algorithm, see Fang2011
 * s: branch-less Badouel SSE4 ray-tracer, see Fang2011
 * g: grid-output using dual-mesh MMC
 */

const char raytracing[]={'p','h','b','s','g','\0'};

/**
 * Output data types
 * x: fluence rate
 * f: fluence
 * e: energy deposit
 * j: jacobian for mua
 * p: scattering counts for computing Jacobians for mus
 */

const char outputtype[]={'x','f','e','j','l','p','\0'};

/**
 * Output file format
 * mc2: binary mc2 format to store fluence volume data
 * nii: output fluence in nii format
 * hdr: output volume in Analyze hdr/img format
 * ubj: output volume in unversal binary json format (not implemented)
 */

const char *outputformat[]={"ascii","bin","json","ubjson",""};

/**
 * Source type specifier
 * User can specify the source type using a string
 */

const char *srctypeid[]={"pencil","isotropic","cone","gaussian","planar",
    "pattern","fourier","arcsine","disk","fourierx","fourierx2d","zgaussian","line","slit",""};

/**
 * @brief Initializing the simulation configuration with default values
 *
 * Constructor of the simulation configuration, initializing all field to default values
 */

void mcx_initcfg(mcconfig *cfg){
     cfg->medianum=0;
     cfg->detnum=0;
     cfg->e0=0;
     cfg->dim.x=0;
     cfg->dim.y=0;
     cfg->dim.z=0;
     cfg->nblocksize=128;
     cfg->nphoton=0;
     cfg->nthread=0;
     cfg->seed=0x623F9A9E;
     cfg->isrowmajor=0;      /* not needed */
     cfg->maxgate=1;
     cfg->isreflect=1;
     cfg->isref3=1;
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
     cfg->seedfile[0]='\0';
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
     cfg->issaveseed=0;
     cfg->issaveexit=0;
     cfg->photonseed=NULL;
     cfg->replaydet=0;
     cfg->replayweight=NULL;
     cfg->replaytime=NULL;
     cfg->isextdet=0;
     cfg->srcdir.w=0.f;
     cfg->isatomic=1;
     cfg->debugphoton=-1;

     cfg->tstart=0.f;
     cfg->tstep=0.f;
     cfg->tend=0.f;

     cfg->mcmethod=mmMCX;

     memset(&(cfg->his),0,sizeof(history));
     cfg->his.version=1;
     cfg->his.unitinmm=1.f;
     cfg->his.normalizer=1.f;
     memcpy(cfg->his.magic,"MCXH",4);

     memset(&(cfg->bary0),0,sizeof(float4));
     memset(&(cfg->srcparam1),0,sizeof(float4));
     memset(&(cfg->srcparam2),0,sizeof(float4));
     cfg->srcpattern=NULL;
     cfg->voidtime=1;
     memset(cfg->checkpt,0,sizeof(unsigned int)*MAX_CHECKPOINT);
	
     memset(&(cfg->detparam1),0,sizeof(float4));
     memset(&(cfg->detparam2),0,sizeof(float4));
     cfg->detpattern=NULL;
}

/**
 * @brief Clearing the simulation configuration data structure
 *
 * Destructor of the simulation configuration, delete all dynamically allocated members
 */

void mcx_clearcfg(mcconfig *cfg){
     if(cfg->medianum)
     	free(cfg->prop);
     if(cfg->detnum)
     	free(cfg->detpos);
     if(cfg->vol)
        free(cfg->vol);
     if(cfg->srcpattern)
        free(cfg->srcpattern);
     if(cfg->detpattern)
        free(cfg->detpattern);
     if(cfg->photonseed)
        free(cfg->photonseed);
     if(cfg->replayweight)
        free(cfg->replayweight);
     if(cfg->replaytime)
        free(cfg->replaytime);
     if(cfg->flog && cfg->flog!=stdout && cfg->flog!=stderr)
        fclose(cfg->flog);
     mcx_initcfg(cfg);
}

/**
 * @brief Save volumetric output (fluence etc) to mc2 format binary file
 *
 * @param[in] dat: volumetric data to be saved
 * @param[in] len: total byte length of the data to be saved
 * @param[in] cfg: simulation configuration
 */

void mcx_savedata(float *dat,int len,mcconfig *cfg){
     FILE *fp;
     char name[MAX_PATH_LENGTH];
     sprintf(name,"%s.mc2",cfg->session);
     fp=fopen(name,"wb");
     fwrite(dat,sizeof(float),len,fp);
     fclose(fp);
}

/**
 * @brief Print a message to the console or a log file
 *
 * @param[in] cfg: simulation configuration
 * @param[in] str: a string to be printed
 */

void mcx_printlog(mcconfig *cfg, char *str){
     if(cfg->flog>0){ /*stdout is 1*/
         MMC_FPRINTF(cfg->flog,"%s\n",str);
     }
}

/**
 * @brief Normalize the solution by multiplying a scaling factor
 *
 * @param[in,out] field: volumetric data before normalization
 * @param[in] scale: the scaling factor (or normalization factor) to be applied
 * @param[in] fieldlen: the length (floating point) of elements in the volume
 */

void mcx_normalize(float field[], float scale, int fieldlen){
     int i;
     for(i=0;i<fieldlen;i++){
         field[i]*=scale;
     }
}

/**
 * @brief Error reporting function
 *
 * @param[in] id: a single integer for the types of the error
 * @param[in] msg: the error message string
 * @param[in] file: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

void mcx_error(const int id,const char *msg,const char *file,const int linenum){
#pragma omp critical
{
#ifdef MCX_CONTAINER
     mmc_throw_exception(id,msg,file,linenum);
#else
     if(id==MMC_INFO)
        MMC_FPRINTF(stdout,"%s\n",msg);
     else
        MMC_FPRINTF(stdout,"\nMMC ERROR(%d):%s in unit %s:%d\n",id,msg,file,linenum);
     exit(id);
#endif
}
}

/**
 * @brief Function to test return value and raise an error
 *
 * @param[in] ret: function return value, non-zero means an error
 * @param[in] file: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

void mcx_assert(const int ret,const char *file,const int linenum){
     if(!ret) mcx_error(ret,"input error",file,linenum);
}

/**
 * @brief Read simulation settings from a configuration file (.inp or .json)
 *
 * @param[in] fname: the name of the input file (.inp or .json)
 * @param[in] cfg: simulation configuration
 */

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

            fclose(fp);
            fp=fopen(fname,"rb");
            fseek (fp, 0, SEEK_END);
            len=ftell(fp)+1;
            jbuf=(char *)malloc(len);
            rewind(fp);
            if(fread(jbuf,len-1,1,fp)!=1)
                MMC_ERROR(-2,"reading input file is terminated");
            jbuf[len-1]='\0';
            if(mcx_loadfromjson(jbuf,cfg)){
	        fclose(fp);
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

/**
 * @brief Load a .json input file into memory
 *
 * This function loads a JSON file using cJSON
 *
 * @param[out] jbuf: json data structure pointer
 * @param[in] cfg: simulation configuration
 */

int mcx_loadfromjson(char *jbuf,mcconfig *cfg){
     cJSON *jroot;
     jroot = cJSON_Parse(jbuf);
     if(jroot){
        mcx_loadjson(jroot,cfg);
        cJSON_Delete(jroot);
     }else{
        char *ptrold, *ptr=(char*)cJSON_GetErrorPtr();
        if(ptr) ptrold=strstr(jbuf,ptr);
        if(ptr && ptrold){
           char *offs=(ptrold-jbuf>=50) ? ptrold-50 : jbuf;
           while(offs<ptrold){
              MMC_FPRINTF(stderr,"%c",*offs);
              offs++;
           }
           MMC_FPRINTF(stderr,"<error>%.50s\n",ptrold);
        }
        return 1;
     }
     return 0;
}

/**
 * @brief Load user inputs from a .json input file
 *
 * This function loads user input from a JSON format in a .json extension
 *
 * @param[out] root: json data structure pointer
 * @param[in] cfg: simulation configuration
 */

int mcx_loadjson(cJSON *root, mcconfig *cfg){
     int i;
     cJSON *Mesh, *Optode, *Forward, *Session, *tmp, *subitem;

     Mesh    = cJSON_GetObjectItem(root,"Mesh");
     if(!Mesh)
         Mesh    = cJSON_GetObjectItem(root,"Domain");
     Optode  = cJSON_GetObjectItem(root,"Optode");
     Session = cJSON_GetObjectItem(root,"Session");
     Forward = cJSON_GetObjectItem(root,"Forward");

     if(Mesh){
        strncpy(cfg->meshtag, FIND_JSON_KEY("MeshID","Mesh.MeshID",Mesh,(MMC_ERROR(-1,"You must specify mesh files"),""),valuestring), MAX_PATH_LENGTH);
        cfg->e0=FIND_JSON_KEY("InitElem","Mesh.InitElem",Mesh,(MMC_ERROR(-1,"InitElem must be given"),0.0),valueint);
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
	      if(subitem->child->next->next->next)
	         cfg->srcdir.w=subitem->child->next->next->next->valuedouble;
           }
           subitem=FIND_JSON_OBJ("Type","Optode.Source.Type",src);
           if(subitem){
              cfg->srctype=mcx_keylookup(subitem->valuestring,srctypeid);
           }
           subitem=FIND_JSON_OBJ("Param1","Optode.Source.Param1",src);
           if(subitem && cJSON_GetArraySize(subitem)==4){
              cfg->srcparam1.x=subitem->child->valuedouble;
              cfg->srcparam1.y=subitem->child->next->valuedouble;
              cfg->srcparam1.z=subitem->child->next->next->valuedouble;
              cfg->srcparam1.w=subitem->child->next->next->next->valuedouble;
           }
           subitem=FIND_JSON_OBJ("Param2","Optode.Source.Param2",src);
           if(subitem && cJSON_GetArraySize(subitem)==4){
              cfg->srcparam2.x=subitem->child->valuedouble;
              cfg->srcparam2.y=subitem->child->next->valuedouble;
              cfg->srcparam2.z=subitem->child->next->next->valuedouble;
              cfg->srcparam2.w=subitem->child->next->next->next->valuedouble;
           }
        }
        dets=FIND_JSON_OBJ("Detector","Optode.Detector",Optode);
        if(dets){
           cJSON *det=dets;
           if(!FIND_JSON_OBJ("Pos","Optode.Detector.Pos",dets))
              det=dets->child;

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
        if(cfg->seed==0x623F9A9E)      cfg->seed=FIND_JSON_KEY("RNGSeed","Session.RNGSeed",Session,-1,valueint);
        if(cfg->nphoton==0)   cfg->nphoton=FIND_JSON_KEY("Photons","Session.Photons",Session,0,valueint);
        if(cfg->session[0]=='\0') strncpy(cfg->session, FIND_JSON_KEY("ID","Session.ID",Session,"default",valuestring), MAX_SESSION_LENGTH);

        if(!cfg->isreflect)   cfg->isreflect=FIND_JSON_KEY("DoMismatch","Session.DoMismatch",Session,cfg->isreflect,valueint);
        if(cfg->issave2pt)    cfg->issave2pt=FIND_JSON_KEY("DoSaveVolume","Session.DoSaveVolume",Session,cfg->issave2pt,valueint);
        if(cfg->isnormalized) cfg->isnormalized=FIND_JSON_KEY("DoNormalize","Session.DoNormalize",Session,cfg->isnormalized,valueint);
        if(!cfg->issavedet)   cfg->issavedet=FIND_JSON_KEY("DoPartialPath","Session.DoPartialPath",Session,cfg->issavedet,valueint);
        if(!cfg->isspecular)  cfg->isspecular=FIND_JSON_KEY("DoSpecular","Session.DoSpecular",Session,cfg->isspecular,valueint);
        if(!cfg->ismomentum)  cfg->ismomentum=FIND_JSON_KEY("DoDCS","Session.DoDCS",Session,cfg->ismomentum,valueint);
        if(!cfg->issaveexit)  cfg->issaveexit=FIND_JSON_KEY("DoSaveExit","Session.DoSaveExit",Session,cfg->issaveexit,valueint);
        if(!cfg->issaveseed)  cfg->issaveseed=FIND_JSON_KEY("DoSaveSeed","Session.DoSaveSeed",Session,cfg->issaveseed,valueint);
        if(cfg->basisorder)   cfg->basisorder=FIND_JSON_KEY("BasisOrder","Session.BasisOrder",Session,cfg->basisorder,valueint);
        if(!cfg->outputformat)  cfg->outputformat=mcx_keylookup((char *)FIND_JSON_KEY("OutputFormat","Session.OutputFormat",Session,"ascii",valuestring),outputformat);
        if(cfg->outputformat<0)
		MMC_ERROR(-2,"the specified output format is not recognized");

        if(cfg->debuglevel==0)
           cfg->debuglevel=mcx_parsedebugopt((char *)FIND_JSON_KEY("DebugFlag","Session.DebugFlag",Session,"",valuestring));
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
     if(cfg->e0==0)
	 MMC_ERROR(-1,"InitElem must be given");
     return 0;
}

/**
 * @brief Write simulation settings to an inp file
 *
 * @param[in] fname: the name of the output file
 * @param[in] cfg: simulation configuration
 */

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

/**
 * @brief Load user inputs from a .inp input file
 *
 * This function loads user input from a simple text input format in a .inp extension
 *
 * @param[in] in: file handle to the .inp file
 * @param[in] cfg: simulation configuration
 */

void mcx_loadconfig(FILE *in, mcconfig *cfg){
     int i,gates,srctype,itmp;
     float dtmp;
     char comment[MAX_PATH_LENGTH],*comm, srctypestr[MAX_SESSION_LENGTH]={'\0'};
     
     if(in==stdin)
     	MMC_FPRINTF(stdout,"Please specify the total number of photons: [1000000]\n\t");
     MMC_ASSERT(fscanf(in,"%d", &(i) )==1); 
     if(cfg->nphoton==0) cfg->nphoton=i;
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     
     if(in==stdin)
     	MMC_FPRINTF(stdout,">> %d\nPlease specify the random number generator seed: [123456789]\n\t",cfg->nphoton);
     if(cfg->seed==0x623F9A9E)
        MMC_ASSERT(fscanf(in,"%d", &(cfg->seed) )==1);
     else
        MMC_ASSERT(fscanf(in,"%d", &itmp )==1);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	MMC_FPRINTF(stdout,">> %d\nPlease specify the position of the source: [10 10 5]\n\t",cfg->seed);
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->srcpos.x),&(cfg->srcpos.y),&(cfg->srcpos.z) )==3);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	MMC_FPRINTF(stdout,">> %f %f %f\nPlease specify the normal direction of the source: [0 0 1]\n\t",
	                            cfg->srcpos.x,cfg->srcpos.y,cfg->srcpos.z);
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->srcdir.x),&(cfg->srcdir.y),&(cfg->srcdir.z)));
     comm=fgets(comment,MAX_PATH_LENGTH,in);
     if(comm!=NULL && sscanf(comm,"%f",&dtmp)==1)
         cfg->srcdir.w=dtmp;

     if(in==stdin)
        MMC_FPRINTF(stdout,">> %f %f %f %f\nPlease specify the time gates in seconds (start end step) [0.0 1e-9 1e-10]\n\t",
	                            cfg->srcdir.x,cfg->srcdir.y,cfg->srcdir.z,cfg->srcdir.w);
     MMC_ASSERT(fscanf(in,"%f %f %f", &(cfg->tstart),&(cfg->tend),&(cfg->tstep) )==3);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
        MMC_FPRINTF(stdout,">> %f %f %f\nPlease specify the mesh file key {node,elem,velem,facenb}_key.dat :\n\t",
                                    cfg->tstart,cfg->tend,cfg->tstep);
     if(cfg->tstart>cfg->tend || cfg->tstep==0.f){
         MMC_ERROR(-9,"incorrect time gate settings");
     }
     if(cfg->tstep>cfg->tend-cfg->tstart){
         cfg->tstep=cfg->tend-cfg->tstart;
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
     	MMC_FPRINTF(stdout,">> %s\nPlease specify the index to the tetrahedral element enclosing the source [start from 1]:\n\t",cfg->meshtag);
     MMC_ASSERT(fscanf(in,"%d", &(cfg->e0))==1);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	MMC_FPRINTF(stdout,">> %d\nPlease specify the total number of detectors and detector diameter (in mm):\n\t",cfg->e0);
     MMC_ASSERT(fscanf(in,"%d %f", &(cfg->detnum), &(cfg->detradius))==2);
     comm=fgets(comment,MAX_PATH_LENGTH,in);

     if(in==stdin)
     	MMC_FPRINTF(stdout,">> %d %f\n",cfg->detnum,cfg->detradius);
     cfg->detpos=(float4*)malloc(sizeof(float4)*cfg->detnum);
     if(cfg->issavedet)
        cfg->issavedet=(cfg->detpos>0);
     for(i=0;i<cfg->detnum;i++){
        if(in==stdin)
		MMC_FPRINTF(stdout,"Please define detector #%d: x,y,z (in mm): [5 5 5 1]\n\t",i);
     	MMC_ASSERT(fscanf(in, "%f %f %f", &(cfg->detpos[i].x),&(cfg->detpos[i].y),&(cfg->detpos[i].z))==3);
        comm=fgets(comment,MAX_PATH_LENGTH,in);
        if(comm!=NULL && sscanf(comm,"%f",&dtmp)==1)
            cfg->detpos[i].w=dtmp;
        else
            cfg->detpos[i].w=cfg->detradius;

        if(in==stdin)
		MMC_FPRINTF(stdout,">> %f %f %f\n",cfg->detpos[i].x,cfg->detpos[i].y,cfg->detpos[i].z);
     }

     if(in==stdin)
        MMC_FPRINTF(stdout,"Please specify the source type [pencil|isotropic|cone|gaussian|planar|pattern|fourier|arcsine|disk|fourierx|fourierx2d|zgaussian|line|slit]:\n\t");
     if(fscanf(in,"%s", srctypestr)==1 && srctypestr[0]){
        srctype=mcx_keylookup(srctypestr,srctypeid);
	if(srctype==-1)
	   MMC_ERROR(-6,"the specified source type is not supported");
        if(srctype>=0){
           comm=fgets(comment,MAX_PATH_LENGTH,in);
	   cfg->srctype=srctype;
	   if(in==stdin)
                MMC_FPRINTF(stdout,">> %d\nPlease specify the source parameters set 1 (4 floating-points):\n\t",cfg->srctype);
           MMC_ASSERT(fscanf(in, "%f %f %f %f", &(cfg->srcparam1.x),&(cfg->srcparam1.y),&(cfg->srcparam1.z),&(cfg->srcparam1.w))==4);
           comm=fgets(comment,MAX_PATH_LENGTH,in);
           if(in==stdin)
		MMC_FPRINTF(stdout,">> %f %f %f %f\nPlease specify the source parameters set 2 (4 floating-points):\n\t",
                       cfg->srcparam1.x,cfg->srcparam1.y,cfg->srcparam1.z,cfg->srcparam1.w);
           if(fscanf(in, "%f %f %f %f", &(cfg->srcparam2.x),&(cfg->srcparam2.y),&(cfg->srcparam2.z),&(cfg->srcparam2.w))==4){
               comm=fgets(comment,MAX_PATH_LENGTH,in);
               if(in==stdin)
	            MMC_FPRINTF(stdout,">> %f %f %f %f\n",cfg->srcparam2.x,cfg->srcparam2.y,cfg->srcparam2.z,cfg->srcparam2.w);
               if(cfg->srctype==stPattern && cfg->srcparam1.w*cfg->srcparam2.w>0){
		    char srcpatternfile[MAX_PATH_LENGTH];
		    FILE *fp;
		    if(in==stdin)
                    	MMC_FPRINTF(stdout,"Please specify the source pattern file name:\n\t");
		    if(cfg->srcpattern) free(cfg->srcpattern);
		    cfg->srcpattern=(float*)calloc((cfg->srcparam1.w*cfg->srcparam2.w),sizeof(float));
		    MMC_ASSERT(fscanf(in, "%s", srcpatternfile)==1);
		    comm=fgets(comment,MAX_PATH_LENGTH,in);
		    fp=fopen(srcpatternfile,"rb");
		    if(fp==NULL)	MMC_ERROR(-6,"source pattern file can not be opened");
		    MMC_ASSERT(fread(cfg->srcpattern,cfg->srcparam1.w*cfg->srcparam2.w,sizeof(float),fp)==sizeof(float));
		    fclose(fp);
		}
	    }
		if(cfg->detnum==1 && cfg->detpos[0].w==0.0){
		// only one detector and its radius is 0, indicates that we are using a wide-field detector
			if(in==stdin)
				MMC_FPRINTF(stdout,">> \nPlease specify the detector parameters set 1 (4 floating-points):\n\t");
			MMC_ASSERT(fscanf(in, "%f %f %f %f", &(cfg->detparam1.x),&(cfg->detparam1.y),&(cfg->detparam1.z),&(cfg->detparam1.w))==4);
			comm=fgets(comment,MAX_PATH_LENGTH,in);
			if(in==stdin)
		            MMC_FPRINTF(stdout,">> %f %f %f %f\nPlease specify the detector parameters set 2 (4 floating-points):\n\t",
		                   cfg->detparam1.x,cfg->detparam1.y,cfg->detparam1.z,cfg->detparam1.w);
			MMC_ASSERT(fscanf(in, "%f %f %f %f", &(cfg->detparam2.x),&(cfg->detparam2.y),&(cfg->detparam2.z),&(cfg->detparam2.w))==4);
		    comm=fgets(comment,MAX_PATH_LENGTH,in);
			if(in==stdin)
			        MMC_FPRINTF(stdout,">> %f %f %f %f\n",cfg->detparam2.x,cfg->detparam2.y,cfg->detparam2.z,cfg->detparam2.w);
			// only load detection pattern under replay mode
			if(cfg->seed==SEED_FROM_FILE && (cfg->outputtype==otWL || cfg->outputtype==otWP) && cfg->detparam1.w*cfg->detparam2.w>0){
				if(in==stdin)	MMC_FPRINTF(stdout,"Please specify the detector pattern file name:\n\t");
				char detpatternfile[MAX_PATH_LENGTH];
		    	FILE *fp;
				if(cfg->detpattern) free(cfg->detpattern);
		    	cfg->detpattern=(float*)calloc((cfg->detparam1.w*cfg->detparam2.w),sizeof(float));
				MMC_ASSERT(fscanf(in, "%s", detpatternfile)==1);
				comm=fgets(comment,MAX_PATH_LENGTH,in);
				fp=fopen(detpatternfile,"rb");
		    	if(fp==NULL)	MMC_ERROR(-6,"detector pattern file can not be opened");
		    	MMC_ASSERT(fread(cfg->detpattern,cfg->detparam1.w*cfg->detparam2.w,sizeof(float),fp)==sizeof(float));
		    	fclose(fp);
			}
		}
	}else
	   return;
     }else
        return;
}

/**
 * @brief Save simulation settings to an inp file
 *
 * @param[in] out: handle to the output file
 * @param[in] cfg: simulation configuration
 */

void mcx_saveconfig(FILE *out, mcconfig *cfg){
     int i;

     MMC_FPRINTF(out,"%d\n", (cfg->nphoton) ); 
     MMC_FPRINTF(out,"%d\n", (cfg->seed) );
     MMC_FPRINTF(out,"%f %f %f\n", (cfg->srcpos.x),(cfg->srcpos.y),(cfg->srcpos.z) );
     MMC_FPRINTF(out,"%f %f %f\n", (cfg->srcdir.x),(cfg->srcdir.y),(cfg->srcdir.z) );
     MMC_FPRINTF(out,"%f %f %f\n", (cfg->tstart),(cfg->tend),(cfg->tstep) );
     MMC_FPRINTF(out,"%f %d %d %d\n", (cfg->steps.x),(cfg->dim.x),(cfg->crop0.x),(cfg->crop1.x));
     MMC_FPRINTF(out,"%f %d %d %d\n", (cfg->steps.y),(cfg->dim.y),(cfg->crop0.y),(cfg->crop1.y));
     MMC_FPRINTF(out,"%f %d %d %d\n", (cfg->steps.z),(cfg->dim.z),(cfg->crop0.z),(cfg->crop1.z));
     MMC_FPRINTF(out,"%d", (cfg->medianum));
     for(i=0;i<cfg->medianum;i++){
     	MMC_FPRINTF(out, "%f %f %f %f\n", (cfg->prop[i].mus),(cfg->prop[i].g),(cfg->prop[i].mua),(cfg->prop[i].n));
     }
     MMC_FPRINTF(out,"%d", (cfg->detnum));
     for(i=0;i<cfg->detnum;i++){
     	MMC_FPRINTF(out, "%f %f %f %f\n", (cfg->detpos[i].x),(cfg->detpos[i].y),(cfg->detpos[i].z),(cfg->detpos[i].w));
     }
}

/**
 * @brief Load media index data volume (.bin or .vol) to the memory (not used in MMC)
 *
 * @param[in] filename: file name to the binary volume data (support 1,2 and 4 bytes per voxel)
 * @param[in] cfg: simulation configuration
 */

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

/**
 * @brief Parse the debug flag in the letter format
 *
 * The debug flag following the -D can be either a string format, or numerical format.
 * This function converts the string debug flags into number format
 *
 * @param[in] debugopt: string following the -D parameter
 * @return debugflag: the numerical format of the debug flag
 */

int mcx_parsedebugopt(char *debugopt){
    char *c=debugopt,*p;
    int debuglevel=0;

    while(*c){
       p=strchr((char *)debugflag, ((*c<='z' && *c>='a') ? *c-'a'+'A' : *c) );
       if(p!=NULL)
	  debuglevel |= (1 << (p-debugflag));
       c++;
    }
    return debuglevel;
}

/**
 * @brief Print a progress bar
 *
 * When -D P is specified, this function prints and update a progress bar.
 *
 * @param[in] n: the number of completed photons
 * @param[in] cfg: simulation configuration
 */

void mcx_progressbar(unsigned int n, mcconfig *cfg){
    unsigned int percentage, j,colwidth=79;
    static unsigned int oldmarker=0xFFFFFFFF;

#ifndef MCX_CONTAINER
  #ifdef TIOCGWINSZ
    struct winsize ttys={0,0,0,0};
    ioctl(0, TIOCGWINSZ, &ttys);
    colwidth=ttys.ws_col;
    if(colwidth==0)
          colwidth=79;
  #endif
#endif

    percentage=(float)n*(colwidth-18)/cfg->nphoton;

    if(percentage != oldmarker){
        oldmarker=percentage;
	for(j=0;j<colwidth;j++)     MMC_FPRINTF(stdout,"\b");
    	MMC_FPRINTF(stdout,"Progress: [");
    	for(j=0;j<percentage;j++)      MMC_FPRINTF(stdout,"=");
    	MMC_FPRINTF(stdout,(percentage<colwidth-18) ? ">" : "=");
    	for(j=percentage;j<colwidth-18;j++) MMC_FPRINTF(stdout," ");
    	MMC_FPRINTF(stdout,"] %3d%%",percentage*100/(colwidth-18));
#ifdef MCX_CONTAINER
        mcx_matlab_flush();
#else
        fflush(stdout);
#endif
    }
}

/**
 * @brief Function to read a single parameter value followed by a command line option
 *
 * This function reads different types of parameter values following a command line option.
 *
 * @param[in] argc: the number of total command line parameters
 * @param[in] argv: the pointer to all command line options
 * @param[in] id: which parameter to be parsed
 * @param[out] output: the pointer to which the parsed value to be written
 * @param[in] type: the type of data support char, int, float, string, bytenumlist, floatlist
 */

int mcx_readarg(int argc, char *argv[], int id, void *output,const char *type){
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

/**
 * @brief Test if a long command line option is supported
 *
 * This function returns 1 if a long option is found, and 0 otherwise
 *
 * @param[in] opt: the long command line option string
 */

int mcx_remap(char *opt){
    int i=0;
    while(shortopt[i]!='\0'){
	if(strcmp(opt,fullopt[i])==0){
		opt[1]=shortopt[i];
		if(shortopt[i]!='-')
		    opt[2]='\0';
		return 0;
	}
	i++;
    }
    return 1;
}

/**
 * @brief Look up a single character in a string
 *
 * @param[in] key: character to be looked up
 * @param[out] index: the dictionary string where the char is searched
 * @return if found, return 0; otherwise, return 1
 */

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

/**
 * @brief Look up a string in a string list and return the index
 *
 * @param[in] key: string to be looked up
 * @param[out] table: the dictionary where the string is searched
 * @return if found, return the index of the string in the dictionary, otherwise -1.
 */

int mcx_keylookup(char *key, const char *table[]){
    int i=0;
    while(key[i]){
        if(key[i]>='A' && key[i]<='Z')
		key[i]+=('a'-'A');
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

/** 
 * @brief Validate all input fields, and warn incompatible inputs
 *
 * Perform self-checking and raise exceptions or warnings when input error is detected
 *
 * @param[in,out] cfg: the simulation configuration structure
 */

void mcx_validatecfg(mcconfig *cfg){
     if(cfg->nphoton<=0){
         MMC_ERROR(-2,"cfg.nphoton must be a positive number");
     }
     if(cfg->tstart>cfg->tend || cfg->tstep==0.f){
         MMC_ERROR(-2,"incorrect time gate settings or missing tstart/tend/tstep fields");
     }
     if(cfg->tstep>cfg->tend-cfg->tstart){
         cfg->tstep=cfg->tend-cfg->tstart;
     }
     if(fabs(cfg->srcdir.x*cfg->srcdir.x+cfg->srcdir.y*cfg->srcdir.y+cfg->srcdir.z*cfg->srcdir.z - 1.f)>1e-4)
         MMC_ERROR(-2,"field 'srcdir' must be a unitary vector (tolerance is 1e-4)");
     if(cfg->tend<=cfg->tstart)
         MMC_ERROR(-2,"field 'tend' must be greater than field 'tstart'");
     cfg->maxgate=(int)((cfg->tend-cfg->tstart)/cfg->tstep+0.5);

     if(cfg->srctype==stPattern && cfg->srcpattern==NULL)
        MMC_ERROR(-2,"the 'srcpattern' field can not be empty when your 'srctype' is 'pattern'");

     if(cfg->seed<0 && cfg->seed!=SEED_FROM_FILE)
        cfg->seed=time(NULL);
     if(cfg->method==rtBLBadouelGrid){
	cfg->basisorder=0;
     }
}

/**
 * @brief Preprocess configuration and set option dependency
 *
 * This function preprocess the user input and set dependent flags
 *
 * @param[in,out] cfg: simulation configuration
 */

void mcx_prep(mcconfig *cfg){
     if(cfg->issavedet && cfg->detnum==0 && cfg->isextdet==0) 
      	cfg->issavedet=0;
     if(cfg->issavedet==0){
        cfg->ismomentum=0;
        cfg->issaveexit=0;
     }
}

/**
 * @brief Main function to read user command line options
 *
 * This function process user command line inputs and parse all short and long options.
 *
 * @param[in] argc: the number of total command line parameters
 * @param[in] argv: the pointer to all command line options
 * @param[in] cfg: simulation configuration
 */

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
                                MMC_FPRINTF(cfg->flog,"option: %s\n",argv[i]);
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
		     case 'q':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->issaveseed),"bool");
		     	        break;
		     case 'g':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->maxgate),"int");
		     	        break;
		     case 'b':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isreflect),"bool");
		     	        break;
		     case 'd':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->issavedet),"bool");
		     	        break;
		     case 'm':
		                i=mcx_readarg(argc,argv,i,&(cfg->mcmethod),"int");
				break;
		     case 'x':
		                i=mcx_readarg(argc,argv,i,&(cfg->issaveexit),"bool");
				if (cfg->issaveexit) cfg->issavedet=1;
				break;
		     case 'C':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->basisorder),"bool");
		     	        break;
		     case 'V':
		     	        i=mcx_readarg(argc,argv,i,&(cfg->isspecular),"bool");
		     	        break;
		     case 'v':
                                mcx_version(cfg);
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
				if(i+1<argc && strstr(argv[i+1],".mch")!=NULL){ /*give an mch file to initialize the seed*/
#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
					MMC_ERROR(-1,"seeding file is not supported in this binary");
#else
                                        i=mcx_readarg(argc,argv,i,cfg->seedfile,"string");
					cfg->seed=SEED_FROM_FILE;
#endif
		     	        }else
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
                     case 'P':
                                i=mcx_readarg(argc,argv,i,&(cfg->replaydet),"int");
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
                     case 'k':
                                i=mcx_readarg(argc,argv,i,&(cfg->voidtime),"int");
                                break;
                     case '-':  /*additional verbose parameters*/
                                if(strcmp(argv[i]+2,"momentum")==0){
		                     i=mcx_readarg(argc,argv,i,&(cfg->ismomentum),"bool");
                                     if (cfg->ismomentum) cfg->issavedet=1;
                                }else if(strcmp(argv[i]+2,"atomic")==0){
		                     i=mcx_readarg(argc,argv,i,&(cfg->isatomic),"bool");
                                }else if(strcmp(argv[i]+2,"debugphoton")==0){
		                     i=mcx_readarg(argc,argv,i,&(cfg->debugphoton),"int");
                                }else
                                     MMC_FPRINTF(cfg->flog,"unknown verbose option: --%s\n",argv[i]+2);
                                break;
                     default:{
                                MMC_FPRINTF(cfg->flog,"option: %s\n",argv[i]);
				MMC_ERROR(-1,"unsupported command line option");
                    }
		}
	    }
	    i++;
     }
     if(issavelog && cfg->session){
          sprintf(logfile,"%s.log",cfg->session);
          cfg->flog=fopen(logfile,"wt");
          if(cfg->flog==NULL){
		cfg->flog=stdout;
		MMC_FPRINTF(cfg->flog,"unable to save to log file, will print from stdout\n");
          }
     }
     if((cfg->outputtype==otJacobian || cfg->outputtype==otWL || cfg->outputtype==otWP) && cfg->seed!=SEED_FROM_FILE)
         MMC_ERROR(-1,"Jacobian output is only valid in the reply mode. Please give an mch file after '-E'.");
     if(cfg->isgpuinfo!=2){ /*print gpu info only*/
       if(isinteractive){
          mcx_readconfig("",cfg);
       }else{
     	  mcx_readconfig(filename,cfg);
       }
     }
     mcx_validatecfg(cfg);
}

/**
 * @brief Print MCX software version
 *
 * @param[in] cfg: simulation configuration
 */

void mcx_version(mcconfig *cfg){
    MMC_ERROR(MMC_INFO,"MMC $Rev::      $");
}

/**
 * @brief Print MCX help information
 *
 * @param[in] exename: path and name of the mcx executable
 */

void mcx_usage(char *exename){
     printf("\
###############################################################################\n\
#                         Mesh-based Monte Carlo (MMC)                        #\n\
#          Copyright (c) 2010-2018 Qianqian Fang <q.fang at neu.edu>          #\n\
#                            http://mcx.space/#mmc                            #\n\
#                                                                             #\n\
#Computational Optics & Translational Imaging (COTI) Lab  [http://fanglab.org]#\n\
#            Department of Bioengineering, Northeastern University            #\n\
#                                                                             #\n\
#                Research funded by NIH/NIGMS grant R01-GM114365              #\n\
###############################################################################\n\
$Rev::       $ Last $Date::                       $ by $Author::              $\n\
###############################################################################\n\
\n\
usage: %s <param1> <param2> ...\n\
where possible parameters include (the first item in [] is the default value)\n\
\n\
== Required option ==\n\
 -f config     (--input)       read an input file in .inp or .json format\n\
\n\
== MC options ==\n\
 -n [0.|float] (--photon)      total photon number, max allowed value is 2^32-1\n\
 -b [0|1]      (--reflect)     1 do reflection at int&ext boundaries, 0 no ref.\n\
 -U [1|0]      (--normalize)   1 to normalize the fluence to unitary,0 save raw\n\
 -m [0|1]      (--mc)          0 use MCX-styled MC method, 1 use MCML style MC\n\
 -C [1|0]      (--basisorder)  1 piece-wise-linear basis for fluence,0 constant\n\
 -u [1.|float] (--unitinmm)    define the mesh data length unit in mm\n\
 -E [1648335518|int|mch](--seed) set random-number-generator seed;\n\
                               if an mch file is followed, MMC \"replays\" \n\
                               the detected photons; the replay mode can be used\n\
                               to calculate the mua/mus Jacobian matrices\n\
 -P [0|int]    (--replaydet)   replay only the detected photons from a given \n\
                               detector (det ID starts from 1), use with -E \n\
 -M [%c|PHBSG] (--method)      choose ray-tracing algorithm (only use 1 letter)\n\
                               P - Plucker-coordinate ray-tracing algorithm\n\
			       H - Havel's SSE4 ray-tracing algorithm\n\
			       B - partial Badouel's method (used by TIM-OS)\n\
			       S - branch-less Badouel's method with SSE\n\
			       G - dual-grid MMC (DMMC) with voxel data output\n\
 -e [1e-6|float](--minenergy)  minimum energy level to trigger Russian roulette\n\
 -V [0|1]      (--specular)    1 source located in the background,0 inside mesh\n\
 -k [1|0]      (--voidtime)    when src is outside, 1 enables timer inside void\n\
 --atomic [1|0]                1 use atomic operations, 0 use non-atomic ones\n\
\n\
== Output options ==\n\
 -O [X|XFEJLP] (--outputtype)  X - output flux, F - fluence, E - energy deposit\n\
                               J - Jacobian, L - weighted path length, P -\n\
                               weighted scattering count (J,L,P: replay mode)\n\
 -s sessionid  (--session)     a string used to tag all output file names\n\
 -S [1|0]      (--save2pt)     1 to save the fluence field, 0 do not save\n\
 -d [0|1]      (--savedet)     1 to save photon info at detectors,0 not to save\n\
 -x [0|1]      (--saveexit)    1 to save photon exit positions and directions\n\
                               setting -x to 1 also implies setting '-d' to 1\n\
 -q [0|1]      (--saveseed)    1 save RNG seeds of detected photons for replay\n\
 -F format     (--outputformat)'ascii', 'bin' (in 'double'), 'json' or 'ubjson'\n\
\n\
== User IO options ==\n\
 -i 	       (--interactive) interactive mode\n\
 -h            (--help)        print this message\n\
 -v            (--version)     print MMC version information\n\
 -l            (--log)         print messages to a log file instead\n\
\n\
== Debug options ==\n\
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
 --debugphoton [-1|int]        to print the debug info specified by -D only for\n\
                               a single photon, followed by its index (start 0)\n\
\n\
== Additional options ==\n\
 --momentum     [0|1]          1 to save photon momentum transfer,0 not to save\n\
\n\
== Example ==\n\
       %s -n 1000000 -f input.json -s test -b 0 -D TP\n",exename,
#ifdef MMC_USE_SSE
'H',
#else
'P',
#endif
exename);
}
