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
\file    simpmesh.c

\brief   Basic vector math and mesh operations
*******************************************************************************/

#include <stdlib.h>
#include "simpmesh.h"
#include <string.h>
#include "highordermesh.h"

#ifdef WIN32
         char pathsep='\\';  /**< path separator on Windows */
#else
         char pathsep='/';   /**< path separator on Linux/Unix/OSX */
#endif

/** 
 * \brief Tetrahedron faces, in counter-clock-wise orders, represented using local node indices
 *
 * node-connectivity, i.e. nc[4] points to the 4 facets of a tetrahedron, with each
 * triangular face made of 3 nodes. The numbers [0-4] are the
 * local node indices (starting from 0). The order of the nodes
 * are in counter-clock-wise orders.
 */

const int out[4][3]={{0,3,1},{3,2,1},{0,2,3},{0,1,2}};

/** 
 * \brief The local index of the node with an opposite face to the i-th face defined in nc[][]
 *
 * nc[i] <-> node[facemap[i]]
 * the 1st face of this tet, i.e. nc[0]={3,0,1}, is opposite to the 3rd node
 * the 2nd face of this tet, i.e. nc[1]={3,1,2}, is opposite to the 1st node
 * etc.
 */

const int facemap[]={2,0,1,3};

/** 
 * \brief Inverse mapping between the local index of the node and the corresponding opposite face in nc[]'s order
 *
 * nc[ifacemap[i]] <-> node[i]
 * the 1st node of this tet is in opposite to the 2nd face, i.e. nc[1]={3,1,2}
 * the 2nd node of this tet is in opposite to the 3rd face, i.e. nc[1]={2,0,3}
 * etc.
 */

const int ifacemap[]={1,2,0,3};

/**
 * @brief Initializing the mesh data structure with default values
 *
 * Constructor of the mesh object, initializing all field to default values
 */


void mesh_init(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->nf=0;
	mesh->prop=0;
	mesh->elemlen=4;
	mesh->node=NULL;
	mesh->elem=NULL;
	mesh->elem2=NULL;
	mesh->srcelemlen=0;
	mesh->srcelem=NULL;
	mesh->detelemlen=0;
	mesh->detelem=NULL;
	mesh->facenb=NULL;
	mesh->type=NULL;
	mesh->med=NULL;
	mesh->atte=NULL;
	mesh->weight=NULL;
	mesh->evol=NULL;
	mesh->nvol=NULL;
	mesh->dref=NULL;
        mesh->nmin.x=VERY_BIG;
	mesh->nmin.y=VERY_BIG;
	mesh->nmin.z=VERY_BIG;
	mesh->nmin.w=1.f;
        mesh->nmax.x=-VERY_BIG;
	mesh->nmax.y=-VERY_BIG;
	mesh->nmax.z=-VERY_BIG;
	mesh->nmax.w=1.f;
}

/**
 * @brief Loading user-specified mesh data
 *
 * Loading node, element etc from files into memory
 */

void mesh_init_from_cfg(tetmesh *mesh,mcconfig *cfg){
        mesh_init(mesh);
        mesh_loadnode(mesh,cfg);
        mesh_loadelem(mesh,cfg);
	if(cfg->basisorder==2)
	  mesh_10nodetet(mesh,cfg);
        mesh_loadfaceneighbor(mesh,cfg);
        mesh_loadmedia(mesh,cfg);
        mesh_loadelemvol(mesh,cfg);
	if(cfg->seed==SEED_FROM_FILE && cfg->seedfile[0]){
          mesh_loadseedfile(mesh,cfg);
        }
}

/**
 * @brief Error-handling in mesh operations
 *
 * @param[in] msg: the error message string
 * @param[in] file: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

void mesh_error(const char *msg,const char *file,const int linenum){
#ifdef MCX_CONTAINER
        mmc_throw_exception(1,msg,file,linenum);
#else
	fprintf(stderr,"Mesh error: %s in unit %s line#%d\n",msg,file,linenum);
        exit(1);
#endif
}

/**
 * @brief Construct a full mesh file name using cfg session and root path
 *
 * @param[in] format: a format string to form the file name
 * @param[in] foutput: pointer to the output string buffer
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_filenames(const char *format,char *foutput,mcconfig *cfg){
	char filename[MAX_PATH_LENGTH];
	sprintf(filename,format,cfg->meshtag);

	if(cfg->rootpath[0]) 
		sprintf(foutput,"%s%c%s",cfg->rootpath,pathsep,filename);
	else
		sprintf(foutput,"%s",filename);
}

/**
 * @brief Load node file and initialize the related mesh properties
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadnode(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int tmp,len,i;
	char fnode[MAX_PATH_LENGTH];
	mesh_filenames("node_%s.dat",fnode,cfg);
	if((fp=fopen(fnode,"rt"))==NULL){
		MESH_ERROR("can not open node file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->nn));
	if(len!=2 || mesh->nn<=0){
		MESH_ERROR("node file has wrong format");
	}
	mesh->node=(float3 *)calloc(sizeof(float3),mesh->nn);
	for(i=0;i<mesh->nn;i++){
		if(fscanf(fp,"%d %f %f %f",&tmp,&(mesh->node[i].x),&(mesh->node[i].y),&(mesh->node[i].z))!=4)
			MESH_ERROR("node file has wrong format");
	}
	fclose(fp);
        if(cfg->method==rtBLBadouelGrid)
	        mesh_createdualmesh(mesh,cfg);
}

void mesh_createdualmesh(tetmesh *mesh,mcconfig *cfg){
        int i;
        mesh->nmin.x=VERY_BIG;
	mesh->nmin.y=VERY_BIG;
	mesh->nmin.z=VERY_BIG;
        mesh->nmax.x=-VERY_BIG;
	mesh->nmax.y=-VERY_BIG;
	mesh->nmax.z=-VERY_BIG;

	for(i=0;i<mesh->nn;i++){
		mesh->nmin.x=MIN(mesh->node[i].x, mesh->nmin.x);
		mesh->nmin.y=MIN(mesh->node[i].y, mesh->nmin.y);
		mesh->nmin.z=MIN(mesh->node[i].z, mesh->nmin.z);
		mesh->nmax.x=MAX(mesh->node[i].x, mesh->nmax.x);
		mesh->nmax.y=MAX(mesh->node[i].y, mesh->nmax.y);
		mesh->nmax.z=MAX(mesh->node[i].z, mesh->nmax.z);
	}
        mesh->nmin.x-=EPS;
	mesh->nmin.y-=EPS;
	mesh->nmin.z-=EPS;
        mesh->nmax.x+=EPS;
	mesh->nmax.y+=EPS;
	mesh->nmax.z+=EPS;

	cfg->dim.x=(int)((mesh->nmax.x-mesh->nmin.x)/cfg->unitinmm)+1;
	cfg->dim.y=(int)((mesh->nmax.y-mesh->nmin.y)/cfg->unitinmm)+1;
	cfg->dim.z=(int)((mesh->nmax.z-mesh->nmin.z)/cfg->unitinmm)+1;

	cfg->crop0.x=cfg->dim.x;
        cfg->crop0.y=cfg->dim.y*cfg->dim.x;
	cfg->crop0.z=cfg->dim.y*cfg->dim.x*cfg->dim.z;
}

/**
 * @brief Load optical property file and initialize the related mesh properties
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadmedia(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int tmp,len,i;
	char fmed[MAX_PATH_LENGTH];
	mesh_filenames("prop_%s.dat",fmed,cfg);
	if((fp=fopen(fmed,"rt"))==NULL){
		MESH_ERROR("can not open media property file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->prop));
	if(len!=2 || mesh->prop<=0){
		MESH_ERROR("property file has wrong format");
	}
        /*when there is an external detector, reindex the property to medianum+1*/
	mesh->med=(medium *)calloc(sizeof(medium),mesh->prop+1+cfg->isextdet);
        mesh->atte=(float *)calloc(sizeof(float),mesh->prop+1+cfg->isextdet);

	mesh->med[0].mua=0.f;
	mesh->med[0].mus=0.f;
	mesh->med[0].n=cfg->nout;
	mesh->med[0].g=1.f;

        /*make medianum+1 the same as medium 0*/
        if(cfg->isextdet){
		memcpy(mesh->med+mesh->prop+1,mesh->med,sizeof(medium));
                for(i=0;i<mesh->ne;i++){
                     if(mesh->type[i]==-2)
                           mesh->type[i]=mesh->prop+1;
                }
        }
	for(i=1;i<=mesh->prop;i++){
		if(fscanf(fp,"%d %f %f %f %f",&tmp,&(mesh->med[i].mua),&(mesh->med[i].mus),
		                                   &(mesh->med[i].g),&(mesh->med[i].n))!=5)
			MESH_ERROR("property file has wrong format");
		/*mesh->atte[i]=expf(-cfg->minstep*mesh->med[i].mua);*/
	}
	fclose(fp);

        if(cfg->method!=rtBLBadouelGrid && cfg->unitinmm!=1.f){
           for(i=1;i<=mesh->prop;i++){
                   mesh->med[i].mus*=cfg->unitinmm;
                   mesh->med[i].mua*=cfg->unitinmm;
           }
        }
	cfg->his.maxmedia=mesh->prop; /*skip media 0*/
}

/**
 * @brief Load element file and initialize the related mesh properties
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadelem(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int tmp,len,i,j,datalen;
	int *pe;
	char felem[MAX_PATH_LENGTH];

	mesh_filenames("elem_%s.dat",felem,cfg);
	if((fp=fopen(felem,"rt"))==NULL){
		MESH_ERROR("can not open element file");
	}
	len=fscanf(fp,"%d %d",&(mesh->elemlen),&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		MESH_ERROR("element file has wrong format");
	}
	if(mesh->elemlen<4)
	    mesh->elemlen=4;

	mesh->elem=(int *)malloc(sizeof(int)*mesh->elemlen*mesh->ne);
	mesh->type=(int *)malloc(sizeof(int )*mesh->ne);
	
	datalen=(cfg->method==rtBLBadouelGrid) ? cfg->crop0.z : ( (cfg->basisorder) ? mesh->nn : mesh->ne);
	mesh->weight=(double *)calloc(sizeof(double)*datalen,cfg->maxgate*cfg->srcnum);

	for(i=0;i<mesh->ne;i++){
		pe=mesh->elem+i*mesh->elemlen;
		if(fscanf(fp,"%d",&tmp)!=1)
			break;
		for(j=0;j<mesh->elemlen;j++)
		        if(fscanf(fp,"%d",pe+j)!=1)
			    break;
		if(fscanf(fp,"%d",mesh->type+i)!=1)
			break;
        }
	fclose(fp);
	if(i<mesh->ne)
	    MESH_ERROR("element file has wrong format");
	mesh_srcdetelem(mesh,cfg);
}

/**
 * @brief Identify wide-field source and detector-related elements (type=-1 for source, type=-2 for det)
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_srcdetelem(tetmesh *mesh,mcconfig *cfg){
	int i;

        mesh->srcelemlen=0;
        mesh->detelemlen=0;
	for(i=0;i<mesh->ne;i++){
		if(mesh->type[i]==-1){	/*number of elements in the initial candidate list*/
			mesh->srcelemlen++;
			cfg->e0=(cfg->e0==0) ? i+1 : cfg->e0;
		}
		if(mesh->type[i]==-2){	/*number of elements in the initial candidate list*/
			mesh->detelemlen++;
			cfg->isextdet=1;
			cfg->detnum=0;	// when detecting wide-field detectors, suppress point detectors
                }
	}
	/*Record the index of inital elements to initiate source search*/
	/*Then change the type of initial elements back to 0 to continue propogation*/
	if(mesh->srcelemlen>0 ||  mesh->detelemlen>0){
             int is=0,id=0;
             mesh->srcelem=(int *)calloc(mesh->srcelemlen,sizeof(int));
             mesh->detelem=(int *)calloc(mesh->detelemlen,sizeof(int));
             for(i=0;i<mesh->ne;i++){
		if(mesh->type[i]<0){
                     if(mesh->type[i]==-1){
			mesh->srcelem[is++]=i+1;
			mesh->type[i]=0;
                     }else if(mesh->type[i]==-2) /*keep -2, will be replaced to medianum+1 in loadmedia*/
			mesh->detelem[id++]=i+1;
		}
             }
	}
}

/**
 * @brief Load tet element volume file and initialize the related mesh properties
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadelemvol(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int tmp,len,i,j,*ee;
	char fvelem[MAX_PATH_LENGTH];
	mesh_filenames("velem_%s.dat",fvelem,cfg);
	if((fp=fopen(fvelem,"rt"))==NULL){
		MESH_ERROR("can not open element volume file");
	}
	len=fscanf(fp,"%d %d",&tmp,&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		MESH_ERROR("mesh file has wrong format");
	}
        mesh->evol=(float *)malloc(sizeof(float)*mesh->ne);
	mesh->nvol=(float *)calloc(sizeof(float),mesh->nn);

	for(i=0;i<mesh->ne;i++){
		if(fscanf(fp,"%d %f",&tmp,mesh->evol+i)!=2)
			MESH_ERROR("mesh file has wrong format");
                if(mesh->type[i]==0)
			continue;
		ee=(int *)(mesh->elem+i*mesh->elemlen);
		for(j=0;j<mesh->elemlen;j++)
			mesh->nvol[ee[j]-1]+=mesh->evol[i]*0.25f;
	}
	fclose(fp);
}

/**
 * @brief Load face-neightbor element list and initialize the related mesh properties
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadfaceneighbor(tetmesh *mesh,mcconfig *cfg){
	FILE *fp;
	int len,i,j;
	int *pe;
	char ffacenb[MAX_PATH_LENGTH];
	mesh_filenames("facenb_%s.dat",ffacenb,cfg);

	if((fp=fopen(ffacenb,"rt"))==NULL){
		MESH_ERROR("can not open face-neighbor list file");
	}
	len=fscanf(fp,"%d %d",&(mesh->elemlen),&(mesh->ne));
	if(len!=2 || mesh->ne<=0){
		MESH_ERROR("mesh file has wrong format");
	}
	if(mesh->elemlen<4)
	    mesh->elemlen=4;
	mesh->facenb=(int *)malloc(sizeof(int)*mesh->elemlen*mesh->ne);

	for(i=0;i<mesh->ne;i++){
		pe=mesh->facenb+i*mesh->elemlen;
		for(j=0;j<mesh->elemlen;j++)
		        if(fscanf(fp,"%d",pe+j)!=1)
			    MESH_ERROR("face-neighbor list file has wrong format");
	}
	fclose(fp);
}

/**
 * @brief Load previously saved photon seeds from an .mch file for replay
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration structure
 */

void mesh_loadseedfile(tetmesh *mesh, mcconfig *cfg){
    history his;
    FILE *fp=fopen(cfg->seedfile,"rb");
    if(fp==NULL)
        MESH_ERROR("can not open the specified history file");
    if(fread(&his,sizeof(history),1,fp)!=1)
        MESH_ERROR("error when reading the history file");
    if(his.savedphoton==0 || his.seedbyte==0){
        fclose(fp);
        return;
    }
    if(his.maxmedia!=mesh->prop)
        MESH_ERROR("the history file was generated with a different media setting");
    if(fseek(fp,his.savedphoton*his.colcount*sizeof(float),SEEK_CUR))
        MESH_ERROR("illegal history file");
    cfg->photonseed=malloc(his.savedphoton*his.seedbyte);
    if(cfg->photonseed==NULL)
        MESH_ERROR("can not allocate memory");
    if(fread(cfg->photonseed,his.seedbyte,his.savedphoton,fp)!=his.savedphoton)
        MESH_ERROR("error when reading the seed data");
    cfg->seed=SEED_FROM_FILE;
    cfg->nphoton=his.savedphoton;

    if(cfg->outputtype==otJacobian || cfg->outputtype==otWL || cfg->outputtype==otWP || cfg->replaydet>0){
       int i,j;
       float *ppath=(float*)malloc(his.savedphoton*his.colcount*sizeof(float));
       cfg->replayweight=(float*)malloc(his.savedphoton*sizeof(float));
       cfg->replaytime=(float*)malloc(his.savedphoton*sizeof(float));
       fseek(fp,sizeof(his),SEEK_SET);
       if(fread(ppath,his.colcount*sizeof(float),his.savedphoton,fp)!=his.savedphoton)
           MESH_ERROR("error when reading the partial path data");

       cfg->nphoton=0;
       for(i=0;i<his.savedphoton;i++)
           if(cfg->replaydet==0 || cfg->replaydet==(int)(ppath[i*his.colcount])){
               memcpy((char *)(cfg->photonseed)+cfg->nphoton*his.seedbyte, (char *)(cfg->photonseed)+i*his.seedbyte, his.seedbyte);
		// replay with wide-field detection pattern, the partial path has to contain photon exit information
		if((cfg->detparam1.w*cfg->detparam2.w>0) && (cfg->detpattern!=NULL)){
		    cfg->replayweight[cfg->nphoton]=mesh_getdetweight(i,his.colcount,ppath,cfg);
		}else{
               	    cfg->replayweight[cfg->nphoton]=ppath[(i+1)*his.colcount-1];
		}
               for(j=2;j<his.maxmedia+2;j++)
                 cfg->replayweight[cfg->nphoton]*=expf(-mesh->med[j-1].mua*ppath[i*his.colcount+j]*his.unitinmm);
               cfg->replaytime[cfg->nphoton]=0.f;
               for(j=2;j<his.maxmedia+2;j++)
                 cfg->replaytime[cfg->nphoton]+=mesh->med[j-1].n*ppath[i*his.colcount+j]*R_C0;
           cfg->nphoton++;
           }
	free(ppath);
        cfg->photonseed=realloc(cfg->photonseed, cfg->nphoton*his.seedbyte);
        cfg->replayweight=(float*)realloc(cfg->replayweight, cfg->nphoton*sizeof(float));
	cfg->replaytime=(float*)realloc(cfg->replaytime, cfg->nphoton*sizeof(float));
	cfg->minenergy=0.f;
    }
    fclose(fp);
}

/**
 * @brief Clearing the mesh data structure
 *
 * Destructor of the mesh data structure, delete all dynamically allocated members
 */


void mesh_clear(tetmesh *mesh){
	mesh->nn=0;
	mesh->ne=0;
	mesh->nf=0;
        mesh->srcelemlen=0;
        mesh->detelemlen=0;
	if(mesh->node){
		free(mesh->node);
		mesh->node=NULL;
	}
	if(mesh->elem){
		free(mesh->elem);
		mesh->elem=NULL;
	}
	if(mesh->elem2){
		free(mesh->elem2);
		mesh->elem2=NULL;
	}
	if(mesh->facenb){
		free(mesh->facenb);
		mesh->facenb=NULL;
	}
	if(mesh->dref){
		free(mesh->dref);
		mesh->dref=NULL;
	}
	if(mesh->type){
		free(mesh->type);
		mesh->type=NULL;
	}
	if(mesh->med){
		free(mesh->med);
		mesh->med=NULL;
	}
        if(mesh->atte){
                free(mesh->atte);
                mesh->atte=NULL;
        }
	if(mesh->weight){
		free(mesh->weight);
		mesh->weight=NULL;
	}
        if(mesh->evol){
                free(mesh->evol);
                mesh->evol=NULL;
        }
	if(mesh->nvol){
		free(mesh->nvol);
		mesh->nvol=NULL;
	}
	if(mesh->srcelem){
		free(mesh->srcelem);
		mesh->srcelem=NULL;
	}
	if(mesh->detelem){
		free(mesh->detelem);
		mesh->detelem=NULL;
	}
}


/**
 * @brief Initialize a data structure storing all pre-computed ray-tracing related data
 *
 * the pre-computed ray-tracing data include 
 * d: the vector of the 6 edges
 * m: the cross-product of the end-nodes of the 6-edges n1 x n2
 * n: the normal vector of the 4 facets, pointing outwards
 *
 * @param[out] tracer: the ray-tracer data structure
 * @param[in] pmesh: the mesh object
 * @param[in] methodid: the ray-tracing algorithm to be used
 */

void tracer_init(raytracer *tracer,tetmesh *pmesh,char methodid){
	tracer->d=NULL;
	tracer->m=NULL;
	tracer->n=NULL;
	tracer->mesh=pmesh;
	tracer->method=methodid;
	tracer_build(tracer);
}

/**
 * @brief Preparing for the ray-tracing calculations
 *
 * This function build the ray-tracing precomputed data and test if the initial
 * element (e0) encloses the photon launch position.
 *
 * @param[out] tracer: the ray-tracer data structure
 * @param[in] cfg: the simulation configuration structure
 */

void tracer_prep(raytracer *tracer,mcconfig *cfg){
        int i, ne=tracer->mesh->ne;
	if(tracer->n==NULL && tracer->m==NULL && tracer->d==NULL){
	    if(tracer->mesh!=NULL)
		tracer_build(tracer);
	    else
	    	MESH_ERROR("tracer is not associated with a mesh");
	}else if( (cfg->srctype==stPencil || cfg->srctype==stIsotropic || cfg->srctype==stCone || cfg->srctype==stArcSin)  && cfg->e0>0){
            int eid=cfg->e0-1;
	    float3 vecS={0.f}, *nodes=tracer->mesh->node, vecAB, vecAC, vecN;
	    int ea,eb,ec;
	    float s=0.f, *bary=&(cfg->bary0.x);
	    int *elems=(int *)(tracer->mesh->elem+eid*tracer->mesh->elemlen); // convert int4* to int*
	    if(eid>=tracer->mesh->ne)
	        MESH_ERROR("initial element index exceeds total element count");
	    for(i=0;i<4;i++){
            	ea=elems[out[i][0]]-1;
            	eb=elems[out[i][1]]-1;
	    	ec=elems[out[i][2]]-1;
            	vec_diff(&nodes[ea],&nodes[eb],&vecAB);
            	vec_diff(&nodes[ea],&nodes[ec],&vecAC);
	    	vec_diff(&nodes[ea],&(cfg->srcpos),&vecS);
            	vec_cross(&vecAB,&vecAC,&vecN);
	    	bary[facemap[i]]=-vec_dot(&vecS,&vecN);
	    }
	    if(cfg->debuglevel&dlWeight)
	       fprintf(cfg->flog,"initial bary-centric volumes [%e %e %e %e]\n",
	           bary[0]/6.,bary[1]/6.,bary[2]/6.,bary[3]/6.);
	    for(i=0;i<4;i++){
	        if(bary[i]<0.f)
		    MESH_ERROR("initial element does not enclose the source!");
	        s+=bary[i];
	    }
	    for(i=0;i<4;i++){
	        bary[i]/=s;
	    }
	}
	ne=tracer->mesh->ne*tracer->mesh->elemlen;
	tracer->mesh->nf=0;
	for(i=0;i<ne;i++){
		if(tracer->mesh->facenb[i]==0)
		    tracer->mesh->facenb[i]=-(++tracer->mesh->nf);
	}
	if(tracer->mesh->dref)
	    free(tracer->mesh->dref);
        if(cfg->issaveref)
            tracer->mesh->dref=(double *)calloc(sizeof(double)*tracer->mesh->nf*cfg->srcnum,cfg->maxgate);
}

/**
 * @brief Pre-computed ray-tracing related data
 *
 * the pre-computed ray-tracing data include 
 * d: the vector of the 6 edges: edge[i]=[ni1,ni2]: d[i]=ni2 - ni1
 * m: the cross-product of the end-nodes of the 6-edges: edge[i]=[ni1,ni2]: m[i]=ni1 x ni2
 * n: the normal vector of the 4 facets, pointing outwards
 *
 * @param[out] tracer: the ray-tracer data structure
 */

void tracer_build(raytracer *tracer){
	int ne,i,j;
	const int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

	float3 *nodes;
	int *elems,ebase;
	int e1,e0;
	float Rn2;

	if(tracer->d || tracer->m || tracer->n || tracer->mesh==NULL) return;

        if(tracer->mesh->node==NULL||tracer->mesh->elem==NULL||
	   tracer->mesh->facenb==NULL||tracer->mesh->med==NULL)
                MESH_ERROR("mesh is missing");

	ne=tracer->mesh->ne;
	nodes=tracer->mesh->node;
	elems=(int *)(tracer->mesh->elem); // convert int4* to int*
	if(tracer->method==rtPlucker){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f};

		tracer->d=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		tracer->m=(float3*)calloc(sizeof(float3),ne*6); // 6 edges/elem
		tracer->n=(float3*)calloc(sizeof(float3),ne*4); // 4 face norms
		for(i=0;i<ne;i++){
			ebase=i<<2;
			for(j=0;j<6;j++){
				e1=elems[ebase+pairs[j][1]]-1;
				e0=elems[ebase+pairs[j][0]]-1;
				vec_diff(&nodes[e0],&nodes[e1],tracer->d+i*6+j);
				vec_cross(&nodes[e0],&nodes[e1],tracer->m+i*6+j);
			}
			for(j=0;j<4;j++){
                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);
				vec_cross(&vecAB,&vecAC,tracer->n+ebase+j);
				Rn2=1.f/sqrt(vec_dot(tracer->n+ebase+j,tracer->n+ebase+j));
				vec_mult(tracer->n+ebase+j,Rn2,tracer->n+ebase+j);
			}
		}
	}else if(tracer->method==rtHavel || tracer->method==rtBadouel){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f};

		tracer->d=NULL;
		tracer->m=(float3*)calloc(sizeof(float3),ne*12);
                for(i=0;i<ne;i++){
                        ebase=i<<2;
			for(j=0;j<4;j++){
				float3 *vecN=tracer->m+3*(ebase+j);

                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);
				vec_cross(&vecAB,&vecAC,vecN); /*N is defined as ACxAB in Jiri's code, but not the paper*/
                                vec_cross(&vecAC,vecN,vecN+1);
                                vec_cross(vecN,&vecAB,vecN+2);

				Rn2=1.f/sqrt(vec_dot(vecN,vecN));

				vec_mult(vecN,Rn2,vecN);
				
				Rn2*=Rn2;
				vec_mult(vecN+1,Rn2,vecN+1);
                                vec_mult(vecN+2,Rn2,vecN+2);
#ifdef MMC_USE_SSE
				vecN->w    = vec_dot(vecN,  &nodes[ea]);
				(vecN+1)->w=-vec_dot(vecN+1,&nodes[ea]);
                                (vecN+2)->w=-vec_dot(vecN+2,&nodes[ea]);
#endif
			}
                }
	}else if(tracer->method==rtBLBadouel || tracer->method==rtBLBadouelGrid){
		int ea,eb,ec;
		float3 vecAB={0.f},vecAC={0.f},vN={0.f};

		tracer->d=NULL;
		tracer->n=(float3*)calloc(sizeof(float3),ne*4);
                for(i=0;i<ne;i++){
                        ebase=i<<2;
			float *vecN=&(tracer->n[ebase].x);
			for(j=0;j<4;j++){
                                ea=elems[ebase+out[j][0]]-1;
                                eb=elems[ebase+out[j][1]]-1;
				ec=elems[ebase+out[j][2]]-1;
                                vec_diff(&nodes[ea],&nodes[eb],&vecAB);
                                vec_diff(&nodes[ea],&nodes[ec],&vecAC);

				vec_cross(&vecAB,&vecAC,&vN); /*N is defined as ACxAB in Jiri's code, but not the paper*/

				Rn2=1.f/sqrt(vec_dot(&vN,&vN));
				vec_mult(&vN,Rn2,&vN);
				vecN[j]=vN.x;
				vecN[j+4]=vN.y;
				vecN[j+8]=vN.z;
#ifdef MMC_USE_SSE
				vecN[j+12]    = vec_dot(&vN, &nodes[ea]);
#endif
			}
                }
	}
}

/**
 * @brief Clear the ray-tracing data structure
 *
 * Deconstructor of the ray-tracing data structure
 *
 * @param[out] tracer: the ray-tracer data structure
 */
 
void tracer_clear(raytracer *tracer){
	if(tracer->d) {
		free(tracer->d);
		tracer->d=NULL;
	}
	if(tracer->m) {
		free(tracer->m);
		tracer->m=NULL;
	}
	if(tracer->n) {
		free(tracer->n);
		tracer->n=NULL;
	}	
	tracer->mesh=NULL;
}

/**
 * @brief Performing one scattering event of the photon
 *
 * This function updates the direction of the photon by performing a scattering calculation
 *
 * @param[in] g: anisotropy g
 * @param[out] dir: current ray direction vector
 * @param[out] ran: random number generator states
 * @param[out] ran0: additional random number generator states
 * @param[out] cfg: the simulation configuration
 * @param[out] pmom: buffer to store momentum transfer data if needed
 */

float mc_next_scatter(float g, float3 *dir,RandType *ran, RandType *ran0, mcconfig *cfg, float *pmom){

    float nextslen;
    float sphi=0.f,cphi=0.f,tmp0,theta,stheta,ctheta,tmp1;
    float3 p;

    rand_need_more(ran,ran0);

    //random scattering length (normalized)
#ifdef MMC_USE_SSE_MATH
    nextslen=rand_next_scatlen_ps(ran);
#else
    nextslen=rand_next_scatlen(ran);
#endif

    //random arimuthal angle
#ifdef MMC_USE_SSE_MATH
    rand_next_aangle_sincos(ran,&sphi,&cphi);
#else
    tmp0=TWO_PI*rand_next_aangle(ran); //next arimuth angle
    mmc_sincosf(tmp0,&sphi,&cphi);
#endif

    //Henyey-Greenstein Phase Function, "Handbook of Optical Biomedical Diagnostics",2002,Chap3,p234
    //see Boas2002

    if(g>EPS){  //if g is too small, the distribution of theta is bad
	tmp0=(1.f-g*g)/(1.f-g+2.f*g*rand_next_zangle(ran));
	tmp0*=tmp0;
	tmp0=(1.f+g*g-tmp0)/(2.f*g);

    	// when ran=1, CUDA will give me 1.000002 for tmp0 which produces nan later
    	if(tmp0> 1.f) tmp0=1.f;
        if(tmp0<-1.f) tmp0=-1.f;

	theta=acosf(tmp0);
	stheta=sqrt(1.f-tmp0*tmp0);
	//stheta=sinf(theta);
	ctheta=tmp0;
    }else{
	theta=acosf(2.f*rand_next_zangle(ran)-1.f);
    	mmc_sincosf(theta,&stheta,&ctheta);
    }

    if( dir->z>-1.f+EPS && dir->z<1.f-EPS ) {
	tmp0=1.f-dir->z*dir->z;   //reuse tmp to minimize registers
	tmp1=1.f/sqrtf(tmp0);
	tmp1=stheta*tmp1;

	p.x=tmp1*(dir->x*dir->z*cphi - dir->y*sphi) + dir->x*ctheta;
	p.y=tmp1*(dir->y*dir->z*cphi + dir->x*sphi) + dir->y*ctheta;
	p.z=-tmp1*tmp0*cphi			    + dir->z*ctheta;
    }else{
	p.x=stheta*cphi;
	p.y=stheta*sphi;
	p.z=(dir->z>0.f)?ctheta:-ctheta;
    }
    if (cfg->ismomentum)
        pmom[0]+=(1.f-ctheta);

    dir->x=p.x;
    dir->y=p.y;
    dir->z=p.z;
    return nextslen;
}

/**
 * @brief Save a snapshot of the simulation output
 *
 * save snapshot of the simulation output during a lengthy simulation
 * the snapshots are specified by a set of "check-points", i.e. the index 
 * of the photons that being completed.
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration
 * @param[in] id: the index of the current photon
 */

void mesh_saveweightat(tetmesh *mesh,mcconfig *cfg,int id){
	char sess[MAX_SESSION_LENGTH];
	int i,found=0;
	for(i=0;i<MAX_CHECKPOINT;i++){
           if(cfg->checkpt[i]==0)  return;
	   if(id==cfg->checkpt[i]) {
		found=1;
		break;
	   }
	}
        if(!found) return;
	memcpy(sess,cfg->session,MAX_SESSION_LENGTH);
	sprintf(cfg->session,"%s_%d",sess,id);
	mesh_saveweight(mesh,cfg,0);
	memcpy(cfg->session,sess,MAX_SESSION_LENGTH);
}


/**
 * @brief Save the fluence output to a file
 *
 * @param[in] mesh: the mesh object
 * @param[in] cfg: the simulation configuration
 */

void mesh_saveweight(tetmesh *mesh,mcconfig *cfg,int isref){
	FILE *fp;
	int i,j, datalen=(cfg->method==rtBLBadouelGrid) ? cfg->crop0.z : ( (cfg->basisorder) ? mesh->nn : mesh->ne);
	char fweight[MAX_PATH_LENGTH];
	double *data=mesh->weight;

	if(isref){
	    data=mesh->dref;
	    datalen=mesh->nf;
	}
        if(cfg->rootpath[0])
                sprintf(fweight,"%s%c%s%s.dat",cfg->rootpath,pathsep,cfg->session,(isref? "_dref": ""));
        else
                sprintf(fweight,"%s%s.dat",cfg->session,(isref? "_dref": ""));

        if(cfg->outputformat>=ofBin && cfg->outputformat<=ofTX3){
	        uint3 dim0=cfg->dim;
	        if(cfg->method!=rtBLBadouelGrid){
		    cfg->dim.x=cfg->srcnum;
		    cfg->dim.y=cfg->maxgate;
		    cfg->dim.z=datalen;
		}
		mcx_savedata(mesh->weight,datalen*cfg->maxgate*cfg->srcnum,cfg,isref);
		cfg->dim=dim0;
		return;
	}
	if((fp=fopen(fweight,"wt"))==NULL){
		MESH_ERROR("can not open weight file to write");
	}
	for(i=0;i<cfg->maxgate;i++){
	    for(j=0;j<datalen;j++){
	    	if(1==cfg->srcnum){
	    	    if(fprintf(fp,"%d\t%e\n",j+1,data[i*datalen+j])==0)
			MESH_ERROR("can not write to weight file");
	    	}else{  // multiple sources for pattern illumination type
	    	    int k, shift;
		    for(k=0;k<cfg->srcnum;k++){
		    	shift = (i*datalen+j)*cfg->srcnum+k;
			    if(fprintf(fp,"%d\t%d\t%e\n",j+1,k+1,data[shift])==0)
				MESH_ERROR("can not write to weight file");
		    }
	    	}
	    }
	}
	fclose(fp);
}

/**
 * @brief Save detected photon data into an .mch history file
 *
 * @param[in] ppath: buffer points to the detected photon data (partial-path, det id, etc)
 * @param[in] seeds: buffer points to the detected photon seeds
 * @param[in] count: how many photons are detected
 * @param[in] seedbyte: how many bytes per detected photon seed
 * @param[in] cfg: the simulation configuration
 */

void mesh_savedetphoton(float *ppath, void *seeds, int count, int seedbyte, mcconfig *cfg){
	FILE *fp;
	char fhistory[MAX_PATH_LENGTH];
        if(cfg->rootpath[0])
                sprintf(fhistory,"%s%c%s.mch",cfg->rootpath,pathsep,cfg->session);
        else
                sprintf(fhistory,"%s.mch",cfg->session);

	if((fp=fopen(fhistory,"wb"))==NULL){
		MESH_ERROR("can not open history file to write");
	}
	cfg->his.totalphoton=cfg->nphoton;
	cfg->his.unitinmm=1.f;
        if(cfg->method!=rtBLBadouelGrid)
	    cfg->his.unitinmm=cfg->unitinmm;
        cfg->his.detected=count;
	cfg->his.savedphoton=count;
	cfg->his.srcnum=cfg->srcnum;
	cfg->his.detnum=cfg->detnum;
	if(cfg->issaveseed && seeds!=NULL){
	   cfg->his.seedbyte=seedbyte;
        }
        cfg->his.colcount=(2+(cfg->ismomentum>0))*cfg->his.maxmedia+(cfg->issaveexit>0)*6+2; /*column count=maxmedia+3*/
	
	if(count>0 && cfg->exportdetected==NULL){
            cfg->detectedcount=count;
            cfg->exportdetected=(float*)malloc(cfg->his.colcount*cfg->detectedcount*sizeof(float));
        }
        memcpy(cfg->exportdetected,ppath,count*cfg->his.colcount*sizeof(float));
	
	fwrite(&(cfg->his),sizeof(history),1,fp);
	fwrite(ppath,sizeof(float),count*cfg->his.colcount,fp);
	if(cfg->issaveseed && seeds!=NULL)
           fwrite(seeds,seedbyte,count,fp);
	fclose(fp);
}


/**
 * @brief Save binned detected photon data over an area-detector as time-resolved 2D images
 *
 * When an area detector is used (such as a CCD), storing all detected photons can generate
 * a huge output file. This can be mitigated by accumulate the data first to a rasterized 
 * detector grid, and then save only the integrated data.
 *
 * @param[out] detmap: buffer points to the output detector data array
 * @param[in] ppath: buffer points to the detected photon data (partial-path, det id, etc)
 * @param[in] count: how many photons are detected
 * @param[in] cfg: the simulation configuration
 * @param[in] mesh: the mesh object 
 */

void mesh_getdetimage(float *detmap, float *ppath, int count, mcconfig *cfg, tetmesh *mesh){
	// cfg->issaveexit is 2 for this mode
	int colcount=(2+(cfg->ismomentum>0))*cfg->his.maxmedia+6+2;
	float x0=cfg->detpos[0].x;
	float y0=cfg->detpos[0].y;
	float xrange=cfg->detparam1.x+cfg->detparam2.x; 
	float yrange=cfg->detparam1.y+cfg->detparam2.y;
	int xsize=cfg->detparam1.w;
	int ysize=cfg->detparam2.w;
	int i,j,xindex,yindex,ntg,offset;
	float unitinmm=(cfg->method!=rtBLBadouelGrid)? cfg->his.unitinmm : 1.f;

	float xloc, yloc, weight, path;
	for(i=0; i<count; i++){
		path = 0;
		weight = ppath[(i+1)*colcount-1];
		for(j=1;j<=cfg->his.maxmedia;j++){
			path += ppath[i*colcount+j+cfg->his.maxmedia]*mesh->med[j].n;
			weight *= expf(-ppath[i*colcount+j+cfg->his.maxmedia]*mesh->med[j].mua*unitinmm);
		}
		ntg = (int) path*R_C0/cfg->tstep;
		if(ntg>cfg->maxgate-1)
			ntg = cfg->maxgate-1;
		xloc = ppath[(i+1)*colcount-7];
		yloc = ppath[(i+1)*colcount-6];
		xindex = (xloc-x0)/xrange*xsize;
		if(xindex<0 || xindex>xsize-1) continue;
		yindex = (yloc-y0)/yrange*ysize;
		if(yindex<0 || yindex>ysize-1) continue;
		offset = ntg*xsize*ysize;
		detmap[offset+yindex*xsize+xindex] += weight;
	}
}

/**
 * @brief Save binned detected photon data over an area-detector
 *
 * function for saving binned detected photon data into time-resolved 2D images
 *
 * @param[in] detmap: buffer points to the output detector data array
 * @param[in] cfg: the simulation configuration
 */

void mesh_savedetimage(float *detmap, mcconfig *cfg){
	
	FILE *fp;
	char fhistory[MAX_PATH_LENGTH];
        if(cfg->rootpath[0])
                sprintf(fhistory,"%s%c%s.img",cfg->rootpath,pathsep,cfg->session);
        else
                sprintf(fhistory,"%s.img",cfg->session);
	if((fp=fopen(fhistory,"wb"))==NULL){
		MESH_ERROR("can not open detector image file to write");
	}
	fwrite(detmap,sizeof(float),cfg->detparam1.w*cfg->detparam2.w*cfg->maxgate,fp);
	fclose(fp);
}

/**
 * @brief Recompute the detected photon weight from the partial-pathlengths
 *
 * This function currently does not consider the final transmission coeff before
 * the photon being detected.
 *
 * @param[in] photonid: index of the detected photon
 * @param[in] colcount: how many 4-byte records per detected photon
 * @param[in] ppath: buffer points to the detected photon data (partial-path, det id, etc)
 * @param[in] cfg: the simulation configuration
 */

float mesh_getdetweight(int photonid, int colcount, float* ppath, mcconfig* cfg){
	
	float x0=cfg->detpos[0].x;
	float y0=cfg->detpos[0].y;
	float xrange=cfg->detparam1.x+cfg->detparam2.x; 
	float yrange=cfg->detparam1.y+cfg->detparam2.y;
	int xsize=cfg->detparam1.w;
	int ysize=cfg->detparam2.w;
	float xloc=ppath[(photonid+1)*colcount-7];
	float yloc=ppath[(photonid+1)*colcount-6];
	int xindex = (xloc-x0)/xrange*xsize;
	int yindex = (yloc-y0)/yrange*ysize;
	if(xindex<0 || xindex>xsize-1 || yindex<0 || yindex>ysize-1)
		MESH_ERROR("photon location not within the detection plane");
	return cfg->detpattern[yindex*xsize+xindex];
}

/**
 * @brief Function to normalize the fluence and remove influence from photon number and volume
 *
 * This function outputs the Green's function from the raw simulation data. This needs
 * to divide the total simulated photon energy, normalize the volumes of each node/elem,
 * and consider the length unit and time-gates
 *
 * @param[in] mesh: the mesh object 
 * @param[in] cfg: the simulation configuration
 * @param[in] Eabsorb: total absorbed energy from ray-tracing accummulation
 * @param[in] Etotal: total launched energy, equal to photon number if not pattern-source
 */

/*see Eq (1) in Fang&Boas, Opt. Express, vol 17, No.22, pp. 20178-20190, Oct 2009*/
float mesh_normalize(tetmesh *mesh,mcconfig *cfg, float Eabsorb, float Etotal, int pair){
        int i,j,k;
	double energydeposit=0.f, energyelem,normalizor;
	int *ee;
        int datalen=(cfg->method==rtBLBadouelGrid) ? cfg->crop0.z : ( (cfg->basisorder) ? mesh->nn : mesh->ne);
	
	if(cfg->issaveref && mesh->dref){
	    float normalizor=1.f/Etotal;
            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<mesh->nf;j++)
                  mesh->dref[i*mesh->nf+j]*=normalizor;
        }

	if(cfg->seed==SEED_FROM_FILE && (cfg->outputtype==otJacobian || cfg->outputtype==otWL || cfg->outputtype==otWP)){
            float normalizor=1.f/(DELTA_MUA*cfg->nphoton);
            if(cfg->outputtype==otWL || cfg->outputtype==otWP)
               normalizor=1.f/Etotal; /*Etotal is total detected photon weight in the replay mode*/

            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<datalen;j++)
                  mesh->weight[(i*datalen+j)*cfg->srcnum+pair]*=normalizor;
           return normalizor;
        }
	if(cfg->outputtype==otEnergy){
            normalizor=1.f/Etotal;
            for(i=0;i<cfg->maxgate;i++)
               for(j=0;j<datalen;j++)
                  mesh->weight[(i*datalen+j)*cfg->srcnum+pair]*=normalizor;
	    return normalizor;
        }
	if(cfg->method==rtBLBadouelGrid){
            normalizor=1.0/(Etotal*cfg->unitinmm*cfg->unitinmm*cfg->unitinmm); /*scaling factor*/
	}else{
	  if(cfg->basisorder){
            for(i=0;i<cfg->maxgate;i++)
              for(j=0;j<datalen;j++)
        	if(mesh->nvol[j]>0.f)
                   mesh->weight[(i*datalen+j)*cfg->srcnum+pair]/=mesh->nvol[j];

            for(i=0;i<mesh->ne;i++){
	      ee=(int *)(mesh->elem+i*mesh->elemlen);
	      energyelem=0.f;
	      for(j=0;j<cfg->maxgate;j++)
		for(k=0;k<4;k++)
		   energyelem+=mesh->weight[(j*mesh->nn+ee[k]-1)*cfg->srcnum+pair]; /*1/4 factor is absorbed two lines below*/
	      energydeposit+=energyelem*mesh->evol[i]*mesh->med[mesh->type[i]].mua; /**mesh->med[mesh->type[i]].n;*/
	    }
	    normalizor=Eabsorb/(Etotal*energydeposit*0.25f); /*scaling factor*/
	  }else{
            for(i=0;i<datalen;i++)
	      for(j=0;j<cfg->maxgate;j++)
	         energydeposit+=mesh->weight[(j*datalen+i)*cfg->srcnum+pair];

            for(i=0;i<datalen;i++){
	      energyelem=mesh->evol[i]*mesh->med[mesh->type[i]].mua;
              for(j=0;j<cfg->maxgate;j++)
        	mesh->weight[(j*datalen+i)*cfg->srcnum+pair]/=energyelem;
	    }
            normalizor=Eabsorb/(Etotal*energydeposit); /*scaling factor*/
	  }
	}
	if(cfg->outputtype==otFlux)
               normalizor/=cfg->tstep;

	for(i=0;i<cfg->maxgate;i++)
	    for(j=0;j<datalen;j++)
	        mesh->weight[(i*datalen+j)*cfg->srcnum+pair]*=normalizor;
	return normalizor;
}
