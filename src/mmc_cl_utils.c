/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Pluker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2009) Qianqian Fang and David A. Boas, 
**          <a href="http://www.opticsinfobase.org/abstract.cfm?uri=oe-17-22-20178">
**          "Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated 
**          by Graphics Processing Units,"</a> Optics Express, 17(22) 20178-20190 (2009).
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmc_utils_cl.c

\brief   << Utilities for MMC OpenCL edition >>
*******************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "mmc_cl_utils.h"
#include "mmc_cl_host.h"

const char *VendorList[]={"Unknown","NVIDIA","AMD","Intel","IntelGPU"};

char *print_cl_errstring(cl_int err) {
    switch (err) {
        case CL_SUCCESS:                          return strdup("Success!");
        case CL_DEVICE_NOT_FOUND:                 return strdup("Device not found.");
        case CL_DEVICE_NOT_AVAILABLE:             return strdup("Device not available");
        case CL_COMPILER_NOT_AVAILABLE:           return strdup("Compiler not available");
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return strdup("Memory object allocation failure");
        case CL_OUT_OF_RESOURCES:                 return strdup("Out of resources");
        case CL_OUT_OF_HOST_MEMORY:               return strdup("Out of host memory");
        case CL_PROFILING_INFO_NOT_AVAILABLE:     return strdup("Profiling information not available");
        case CL_MEM_COPY_OVERLAP:                 return strdup("Memory copy overlap");
        case CL_IMAGE_FORMAT_MISMATCH:            return strdup("Image format mismatch");
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return strdup("Image format not supported");
        case CL_BUILD_PROGRAM_FAILURE:            return strdup("Program build failure");
        case CL_MAP_FAILURE:                      return strdup("Map failure");
        case CL_INVALID_VALUE:                    return strdup("Invalid value");
        case CL_INVALID_DEVICE_TYPE:              return strdup("Invalid device type");
        case CL_INVALID_PLATFORM:                 return strdup("Invalid platform");
        case CL_INVALID_DEVICE:                   return strdup("Invalid device");
        case CL_INVALID_CONTEXT:                  return strdup("Invalid context");
        case CL_INVALID_QUEUE_PROPERTIES:         return strdup("Invalid queue properties");
        case CL_INVALID_COMMAND_QUEUE:            return strdup("Invalid command queue");
        case CL_INVALID_HOST_PTR:                 return strdup("Invalid host pointer");
        case CL_INVALID_MEM_OBJECT:               return strdup("Invalid memory object");
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return strdup("Invalid image format descriptor");
        case CL_INVALID_IMAGE_SIZE:               return strdup("Invalid image size");
        case CL_INVALID_SAMPLER:                  return strdup("Invalid sampler");
        case CL_INVALID_BINARY:                   return strdup("Invalid binary");
        case CL_INVALID_BUILD_OPTIONS:            return strdup("Invalid build options");
        case CL_INVALID_PROGRAM:                  return strdup("Invalid program");
        case CL_INVALID_PROGRAM_EXECUTABLE:       return strdup("Invalid program executable");
        case CL_INVALID_KERNEL_NAME:              return strdup("Invalid kernel name");
        case CL_INVALID_KERNEL_DEFINITION:        return strdup("Invalid kernel definition");
        case CL_INVALID_KERNEL:                   return strdup("Invalid kernel");
        case CL_INVALID_ARG_INDEX:                return strdup("Invalid argument index");
        case CL_INVALID_ARG_VALUE:                return strdup("Invalid argument value");
        case CL_INVALID_ARG_SIZE:                 return strdup("Invalid argument size");
        case CL_INVALID_KERNEL_ARGS:              return strdup("Invalid kernel arguments");
        case CL_INVALID_WORK_DIMENSION:           return strdup("Invalid work dimension");
        case CL_INVALID_WORK_GROUP_SIZE:          return strdup("Invalid work group size");
        case CL_INVALID_WORK_ITEM_SIZE:           return strdup("Invalid work item size");
        case CL_INVALID_GLOBAL_OFFSET:            return strdup("Invalid global offset");
        case CL_INVALID_EVENT_WAIT_LIST:          return strdup("Invalid event wait list");
        case CL_INVALID_EVENT:                    return strdup("Invalid event");
        case CL_INVALID_OPERATION:                return strdup("Invalid operation");
        case CL_INVALID_GL_OBJECT:                return strdup("Invalid OpenGL object");
        case CL_INVALID_BUFFER_SIZE:              return strdup("Invalid buffer size");
        case CL_INVALID_MIP_LEVEL:                return strdup("Invalid mip-map level");
        default:                                  return strdup("Unknown");
    }
}


/*
   assert cuda memory allocation result
*/
void ocl_assess(int cuerr,const char *file,const int linenum){
     if(cuerr!=CL_SUCCESS){
         mcx_error(-(int)cuerr,print_cl_errstring(cuerr),file,linenum);
     }
}


/**
  obtain GPU core number per MP, this replaces 
  ConvertSMVer2Cores() in libcudautils to avoid 
  extra dependency.
*/
int mcx_nv_corecount(int v1, int v2){  
     int v=v1*10+v2;
     if(v<20)      return 8;
     else if(v<21) return 32;
     else if(v<30) return 48;
     else if(v<50) return 192;
     else if(v<60) return 128;
     else if(v<61) return 64;
     else          return 128;
}

/**
 * @brief Utility function to query GPU info and set active GPU
 *
 * This function query and list all available GPUs on the system and print
 * their parameters. This is used when -L or -I is used.
 *
 * @param[in,out] cfg: the simulation configuration structure
 * @param[out] info: the GPU information structure
 */

cl_platform_id mcx_list_gpu(mcconfig *cfg,unsigned int *activedev,cl_device_id *activedevlist,GPUInfo **info){

    uint i,j,k,cuid=0,devnum;
    cl_uint numPlatforms;
    cl_platform_id platform = NULL, activeplatform=NULL;
    cl_device_type devtype[]={CL_DEVICE_TYPE_GPU,CL_DEVICE_TYPE_CPU};
    cl_context context;                 // compute context
    const char *devname[]={"GPU","CPU"};
    char pbuf[100];
    cl_context_properties cps[3]={CL_CONTEXT_PLATFORM, 0, 0};
    cl_int status = 0;
    size_t deviceListSize;
    int totaldevice=0;

    OCL_ASSERT((clGetPlatformIDs(0, NULL, &numPlatforms)));
    if(activedev) *activedev=0;

    *info=(GPUInfo *)calloc(MAX_DEVICE,sizeof(GPUInfo));

    if (numPlatforms>0) {
        cl_platform_id* platforms =(cl_platform_id*)malloc(sizeof(cl_platform_id)*numPlatforms);
        OCL_ASSERT((clGetPlatformIDs(numPlatforms, platforms, NULL)));
        for (i = 0; i < numPlatforms; ++i) {

            platform = platforms[i];
	    if(1){
                OCL_ASSERT((clGetPlatformInfo(platforms[i],
                          CL_PLATFORM_NAME,sizeof(pbuf),pbuf,NULL)));
	        if(cfg->isgpuinfo) printf("Platform [%d] Name %s\n",i,pbuf);
                cps[1]=(cl_context_properties)platform;

        	for(j=0; j<2; j++){
		    cl_device_id * devices=NULL;
		    context=clCreateContextFromType(cps,devtype[j],NULL,NULL,&status);
		    if(status!=CL_SUCCESS){
		            clReleaseContext(context);
			    continue;
		    }
		    OCL_ASSERT((clGetContextInfo(context, CL_CONTEXT_DEVICES,0,NULL,&deviceListSize)));
                    devices = (cl_device_id*)malloc(deviceListSize);
                    OCL_ASSERT((clGetContextInfo(context,CL_CONTEXT_DEVICES,deviceListSize,devices,NULL)));
		    devnum=deviceListSize/sizeof(cl_device_id);

                    totaldevice+=devnum;

                    for(k=0;k<devnum;k++){
		          GPUInfo cuinfo={0};
                          cuinfo.platformid=i;
        		  cuinfo.id=cuid+1;
			  cuinfo.devcount=devnum;

        	          OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_NAME,100,(void*)cuinfo.name,NULL)));
			  OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),(void*)&cuinfo.sm,NULL)));
                	  OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),(void*)&cuinfo.globalmem,NULL)));
                          OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),(void*)&cuinfo.sharedmem,NULL)));
                	  OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,sizeof(cl_ulong),(void*)&cuinfo.constmem,NULL)));
                	  OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint),(void*)&cuinfo.clock,NULL)));
        		  cuinfo.maxgate=cfg->maxgate;
        		  cuinfo.autoblock=64;
			  cuinfo.core=cuinfo.sm;
                          if(strstr(pbuf,"NVIDIA") && j==0){
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV,sizeof(cl_uint),(void*)&cuinfo.major,NULL)));
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV,sizeof(cl_uint),(void*)&cuinfo.minor,NULL)));
                               cuinfo.core=cuinfo.sm*mcx_nv_corecount(cuinfo.major,cuinfo.minor);
			       cuinfo.vendor=dvNVIDIA;
                          }else if(strstr(pbuf,"AMD") && j==0){
                               int corepersm=0;
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_GFXIP_MAJOR_AMD,sizeof(cl_uint),(void*)&cuinfo.major,NULL)));
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_GFXIP_MINOR_AMD,sizeof(cl_uint),(void*)&cuinfo.minor,NULL)));
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_BOARD_NAME_AMD,100,(void*)cuinfo.name,NULL)));
                               OCL_ASSERT((clGetDeviceInfo(devices[k],CL_DEVICE_SIMD_PER_COMPUTE_UNIT_AMD,sizeof(cl_uint),(void*)&corepersm,NULL)));
                               corepersm=(corepersm==0) ? 2 : corepersm;
                               cuinfo.core=(cuinfo.sm*corepersm<<4);
                               cuinfo.autoblock=64;
                               cuinfo.vendor=dvAMD;
                          }else if(strstr(pbuf,"Intel") && strstr(cuinfo.name,"Graphics") && j==0){
                               cuinfo.autoblock=64;
			       cuinfo.vendor=dvIntelGPU;
                          }else if(strstr(pbuf,"Intel")){
			       cuinfo.vendor=dvIntel;
			  }
        		  cuinfo.autothread=cuinfo.autoblock * cuinfo.core;

			  if(cfg->isgpuinfo){
                		printf("============ %s device ID %d [%d of %d]: %s  ============\n",devname[j],cuid,k+1,devnum,cuinfo.name);
                		printf(" Device %d of %d:\t\t%s\n",cuid+1,devnum,cuinfo.name);
                		printf(" Compute units   :\t%d core(s)\n",(uint)cuinfo.sm);
                		printf(" Global memory   :\t%ld B\n",(unsigned long)cuinfo.globalmem);
                		printf(" Local memory    :\t%ld B\n",(unsigned long)cuinfo.sharedmem);
                		printf(" Constant memory :\t%ld B\n",(unsigned long)cuinfo.constmem);
                		printf(" Clock speed     :\t%d MHz\n",cuinfo.clock);

                                if(strstr(pbuf,"NVIDIA")){
                                     printf(" Compute Capacity:\t%d.%d\n",cuinfo.major,cuinfo.minor);
                                     printf(" Stream Processor:\t%d\n",cuinfo.core);
                                }else if(strstr(pbuf,"AMD") && j==0){
                                     printf(" GFXIP version:   \t%d.%d\n",cuinfo.major,cuinfo.minor);
                                     printf(" Stream Processor:\t%d\n",cuinfo.core);
                                }
                		printf(" Vendor name    :\t%s\n",VendorList[cuinfo.vendor]);
                		printf(" Auto-thread    :\t%ld\n",cuinfo.autothread);
                		printf(" Auto-block     :\t%ld\n",cuinfo.autoblock);
                      	  }
		          if(activedevlist!=NULL){
			     if(cfg->deviceid[cuid++]=='1'){
			        memcpy((*info)+(*activedev),&cuinfo,sizeof(GPUInfo));
				activedevlist[(*activedev)++]=devices[k];
				if(activeplatform && activeplatform!=platform){
					MMC_FPRINTF(stderr,"Error: one can not mix devices between different platforms\n");fflush(cfg->flog);
					exit(-1);
				}
				activeplatform=platform;
			     }
                          }else{
			        cuid++;
			        memcpy((*info)+((*activedev)++),&cuinfo,sizeof(GPUInfo));
			  }
                    }
                    if(devices) free(devices);
                    clReleaseContext(context);
               }
	    }
        }
	*info=(GPUInfo *)realloc(*info,(*activedev)*sizeof(GPUInfo));
	for (i=0;i<*activedev;i++)
	     (*info)[i].devcount=totaldevice;
        free(platforms);
    }
    if(cfg->isgpuinfo==2 && cfg->parentid==mpStandalone)
        exit(0);

    return activeplatform;
}
