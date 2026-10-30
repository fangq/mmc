/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2026
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing
**          in Plucker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli,
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a>
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**  \li \c (\b Fang2019) Qianqian Fang and Shijie Yan,
**          <a href="http://dx.doi.org/10.1117/1.JBO.24.11.115002">
**          "Graphics processing unit-accelerated mesh-based Monte Carlo photon transport
**           simulations,"</a> J. of Biomedical Optics, 24(11), 115002 (2019)
**  \li \c (\b Yuan2021) Yaoshen Yuan, Shijie Yan, and Qianqian Fang,
**          <a href="https://www.osapublishing.org/boe/fulltext.cfm?uri=boe-12-1-147">
**          "Light transport modeling in highly complex tissues using the implicit
**           mesh-based Monte Carlo algorithm,"</a> Biomed. Optics Express, 12(1) 147-161 (2021)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    mmclab.cpp

@brief   mex function for MMCLAB
*******************************************************************************/

#include <cctype>
#include <cerrno>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <exception>
#include <stdexcept>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "mex.h"

/**
 * mexErrMsgTxt raises an exception from MATLAB's own C++ runtime; the released mex files link a
 * static libstdc++, which aborts MATLAB when that exception unwinds into the try/catch blocks of
 * mmclab (fangq/mcx#278). Inside mmclab, errors are thrown as std::runtime_error instead, and
 * mexFunction passes the message to MATLAB's mexErrMsgTxt once, outside every try block.
 */
static void (*const mcxlab_mexerror)(const char*) = mexErrMsgTxt;

static void mcxlab_throw(const char* msg) {
    throw std::runtime_error(msg);
}

#undef mexErrMsgTxt  /* a macro in the mex.h of newer MATLAB releases */
#define mexErrMsgTxt mcxlab_throw
#include "mmc_const.h"
#include "mmc_mesh.h"
#include "mmc_host.h"
#ifdef USE_OPENCL
    #include "mmc_cl_host.h"
#endif
#ifdef USE_CUDA
    #include "mmc_cu_host.h"
#endif
#include "mmc_tictoc.h"
#include "mmc_raytrace.h"
#if USE_OPTIX
    #include "mmc_optix_host.h"
#endif

//! Macro to read the 1st scalar cfg member
#define GET_1ST_FIELD(x,y)  if(strcmp(name,#y)==0) {double *val=mxGetPr(item);x->y=val[0];printf("mmc.%s=%g;\n",#y,(float)(x->y));}

//! Macro to read one scalar cfg member
#define GET_ONE_FIELD(x,y)  else GET_1ST_FIELD(x,y)

//! Macro to read one 3-element vector member of cfg
#define GET_VEC3_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];\
        printf("mmc.%s=[%g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z));}

//! Macro to read one 3- or 4-element vector member of cfg
#define GET_VEC34_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];if(mxGetNumberOfElements(item)==4) u->v.w=val[3];\
        printf("mmc.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w));}

//! Macro to read one 4-element vector member of cfg
#define GET_VEC4_FIELD(u,v) else if(strcmp(name,#v)==0) {double *val=mxGetPr(item);u->v.x=val[0];u->v.y=val[1];u->v.z=val[2];u->v.w=val[3];\
        printf("mmc.%s=[%g %g %g %g];\n",#v,(float)(u->v.x),(float)(u->v.y),(float)(u->v.z),(float)(u->v.w));}
/**<  Macro to output GPU parameters as field */
#define SET_GPU_INFO(output,id,v)  mxSetField(output,id,#v,mxCreateDoubleScalar(gpuinfo[i].v));

#define MEXERROR(a)  mcx_error(999,a,__FILE__,__LINE__)   //! Macro to add unit name and line number in error printing

typedef mwSize dimtype;                                   //! MATLAB type alias for integer type to use for array sizes and dimensions

#include <stdarg.h>

/**
 * The messages mmclab prints while it reads cfg (the echo of the fields, the progress lines and
 * the warnings) go where cfg.flog sends the simulation log: the MATLAB/Octave console (mexPrintf)
 * for 2 (default), the stdout of the process for 1 and 0, or the log file; 'none'
 * silences everything. The errors and the usage are always printed to the console.
 */
/** the number of dimensions and the 1st dimension of cfg.srcpattern, to set srcnum after srctype is read */
static int mmclab_patndim = 0, mmclab_patdim0 = 1;

/**
 * cfg.flog='none' (any case) silences the log; '/dev/null' and 'NUL', the null devices of
 * Linux/macOS and of Windows, are accepted on every OS too: cfg.flog='NUL' on Linux or macOS
 * wrote, or failed to write, a file of that name (fangq/mcx#280)
 */
static const char* mcx_logpath(const char* path) {
#ifdef _WIN32
    const char* nulldev = "NUL";
#else
    const char* nulldev = "/dev/null";
#endif
    const char* names[] = {"none", "nul"};
    int i, j;

    if (strcmp(path, "/dev/null") == 0) {
        return nulldev;
    }

    for (i = 0; i < 2; i++) {
        for (j = 0; names[i][j] && tolower((unsigned char)path[j]) == names[i][j]; j++);

        if (names[i][j] == '\0' && path[j] == '\0') {
            return nulldev;
        }
    }

    return path;
}

static int mcxlab_lognull = 0;     /**< 1: cfg.flog is the null device, which also silences the progress bar */
static FILE* mcxlab_logfp = NULL;   /**< where printf() of this unit goes; NULL or stderr: the console */

static int mcxlab_printf(const char* format, ...) {
    va_list args;
    int len;

    va_start(args, format);

    if (mcxlab_logfp == NULL || mcxlab_logfp == stderr) {
        va_list args2;
        va_copy(args2, args);
        len = vsnprintf(NULL, 0, format, args2);
        va_end(args2);

        if (len >= 0) {
            char* buf = (char*)malloc(len + 1);

            if (buf) {
                vsnprintf(buf, len + 1, format, args);
                mexPrintf("%s", buf);
                free(buf);
            }
        }
    } else {
        len = vfprintf(mcxlab_logfp, format, args);
    }

    va_end(args);
    return len;
}

#ifdef printf
    #undef printf
#endif
#define printf mcxlab_printf

void mmc_set_field(const mxArray* root, const mxArray* item, int idx, int jstruct, mcconfig* cfg, tetmesh* mesh);
static void mcxlab_setflog(const mxArray* item, mcconfig* cfg);
void mmclab_usage();
static int mmclab_task_guard(void (*task)(void*), void* arg);
mxArray* mmclab_assert(mxArray* arr);

extern const char debugflag[];

float* detps = NULL;       //! buffer to receive data from cfg.detphotons field
int    dimdetps[2] = {0, 0}; //! dimensions of the cfg.detphotons array
int    seedbyte = 0;

/** @brief Mex function for the MMC host function for MATLAB/Octave
 *  This is the master function to interface all MMC features inside MATLAB.
 *  In MMCLAB, all inputs are read from the cfg structure, which contains all
 *  simuation parameters and data.
 */

static void mcxlab_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    mcconfig cfg;
    tetmesh mesh;
    raytracer tracer = {NULL, 0, NULL, NULL, NULL};
    unsigned int threadid = 0, t0, dt;
    GPUInfo* gpuinfo = NULL;

    mxArray*    tmp;
    int        ifield, jstruct;
    int        ncfg, nfields;
    char       cfgerr[MAX_PATH_LENGTH] = "";  /** the error of a simulation, raised after its clean-up */
    dimtype     fielddim[6];
    int        errorflag = 0;
    char       errmsg[MAX_FULL_PATH + 64] = "MMCLAB terminated due to an exception";  /* copied by the exception mexErrMsgTxt throws */
    cl_uint    workdev;

    const char*       outputtag[] = {"data"};
    const char*       datastruct[] = {"data", "dref"};
    const char*       gpuinfotag[] = {"name", "id", "devcount", "major", "minor", "globalmem",
                                      "constmem", "sharedmem", "regcount", "clock", "sm", "core",
                                      "autoblock", "autothread", "maxgate"
                                     };

    /**
     * If no input is given for this function, it prints help information and return.
     */
    if (nrhs == 0) {
        mmclab_usage();
        return;
    }

    /**
     * If a single string is passed, and if this string is 'gpuinfo', this function
     * returns the list of GPUs on this host and return.
     */
    if (nrhs == 1 && mxIsChar(prhs[0])) {
        char shortcmd[MAX_SESSION_LENGTH];
        mxGetString(prhs[0], shortcmd, MAX_SESSION_LENGTH);
        shortcmd[MAX_SESSION_LENGTH - 1] = '\0';

        if (strcmp(shortcmd, "gpuinfo") == 0) {
            mcx_initcfg(&cfg);
            cfg.isgpuinfo = 3;

            try {
                mcx_list_cl_gpu(&cfg, &workdev, NULL, &gpuinfo);
            } catch (...) {
                mexErrMsgTxt("OpenCL is not supported or not fully installed on your system");
            }

            if (!workdev) {
                mexErrMsgTxt("no active GPU device found");
            }

            if (workdev > MAX_DEVICE) {
                workdev = MAX_DEVICE;
            }

            plhs[0] = mxCreateStructMatrix(gpuinfo[0].devcount, 1, 15, gpuinfotag);

            for (cl_uint i = 0; i < workdev; i++) {
                mxSetField(plhs[0], i, "name", mxCreateString(gpuinfo[i].name));
                SET_GPU_INFO(plhs[0], i, id);
                SET_GPU_INFO(plhs[0], i, devcount);
                SET_GPU_INFO(plhs[0], i, major);
                SET_GPU_INFO(plhs[0], i, minor);
                SET_GPU_INFO(plhs[0], i, globalmem);
                SET_GPU_INFO(plhs[0], i, constmem);
                SET_GPU_INFO(plhs[0], i, sharedmem);
                SET_GPU_INFO(plhs[0], i, regcount);
                SET_GPU_INFO(plhs[0], i, clock);
                SET_GPU_INFO(plhs[0], i, sm);
                SET_GPU_INFO(plhs[0], i, core);
                SET_GPU_INFO(plhs[0], i, autoblock);
                SET_GPU_INFO(plhs[0], i, autothread);
                SET_GPU_INFO(plhs[0], i, maxgate);
            }

            mcx_cleargpuinfo(&gpuinfo);
            mcx_clearcfg(&cfg);
        }

        return;
    }

    /**
     * If a structure is passed to this function, a simulation will be launched.
     */
    if (!mxIsStruct(prhs[0])) {
        MEXERROR("Input must be a structure.");
    }

    /**
     * Find out information about input and output.
     */
    nfields = mxGetNumberOfFields(prhs[0]);
    ncfg = mxGetNumberOfElements(prhs[0]);

    /**
     * The function can return 1-3 outputs (i.e. the LHS)
     */
    if (nlhs >= 1) {
        plhs[0] = mxCreateStructMatrix(ncfg, 1, 2, datastruct);
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateStructMatrix(ncfg, 1, 1, outputtag);
    }

    if (nlhs >= 3) {
        plhs[2] = mxCreateStructMatrix(ncfg, 1, 1, outputtag);
    }

    if (nlhs >= 4) {
        plhs[3] = mxCreateStructMatrix(ncfg, 1, 1, outputtag);
    }

    /**
     * Loop over each element of the struct if it is an array, each element is a simulation
     */
    for (jstruct = 0; jstruct < ncfg; jstruct++) {  /* how many configs */

        /** Enclose all simulation calls inside a try/catch construct for exception handling */
        try {

            /** Initialize cfg with default values first */
            t0 = StartTimer();
            mcx_initcfg(&cfg);
            cfg.flog = stderr;  /* the MATLAB/Octave console, unless cfg.flog says otherwise */
            mcxlab_logfp = stderr;

            /** cfg.flog first: the messages below go where it sends the log */
            mcxlab_lognull = 0;

            if (mxGetField(prhs[0], jstruct, "flog") && !mxIsEmpty(mxGetField(prhs[0], jstruct, "flog"))) {
                mcxlab_setflog(mxGetField(prhs[0], jstruct, "flog"), &cfg);
            }

            if (jstruct == 0) {
                printf("Launching MMCLAB - Mesh-based Monte Carlo for MATLAB & GNU Octave ...\n");
            }

            printf("Running simulations for configuration #%d ...\n", jstruct + 1);
            MMCDEBUG(&cfg, dlTime, (cfg.flog, "initializing ... "));
            mesh_init(&mesh);

            mmclab_patndim = 0;
            mmclab_patdim0 = 1;

            /** Read each struct element from input and set value to the cfg configuration */
            for (ifield = 0; ifield < nfields; ifield++) { /* how many input struct fields */
                tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield);

                /** an empty field, as a struct array gives the elements that do not set it, keeps
                    the default; it was read as 0 (e.g. isreflect), or past its end for vectors */
                if (tmp == NULL || mxIsEmpty(tmp)) {
                    continue;
                }

                mmc_set_field(prhs[0], tmp, ifield, jstruct, &cfg, &mesh);
            }

            /** the progress bar is written to the console, not to the log: drop it when the log is the
                null device, so that cfg.flog='none' silences everything */
            if (mcxlab_lognull) {
                cfg.debuglevel &= ~dlProgress;
            }

            /** a stack of patterns for photon sharing is [srcnum Nx Ny] for 'pattern' and
                [srcnum Nx Ny Nz] for 'pattern3d', whose 3-D srcpattern is one volume; srcnum was
                taken from any 3-D srcpattern, which broke pattern3d */
            if (mmclab_patndim == (cfg.srctype == stPattern3D ? 4 : 3)) {
                cfg.srcnum = mmclab_patdim0;
            }

            mexEvalString("pause(.001);");

            /** Overwite the output flags using the number of output present */
            if (nlhs < 1) {
                cfg.issave2pt = 0;    /** issave2pt defaults to 1, users can disable it; auto disable only if there is no output */
            }

            cfg.issavedet = (nlhs >= 2); /** save detected photon data to the 2nd output if present */
            cfg.issaveseed = (nlhs >= 3); /** save detected photon seeds to the 3rd output if present */

            if (nlhs >= 4) {
                cfg.exportdebugdata = (float*)malloc(cfg.maxjumpdebug * sizeof(float) * MCX_DEBUG_REC_LEN);
                cfg.debuglevel |= dlTraj;
            }

#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
            cfg.issaveseed = 0;
#endif
            mesh_srcdetelem(&mesh, &cfg);

            /** Build srcdata from detectors for adjoint / srcid=-2 mode.
             *  The field buffer uses nsrcslots = Ns + detnum slots:
             *    Slots 0..Ns-1                : forward source photons (phi_src)
             *    Slots Ns..Ns+detnum-1        : detector-as-source photons (phi_det)
             *  The adjoint Jacobian kernel then computes J[vox,s,d] = phi_src[s] * phi_det[d].
             *  srcid==-2 builds the same slot layout for forward fluence only (no Jacobian).
             *
             *  Two routes feed Ns:
             *    (a) cfg.srcpos was an Mx{3,4} matrix -> the multi-source srcpos parser
             *        already populated cfg.srcdata with M forward slots and set
             *        cfg.extrasrclen = M. Keep those slots and only append Nd detectors.
             *    (b) cfg.srcpos was a single row -> srcdata is unpopulated. Replicate
             *        the single source cfg.srcnum times as before (single-source-with-
             *        photon-sharing convention).
             */
            if ((MCX_IS_ADJOINT_TYPE(cfg.outputtype) || cfg.srcid == -2)
                    && cfg.detnum > 0 && cfg.detdir != NULL) {
                int Nd = cfg.detnum;
                int Ns;
                int already_populated = (cfg.srcdata != NULL && cfg.extrasrclen > 0);

                if (already_populated) {
                    Ns = cfg.extrasrclen;
                    cfg.srcdata = (ExtraSrc*)realloc(cfg.srcdata, (Ns + Nd) * sizeof(ExtraSrc));
                    memset(cfg.srcdata + Ns, 0, Nd * sizeof(ExtraSrc));

                    /* Uniform per-slot launch weight: every slot (forward source and
                     * detector-as-adjoint source) launches photons with srcpos.w = 1,
                     * matching MCX's convention (mcx_utils.c:1822, :1985-1986). With
                     * uniform random slot picking, each slot receives N_photon/(Ns+Nd)
                     * photons of unit weight, so the total simulated energy equals
                     * N_photon and per-slot Eabsorb/Etotal ratios are comparable. */
                    for (int is = 0; is < Ns; is++) {
                        if (cfg.srcdata[is].srcpos.w == 0.f) {
                            cfg.srcdata[is].srcpos.w = 1.f;
                        }
                    }
                } else {
                    if (cfg.srcdata) {
                        free(cfg.srcdata);
                    }

                    Ns = cfg.srcnum;   /* photon-sharing replication count */
                    cfg.srcdata = (ExtraSrc*)calloc(Ns + Nd, sizeof(ExtraSrc));

                    /* Slots 0..Ns-1: forward sources (replicate the single main source).
                     * Uniform unit weight per slot (MCX parity); preserve cfg.srcdir.w
                     * (focal length / lens parameter) so users can control beam profiles
                     * (e.g., planar wavefront via srcdir(4)=-inf) for the source slots;
                     * the detector slots below honor cfg.detdir.w. */
                    for (int is = 0; is < Ns; is++) {
                        cfg.srcdata[is].srcpos    = {cfg.srcpos.x, cfg.srcpos.y,
                                                     cfg.srcpos.z, 1.f
                                                    };
                        cfg.srcdata[is].srcdir    = {cfg.srcdir.x, cfg.srcdir.y,
                                                     cfg.srcdir.z, cfg.srcdir.w
                                                    };
                        cfg.srcdata[is].srcparam1 = cfg.srcparam1;
                        cfg.srcdata[is].srcparam2 = cfg.srcparam2;
                    }
                }

                cfg.extrasrclen = Ns + Nd;

                /* Slots Ns..Ns+Nd-1: detector-as-reversed-source, unit weight per slot. */
                for (int id = 0; id < Nd; id++) {
                    cfg.srcdata[Ns + id].srcpos   = {cfg.detpos[id].x, cfg.detpos[id].y,
                                                     cfg.detpos[id].z, 1.f
                                                    };
                    cfg.srcdata[Ns + id].srcdir   = {cfg.detdir[id].x, cfg.detdir[id].y,
                                                     cfg.detdir[id].z, cfg.detdir[id].w
                                                    };
                    cfg.srcdata[Ns + id].srcparam1 = {cfg.detpos[id].w, 0.f, 0.f, 0.f};
                    cfg.srcdata[Ns + id].srcparam2 = {0.f, 0.f, 0.f, 0.f};
                }

                if (MCX_IS_ADJOINT_TYPE(cfg.outputtype)) {
                    cfg.srcid = -1;
                }
            }

            /** Validate all input fields, and warn incompatible inputs */
            mmc_validate_config(&cfg, detps, dimdetps, seedbyte);
            mesh_validate(&mesh, &cfg);
            /* mmc_prep -> mesh_init_srcdata_eid fills per-slot e0 into
             * cfg.srcdata[].srcparam2.w so the kernel's launchnewphoton can
             * use the correct initial tet for each slot in multi-source mode. */

            if (cfg.isgpuinfo == 0) {
                mmc_prep(&cfg, &mesh, &tracer);
            }

            dt = GetTimeMillis();
            MMCDEBUG(&cfg, dlTime, (cfg.flog, "\tdone\t%d\nsimulating ... \n", dt - t0));


            /** \subsection ssimu Parallel photon transport simulation */

            mmc_task_guard = mmclab_task_guard;  /* catches the errors of the photons of mmc_run_mp */

            try {
                if (cfg.compute == cbSSE || cfg.gpuid > MAX_DEVICE) {
                    mmc_run_mp(&cfg, &mesh, &tracer);
                }

#ifdef USE_CUDA
                else if (cfg.compute == cbCUDA) {
                    mmc_run_cu(&cfg, &mesh, &tracer);
                }

#endif
#ifdef USE_OPTIX
                else if (cfg.compute == cbOptiX) {
                    mmc_run_optix(&cfg, &mesh, &tracer);
                }

#endif
#ifdef USE_OPENCL
                else {
                    mmc_run_cl(&cfg, &mesh, &tracer);
                }

#endif
            } catch (const char* err) {
                mexPrintf("Error from thread (%d): %s\n", threadid, err);
                snprintf(errmsg, sizeof(errmsg), "MMCLAB terminated due to an exception: %s", err);
                errorflag++;
            } catch (const std::exception& err) {
                mexPrintf("C++ Error from thread (%d): %s\n", threadid, err.what());
                snprintf(errmsg, sizeof(errmsg), "MMCLAB terminated due to an exception: %s", err.what());
                errorflag++;
            } catch (...) {
                mexPrintf("Unknown Exception from thread (%d)", threadid);
                errorflag++;
            }

            /** \subsection sreport Post simulation */

            dt = GetTimeMillis() - dt;

            /** Clear up simulation data structures by calling the destructors */

            tracer_clear(&tracer);
            MMCDEBUG(&cfg, dlTime, (cfg.flog, "\tdone\t%d\n", GetTimeMillis() - t0));

            /** if 5th output presents, output the photon trajectory data */
            if (nlhs >= 4) {
                int outputidx = 3;
                fielddim[0] = MCX_DEBUG_REC_LEN;
                fielddim[1] = cfg.debugdatalen; // his.savedphoton is for one repetition, should correct
                fielddim[2] = 0;
                fielddim[3] = 0;
                mxSetFieldByNumber(plhs[outputidx], jstruct, 0, mmclab_assert(mxCreateNumericArray(2, fielddim, mxSINGLE_CLASS, mxREAL)));

                if ((cfg.debuglevel & dlTraj) && cfg.exportdebugdata) {
                    memcpy((float*)mxGetPr(mxGetFieldByNumber(plhs[outputidx], jstruct, 0)), cfg.exportdebugdata, fielddim[0]*fielddim[1]*sizeof(float));
                }
            }

            if (cfg.exportdebugdata) {
                free(cfg.exportdebugdata);
                cfg.exportdebugdata = NULL;
            }

            if (nlhs >= 3) {
                fielddim[0] = (sizeof(RandType) * RAND_BUF_LEN);
                fielddim[1] = cfg.detectedcount;
                fielddim[2] = 0;
                fielddim[3] = 0;
                mxSetFieldByNumber(plhs[2], jstruct, 0, mmclab_assert(mxCreateNumericArray(2, fielddim, mxUINT8_CLASS, mxREAL)));
                memcpy((unsigned char*)mxGetPr(mxGetFieldByNumber(plhs[2], jstruct, 0)), cfg.exportseed, fielddim[0]*fielddim[1]);
            }

            if (cfg.exportseed) {
                free(cfg.exportseed);
                cfg.exportseed = NULL;
            }

            /** if the 2nd output presents, output the detected photon partialpath data */
            if (nlhs >= 2) {
                if (cfg.issaveexit != 2) {
                    int hostdetreclen = (2 + ((cfg.ismomentum) > 0)) * mesh.prop + (cfg.issaveexit > 0) * 6 + 2;
                    fielddim[0] = hostdetreclen;
                    fielddim[1] = cfg.detectedcount;
                    fielddim[2] = 0;
                    fielddim[3] = 0;

                    if (cfg.detectedcount > 0) {
                        mxSetFieldByNumber(plhs[1], jstruct, 0, mmclab_assert(mxCreateNumericArray(2, fielddim, mxSINGLE_CLASS, mxREAL)));
                        memcpy((float*)mxGetPr(mxGetFieldByNumber(plhs[1], jstruct, 0)), cfg.exportdetected,
                               fielddim[0]*fielddim[1]*sizeof(float));
                    }
                } else {
                    fielddim[0] = cfg.detparam1.w;
                    fielddim[1] = cfg.detparam2.w;
                    fielddim[2] = cfg.maxgate;
                    fielddim[3] = 0;
                    mxSetFieldByNumber(plhs[1], jstruct, 0, mmclab_assert(mxCreateNumericArray(3, fielddim, mxSINGLE_CLASS, mxREAL)));
                    float* detmap = (float*)mxGetPr(mxGetFieldByNumber(plhs[1], jstruct, 0));
                    memset(detmap, cfg.detparam1.w * cfg.detparam2.w * cfg.maxgate, sizeof(float));
                    mesh_getdetimage(detmap, cfg.exportdetected, cfg.detectedcount, &cfg, &mesh);
                }
            }

            if (cfg.exportdetected) {
                free(cfg.exportdetected);
                cfg.exportdetected = NULL;
            }

            if (nlhs >= 1) {
                int isrfforward = (cfg.omega > 0.f && cfg.seed != SEED_FROM_FILE);
                int isadjoint = MCX_IS_ADJOINT_TYPE(cfg.outputtype);

                /* When extrasrclen exceeds srcnum, the kernel runs nsrcslots = extrasrclen
                 * slots in the weight buffer regardless of whether outputtype is adjoint
                 * (e.g. cfg.srcid = -2 forward-only multi-source mode). */
                int nsrcslots = mesh_outputslots(&cfg);  /* one per source only for srcid<0 (also set by the adjoint outputs) */

                int datalen = (cfg.method == rtBLBadouelGrid) ? cfg.crop0.z : ( (cfg.basisorder) ? mesh.nn : mesh.ne);
                fielddim[0] = nsrcslots;
                fielddim[1] = datalen;
                fielddim[2] = cfg.maxgate;
                fielddim[3] = 0;
                fielddim[4] = 0;

                /* Always output forward fluence in flux.data.
                 * In adjoint mode, flux.data has nsrcslots = srcnum+detnum slices (sources then
                 * detectors-as-sources), and the Jacobian is returned separately in flux.jmua/jd/etc. */
                /* multi-slot layouts use either:
                 *   detector-adjoint slots (nsrcslots > srcnum; adjoint or srcid==-2 forward):
                 *       voxel-fastest, gate-middle, slot-slowest  ->  [..., maxgate, nsrcslots]
                 *   pattern source (nsrcslots == srcnum > 1):
                 *       pidx-fastest                             ->  [srcnum, ..., maxgate] */
                int ismultislot = (nsrcslots > cfg.srcnum);

                if (cfg.issave2pt) {   /** cfg.issave2pt=0 leaves the 1st output's data empty, as mcxlab does */
                    if (cfg.method == rtBLBadouelGrid) {
                        if (ismultislot) {
                            fielddim[0] = cfg.dim.x;
                            fielddim[1] = cfg.dim.y;
                            fielddim[2] = cfg.dim.z;
                            fielddim[3] = cfg.maxgate;
                            fielddim[4] = nsrcslots;
                        } else {
                            fielddim[0] = cfg.srcnum;
                            fielddim[1] = cfg.dim.x;
                            fielddim[2] = cfg.dim.y;
                            fielddim[3] = cfg.dim.z;
                            fielddim[4] = cfg.maxgate;
                        }

                        mxComplexity cplx = (isrfforward && cfg.exportadjoint) ? mxCOMPLEX : mxREAL;

                        if (ismultislot || cfg.srcnum > 1) {
                            mxSetFieldByNumber(plhs[0], jstruct, 0, mmclab_assert(mxCreateNumericArray(5, fielddim, mxDOUBLE_CLASS, cplx)));
                        } else {
                            mxSetFieldByNumber(plhs[0], jstruct, 0, mmclab_assert(mxCreateNumericArray(4, &fielddim[1], mxDOUBLE_CLASS, cplx)));
                        }
                    } else {
                        /* mesh-mode layout. Kernel write order:
                         *   detector-adjoint slots (nsrcslots > srcnum, srcnum==1):
                         *       field[node + gate*nn + slot*nn*maxgate]  -> dims [datalen, maxgate, nsrcslots]
                         *   pattern source (nsrcslots == srcnum > 1):
                         *       field[(gate*ne + eid)*srcnum + pidx]     -> dims [nsrcslots, datalen, maxgate]
                         *   single source:
                         *       field[node + gate*nn]                    -> dims [datalen, maxgate] */
                        mxComplexity mcplx = (isrfforward && cfg.exportadjoint) ? mxCOMPLEX : mxREAL;

                        if (nsrcslots > 1) {
                            dimtype meshfielddim[3];

                            if (ismultislot) {
                                meshfielddim[0] = datalen;
                                meshfielddim[1] = cfg.maxgate;
                                meshfielddim[2] = nsrcslots;
                            } else {
                                meshfielddim[0] = nsrcslots;
                                meshfielddim[1] = datalen;
                                meshfielddim[2] = cfg.maxgate;
                            }

                            mxSetFieldByNumber(plhs[0], jstruct, 0, mmclab_assert(mxCreateNumericArray(3, meshfielddim, mxDOUBLE_CLASS, mcplx)));
                        } else {
                            mxSetFieldByNumber(plhs[0], jstruct, 0, mmclab_assert(mxCreateNumericArray(2, &fielddim[1], mxDOUBLE_CLASS, mcplx)));
                        }
                    }

                    mxArray* outfield = mxGetFieldByNumber(plhs[0], jstruct, 0);
                    size_t nelem = (size_t)nsrcslots * datalen * cfg.maxgate;

                    if (isrfforward && cfg.exportadjoint) {
#if MX_HAS_INTERLEAVED_COMPLEX
                        mxComplexDouble* cptr = mxGetComplexDoubles(outfield);

                        if (cptr) {
                            for (size_t ii = 0; ii < nelem; ii++) {
                                cptr[ii].real = mesh.weight[ii];
                                cptr[ii].imag = (double)cfg.exportadjoint[ii];
                            }
                        }

#else
                        double* output = (double*)mxGetPr(outfield);
                        double* imoutput = (double*)mxGetPi(outfield);

                        if (output == NULL) {
                            mexErrMsgTxt("MMCLAB failed to allocate memory!");
                        }

                        memcpy(output, mesh.weight, nelem * sizeof(double));

                        if (imoutput) {
                            for (size_t ii = 0; ii < nelem; ii++) {
                                imoutput[ii] = (double)cfg.exportadjoint[ii];
                            }
                        }

#endif
                    } else {
                        double* output = (double*)mxGetPr(outfield);

                        if (output == NULL) {
                            mexErrMsgTxt("MMCLAB failed to allocate memory!");
                        }

                        memcpy(output, mesh.weight, nelem * sizeof(double));
                    }
                }

                if (cfg.issaveref) {      /** save diffuse reflectance */
                    fielddim[1] = mesh.nf;
                    fielddim[2] = cfg.maxgate;
                    mxSetFieldByNumber(plhs[0], jstruct, 1, mmclab_assert(mxCreateNumericArray(2, &fielddim[1], mxDOUBLE_CLASS, mxREAL)));
                    memcpy((double*)mxGetPr(mxGetFieldByNumber(plhs[0], jstruct, 1)), mesh.dref, fielddim[1]*fielddim[2]*sizeof(double));
                }

                /* Output adjoint Jacobian in separate struct fields (flux.jmua, flux.jd, etc.).
                 * Grid mode: shape [Nx, Ny, Nz, maxgate, Nch].
                 * Mesh mode (basisorder=1): shape [nn, Nch] (CW, no maxgate dim).
                 * Nch: Ns*Nd, or the rows of cfg.sdmap */
                if (isadjoint && cfg.exportjacob) {
                    int nsrcpairs = (int)mcx_jacchannels(&cfg);
                    int isdual = MCX_IS_DUAL_ADJOINT_TYPE(cfg.outputtype);
                    mxComplexity jcplx = isrfforward ? mxCOMPLEX : mxREAL;
                    int jndim;
                    size_t adjlen;

                    if (cfg.method == rtBLBadouelGrid) {
                        adjlen = (size_t)cfg.dim.x * cfg.dim.y * cfg.dim.z * cfg.maxgate * nsrcpairs;
                        fielddim[0] = cfg.dim.x;
                        fielddim[1] = cfg.dim.y;
                        fielddim[2] = cfg.dim.z;
                        fielddim[3] = cfg.maxgate;
                        fielddim[4] = nsrcpairs;
                        jndim = 5;
                    } else {
                        /* mesh mode: nodal output, CW only */
                        adjlen = (size_t)mesh.nn * nsrcpairs;
                        fielddim[0] = mesh.nn;
                        fielddim[1] = nsrcpairs;
                        fielddim[2] = 0;
                        fielddim[3] = 0;
                        fielddim[4] = 0;
                        jndim = 2;
                    }

                    /* Determine Jacobian field name(s) based on output type */
                    const char* jname1 = "jmua";
                    const char* jname2 = "jd";

                    switch (cfg.outputtype) {
                        case otAdjoint:
                            jname1 = "jmua";
                            break;

                        case otAdjointDcoeff:
                            jname1 = "jd";
                            break;

                        case otAdjointMus:
                            jname1 = "jmus";
                            break;

                        case otAdjointMusp:
                            jname1 = "jmusp";
                            break;

                        case otAdjointMuaD:
                            jname1 = "jmua";
                            jname2 = "jd";
                            break;

                        case otAdjointMuaMusp:
                            jname1 = "jmua";
                            jname2 = "jmusp";
                            break;

                        default:
                            break;
                    }

                    /* Helper: create and populate a Jacobian field (float32, real or complex) */
                    auto add_jac_field = [&](const char* fname, float * re_data, float * im_data) {
                        if (mxGetFieldNumber(plhs[0], fname) < 0) {
                            mxAddField(plhs[0], fname);
                        }

                        mxArray* jfield = mmclab_assert(mxCreateNumericArray(jndim, fielddim, mxSINGLE_CLASS, jcplx));

                        if (!isrfforward) {
                            memcpy((float*)mxGetData(jfield), re_data, adjlen * sizeof(float));
                        } else {
#if MX_HAS_INTERLEAVED_COMPLEX
                            mxComplexSingle* cptr = mxGetComplexSingles(jfield);

                            if (cptr) {
                                for (size_t ii = 0; ii < adjlen; ii++) {
                                    cptr[ii].real = re_data[ii];
                                    cptr[ii].imag = im_data[ii];
                                }
                            }

#else
                            memcpy((float*)mxGetData(jfield), re_data, adjlen * sizeof(float));
                            float* jpi = (float*)mxGetImagData(jfield);

                            if (jpi) {
                                memcpy(jpi, im_data, adjlen * sizeof(float));
                            }

#endif
                        }

                        mxSetField(plhs[0], jstruct, fname, jfield);
                    };

                    /* Pointer arithmetic into exportjacob:
                     * CW non-dual:  [re1(adjlen)]
                     * CW dual:      [re1(adjlen), re2(adjlen)]
                     * RF non-dual:  [re1(adjlen), im1(adjlen)]
                     * RF dual:      [re1(adjlen), re2(adjlen), im1(adjlen), im2(adjlen)] */
                    float* re1 = cfg.exportjacob;
                    float* re2 = isdual               ? cfg.exportjacob + adjlen         : NULL;
                    float* im1 = isrfforward          ? cfg.exportjacob + (isdual ? 2 : 1) * adjlen : NULL;
                    float* im2 = (isrfforward && isdual) ? cfg.exportjacob + 3 * adjlen : NULL;

                    add_jac_field(jname1, re1, im1);

                    if (isdual) {
                        add_jac_field(jname2, re2, im2);
                    }
                }
            }

            if (errorflag) {
                mexErrMsgTxt(errmsg);  /* with the message: it was only printed, where flog could hide it */
            }
        } catch (const char* err) {
            snprintf(cfgerr, sizeof(cfgerr), "%s", err);
        } catch (const std::exception& err) {
            snprintf(cfgerr, sizeof(cfgerr), "%s", err.what());
        } catch (...) {
            snprintf(cfgerr, sizeof(cfgerr), "unknown exception");
        }

        /** \subsection sclean End the simulation */
        mesh_clear(&mesh, &cfg);
        mcx_clearcfg(&cfg);  /* closes the log file */
        mcxlab_logfp = NULL;

        /** report the error to MATLAB after the clean-up */
        if (cfgerr[0]) {
            throw std::runtime_error(cfgerr);
        }
    }

    return;
}

/**
 * @brief Set cfg.flog (2 console, 1 stdout, 0 stdout without the banner, or a file name), and send
 * the messages of mmclab there too; called for each configuration before its other fields
 */

static void mcxlab_setflog(const mxArray* item, mcconfig* cfg) {
    int len = mxGetNumberOfElements(item);
    char logfile[MAX_PATH_LENGTH] = {'\0'};

    if (mxIsChar(item)) {
        if (len > 0) {
            mxGetString(item, logfile, MAX_PATH_LENGTH);
            mcxlab_lognull = (mcx_logpath(logfile) != logfile);
            cfg->flog = fopen(mcx_logpath(logfile), "a+");

            if (cfg->flog == NULL) {
                char msg[MAX_PATH_LENGTH + 128];
                snprintf(msg, sizeof(msg), "the log file '%s' of cfg.flog can not be written: %s", logfile, strerror(errno));
                mexErrMsgTxt(msg);
            }
        } else {
            cfg->flog = stdout;
        }
    } else {
        double* val = mxGetPr(item);

        if (len > 0 && val[0] <= 2) {
            cfg->flog = ((int)val[0] == 2 ? stderr : ((int)val[0] == 1 ? stdout : (cfg->printnum = -1, stdout)));
        }
    }

    mcxlab_logfp = cfg->flog;
}

/**
 * @brief The mex entry: runs mcxlab_main and raises its errors in MATLAB/Octave outside any try block
 */

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    char errmsg[MAX_PATH_LENGTH] = "";

    try {
        mcxlab_main(nlhs, plhs, nrhs, prhs);
    } catch (const std::exception& err) {
        snprintf(errmsg, sizeof(errmsg), "%s", err.what());
    } catch (const char* err) {
        snprintf(errmsg, sizeof(errmsg), "%s", err);
    } catch (...) {
        snprintf(errmsg, sizeof(errmsg), "unknown exception");
    }

    if (errmsg[0]) {
        mcxlab_mexerror(errmsg);
    }
}

/**
 * @brief Function to parse one subfield of the input structure
 *
 * This function reads in all necessary information from the cfg input structure.
 * it can handle single scalar inputs, short vectors (3-4 elem), strings and arrays.
 *
 * @param[in] root: the cfg input data structure
 * @param[in] item: the current element of the cfg input data structure
 * @param[in] idx: the index of the current element (starting from 0)
 * @param[out] cfg: the simulation configuration structure to store all input read from the parameters
 * @param[out] mesh: the mesh data structure
 */

void mmc_set_field(const mxArray* root, const mxArray* item, int idx, int jstruct, mcconfig* cfg, tetmesh* mesh) {
    const char* name = mxGetFieldNameByNumber(root, idx);
    const dimtype* arraydim;
    char* jsonshapes = NULL;
    int i, j;

    if (strcmp(name, "nphoton") == 0 && cfg->photonseed != NULL) {
        return;
    }

    GET_1ST_FIELD(cfg, nphoton)
    GET_ONE_FIELD(cfg, nblocksize)
    GET_ONE_FIELD(cfg, nthread)
    GET_ONE_FIELD(cfg, tstart)
    GET_ONE_FIELD(cfg, tstep)
    GET_ONE_FIELD(cfg, tend)
    GET_ONE_FIELD(cfg, isreflect)
    GET_ONE_FIELD(cfg, isspecular)
    GET_ONE_FIELD(cfg, ismomentum)
    GET_ONE_FIELD(cfg, issaveexit)
    GET_ONE_FIELD(cfg, issaveseed)
    GET_ONE_FIELD(cfg, optlevel)
    else if (strcmp(name, "isatomic") == 0) {
        if (mxGetScalar(item) == 0.0) {
            printf("WARNING: cfg.isatomic=0 is deprecated and ignored; non-atomic accumulation raced and lost part of the deposited energy, at no speed benefit\n");
        }
    }

    GET_ONE_FIELD(cfg, basisorder)
    GET_ONE_FIELD(cfg, outputformat)
    GET_ONE_FIELD(cfg, roulettesize)
    GET_ONE_FIELD(cfg, nout)
    GET_ONE_FIELD(cfg, isref3)
    GET_ONE_FIELD(cfg, isnormalized)
    GET_ONE_FIELD(cfg, issave2pt)
    GET_ONE_FIELD(cfg, issaveref)
    GET_ONE_FIELD(cfg, debugphoton)
    GET_ONE_FIELD(cfg, minenergy)
    GET_ONE_FIELD(cfg, replaydet)
    GET_ONE_FIELD(cfg, unitinmm)
    GET_ONE_FIELD(cfg, voidtime)
    GET_ONE_FIELD(cfg, mcmethod)
    GET_ONE_FIELD(cfg, maxdetphoton)
    GET_ONE_FIELD(cfg, maxjumpdebug)
    GET_ONE_FIELD(cfg, srcid)
    GET_ONE_FIELD(cfg, adjointmode)
    GET_ONE_FIELD(cfg, isnodalmua)
    GET_ONE_FIELD(cfg, isnodalmusp)
    GET_VEC3_FIELD(cfg, steps)
    GET_VEC4_FIELD(cfg, detparam1)
    GET_VEC4_FIELD(cfg, detparam2)
    /* srcpos / srcdir / srcparam1 / srcparam2 can be a 1xN row (single source)
     * or an Mx{3,4} matrix (multi-source); first row populates cfg->srcpos etc,
     * remaining rows populate cfg->srcdata[0..M-2]. Mirrors mcxlab.cpp. */
    else if (strcmp(name, "srcpos") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] == 0 || arraydim[1] < 3 || arraydim[1] > 4) {
            MEXERROR("the 'srcpos' field must have 3 or 4 columns (x,y,z,w0)");
        }

        double* val = mxGetPr(item);

        for (int kk = 0; kk < (int)arraydim[1]; kk++) {
            ((float*)(&cfg->srcpos.x))[kk] = val[kk * arraydim[0]];
        }

        printf("mmc.srcpos=[%g %g %g %g];\n", cfg->srcpos.x, cfg->srcpos.y, cfg->srcpos.z, cfg->srcpos.w);

        if (arraydim[0] == 1 && cfg->extrasrclen == 0) {
            return;
        }

        /* mmc convention: srcdata[] holds EVERY source slot, including row 0
         * (the kernel's multi-source dispatch never reads cfg->srcpos). */
        int nrows = (int)arraydim[0];

        if (cfg->extrasrclen && cfg->extrasrclen != nrows) {
            MEXERROR("Length of sub-elements of srcpos/srcdir/srcparam1/srcparam2 must match");
        } else {
            cfg->extrasrclen = nrows;
        }

        if (cfg->srcdata == NULL) {
            cfg->srcdata = (ExtraSrc*)calloc(sizeof(ExtraSrc), cfg->extrasrclen);
        }

        for (int jj = 0; jj < (int)arraydim[1]; jj++)
            for (int ii = 0; ii < cfg->extrasrclen; ii++) {
                ((float*)(&cfg->srcdata[ii].srcpos.x))[jj] = val[jj * arraydim[0] + ii];
            }

        printf("mmc.extrasrclen=%d;\n", cfg->extrasrclen);
    } else if (strcmp(name, "srcdir") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] == 0 || arraydim[1] < 3 || arraydim[1] > 4) {
            MEXERROR("the 'srcdir' field must have 3 or 4 columns (vx,vy,vz,focallength)");
        }

        double* val = mxGetPr(item);

        for (int kk = 0; kk < (int)arraydim[1]; kk++) {
            ((float*)(&cfg->srcdir.x))[kk] = val[kk * arraydim[0]];
        }

        printf("mmc.srcdir=[%g %g %g %g];\n", cfg->srcdir.x, cfg->srcdir.y, cfg->srcdir.z, cfg->srcdir.w);

        if (arraydim[0] == 1 && cfg->extrasrclen == 0) {
            return;
        }

        int nrows = (int)arraydim[0];

        if (cfg->extrasrclen && cfg->extrasrclen != nrows) {
            MEXERROR("Length of sub-elements of srcpos/srcdir/srcparam1/srcparam2 must match");
        } else {
            cfg->extrasrclen = nrows;
        }

        if (cfg->srcdata == NULL) {
            cfg->srcdata = (ExtraSrc*)calloc(sizeof(ExtraSrc), cfg->extrasrclen);
        }

        for (int jj = 0; jj < (int)arraydim[1]; jj++)
            for (int ii = 0; ii < cfg->extrasrclen; ii++) {
                ((float*)(&cfg->srcdata[ii].srcdir.x))[jj] = val[jj * arraydim[0] + ii];
            }

        printf("mmc.extrasrclen=%d;\n", cfg->extrasrclen);
    } else if (strcmp(name, "srcparam1") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] == 0 || (int)arraydim[1] < 1 || (int)arraydim[1] > 4) {
            MEXERROR("the 'srcparam1' field must have 1-4 columns");
        }

        double* val = mxGetPr(item);

        for (int kk = 0; kk < (int)arraydim[1]; kk++) {
            ((float*)(&cfg->srcparam1.x))[kk] = val[kk * arraydim[0]];
        }

        if (arraydim[0] == 1 && cfg->extrasrclen == 0) {
            return;
        }

        int nrows = (int)arraydim[0];

        if (cfg->extrasrclen && cfg->extrasrclen != nrows) {
            MEXERROR("Length of sub-elements of srcpos/srcdir/srcparam1/srcparam2 must match");
        } else {
            cfg->extrasrclen = nrows;
        }

        if (cfg->srcdata == NULL) {
            cfg->srcdata = (ExtraSrc*)calloc(sizeof(ExtraSrc), cfg->extrasrclen);
        }

        for (int jj = 0; jj < (int)arraydim[1]; jj++)
            for (int ii = 0; ii < cfg->extrasrclen; ii++) {
                ((float*)(&cfg->srcdata[ii].srcparam1.x))[jj] = val[jj * arraydim[0] + ii];
            }
    } else if (strcmp(name, "srcparam2") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] == 0 || (int)arraydim[1] < 1 || (int)arraydim[1] > 4) {
            MEXERROR("the 'srcparam2' field must have 1-4 columns");
        }

        double* val = mxGetPr(item);

        for (int kk = 0; kk < (int)arraydim[1]; kk++) {
            ((float*)(&cfg->srcparam2.x))[kk] = val[kk * arraydim[0]];
        }

        if (arraydim[0] == 1 && cfg->extrasrclen == 0) {
            return;
        }

        int nrows = (int)arraydim[0];

        if (cfg->extrasrclen && cfg->extrasrclen != nrows) {
            MEXERROR("Length of sub-elements of srcpos/srcdir/srcparam1/srcparam2 must match");
        } else {
            cfg->extrasrclen = nrows;
        }

        if (cfg->srcdata == NULL) {
            cfg->srcdata = (ExtraSrc*)calloc(sizeof(ExtraSrc), cfg->extrasrclen);
        }

        for (int jj = 0; jj < (int)arraydim[1]; jj++)
            for (int ii = 0; ii < cfg->extrasrclen; ii++) {
                ((float*)(&cfg->srcdata[ii].srcparam2.x))[jj] = val[jj * arraydim[0] + ii];
            }
    } else if (strcmp(name, "e0") == 0) {
        double* val = mxGetPr(item);
        cfg->e0 = val[0];
        printf("mmc.e0=%d;\n", cfg->e0);
    } else if (strcmp(name, "node") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 3) {
            MEXERROR("the 'node' field must have 3 columns (x,y,z)");
        }

        double* val = mxGetPr(item);
        mesh->nn = arraydim[0];

        if (mesh->node) {
            free(mesh->node);
        }

        mesh->node = (FLOAT3*)calloc(sizeof(FLOAT3), mesh->nn);
        mesh->fnode = (FLOAT3*)calloc(sizeof(FLOAT3), mesh->nn);

        for (j = 0; j < 3; j++)
            for (i = 0; i < mesh->nn; i++) {
                ((float*)(&mesh->node[i]))[j] = val[j * mesh->nn + i];
                ((float*)(&mesh->fnode[i]))[j] = val[j * mesh->nn + i];
            }

        printf("mmc.nn=%d;\n", mesh->nn);
    } else if (strcmp(name, "elem") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] < 4) {
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        }

        double* val = mxGetPr(item);
        mesh->ne = arraydim[0];
        mesh->elemlen = arraydim[1];

        if (mesh->elem) {
            free(mesh->elem);
        }

        mesh->elem = (int*)calloc(sizeof(int) * arraydim[1], mesh->ne);

        for (j = 0; j < mesh->elemlen; j++)
            for (i = 0; i < mesh->ne; i++) {
                mesh->elem[i * mesh->elemlen + j] = val[j * mesh->ne + i];
            }

        printf("mmc.elem=[%d,%d];\n", mesh->ne, mesh->elemlen);
    } else if (strcmp(name, "face") == 0) {
        /* a surface-only mesh for OptiX: 3 node indices (1-based), the inner and the outer labels,
           given instead of cfg.elem.  With cfg.elem, cfg.face is not mmc's -- scripts often keep the
           boundary triangles iso2mesh returns there (redbird does) -- and is ignored, as it was
           before surface-only meshes were supported */
        const mxArray* elem = mxGetField(root, jstruct, "elem");

        if (elem != NULL && !mxIsEmpty(elem)) {
            return;
        }

        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 5) {
            MEXERROR("the 'face' field (a surface-only mesh for OptiX) must have 5 columns: 3 node indices, the inner and the outer labels");
        }

        double* val = mxGetPr(item);
        cfg->surfnum = arraydim[0];

        if (cfg->surf) {
            free(cfg->surf);
        }

        cfg->surf = (int*)malloc(sizeof(int) * 5 * cfg->surfnum);

        for (j = 0; j < 5; j++)
            for (i = 0; i < (int)cfg->surfnum; i++) {
                cfg->surf[i * 5 + j] = (int)val[j * cfg->surfnum + i];
            }

        printf("mmc.face=[%d,5];\n", cfg->surfnum);
    } else if (strcmp(name, "noderoi") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'noderoi' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->nn = MAX(arraydim[0], arraydim[1]);

        if (mesh->noderoi) {
            free(mesh->noderoi);
        }

        mesh->noderoi = (float*)malloc(sizeof(float) * mesh->nn);

        for (i = 0; i < mesh->nn; i++) {
            mesh->noderoi[i] = val[i];
        }

        cfg->implicit = 1;
        printf("mmc.noderoi=%d;\n", mesh->nn);
    } else if (strcmp(name, "edgeroi") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 6) {
            MEXERROR("the 'edgeroi' field must have 6 columns (e1,e2)");
        }

        double* val = mxGetPr(item);
        mesh->ne = arraydim[0];

        if (mesh->edgeroi) {
            free(mesh->edgeroi);
        }

        mesh->edgeroi = (float*)calloc(sizeof(float) * arraydim[1], mesh->ne);

        for (j = 0; j < 6; j++)
            for (i = 0; i < mesh->ne; i++) {
                mesh->edgeroi[i * 6 + j] = val[j * mesh->ne + i];
            }

        cfg->implicit = 1;
        printf("mmc.edgeroi=[%d,%d];\n", mesh->ne, 6);
    } else if (strcmp(name, "faceroi") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] != 4) {
            MEXERROR("the 'faceroi' field must have 4 columns (f1,f2,f3,f4)");
        }

        double* val = mxGetPr(item);
        mesh->ne = arraydim[0];

        if (mesh->faceroi) {
            free(mesh->faceroi);
        }

        mesh->faceroi = (float*)calloc(sizeof(float) * arraydim[1], mesh->ne);

        for (j = 0; j < 4; j++)
            for (i = 0; i < mesh->ne; i++) {
                mesh->faceroi[i * 4 + j] = val[j * mesh->ne + i];
            }

        cfg->implicit = 2;
        printf("mmc.faceroi=[%d,%d];\n", mesh->ne, 4);
    } else if (strcmp(name, "elemprop") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'elemprop' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ne = MAX(arraydim[0], arraydim[1]);

        if (mesh->type) {
            free(mesh->type);
        }

        mesh->type = (int*)malloc(sizeof(int ) * mesh->ne);

        for (i = 0; i < mesh->ne; i++) {
            mesh->type[i] = val[i];
        }

        printf("mmc.ne=%d;\n", mesh->ne);
    } else if (strcmp(name, "facenb") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] <= 0 || arraydim[1] < 4) {
            MEXERROR("the 'elem' field must have 4 columns (e1,e2,e3,e4)");
        }

        double* val = mxGetPr(item);
        mesh->ne = arraydim[0];
        mesh->elemlen = arraydim[1];

        if (mesh->facenb) {
            free(mesh->facenb);
        }

        mesh->facenb = (int*)malloc(sizeof(int) * arraydim[1] * mesh->ne);

        for (dimtype j = 0; j < arraydim[1]; j++)
            for (i = 0; i < mesh->ne; i++) {
                mesh->facenb[i * arraydim[1] + j] = val[j * mesh->ne + i];
            }

        printf("mmc.facenb=[%d,%d];\n", mesh->ne, mesh->elemlen);
    } else if (strcmp(name, "evol") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'evol' field can not be empty");
        }

        double* val = mxGetPr(item);
        mesh->ne = MAX(arraydim[0], arraydim[1]);

        if (mesh->evol) {
            free(mesh->evol);
        }

        mesh->evol = (float*)malloc(sizeof(float) * mesh->ne);

        for (i = 0; i < mesh->ne; i++) {
            mesh->evol[i] = val[i];
        }

        printf("mmc.evol=%d;\n", mesh->ne);
    } else if (strcmp(name, "detpos") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] > 0 && arraydim[1] != 4) {
            MEXERROR("the 'detpos' field must have 4 columns (x,y,z,radius)");
        }

        double* val = mxGetPr(item);
        cfg->detnum = arraydim[0];

        if (cfg->detpos) {
            free(cfg->detpos);
        }

        cfg->detpos = (float4*)malloc(cfg->detnum * sizeof(float4));

        for (j = 0; j < 4; j++)
            for (i = 0; i < cfg->detnum; i++) {
                ((float*)(&cfg->detpos[i]))[j] = val[j * cfg->detnum + i];
            }

        printf("mmc.detnum=%d;\n", cfg->detnum);
    } else if (strcmp(name, "prop") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] > 0 && arraydim[1] != 4) {
            MEXERROR("the 'prop' field must have 4 columns (mua,mus,g,n)");
        }

        double* val = mxGetPr(item);
        mesh->prop = arraydim[0] - 1;

        if (mesh->med) {
            free(mesh->med);
        }

        mesh->med = (medium*)calloc(sizeof(medium), mesh->prop + 1);

        for (j = 0; j < 4; j++)
            for (i = 0; i <= mesh->prop; i++) {
                ((float*)(&mesh->med[i]))[j] = val[j * (mesh->prop + 1) + i];
            }

        cfg->his.maxmedia = mesh->prop;
        printf("mmc.prop=%d;\n", mesh->prop);
    } else if (strcmp(name, "debuglevel") == 0) {
        int len = mxGetNumberOfElements(item);
        char buf[MAX_SESSION_LENGTH];

        if (!mxIsChar(item) || len == 0) {
            MEXERROR("the 'debuglevel' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            MEXERROR("the 'debuglevel' field is too long");
        }

        int status = mxGetString(item, buf, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        cfg->debuglevel = mcx_parsedebugopt(buf, debugflag);
        printf("mmc.debuglevel='%s';\n", buf);
    } else if (strcmp(name, "compileropt") == 0) {
        int len = mxGetNumberOfElements(item);

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'compileropt' field must be a non-empty string");
        }

        if (len > MAX_PATH_LENGTH) {
            mexErrMsgTxt("the 'compileropt' field is too long");
        }

        int status = mxGetString(item, cfg->compileropt, MAX_PATH_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        printf("mcx.compileropt='%s';\n", cfg->compileropt);
    } else if (strcmp(name, "kernelfile") == 0) {
        int len = mxGetNumberOfElements(item);

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'kernelfile' field must be a non-empty string");
        }

        if (len >= MAX_PATH_LENGTH) {
            mexErrMsgTxt("the 'kernelfile' field is too long");
        }

        int status = mxGetString(item, cfg->kernelfile, MAX_PATH_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        printf("mcx.kernelfile=string of %d;\n", len);
    } else if (strcmp(name, "srctype") == 0) {
        int len = mxGetNumberOfElements(item);
        const char* srctypeid[] = {"pencil", "isotropic", "cone", "gaussian", "planar", "pattern", "fourier", "arcsine", "disk", "fourierx", "fourierx2d", "zgaussian", "line", "slit", "pencilarray", "pattern3d", "hyperboloid", "ring", ""};
        char strtypestr[MAX_SESSION_LENGTH] = {'\0'};

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'srctype' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            mexErrMsgTxt("the 'srctype' field is too long");
        }

        int status = mxGetString(item, strtypestr, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        cfg->srctype = mcx_keylookup(strtypestr, srctypeid);

        if (cfg->srctype == -1) {
            mexErrMsgTxt("the specified source type is not supported");
        }

        printf("mmc.srctype='%s';\n", strtypestr);
    } else if (strcmp(name, "session") == 0) {
        int len = mxGetNumberOfElements(item);

        if (!mxIsChar(item) || len == 0) {
            MEXERROR("the 'session' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            MEXERROR("the 'session' field is too long");
        }

        int status = mxGetString(item, cfg->session, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        printf("mmc.session='%s';\n", cfg->session);
    } else if (strcmp(name, "srcpattern") == 0) {
        arraydim = mxGetDimensions(item);
        dimtype k, len = mxGetNumberOfElements(item);

        /* srcnum is set once srctype is known, see mmclab_patndim */
        mmclab_patndim = (int)mxGetNumberOfDimensions(item);
        mmclab_patdim0 = (int)arraydim[0];

        double* val = mxGetPr(item);

        if (cfg->srcpattern) {
            free(cfg->srcpattern);
        }

        cfg->srcpattern = (float*)malloc(len * sizeof(float));

        for (k = 0; k < len; k++) {
            cfg->srcpattern[k] = val[k];
        }

        printf("mmc.srcpattern=[%lld %lld %lld];\n", (long long)arraydim[0], (long long)arraydim[1], (long long)(len / (arraydim[0] * arraydim[1])));
    } else if (strcmp(name, "method") == 0) {
        int len = mxGetNumberOfElements(item);
        const char* methods[] = {"plucker", "havel", "badouel", "elem", "grid", ""};
        char methodstr[MAX_SESSION_LENGTH] = {'\0'};

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'method' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            mexErrMsgTxt("the 'method' field is too long");
        }

        int status = mxGetString(item, methodstr, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        cfg->method = mcx_keylookup(methodstr, methods);

        if (cfg->method == -1) {
            mexErrMsgTxt("the specified method is not supported");
        }

        printf("mmc.method='%s';\n", methodstr);
    } else if (strcmp(name, "outputtype") == 0) {
        int len = mxGetNumberOfElements(item);
        const char* outputtype[] = {"flux", "fluence", "energy", "jacobian", "wl", "wp",
                                    "rf", "rfmus", "adjoint", "adjoint_dcoeff", "adjoint_mus", "adjoint_musp", "adjoint_mua_d", "adjoint_mua_musp", ""
                                   };
        char outputstr[MAX_SESSION_LENGTH] = {'\0'};

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'outputtype' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            mexErrMsgTxt("the 'outputtype' field is too long");
        }

        int status = mxGetString(item, outputstr, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        cfg->outputtype = mcx_keylookup(outputstr, outputtype);

        if (cfg->outputtype == -1) {
            mexErrMsgTxt("the specified output type is not supported");
        }

        printf("mmc.outputtype='%s';\n", outputstr);
    } else if (strcmp(name, "omega") == 0) {
        double* val = mxGetPr(item);
        cfg->omega = (float)val[0];
        printf("mmc.omega=%g;\n", cfg->omega);
    } else if (strcmp(name, "detdir") == 0) {
        arraydim = mxGetDimensions(item);

        if (arraydim[0] > 0 && arraydim[1] != 4) {
            MEXERROR("the 'detdir' field must have 4 columns (dx,dy,dz,focal)");
        }

        double* val = mxGetPr(item);
        int nd = arraydim[0];

        if (cfg->detdir) {
            free(cfg->detdir);
        }

        cfg->detdir = (float4*)malloc(nd * sizeof(float4));

        for (j = 0; j < 4; j++)
            for (i = 0; i < nd; i++) {
                ((float*)(&cfg->detdir[i]))[j] = val[j * nd + i];
            }

        printf("mmc.detdir=[%d,4];\n", nd);
    } else if (strcmp(name, "sdmap") == 0) {
        arraydim = mxGetDimensions(item);

        if (cfg->sdmap) {
            free(cfg->sdmap);
            cfg->sdmap = NULL;
        }

        cfg->sdmaplen = 0;

        if (mxGetNumberOfElements(item) == 0) {
            return;
        }

        if (!mxIsDouble(item) || mxIsComplex(item) || arraydim[1] < 2) {
            MEXERROR("the 'sdmap' field must be a real matrix with 2 or more columns (srcid,detid,...)");
        }

        double* val = mxGetPr(item);
        cfg->sdmaplen = arraydim[0];
        cfg->sdmap = (unsigned int*)malloc(sizeof(unsigned int) * 2 * cfg->sdmaplen);

        for (i = 0; i < (int)cfg->sdmaplen; i++) {
            if (val[i] < 1 || val[arraydim[0] + i] < 1) {
                MEXERROR("the source and detector indices in 'sdmap' start from 1");
            }

            cfg->sdmap[i << 1] = (unsigned int)val[i] - 1;
            cfg->sdmap[(i << 1) + 1] = (unsigned int)val[arraydim[0] + i] - 1;
        }

        printf("mmc.sdmap=[%d channels];\n", cfg->sdmaplen);
    } else if (strcmp(name, "nodemua") == 0) {
        /* per-node absorption coefficient array, length mesh.nn; auto-sets isnodalmua.
         * Accepts double or single precision input (matlab default vs. explicit single()). */
        arraydim = mxGetDimensions(item);
        int nn = (int)(arraydim[0] * arraydim[1]);

        if (cfg->nodemua) {
            free(cfg->nodemua);
        }

        cfg->nodemua = (float*)malloc(nn * sizeof(float));

        if (mxIsDouble(item)) {
            double* val = mxGetPr(item);

            for (i = 0; i < nn; i++) {
                cfg->nodemua[i] = (float)val[i];
            }
        } else if (mxIsSingle(item)) {
            float* val = (float*)mxGetData(item);
            memcpy(cfg->nodemua, val, nn * sizeof(float));
        } else {
            MEXERROR("cfg.nodemua must be a double or single-precision array");
        }

        cfg->isnodalmua = 1;
        printf("mmc.nodemua=[%d,1];\n", nn);
    } else if (strcmp(name, "nodemusp") == 0) {
        /* per-node reduced scattering coefficient, length mesh.nn; auto-sets isnodalmusp. */
        arraydim = mxGetDimensions(item);
        int nn = (int)(arraydim[0] * arraydim[1]);

        if (cfg->nodemusp) {
            free(cfg->nodemusp);
        }

        cfg->nodemusp = (float*)malloc(nn * sizeof(float));

        if (mxIsDouble(item)) {
            double* val = mxGetPr(item);

            for (i = 0; i < nn; i++) {
                cfg->nodemusp[i] = (float)val[i];
            }
        } else if (mxIsSingle(item)) {
            float* val = (float*)mxGetData(item);
            memcpy(cfg->nodemusp, val, nn * sizeof(float));
        } else {
            MEXERROR("cfg.nodemusp must be a double or single-precision array");
        }

        cfg->isnodalmusp = 1;
        printf("mmc.nodemusp=[%d,1];\n", nn);
    } else if (strcmp(name, "compute") == 0) {
        int len = mxGetNumberOfElements(item);
        const char* computebackend[] = {"sse", "opencl", "cuda", "optix", ""};
        char computestr[MAX_SESSION_LENGTH] = {'\0'};

        if (!mxIsChar(item) || len == 0) {
            mexErrMsgTxt("the 'compute' field must be a non-empty string");
        }

        if (len > MAX_SESSION_LENGTH) {
            mexErrMsgTxt("the 'compute' field is too long");
        }

        int status = mxGetString(item, computestr, MAX_SESSION_LENGTH);

        if (status != 0) {
            mexWarnMsgTxt("not enough space. string is truncated.");
        }

        cfg->compute = mcx_keylookup(computestr, computebackend);

        if (cfg->compute == -1) {
            mexErrMsgTxt("the specified compute is not supported");
        }

#ifndef USE_OPENCL

        if (cfg->compute == cbOpenCL) {
            mexErrMsgTxt("this mmc mex file was built without OpenCL; rebuild it with 'make mex' or 'make oct'");
        }

#endif
#ifndef USE_CUDA

        if (cfg->compute == cbCUDA) {
            mexErrMsgTxt("this mmc mex file was built without CUDA; rebuild it with 'make cudamex' or 'make cudaoct'");
        }

#endif
#ifndef USE_OPTIX

        if (cfg->compute == cbOptiX) {
            mexErrMsgTxt("this mmc mex file was built without OptiX; rebuild it with 'make optixmex' or 'make optixoct'");
        }

#endif

        printf("mmc.compute='%s';\n", computestr);
    } else if (strcmp(name, "shapes") == 0) {
        int len = mxGetNumberOfElements(item);

        if (!mxIsChar(item) || len == 0) {
            MEXERROR("the 'shapes' field must be a non-empty string");
        }

        jsonshapes = new char[len + 1];
        mxGetString(item, jsonshapes, len + 1);
        jsonshapes[len] = '\0';
    } else if (strcmp(name, "detphotons") == 0) {
        arraydim = mxGetDimensions(item);
        dimdetps[0] = arraydim[0];
        dimdetps[1] = arraydim[1];
        detps = (float*)malloc(arraydim[0] * arraydim[1] * sizeof(float));
        memcpy(detps, mxGetData(item), arraydim[0]*arraydim[1]*sizeof(float));
        printf("mmc.detphotons=[%lld %lld];\n", (long long)arraydim[0], (long long)arraydim[1]);
    } else if (strcmp(name, "seed") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'seed' field can not be empty");
        }

        if (!mxIsUint8(item)) {
            double* val = mxGetPr(item);

            cfg->seed = val[0];
            printf("mmc.seed=%d;\n", cfg->seed);
        } else {
            seedbyte = arraydim[0];
            cfg->photonseed = malloc(arraydim[0] * arraydim[1]);

            if (arraydim[0] != (sizeof(RandType)*RAND_BUF_LEN)) {
                MEXERROR("the row number of cfg.seed does not match RNG seed byte-length");
            }

            memcpy(cfg->photonseed, mxGetData(item), arraydim[0]*arraydim[1]);
            cfg->seed = SEED_FROM_FILE;
            cfg->nphoton = arraydim[1];
            printf("mmc.nphoton=%zu;\n", cfg->nphoton);
        }
    } else if (strcmp(name, "replayweight") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'replayweight' field can not be empty");
        }

        cfg->his.detected = arraydim[0] * arraydim[1];
        cfg->replayweight = (float*)malloc(cfg->his.detected * sizeof(float));
        memcpy(cfg->replayweight, mxGetData(item), cfg->his.detected * sizeof(float));
        printf("mmc.replayweight=%d;\n", cfg->his.detected);
    } else if (strcmp(name, "replaytime") == 0) {
        arraydim = mxGetDimensions(item);

        if (MAX(arraydim[0], arraydim[1]) == 0) {
            MEXERROR("the 'replaytime' field can not be empty");
        }

        cfg->his.detected = arraydim[0] * arraydim[1];
        cfg->replaytime = (float*)malloc(cfg->his.detected * sizeof(float));
        memcpy(cfg->replaytime, mxGetData(item), cfg->his.detected * sizeof(float));
        printf("mmc.replaytime=%d;\n", cfg->his.detected);
    } else if (strcmp(name, "gpuid") == 0) {
        int len = mxGetNumberOfElements(item);

        if (mxIsChar(item)) {
            if (len == 0) {
                mexErrMsgTxt("the 'gpuid' field must be an integer or non-empty string");
            }

            if (len > MAX_DEVICE) {
                mexErrMsgTxt("the 'gpuid' field is too long");
            }

            int status = mxGetString(item, cfg->deviceid, MAX_DEVICE);

            if (status != 0) {
                mexWarnMsgTxt("not enough space. string is truncated.");
            }

            printf("mmc.gpuid='%s';\n", cfg->deviceid);
        } else {
            double* val = mxGetPr(item);
            cfg->gpuid = val[0];
            memset(cfg->deviceid, 0, MAX_DEVICE);

            if (cfg->gpuid > 0 && cfg->gpuid < MAX_DEVICE) {
                memset(cfg->deviceid, '0', cfg->gpuid - 1);
                cfg->deviceid[cfg->gpuid - 1] = '1';
            }

            printf("mmc.gpuid=%d;\n", cfg->gpuid);
        }

        /* a 0/1 mask: '\0' marks the unused devices; a device name pattern ("4070") stays as is */
        if (!mcx_isdevicename(cfg->deviceid)) {
            for (int i = 0; i < MAX_DEVICE; i++)
                if (cfg->deviceid[i] == '0') {
                    cfg->deviceid[i] = '\0';
                }
        }
    } else if (strcmp(name, "workload") == 0) {
        double* val = mxGetPr(item);
        arraydim = mxGetDimensions(item);

        if (arraydim[0]*arraydim[1] > MAX_DEVICE) {
            mexErrMsgTxt("the workload list can not be longer than 256");
        }

        for (dimtype i = 0; i < arraydim[0]*arraydim[1]; i++) {
            cfg->workload[i] = val[i];
        }

        printf("mmc.workload=<<%.0f>>;\n", (double)arraydim[0]*arraydim[1]);
    } else if (strcmp(name, "flog") == 0) {
        /* read by mcxlab_setflog before the other fields, so that their echo follows it */
        printf("mmc.flog=%d;\n", (cfg->flog == stderr) ? 2 : (cfg->printnum < 0 ? 0 : 1));
    } else if (strcmp(name, "isreoriented") == 0) {
        /*internal flag, don't need to do anything*/
    } else {
        printf("WARNING: redundant field '%s'\n", name);
    }

    if (jsonshapes) {
        delete [] jsonshapes;
    }
}

/**
 * @brief Print a brief help information if nothing is provided
 */

mxArray* mmclab_assert(mxArray* arr) {
    if (arr == NULL) {
        mexErrMsgTxt("MMCLAB failed to allocate memory!");
    }

    return arr;
}


/**
 * @brief Run one photon of mmc_run_mp in a try block, see mmc_task_guard in mmc_utils.h
 *
 * An exception can not leave the OpenMP parallel region of mmc_run_mp: the error of a photon
 * is recorded by mcx_error/mesh_error and raised again after the region; the exception is
 * caught here, in the thread that threw it.
 */

static int mmclab_task_guard(void (*task)(void*), void* arg) {
    try {
        task(arg);
    } catch (...) {
        mmc_defer_error(-1, "unknown exception in a simulation thread", __FILE__, __LINE__);
        return 1;
    }

    return 0;
}

/**
 * @brief Error reporting function in the mex function, equivallent to mcx_error in binary mode
 *
 * @param[in] id: a single integer for the types of the error
 * @param[in] msg: the error message string
 * @param[in] filename: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

extern "C" int mmc_throw_exception(const int id, const char* msg, const char* filename, const int linenum) {
    mexPrintf("MMCLAB ERROR (%d): %s in unit %s:%d\n", id, msg, filename, linenum);
    throw std::runtime_error(msg);  /* a copy: msg may be on the stack of the caller */
    return id;
}

/**
 * @brief Print a brief help information if nothing is provided
 */

void mmclab_usage() {
    mexPrintf("MMCLAB " MMC_VERSION "\nUsage:\n    [flux,detphoton]=mmclab(cfg);\n\nPlease run 'help mmclab' for more details.\n");
}

/**
 * @brief Force matlab refresh the command window to print all buffered messages
 */

extern "C" void mcx_matlab_flush() {
#if defined(MATLAB_MEX_FILE)
    mexEvalString("pause(.0001);");
#else
    mexEvalString("fflush(stdout);");
#endif
}

#if defined(__APPLE__)
/**
 * @brief Phantom main function to let macos to build mmclab
 */

int main(void) {
    return 1;
}

#endif
