/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2021
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
\file    mmc.c

\brief   << Main program of MMC >>
*******************************************************************************/

#include "mmc_host.h"

#ifdef USE_OPENCL
    #include "mmc_cl_host.h"
#endif

#ifdef USE_CUDA
    #include "mmc_cu_host.h"
#endif

#ifdef USE_OPTIX
    #include "mmc_optix_host.h"
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

int main(int argc, char** argv) {
    mcconfig cfg;          /** cfg: structure to store all simulation parameters */
    tetmesh mesh;          /** mesh: structure to store mesh information */
    raytracer tracer;      /** tracer: structure to store  */

    /**
       To start an MMC simulation, we first create a simulation configuration,
    initialize all elements to its default settings, and then set user-specified
    settings via command line or input files.
     */
    mmc_init_from_cmd(&cfg, &mesh, &tracer, argc, argv);

    /**
           In the second step, we pre-compute all needed mesh and ray-tracing data
       and get ready for launching photon simulations.
        */
    if (cfg.isgpuinfo == 0) {
        mmc_prep(&cfg, &mesh, &tracer);
    }

    /**
           The core simulation loop is executed in the mmc_run_mp() function where
       multiple threads are executed to simulate all photons.
         */
    if (cfg.compute == cbSSE || cfg.gpuid > MAX_DEVICE) {
        mmc_run_mp(&cfg, &mesh, &tracer, mcx_progressbar, &cfg);
    }

#ifdef USE_CUDA
    else if (cfg.compute == cbCUDA) {
        mmc_run_cu(&cfg, &mesh, &tracer, mcx_progressbar, &cfg);
    }

#endif
#ifdef USE_OPTIX
    else if (cfg.compute == cbOptiX) {
        mmc_run_optix(&cfg, &mesh, &tracer, mcx_progressbar, &cfg);
    }

#endif
#ifdef USE_OPENCL
    else {
        mmc_run_cl(&cfg, &mesh, &tracer, mcx_progressbar, &cfg);
    }

#endif

    /**
           Once all photon simulations are complete, we clean up all allocated memory
       and finish the execution.
         */
    mmc_cleanup(&cfg, &mesh, &tracer);

    return 0;
}
