#include <cstdlib>
#ifdef _OPENMP
    #include <omp.h>
#endif

#include "mmc_optix_utils.h"
#include "mmc_optix_host.h"
#include "mmc_cuda_query_gpu.h"

/************************************************************************** In
this unit, we first launch a master thread and initialize the necessary data
structures.This include the command line options(cfg), tetrahedral mesh(mesh)
and the ray tracer precomputed data (tracer).
******************************************************************************/
void mmc_run_optix(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {
    GPUInfo* gpuinfo = NULL;      /** gpuinfo: structure to store GPU information */
    unsigned int activedev = 0;   /** activedev: count of total active GPUs to be used */

    if (!(activedev = mcx_list_cu_gpu(cfg, &gpuinfo))) {
        mcx_error(-1, "No GPU device found\n", __FILE__, __LINE__);
    }

    mmc_error_reset();

#ifdef _OPENMP
    /**
        Now we are ready to launch one thread for each involked GPU to run the simulation
     */
    omp_set_num_threads(activedev);
    #pragma omp parallel
    {
        /* an exception can not leave the parallel region: the error of a GPU thread is recorded
           (see mmc_utils.h) and raised again once all threads have finished */
        try {
            optix_run_simulation(cfg, mesh, tracer, gpuinfo);
        } catch (...) {
            mmc_defer_error(-1, "unknown exception in a GPU thread", __FILE__, __LINE__);
        }
    }

    mmc_raise_deferred();
#else
    /**
        This line runs the main MCX simulation for each GPU inside each thread
     */
    optix_run_simulation(cfg, mesh, tracer, gpuinfo);
#endif
}
