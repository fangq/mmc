#ifndef _MMC_OPTIX_HOST_H
#define _MMC_OPTIX_HOST_H

#include "mmc_utils.h"
#include "mmc_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

void mmc_run_optix(mcconfig* cfg, tetmesh* mesh, raytracer* tracer);

#ifdef __cplusplus
}
#endif

#endif
