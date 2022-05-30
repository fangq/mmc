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
\file    mmc_host.c

\brief   << Driver program of MMC >>
*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "mmc_host.h"
#include "mmc_tictoc.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

/***************************************************************************//**
In this unit, we first launch a master thread and initialize the
necessary data structures. This include the command line options (cfg),
tetrahedral mesh (mesh) and the ray-tracer precomputed data (tracer).
*******************************************************************************/

/**
 * \brief Initialize simulation configuration structure cfg using command line options
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[in] mesh: the mesh data structure
 * \param[in] tracer: the ray-tracer data structure
 * \param[in] argc: total number of command line option strings
 * \param[in] argv: command line option string array
 */

int mmc_init_from_cmd(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, int argc, char** argv) {
    mcx_initcfg(cfg);
    mcx_parsecmd(argc, argv, cfg);

    if (cfg->isgpuinfo == 0) {
        mesh_init_from_cfg(mesh, cfg);
    }

    return 0;
}

/**
 * \brief Initialize simulation configuration structure cfg using a JSON input file
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[in] mesh: the mesh data structure
 * \param[in] tracer: the ray-tracer data structure
 * \param[in] jcfg: JSON data structure parsed from the input file
 * \param[in] jmesh: JSON data structure parsed from the mesh data
 */

int mmc_init_from_json(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, char* jcfg, char* jmesh) {
    mcx_initcfg(cfg);
    mcx_loadfromjson(jcfg, cfg);
    mesh_init_from_cfg(mesh, cfg); /*need to define mesh load from json*/
    tracer_init(tracer, mesh, cfg->method);
    return 0;
}

/**
 * \brief Rest simulation related data structures
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[out] mesh: the mesh data structure
 * \param[out] tracer: the ray-tracer data structure
 */

int mmc_reset(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {
    mmc_cleanup(cfg, mesh, tracer);
    mcx_initcfg(cfg);
    mesh_init(mesh);
    return 0;
}

/**
 * \brief Clear simulation related data structures
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[out] mesh: the mesh data structure
 * \param[out] tracer: the ray-tracer data structure
 */

int mmc_cleanup(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {
    tracer_clear(tracer);
    mesh_clear(mesh);
    mcx_clearcfg(cfg);
    return 0;
}

/**
 * \brief Peprocessing simulation settings and data for simulation
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[out] mesh: the mesh data structure
 * \param[out] tracer: the ray-tracer data structure
 */

int mmc_prep(mcconfig* cfg, tetmesh* mesh, raytracer* tracer) {
    mcx_prep(cfg);
    tracer_init(tracer, mesh, cfg->method);
    tracer_prep(tracer, cfg);
    return 0;
}

/**
 * \brief Main function to launch MMC photon simulation
 *
 * This is the main loop of the Monte Carlo photon simulation. This function
 * run a complete photon simulation session based on one set of user input.
 *
 * \param[out] cfg: the simulation configuration structure
 * \param[out] mesh: the mesh data structure
 * \param[out] tracer: the ray-tracer data structure
 */

int mmc_run_mp(mcconfig* cfg, tetmesh* mesh, raytracer* tracer, void (*progressfun)(float, void*), void* handle) {
    RandType ran0[RAND_BUF_LEN] __attribute__ ((aligned(16)));
    RandType ran1[RAND_BUF_LEN] __attribute__ ((aligned(16)));
    unsigned int i, j;
    float raytri = 0.f, raytri0 = 0.f;
    unsigned int threadid = 0, ncomplete = 0, t0, dt, debuglevel = 0;
    visitor master = {0.f, 0.f, 0.f, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL};
    visitor_init(cfg, &master);

    if (progressfun == NULL) {
        cfg->debuglevel = cfg->debuglevel & (~dlProgress);
    }

    t0 = StartTimer();

    mcx_printheader(cfg);

#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
    cfg->issaveseed = 0;
#endif
    dt = GetTimeMillis();
    MMCDEBUG(cfg, dlTime, (cfg->flog, "seed=%u\nsimulating ... \n", cfg->seed));

    if (cfg->debugphoton >= 0) {
        debuglevel = cfg->debuglevel;
        cfg->debuglevel &= 0xFFFFEA00;
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
    }

    unsigned int* seeds = NULL;

    /***************************************************************************//**
    The master thread then spawn multiple work-threads depending on your
    OpenMP settings. By default, the total thread number (master + work) is
    your total CPU core number. For example, if you have a dual-core CPU,
    the total thread number is 2; if you have two quad-core CPUs, the total
    thread number is 8. If you want to set the total thread number manually,
    you need to set the OMP_NUM_THREADS environment variable. For example,
    \c OMP_NUM_THREADS=3 sets the total thread number to 3.
    *******************************************************************************/

    /** \subsection ssimu Parallel photon transport simulation */

    #pragma omp parallel private(ran0,ran1,threadid,j)
    {
        visitor visit = {0.f, 0.f, 1.f / cfg->tstep, DET_PHOTON_BUF, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL};
        size_t id;

#ifdef _OPENMP
        unsigned int threadnum = omp_get_num_threads();
#else
        unsigned int threadnum = 1;
#endif

        #pragma omp master
        {
            seeds = (unsigned int*)malloc(sizeof(int) * threadnum * RAND_SEED_WORD_LEN);
            srand(cfg->seed);

            for (i = 0; i < threadnum * RAND_SEED_WORD_LEN; i++) {
                seeds[i] = rand();
            }
        }
        #pragma omp barrier
        visit.reclen = (2 + ((cfg->ismomentum) > 0)) * mesh->prop + (cfg->issaveexit > 0) * 6 + 2;
        visitor_init(cfg, &visit);

        if (cfg->issavedet) {
            if (cfg->issaveseed) {
                visit.photonseed = calloc(visit.detcount, (sizeof(RandType) * RAND_BUF_LEN));
            }

            visit.partialpath = (float*)calloc(visit.detcount * visit.reclen, sizeof(float));
        }

#ifdef _OPENMP
        threadid = omp_get_thread_num();
#endif

        rng_init(ran0, ran1, seeds, threadid);

        if ((cfg->debuglevel & dlProgress) && threadid == 0) {
            progressfun(-0.f, handle);
        }

        /*launch photons*/
        #pragma omp for reduction(+:raytri,raytri0)

        for (id = 0; id < cfg->nphoton; id++) {
            visit.raytet = 0.f;
            visit.raytet0 = 0.f;

            if (id == cfg->debugphoton) {
                cfg->debuglevel = debuglevel;
            }

            if (cfg->seed == SEED_FROM_FILE) {
                onephoton(id, tracer, mesh, cfg, ((RandType*)cfg->photonseed) + id * RAND_BUF_LEN, ran1, &visit);
            } else {
                onephoton(id, tracer, mesh, cfg, ran0, ran1, &visit);
            }

            raytri += visit.raytet;
            raytri0 += visit.raytet0;

            if (id == cfg->debugphoton) {
                cfg->debuglevel &= 0xFFFFEA00;
            }

            #pragma omp atomic
            ncomplete++;

            if ((cfg->debuglevel & dlProgress) && threadid == 0) {
                progressfun((float)ncomplete / cfg->nphoton, handle);
            }
        }

        for (j = 0; j < cfg->srcnum; j++) {
            #pragma omp atomic
            master.launchweight[j] += visit.launchweight[j];
            #pragma omp atomic
            master.absorbweight[j] += visit.absorbweight[j];
        }

        if (cfg->issavedet) {
            #pragma omp atomic
            master.detcount += visit.bufpos;
            #pragma omp barrier

            if (threadid == 0) {
                master.partialpath = (float*)calloc(master.detcount * visit.reclen, sizeof(float));

                if (cfg->issaveseed) {
                    master.photonseed = calloc(master.detcount, (sizeof(RandType) * RAND_BUF_LEN));
                }
            }

            #pragma omp barrier
            #pragma omp critical
            {
                memcpy(master.partialpath + master.bufpos * visit.reclen,
                       visit.partialpath, visit.bufpos * visit.reclen * sizeof(float));

                if (cfg->issaveseed)
                    memcpy((unsigned char*)master.photonseed + master.bufpos * (sizeof(RandType)*RAND_BUF_LEN),
                           visit.photonseed, visit.bufpos * (sizeof(RandType)*RAND_BUF_LEN));

                master.bufpos += visit.bufpos;
            }
        }

        #pragma omp barrier
        visitor_clear(&visit);

        if (visit.partialpath) {
            free(visit.partialpath);
        }

        if (cfg->issaveseed && visit.photonseed) {
            free(visit.photonseed);
        }
    }

    if (seeds) {
        free(seeds);
    }

    /** \subsection sreport Post simulation */

    if ((cfg->debuglevel & dlProgress)) {
        progressfun(1.f, handle);
    }

    dt = GetTimeMillis() - dt;
    MMCDEBUG(cfg, dlProgress, (cfg->flog, "\n"));
    MMCDEBUG(cfg, dlTime, (cfg->flog, "\tdone\t%d\n", dt));
    MMCDEBUG(cfg, dlTime, (cfg->flog, "speed ...\t"S_BOLD""S_BLUE"%.2f photon/ms"S_RESET", %.0f ray-tetrahedron tests (%.0f overhead, %.2f test/ms)\n", (double)cfg->nphoton / dt, raytri, raytri0, raytri / dt));

    if (cfg->issavedet) {
        MMC_FPRINTF(cfg->flog, "detected %d photons\n", master.detcount);
    }

    if (cfg->isnormalized) {
        double cur_normalizer, sum_normalizer = 0;

        for (j = 0; j < cfg->srcnum; j++) {
            cur_normalizer = mesh_normalize(mesh, cfg, master.absorbweight[j], master.launchweight[j], j);
            sum_normalizer += cur_normalizer;
            MMCDEBUG(cfg, dlTime, (cfg->flog, "source %d\ttotal simulated energy: %f\tabsorbed: "S_BOLD""S_BLUE"%5.5f%%"S_RESET"\tnormalizor=%g\n",
                                   j + 1, master.launchweight[j], 100.f * master.absorbweight[j] / master.launchweight[j], cur_normalizer));
        }

        cfg->his.normalizer = sum_normalizer / cfg->srcnum; // average normalizer value for all simulated sources
    }

    if (cfg->issave2pt) {
        switch (cfg->outputtype) {
            case otFlux:
                MMCDEBUG(cfg, dlTime, (cfg->flog, "saving flux ..."));
                break;

            case otFluence:
                MMCDEBUG(cfg, dlTime, (cfg->flog, "saving fluence ..."));
                break;

            case otEnergy:
                MMCDEBUG(cfg, dlTime, (cfg->flog, "saving energy deposit ..."));
                break;
        }

        mesh_saveweight(mesh, cfg, 0);
    }

    if (cfg->issavedet) {
        MMCDEBUG(cfg, dlTime, (cfg->flog, "saving detected photons ..."));

        if (cfg->issaveexit) {
            mesh_savedetphoton(master.partialpath, master.photonseed, master.bufpos, (sizeof(RandType)*RAND_BUF_LEN), cfg);
        }

        if (cfg->issaveexit == 2) {
            float* detimage = (float*)calloc(cfg->detparam1.w * cfg->detparam2.w * cfg->maxgate, sizeof(float));
            mesh_getdetimage(detimage, master.partialpath, master.bufpos, cfg, mesh);
            mesh_savedetimage(detimage, cfg);
            free(detimage);
        }

        free(master.partialpath);

        if (cfg->issaveseed && master.photonseed) {
            cfg->exportseed = (unsigned char*)malloc(cfg->detectedcount * sizeof(RandType) * RAND_BUF_LEN);
            memcpy(cfg->exportseed, master.photonseed, cfg->detectedcount * sizeof(RandType)*RAND_BUF_LEN);
            free(master.photonseed);
        }
    }

    if (cfg->issaveref) {
        MMCDEBUG(cfg, dlTime, (cfg->flog, "saving surface diffuse reflectance ..."));
        mesh_saveweight(mesh, cfg, 1);
    }

    MMCDEBUG(cfg, dlTime, (cfg->flog, "\tdone\t%d\n", GetTimeMillis() - t0));
    visitor_clear(&master);

    return 0;
}
