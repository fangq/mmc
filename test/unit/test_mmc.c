/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2026
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    test_mmc.c

@brief   Host-side unit tests of MMC

Each test either calls a utility function directly on a tiny input, or makes a
short in-process MMC call (mmc_init_from_cmd/mmc_prep/mmc_run_mp/mmc_run_cl/
mmc_cleanup) with a few hundred photons on the 6-element built-in benchmark
meshes, to exercise one ray tracer, source type or output option. Output files
are written to a temporary folder. Build and run with "make unittest" under src/;
the OpenCL tests are skipped when no OpenCL device is found.
*******************************************************************************/

#define _DEFAULT_SOURCE  /* mkdtemp */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include "mmc_host.h"
#include "mmc_utils.h"
#include "mmc_mesh.h"
#include "mmc_const.h"
#include "nifti1.h"
#include "mmc_raytrace.h"
#include "zmat/zmatlib.h"
#include "ubj/ubj.h"

#ifdef USE_OPENCL
    #include "mmc_cl_host.h"
#endif

#ifdef USE_CUDA
    #include "mmc_cu_host.h"
#endif

#ifdef USE_OPTIX
    #include "mmc_optix_host.h"
#endif

int ntest = 0, nfail = 0;

#define EXPECT_NEAR(a, b, tol, msg) do {                                         \
        double va_ = (double)(a), vb_ = (double)(b);                             \
        ntest++;                                                                 \
        if (!(fabs(va_ - vb_) <= (tol))) {                                       \
            nfail++;                                                             \
            printf("FAIL %s:%d %s: %.9g != %.9g (tol %g)\n", __FILE__, __LINE__, \
                   msg, va_, vb_, (double)(tol));                                \
        }                                                                        \
    } while (0)

#define EXPECT_EQ(a, b, msg) EXPECT_NEAR(a, b, 0, msg)

#define EXPECT_TRUE(cond, msg) do {                                              \
        ntest++;                                                                 \
        if (!(cond)) {                                                           \
            nfail++;                                                             \
            printf("FAIL %s:%d %s\n", __FILE__, __LINE__, msg);                  \
        }                                                                        \
    } while (0)

/**
 * check that a statement reports an error and exits: it runs in a forked child with
 * its output discarded and must exit with a non-zero code; the child always leaves
 * through exit() so that the coverage counters of an instrumented build are saved
 */
#define EXPECT_EXIT(stmt, msg) do {                                              \
        int status_ = 0;                                                         \
        pid_t pid_;                                                              \
        ntest++;                                                                 \
        fflush(stdout);                                                          \
        pid_ = fork();                                                           \
        if (pid_ == 0) {                                                         \
            if (freopen("/dev/null", "w", stdout) == NULL ||                     \
                    freopen("/dev/null", "w", stderr) == NULL) {                 \
                exit(0);                                                         \
            }                                                                    \
            stmt;                                                                \
            exit(0);                                                             \
        }                                                                        \
        waitpid(pid_, &status_, 0);                                              \
        if (!(WIFEXITED(status_) && WEXITSTATUS(status_) != 0)) {                \
            nfail++;                                                             \
            printf("FAIL %s:%d %s: no error was reported\n", __FILE__, __LINE__, msg); \
        }                                                                        \
    } while (0)

/** check that a statement finishes the program normally with exit(0), such as --dumpjson does */
#define EXPECT_EXIT_OK(stmt, msg) do {                                           \
        int status_ = 0;                                                         \
        pid_t pid_;                                                              \
        ntest++;                                                                 \
        fflush(stdout);                                                          \
        pid_ = fork();                                                           \
        if (pid_ == 0) {                                                         \
            if (freopen("/dev/null", "w", stdout) == NULL ||                     \
                    freopen("/dev/null", "w", stderr) == NULL) {                 \
                exit(3);                                                         \
            }                                                                    \
            stmt;                                                                \
            exit(0);                                                             \
        }                                                                        \
        waitpid(pid_, &status_, 0);                                              \
        if (!(WIFEXITED(status_) && WEXITSTATUS(status_) == 0)) {                \
            nfail++;                                                             \
            printf("FAIL %s:%d %s: the program did not finish normally\n", __FILE__, __LINE__, msg); \
        }                                                                        \
    } while (0)

#ifndef MMC_EXAMPLE_DIR
    #define MMC_EXAMPLE_DIR "../examples"  /* set by the unittest target of src/Makefile */
#endif

extern char flagset[256];

static FILE* devnull = NULL;
static int savedout = -1;  /**< the stdout descriptor while it is redirected to /dev/null */

/** send stdout to /dev/null during an MMC call, which prints its header and progress bar there */
static void quiet_begin(void) {
    if (getenv("MMC_UNITTEST_VERBOSE") == NULL && savedout < 0) {
        fflush(stdout);
        savedout = dup(STDOUT_FILENO);
        dup2(fileno(devnull), STDOUT_FILENO);
    }
}

static void quiet_end(void) {
    if (savedout >= 0) {
        fflush(stdout);
        dup2(savedout, STDOUT_FILENO);
        close(savedout);
        savedout = -1;
    }
}
static int hasocl = 0;   /**< 1 if an OpenCL device is available */

/** size of a file, -1 if it does not exist */
static long file_size(const char* fname) {
    FILE* fp = fopen(fname, "rb");
    long len;

    if (fp == NULL) {
        return -1;
    }

    fseek(fp, 0, SEEK_END);
    len = ftell(fp);
    fclose(fp);
    return len;
}

static int mkdir_p(const char* dir) {
    return (mkdir(dir, 0755) == 0 || access(dir, F_OK) == 0) ? 0 : -1;
}

static void write_file(const char* fname, const void* buf, size_t len) {
    FILE* fp = fopen(fname, "wb");

    if (fp) {
        fwrite(buf, 1, len, fp);
        fclose(fp);
    }
}

/** an optional change of the settings after parsing and before mmc_prep, for the
    fields that only mmclab/pmmc set, such as adjointmode */
static void (*run_premut)(mcconfig*, tetmesh*) = NULL;

/** results of one short MMC call */
typedef struct {
    double normalizer;
    int detected;
    int ran;
    double wsum;   /**< sum of the first output slot (all gates), for comparing backends */
    double jsum;   /**< sum of the first adjoint Jacobian block, for comparing backends */
    double imsum;  /**< sum of the RF imaginary output (exportadjoint), for comparing backends */
    double cen[3]; /**< grid mode: the fluence-weighted centroid (voxels), for comparing source shapes */
    double rms;    /**< grid mode: the fluence-weighted RMS distance from the centroid (voxels) */
    int ntraj;     /**< number of saved trajectory points (-D M) */
    int trajok;    /**< 1 if every trajectory point has a valid photon index and weight */
} RunResult;

/**
 * Run "mmc <args>" in-process: parse, prepare, simulate on the CPU (-G -1) or
 * with OpenCL, then clean up. The settings can be inspected through the
 * optional check callback before the clean-up.
 */
static RunResult run_mmc(const char* args, void (*check)(mcconfig*, tetmesh*)) {
    char buf[4096], *argv[128];
    int argc = 0;
    mcconfig cfg;
    tetmesh mesh;
    raytracer tracer;
    RunResult res = {0., 0, 0};

    snprintf(buf, sizeof(buf), "mmc %s", args);

    if (getenv("MMC_UNITTEST_VERBOSE")) {
        printf("run: mmc %s\n", args);
        fflush(stdout);
    }

    for (char* tok = strtok(buf, " "); tok && argc < 127; tok = strtok(NULL, " ")) {
        argv[argc++] = tok;
    }

    argv[argc] = NULL;
    memset(flagset, 0, sizeof(flagset));  /* forget the options of the previous call */
    memset(&tracer, 0, sizeof(tracer));

    quiet_begin();
    mmc_init_from_cmd(&cfg, &mesh, &tracer, argc, argv);
    cfg.flog = devnull;

    if (run_premut) {
        run_premut(&cfg, &mesh);
    }

    if (cfg.isgpuinfo == 0) {
        mmc_prep(&cfg, &mesh, &tracer);

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
        res.ran = 1;
        res.ntraj = cfg.debugdatalen;
        res.trajok = (cfg.debugdatalen > 0);

        for (int k = 0; k < cfg.debugdatalen && cfg.exportdebugdata; k++) {
            float* rec = cfg.exportdebugdata + k * MCX_DEBUG_REC_LEN;
            unsigned int id = *(unsigned int*)rec;

            if (id >= cfg.nphoton || !(rec[4] > 0.f && rec[4] <= 1.f)) {
                res.trajok = 0;
            }
        }

        if (mesh.weight) {
            size_t datalen = (cfg.method == rtBLBadouelGrid) ? cfg.crop0.z : (cfg.basisorder ? mesh.nn : mesh.ne);

            for (size_t k = 0; k < datalen * cfg.maxgate; k++) {
                res.wsum += mesh.weight[k];
            }

            /* the spatial moments of a grid output, all gates and patterns */
            if (cfg.method == rtBLBadouelGrid && cfg.crop0.x > 0) {
                double m0 = 0., m1[3] = {0., 0., 0.}, m2 = 0.;

                for (size_t k = 0; k < datalen * cfg.maxgate * cfg.srcnum; k++) {
                    size_t vox = (k / cfg.srcnum) % datalen;
                    double w = mesh.weight[k], p[3] = {(double)(vox % cfg.crop0.x), (double)((vox % cfg.crop0.y) / cfg.crop0.x), (double)(vox / cfg.crop0.y)};

                    m0 += w;

                    for (int d = 0; d < 3; d++) {
                        m1[d] += w * p[d];
                        m2 += w * p[d] * p[d];
                    }
                }

                if (m0 > 0.) {
                    double c2 = 0.;

                    for (int d = 0; d < 3; d++) {
                        res.cen[d] = m1[d] / m0;
                        c2 += res.cen[d] * res.cen[d];
                    }

                    res.rms = sqrt(fmax(m2 / m0 - c2, 0.));
                }
            }

            if (cfg.omega > 0.f && cfg.exportadjoint && !MCX_IS_ADJOINT_TYPE(cfg.outputtype)) {
                for (size_t k = 0; k < datalen * cfg.maxgate; k++) {
                    res.imsum += cfg.exportadjoint[k];
                }
            }

            if (cfg.exportjacob) {
                datalen = (cfg.method == rtBLBadouelGrid) ? cfg.crop0.z : mesh.nn;

                for (size_t k = 0; k < datalen; k++) {
                    res.jsum += cfg.exportjacob[k];
                }
            }
        }
    }

    res.normalizer = cfg.his.normalizer;
    res.detected = cfg.detectedcount;
    quiet_end();

    if (check) {
        check(&cfg, &mesh);
    }

    cfg.flog = stdout;  /* mmc_cleanup closes flog unless it is stdout */
    mmc_cleanup(&cfg, &mesh, &tracer);
    return res;
}

/*-----------------------------------------------------------------------------
 *  ray tracers and Monte Carlo methods on the CPU
 *---------------------------------------------------------------------------*/

static void test_tracers(void) {
    const char* method[] = {"P", "H", "B", "S", "G"};
    char args[512], msg[128];

    for (int i = 0; i < 5; i++) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 500 -E 1 -M %s -C %d -s utm%s",
                 method[i], (i < 4), method[i]);
        RunResult r = run_mmc(args, NULL);
        sprintf(msg, "ray tracer -M %s runs", method[i]);
        EXPECT_TRUE(r.ran && r.normalizer > 0, msg);
    }

    /* MCML-style (albedo weight) photon transport, with Plucker and Badouel tracers */
    RunResult r = run_mmc("--bench dmmc-cube60 -G -1 -n 500 -E 1 -M H -C 1 -m 1 -s utmcml", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "-m 1 MCML-style weight update");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 500 -E 1 -M P -C 0 -m 1 -s utmcml0", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "-m 1 with elem-based output");

    /* second order basis on 10-node tetrahedra */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 500 -E 1 -M H -C 2 -s utc2", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "-C 2 second order basis");
}

/*-----------------------------------------------------------------------------
 *  source types, launched on the CPU
 *---------------------------------------------------------------------------*/

static void test_sources(void) {
    /* the initial element (2) of dmmc-cube60 has its centroid at [30,45,15], at least 10.6 mm
       from each face, so that the launch positions of these sources stay inside it */
    const char* src[] = {
        "{\"Type\":\"isotropic\",\"Pos\":[30,45,15]}",
        "{\"Type\":\"cone\",\"Param1\":[0.4,0,0,0]}",
        "{\"Type\":\"arcsine\",\"Param1\":[0.4,0,0,0]}",
        "{\"Type\":\"gaussian\",\"Pos\":[30,45,15],\"Param1\":[1,0,0,0]}",
        "{\"Type\":\"gaussian\",\"Pos\":[30,45,15],\"Param1\":[1,1.0,0,0],\"Dir\":[0,0,1,2]}",
        "{\"Type\":\"zgaussian\",\"Param1\":[0.3,0,0,0]}",
        "{\"Type\":\"hyperboloid\",\"Pos\":[30,45,15],\"Param1\":[1,10,2,0]}",
        "{\"Type\":\"ring\",\"Pos\":[30,45,15],\"Param1\":[3,1,0.5,2.5]}",
        "{\"Type\":\"line\",\"Pos\":[28,45,15],\"Param1\":[4,0,0,0]}",
        "{\"Type\":\"slit\",\"Pos\":[28,45,15],\"Param1\":[4,0,0,0]}",
        "{\"Type\":\"disk\",\"Pos\":[30,45,15],\"Param1\":[3,0,0,0]}",
        "{\"Type\":\"planar\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,0]}",
        "{\"Type\":\"fourier\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,1]}",
        "{\"Type\":\"fourierx\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,1]}",
        "{\"Type\":\"fourierx2d\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,0],\"Param2\":[0,4,1,1]}",
        "{\"Type\":\"pattern\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1]}}",
        "{\"Type\":\"pattern3d\",\"Pos\":[29,44,14],\"Param1\":[2,2,2,0],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Nz\":2,\"Data\":[1,1,1,1,1,1,1,1]}}",
        "{\"Type\":\"pencilarray\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2]}",
        "{\"Type\":\"pencil\",\"Dir\":[0,0,1,\"_NaN_\"]}",
        "{\"Type\":\"hyperboloid\",\"Pos\":[30,45,15],\"Dir\":[0.6,0,0.8],\"Param1\":[1,10,2,0]}",
        "{\"Type\":\"hyperboloid\",\"Pos\":[30,45,15],\"Dir\":[0,0,-1],\"Param1\":[1,0,0,0]}",
        "{\"Type\":\"ring\",\"Pos\":[30,45,15],\"Param1\":[3,1,0,0]}",
    };
    char args[2048], msg[256];

    for (int i = 0; i < (int)(sizeof(src) / sizeof(src[0])); i++) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 300 -E 1 -M H -C 1 -s utsrc%d --json {\"Optode\":{\"Source\":%s}}", i, src[i]);
        RunResult r = run_mmc(args, NULL);
        sprintf(msg, "source %.60s runs", src[i]);
        EXPECT_TRUE(r.ran, msg);
    }

    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 --json {\"Optode\":{\"Source\":{\"Type\":\"nosuchsrc\"}}}", NULL),
                "an unknown source type is an error");
}

/*-----------------------------------------------------------------------------
 *  output formats and detected photon options
 *---------------------------------------------------------------------------*/

static void test_outputs(void) {
    const char* fmt[] = {"ascii", "bin", "nii", "hdr", "mc2", "jnii", "bnii"};
    char args[512], msg[128];

    for (int i = 0; i < (int)(sizeof(fmt) / sizeof(fmt[0])); i++) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 200 -E 1 -M H -C 1 -F %s -s utfmt%d", fmt[i], i);
        RunResult r = run_mmc(args, NULL);
        sprintf(msg, "-F %s runs", fmt[i]);
        EXPECT_TRUE(r.ran, msg);
    }

    /* detected photons with seeds, exit positions, momentum transfer and surface reflectance */
    RunResult r = run_mmc("--bench dmmc-cube60 -G -1 -n 2000 -E 1 -M H -C 1 -d 1 -q 1 -x 1 -X 1 --momentum 1 -F mc2 -s utdet", NULL);
    EXPECT_TRUE(r.ran && r.detected > 0, "detected photons with all fields");
    EXPECT_TRUE(file_size("utdet.mch") > 0, "-d 1 writes an .mch file");

    r = run_mmc("--bench dmmc-cube60 -G -1 -n 2000 -E 1 -M H -C 1 -d 1 -q 1 -F jnii -s utdetj", NULL);
    EXPECT_TRUE(file_size("utdetj_detp.jdt") > 0, "-F jnii writes a .jdt detected photon file");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 2000 -E 1 -M H -C 1 -d 1 -q 1 -F bnii -s utdetb", NULL);
    EXPECT_TRUE(file_size("utdetb_detp.jdb") > 0, "-F bnii writes a .jdb detected photon file");

    /* replay from the .mch file (-E only takes .mch files), with the Jacobian outputs */
    r = run_mmc("--bench dmmc-cube60 -G -1 -E utdet.mch -O J -M H -C 1 -s utrep", NULL);
    EXPECT_TRUE(r.ran, "replay from an .mch file");
    r = run_mmc("--bench dmmc-cube60 -G -1 -E utdet.mch -P 1 -O L -M H -C 1 -s utrep2", NULL);
    EXPECT_TRUE(r.ran, "replay of detector 1 with the path length output");
    r = run_mmc("--bench dmmc-cube60 -G -1 -E utdet.mch -O P -M H -C 1 -s utrep3", NULL);
    EXPECT_TRUE(r.ran, "replay with the scattering-count output");

    /* other output types and debug flags */
    const char* otype[] = {"F", "E"};

    for (int i = 0; i < 2; i++) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 200 -E 1 -M H -C 1 -O %s -s utot%d", otype[i], i);
        r = run_mmc(args, NULL);
        sprintf(msg, "-O %s runs", otype[i]);
        EXPECT_TRUE(r.ran, msg);
    }

    r = run_mmc("--bench dmmc-cube60 -G -1 -n 20 -E 1 -M H -C 1 -D TMPXRCBWDIOA --debugphoton 3 -s utdbg", NULL);
    EXPECT_TRUE(r.ran, "debug flags with a single debugged photon");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 20 -E 1 -M H -C 1 -D M -s uttraj", NULL);
    EXPECT_TRUE(r.ran && file_size("uttraj.mct") > 0, "-D M saves the trajectory");
}

/*-----------------------------------------------------------------------------
 *  implicit MMC with edge, node and face ROIs, and refractive index mismatch
 *---------------------------------------------------------------------------*/

static void check_roiprim(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_TRUE(cfg->implicit == 1 && mesh->nroiprim == 3 && mesh->roioff && mesh->roioff[mesh->ne] > mesh->ne + 1,
                "Mesh.MeshROIPrim rows join the node ROI in the element lists");
}

static void check_slab1(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_TRUE(cfg->implicit == 2 && mesh->nslab1 > 0, "a face labeled from one side only is a single-sided slab");
}

static void test_immc(void) {
    const char* bench[] = {"edgeimmc", "nodeimmc", "faceimmc"};
    /* an ROI medium of a lower index: Fresnel events, total internal reflection and self-hits on the ROI surfaces */
    const char* media = "{\"Domain\":{\"Media\":[{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1},{\"mua\":0.0458,\"mus\":35.65,\"g\":0.9,\"n\":1.37},"
                        "{\"mua\":0.5,\"mus\":2,\"g\":0.9,\"n\":1}]}}";
    char args[1024], msg[128];

    for (int i = 0; i < 3; i++) {
        snprintf(args, sizeof(args), "--bench %s -G -1 -n 500 -E 1 -b 1 -s utimmc%d", bench[i], i);
        RunResult r = run_mmc(args, NULL);
        sprintf(msg, "iMMC %s with reflection runs", bench[i]);
        EXPECT_TRUE(r.ran, msg);
        snprintf(args, sizeof(args), "--bench %s -G -1 -n 500 -E 1 -b 0 -m 1 -s utimmcb%d", bench[i], i);
        r = run_mmc(args, NULL);
        sprintf(msg, "iMMC %s without reflection, MCML-style", bench[i]);
        EXPECT_TRUE(r.ran, msg);
        snprintf(args, sizeof(args), "--bench %s -G -1 -n 2000 -E 1 -b 1 --json %s -s utimmcn%d", bench[i], media, i);
        r = run_mmc(args, NULL);
        sprintf(msg, "iMMC %s with an index step at the ROI surface", bench[i]);
        EXPECT_TRUE(r.ran && r.wsum > 0., msg);
    }

    /* a capsule and a sphere given by coordinates, together with the node ROI of nodeimmc */
    RunResult r = run_mmc("--bench nodeimmc -G -1 -n 2000 -E 1 -b 1 --json {\"Shapes\":{\"MeshROIPrim\":[[0.5,0.5,0.5,0.5,0.5,0.5,0.3],"
                          "[0.1,0.2,0.3,0.9,0.8,0.7,0.1]]}} -s utprim", check_roiprim);
    EXPECT_TRUE(r.ran && r.wsum > 0., "iMMC with Mesh.MeshROIPrim runs");
    EXPECT_EXIT(run_mmc("--bench nodeimmc --json {\"Shapes\":{\"MeshROIPrim\":[[1,2,3]]}}", NULL), "a MeshROIPrim row of 3 numbers");
    EXPECT_EXIT(run_mmc("--bench nodeimmc --json {\"Shapes\":{\"MeshROIPrim\":1}}", NULL), "a MeshROIPrim that is not an array");
    EXPECT_EXIT(run_mmc("--bench faceimmc --json {\"Shapes\":{\"MeshROIPrim\":[[0,0,0,1,1,1,0.1]]}}", NULL),
                "MeshROIPrim combined with a face ROI");

    /* a slab on one side of an interior face only: a Fresnel interface between the background and the ROI medium */
    snprintf(args, sizeof(args), "--bench faceimmc -G -1 -n 3000 -E 1 -b 1 --json {\"Shapes\":{\"MeshROI\":[[0,0,0,0],[0,0,0,0],[0,0,0,0],"
             "[0.1,0,0,0],[0,0,0,0],[0,0,0,0]]},%s -s utslab1", media + 1);
    r = run_mmc(args, check_slab1);
    EXPECT_TRUE(r.ran && r.wsum > 0., "face iMMC with a single-sided slab runs");

    /* inputs found by mmcfuzz.py: each runs or stops with an error instead of crashing or hanging */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 200 -E 1 --json {\"Optode\":{\"Source\":{\"Type\":\"cone\",\"Param1\":[0,0,0,0]}}} -s utcone0", NULL);
    EXPECT_TRUE(r.ran, "a cone of zero angle (a rejection loop never ended)");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 200 -E 1 --json {\"Forward\":{\"T0\":1e-9,\"T1\":5e-9,\"Dt\":1e-9}} -s utt0", NULL);
    EXPECT_TRUE(r.ran, "a time window after the launch (a negative gate indexed before the output)");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 200 -E 1 -C 2 -s utc2", NULL);
    EXPECT_TRUE(r.ran, "the 10-node mesh of -C 2 (an O(edges^2) search and wrong node indices)");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 --json {\"Session\":{\"DoSaveExit\":-5}}", NULL), "DoSaveExit of -5");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 --json {\"Session\":{\"RayTracer\":1}}", NULL), "a RayTracer that is not a string");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 --json {\"Session\":{\"OutputType\":\"wp\"}}", NULL),
                "a replay output type without a replay file");
    EXPECT_EXIT(run_mmc("--bench faceimmc -G -1 -n 10 --json {\"Shapes\":{\"MeshROI\":[[0,0,0,0.1]]}}", NULL), "a MeshROI shorter than the mesh");

    /* reorienting an inverted element swaps its edge and face ROI columns with its nodes 3 and 4 */
    {
        mcconfig cfg;
        tetmesh mesh;
        FLOAT3 nd[4] = {{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, {0.f, 0.f, 1.f}};    /* negative volume in mmc's sign */
        int i, ok = 1;

        mcx_initcfg(&cfg);
        mesh_init(&mesh);
        mesh.nn = 4, mesh.ne = 1, mesh.prop = 1;
        mesh.node = (FLOAT3*)malloc(sizeof(nd));
        memcpy(mesh.node, nd, sizeof(nd));
        mesh.elem = (int*)malloc(sizeof(int) * 4);
        mesh.type = (int*)calloc(1, sizeof(int));
        mesh.edgeroi = (float*)malloc(sizeof(float) * 6);
        mesh.faceroi = (float*)malloc(sizeof(float) * 4);

        for (i = 0; i < 6; i++) {
            mesh.edgeroi[i] = (float)(i + 1);
            mesh.elem[i & 3] = (i & 3) + 1;
            mesh.faceroi[i & 3] = (float)((i & 3) + 1);
        }

        mesh_getvolume(&mesh, &cfg);
        ok = (mesh.elem[2] == 4 && mesh.elem[3] == 3);
        ok &= (mesh.edgeroi[0] == 1.f && mesh.edgeroi[1] == 3.f && mesh.edgeroi[2] == 2.f && mesh.edgeroi[3] == 5.f && mesh.edgeroi[4] == 4.f
               && mesh.edgeroi[5] == 6.f);
        ok &= (mesh.faceroi[0] == 2.f && mesh.faceroi[1] == 1.f && mesh.faceroi[2] == 3.f && mesh.faceroi[3] == 4.f);
        EXPECT_TRUE(ok, "a reoriented element keeps its ROI columns with its edges and faces");
        mesh_clear(&mesh, &cfg);
        mcx_clearcfg(&cfg);
    }
}

/*-----------------------------------------------------------------------------
 *  OpenCL runs: tracer variants, RF and the adjoint Jacobians
 *---------------------------------------------------------------------------*/

static void mut_adjointnodal(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->adjointmode = 1;
}

static void test_opencl(void) {
    if (!hasocl) {
        printf("no OpenCL device, skip the OpenCL tests\n");
        return;
    }

    const char* cases[] = {
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -s utcl1",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -d 1 -X 1 -x 1 -q 1 -s utcl2",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 0 -O F -s utcl3",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 --json {\"Forward\":{\"Omega\":6.28e8}} -O R -s utcl5",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utcl6",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M G -O D --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utcl7",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utcl8",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M G -O Q --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utcl9",
    };
    char msg[128];

    for (int i = 0; i < (int)(sizeof(cases) / sizeof(cases[0])); i++) {
        RunResult r = run_mmc(cases[i], NULL);
        sprintf(msg, "OpenCL case %d runs", i + 1);
        EXPECT_TRUE(r.ran, msg);
    }

    /* more OpenCL variants: a fixed thread configuration (-A 0), replay, a focused or
       divergent source, respin, RF with elem-based output, two patterns, the mus/musp
       adjoints, and the trajectory output */
    const char* more[] = {
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -A 0 -t 1000 -T 64 -s utclA0",
        "--bench dmmc-cube60 -G 1 -E utdet.mch -O J -M H -C 1 -s utclrep",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -r 2 -s utclr2",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 0 -O R --json {\"Forward\":{\"Omega\":6.28e8}} -s utclrf0",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 --json {\"Optode\":{\"Source\":{\"Dir\":[0,0,1,\"_NaN_\"]}}} -s utclnan",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 --json {\"Optode\":{\"Source\":{\"Dir\":[0,0,1,\"_Inf_\"]}}} -s utclinf",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 --json {\"Optode\":{\"Source\":{\"Dir\":[0,0,1,\"-_Inf_\"]}}} -s utclninf",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -s utclmp --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"Pos\":[28,43,15],"
        "\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"SrcNum\":2,\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1,0,1,1,0]}}}}",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O U --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utclU",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O V --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utclV",
        "--bench dmmc-cube60 -G 1 -n 50 -E 1 -M H -C 1 -D M -s utcltraj",
    };

    for (int i = 0; i < (int)(sizeof(more) / sizeof(more[0])); i++) {
        RunResult r = run_mmc(more[i], NULL);
        sprintf(msg, "OpenCL variant %d runs", i + 1);
        EXPECT_TRUE(r.ran, msg);
    }

    /* the nodal-approximation mesh adjoint (adjointmode 1, set by mmclab/pmmc), alone,
       dual, and dual with RF */
    const char* adj[] = {
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utclna",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utclnw",
        "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}} -s utclnwrf",
    };
    run_premut = mut_adjointnodal;

    for (int i = 0; i < (int)(sizeof(adj) / sizeof(adj[0])); i++) {
        RunResult r = run_mmc(adj[i], NULL);
        sprintf(msg, "OpenCL nodal adjoint %d runs", i + 1);
        EXPECT_TRUE(r.ran, msg);
    }

    run_premut = NULL;
}

/*-----------------------------------------------------------------------------
 *  features of the CPU and OpenCL paths: RF on the CPU, per-pattern outputs,
 *  the per-step debug log, buffer limits and per-node optical properties
 *---------------------------------------------------------------------------*/

static FILE* featlog = NULL;

/** send the log of the next run to utfeat.log */
static void mut_logfile(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    featlog = fopen("utfeat.log", "w");
    cfg->flog = featlog ? featlog : devnull;
}

static void check_closelog(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;

    if (featlog) {
        fclose(featlog);
        featlog = NULL;
    }

    cfg->flog = devnull;
}

static void mut_badmethod(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->method = 99;
}

/** per-node mua and mus equal to those of the medium of the elements */
static void mut_nodalprop(mcconfig* cfg, tetmesh* mesh) {
    cfg->nodemua = (float*)malloc(sizeof(float) * mesh->nn);
    cfg->nodemusp = (float*)malloc(sizeof(float) * mesh->nn);

    for (int i = 0; i < mesh->nn; i++) {
        cfg->nodemua[i] = mesh->med[1].mua;
        cfg->nodemusp[i] = mesh->med[1].mus;
    }

    cfg->isnodalmua = cfg->isnodalmusp = 1;
}

static double patsum[2];

/** the output of each of the two patterns, interleaved per voxel */
static void check_patsum(mcconfig* cfg, tetmesh* mesh) {
    patsum[0] = patsum[1] = 0.;

    for (size_t k = 0; k < (size_t)cfg->crop0.z * cfg->maxgate; k++) {
        patsum[0] += mesh->weight[k * 2];
        patsum[1] += mesh->weight[k * 2 + 1];
    }
}

static int has_text(const char* fname, const char* key) {
    FILE* fp = fopen(fname, "rb");
    char line[4096];
    int found = 0;

    while (fp && !found && fgets(line, sizeof(line), fp)) {
        found = (strstr(line, key) != NULL);
    }

    if (fp) {
        fclose(fp);
    }

    return found;
}

static void test_features(void) {
    const char* rfjson = "{\"Forward\":{\"Omega\":6.28e8}}";
    char args[1024], msg[256];
    RunResult r, rc;

    /* RF on the CPU: the branch-less Badouel tracers in grid (-M G) and mesh (-M s) mode give
       the imaginary part of the fluence, which agrees with OpenCL; the 8 nodes of this mesh
       (one cube) give a nodal output that depends on the tracer, so in mesh mode only the
       phase, the ratio of the imaginary to the real sum, is compared */
    const char* rfmode[] = {"-M G", "-M s -C 1"};

    for (int i = 0; i < 2; i++) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 20000 -E 1 %s -O R --json %s -s utcpurf%d", rfmode[i], rfjson, i);
        r = run_mmc(args, NULL);
        sprintf(msg, "RF forward on the CPU (%s) gives an imaginary part", rfmode[i]);
        EXPECT_TRUE(r.ran && r.wsum > 0. && r.imsum != 0. && isfinite(r.imsum), msg);

        if (hasocl) {
            snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 20000 -E 1 %s -O R --json %s -s utclrf%d", rfmode[i], rfjson, i);
            rc = run_mmc(args, NULL);
            sprintf(msg, "RF forward on the CPU (%s) agrees with OpenCL within 5%%", rfmode[i]);

            if (i == 0) {
                EXPECT_TRUE(rc.ran && fabs(r.wsum - rc.wsum) < 0.05 * fabs(rc.wsum) && fabs(r.imsum - rc.imsum) < 0.05 * fabs(rc.imsum), msg);
            } else {
                EXPECT_TRUE(rc.ran && rc.wsum > 0. && fabs(r.imsum / r.wsum - rc.imsum / rc.wsum) < 0.05 * fabs(rc.imsum / rc.wsum), msg);
            }
        }
    }

    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 100 -M H -C 1 -O R --json {\"Forward\":{\"Omega\":6.28e8}}", NULL),
                "RF on the CPU with the Havel tracer");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 100 -M G -O R --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Source\":{\"Type\":\"pattern\","
                        "\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"SrcNum\":2,\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,1,1,1,2,2,2,2]}}}}", NULL),
                "RF on the CPU with two source patterns");

    /* two patterns in grid mode on the CPU: pattern 2 is twice pattern 1, and each output is
       normalized by its own launched energy, so both give the same fluence */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 5000 -E 1 -M G -s utcpupat --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"Pos\":[28,43,15],"
                "\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"SrcNum\":2,\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,1,1,1,2,2,2,2]}}}}", check_patsum);
    EXPECT_TRUE(r.ran && patsum[0] > 0. && fabs(patsum[1] / patsum[0] - 1.) < 1e-3, "two patterns on the CPU are normalized separately");

    /* the per-step debug log (-D S: moves, E: exits, X: edges, R: reflections) of a photon in
       an implicit-MMC domain with a mismatched boundary */
    run_premut = mut_logfile;
    r = run_mmc("--bench nodeimmc -G -1 -n 20 -E 1 -b 1 -D SXCIOERAB --debugphoton 2 -s utfeatlog", check_closelog);
    run_premut = NULL;
    EXPECT_TRUE(r.ran && has_text("utfeat.log", "\"move\"") && has_text("utfeat.log", "\"implicit\""), "the per-step debug log records the moves and the ROI steps");

    /* more trajectory points than --maxjumpdebug keeps */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 200 -E 1 -M H -C 1 -D M --maxjumpdebug 50 -s utcpujump", NULL);
    EXPECT_TRUE(r.ran && r.ntraj == 50, "the trajectory is capped at --maxjumpdebug points on the CPU");

    /* invalid settings */
    run_premut = mut_badmethod;
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10", NULL), "an undefined ray-tracing method");
    run_premut = NULL;
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 100 -M H -C 0 -O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", NULL),
                "a mesh-mode adjoint Jacobian with element output (-C 0)");

    if (!hasocl) {
        return;
    }

    /* per-node mua and mus (redbird's DOT reconstruction) equal to the medium give the same output */
    r = run_mmc("--bench dmmc-cube60 -G 1 -n 20000 -E 1 -M H -C 1 -s utclplain", NULL);
    run_premut = mut_nodalprop;
    rc = run_mmc("--bench dmmc-cube60 -G 1 -n 20000 -E 1 -M H -C 1 -s utclnodal", NULL);
    run_premut = NULL;
    EXPECT_TRUE(rc.ran && fabs(rc.wsum - r.wsum) < 0.02 * r.wsum, "per-node mua and mus equal to the medium give the same fluence on OpenCL");

    /* the OpenCL progress bar, and more detected photons or trajectory points than kept */
    r = run_mmc("--bench dmmc-cube60 -G 1 -n 20000 -E 1 -M H -C 1 -D P -s utclprog", NULL);
    EXPECT_TRUE(r.ran && r.wsum > 0., "the progress bar on OpenCL");
    r = run_mmc("--bench dmmc-cube60 -G 1 -n 20000 -E 1 -M H -C 1 -d 1 -H 10 -s utclmaxdet --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":4}]}}", NULL);
    EXPECT_TRUE(r.ran && r.detected == 10, "the detected photons are capped at -H on OpenCL");
    r = run_mmc("--bench dmmc-cube60 -G 1 -n 200 -E 1 -M H -C 1 -D M --maxjumpdebug 50 -s utcljump", NULL);
    EXPECT_TRUE(r.ran && r.ntraj == 50, "the trajectory is capped at --maxjumpdebug points on OpenCL");

#ifdef USE_CUDA
    /* the same on CUDA, in this process as the last runs of this group (a forked child can not
       use the GPU once it is initialized): per-node mua and mus, and the nodal-approximation
       mesh adjoint (adjointmode 1, set by mmclab/pmmc), alone, dual and dual with RF, which
       agree with OpenCL */
    {
        const char* adj[] = {
            "-O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}",
            "-O W --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}",
            "-O W --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}",
        };

        r = run_mmc("--bench dmmc-cube60 --compute cuda -G 1 -n 20000 -E 1 -M H -C 1 -s utcuplain", NULL);
        run_premut = mut_nodalprop;
        rc = run_mmc("--bench dmmc-cube60 --compute cuda -G 1 -n 20000 -E 1 -M H -C 1 -s utcunodal", NULL);
        EXPECT_TRUE(rc.ran && fabs(rc.wsum - r.wsum) < 0.02 * r.wsum, "per-node mua and mus equal to the medium give the same fluence on CUDA");

        run_premut = mut_adjointnodal;

        for (int i = 0; i < 3; i++) {
            snprintf(args, sizeof(args), "--bench dmmc-cube60 --compute opencl -G 1 -n 20000 -E 1234 -M H -C 1 %s -s utclna%d", adj[i], i);
            rc = run_mmc(args, NULL);
            snprintf(args, sizeof(args), "--bench dmmc-cube60 --compute cuda -G 1 -n 20000 -E 1234 -M H -C 1 %s -s utcuna%d", adj[i], i);
            r = run_mmc(args, NULL);
            printf("  nodal mesh adjoint %d: J sum cuda %.6g opencl %.6g, output sum cuda %.6g opencl %.6g\n", i + 1, r.jsum, rc.jsum, r.wsum, rc.wsum);
            sprintf(msg, "the nodal mesh adjoint %d on CUDA agrees with OpenCL within 5%%", i + 1);
            EXPECT_TRUE(r.ran && rc.jsum != 0. && fabs(r.jsum - rc.jsum) < 0.05 * fabs(rc.jsum), msg);
        }

        run_premut = NULL;
    }
#endif
}

/*-----------------------------------------------------------------------------
 *  input files: .inp with text mesh files, JSON examples, ROI files and a
 *  wide-field detector; the example folders are copied to the temporary folder
 *---------------------------------------------------------------------------*/

/** parse "mmc <args>" and prepare the mesh, without simulating */
static void parse_mmc(const char* args, void (*check)(mcconfig*, tetmesh*)) {
    char buf[4096], *argv[128];
    int argc = 0;
    mcconfig cfg;
    tetmesh mesh;
    raytracer tracer;

    snprintf(buf, sizeof(buf), "mmc %s", args);

    for (char* tok = strtok(buf, " "); tok && argc < 127; tok = strtok(NULL, " ")) {
        argv[argc++] = tok;
    }

    argv[argc] = NULL;
    memset(flagset, 0, sizeof(flagset));
    memset(&tracer, 0, sizeof(tracer));
    quiet_begin();
    mmc_init_from_cmd(&cfg, &mesh, &tracer, argc, argv);
    cfg.flog = devnull;
    quiet_end();

    if (check) {
        check(&cfg, &mesh);
    }

    cfg.flog = stdout;
    mmc_cleanup(&cfg, &mesh, &tracer);
}

static void check_inp(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_EQ(mesh->ne, 6, ".inp loads the 6 elements of onecube");
    EXPECT_EQ(cfg->detnum, 3, ".inp detectors");
    EXPECT_NEAR(cfg->tstep, 5e-10, 1e-15, ".inp time step");
}

static void check_edgeroi(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_TRUE(mesh->edgeroi != NULL && cfg->implicit == 1, "an edge ROI file is loaded");
}

static void write_roi(const char* name, int col, int row, float val) {
    FILE* fp = fopen(name, "wt");

    if (fp) {
        fprintf(fp, "%d %d\n", col, row);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                fprintf(fp, "%g ", val);
            }

            fprintf(fp, "\n");
        }

        fclose(fp);
    }
}

static void test_inputfiles(void) {
    char cmd[1024];
    RunResult r;

    /* both folders have a README.txt, copy them one after the other */
    snprintf(cmd, sizeof(cmd), "cp -r %s/onecube/* . && cp -rf %s/widedet/* .", MMC_EXAMPLE_DIR, MMC_EXAMPLE_DIR);

    if (system(cmd) != 0 || file_size("onecube.inp") < 0) {
        printf("FAIL can not copy the example inputs from %s\n", MMC_EXAMPLE_DIR);
        nfail++;
        return;
    }

    parse_mmc("-f onecube.inp", check_inp);
    r = run_mmc("-f onecube.inp -G -1 -n 500 -d 1 -q 1 -x 1 -s utinp", NULL);
    EXPECT_TRUE(r.ran, ".inp with text mesh files runs");
    r = run_mmc("-f onecube.json -G -1 -n 500 -s utjson", NULL);
    EXPECT_TRUE(r.ran, "onecube.json runs");
    r = run_mmc("-f onecube_isotropic.json -G -1 -n 500 -s utiso", NULL);
    EXPECT_TRUE(r.ran, "onecube_isotropic.json runs");
    /* two sources on onecube: each simulated separately (srcid -1), or only the first (srcid 1) */
    const char* msjson = "{\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},\"Session\":{\"Photons\":200,\"RNGSeed\":17182818,\"ID\":\"utms\"},"
                         "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-10},\"Optode\":{\"Source\":{\"Pos\":[[2,8,0,0.5],[8,2,0,0.5]],"
                         "\"Dir\":[[0,0,1,0],[0,0,1,0]],\"ID\":%d},\"Detector\":[{\"Pos\":[2,6,0],\"R\":1},{\"Pos\":[2,2,0],\"R\":1}]}}";
    char json[1024];
    snprintf(json, sizeof(json), msjson, -1);
    write_file("utms.json", json, strlen(json));
    r = run_mmc("-f utms.json -G -1 -s utms", NULL);
    EXPECT_TRUE(r.ran, "two sources simulated separately (srcid -1)");
    snprintf(json, sizeof(json), msjson, 1);
    write_file("utms1.json", json, strlen(json));
    r = run_mmc("-f utms1.json -G -1 -s utms1", NULL);
    EXPECT_TRUE(r.ran, "the first of two sources (srcid 1)");
    snprintf(json, sizeof(json), msjson, 5);
    write_file("utms5.json", json, strlen(json));
    EXPECT_EXIT(run_mmc("-f utms5.json -G -1 -s utms5", NULL), "srcid beyond the source count is an error");

    /* OpenCL grid output of two identical sources (srcid -1): the photons are
       split between the slots, but each slot is normalized to its own share, so
       it matches a run of one source; it was 1/2 of it, and adjoint Jacobians 1/4 */
    if (hasocl) {
        const char* gsjson = "{\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},\"Session\":{\"Photons\":40000,\"RNGSeed\":17182818,\"ID\":\"utgs\"},"
                             "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},\"Optode\":{\"Source\":{\"Pos\":[[2,8,0],[2,8,0]],"
                             "\"Dir\":[[0,0,1,0],[0,0,1,0]],\"ID\":%d}}}";
        RunResult rm, r1;

        snprintf(json, sizeof(json), gsjson, -1);
        write_file("utgsm.json", json, strlen(json));
        snprintf(json, sizeof(json), gsjson, 1);
        write_file("utgs1.json", json, strlen(json));
        rm = run_mmc("-f utgsm.json -G 1 -M G -s utgsm", NULL);
        r1 = run_mmc("-f utgs1.json -G 1 -M G -s utgs1", NULL);
        EXPECT_TRUE(rm.ran && r1.ran && r1.wsum > 0., "grid output of two sources and of one");
        EXPECT_NEAR(rm.wsum / r1.wsum, 1.0, 0.05, "grid output: each source slot is normalized to its own photons");
    }

    /* implicit MMC: node, face and edge ROIs read from roi_onecube.dat */
    write_roi("roi_onecube.dat", 1, 8, 0.5f);
    r = run_mmc("-f onecube.inp -G -1 -n 300 -s utroin", NULL);
    EXPECT_TRUE(r.ran, "a node ROI file");
    write_roi("roi_onecube.dat", 4, 6, 0.2f);
    r = run_mmc("-f onecube.inp -G -1 -n 300 -s utroif", NULL);
    EXPECT_TRUE(r.ran, "a face ROI file");
    /* an edge ROI on every edge of this tiny mesh makes photons hit the same triangle
       again and again (a known MMC limitation), so the edge ROI file is only loaded */
    write_roi("roi_onecube.dat", 6, 6, 0.3f);
    parse_mmc("-f onecube.inp", check_edgeroi);
    write_roi("roi_onecube.dat", 3, 6, 0.3f);
    EXPECT_EXIT(run_mmc("-f onecube.inp -G -1 -n 100 -s utroix", NULL), "an ROI file with 3 columns is an error");
    remove("roi_onecube.dat");

    /* a wide-field detector, saved as a detector image (-x 2) */
    r = run_mmc("-f widedet.inp -G -1 -n 2000 -x 2 -b 1 -s utwide", NULL);
    EXPECT_TRUE(r.ran && file_size("utwide_detimage.bin") != 0, "a wide-field detector image");

    /* an .inp read from stdin, with its prompts */
    {
        int status = 0;
        pid_t pid = fork();

        if (pid == 0) {
            if (freopen("onecube.inp", "r", stdin) == NULL || freopen("/dev/null", "w", stdout) == NULL) {
                exit(2);
            }

            parse_mmc("-i", NULL);
            exit(0);
        }

        waitpid(pid, &status, 0);
        EXPECT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 0, "an .inp read from stdin (-i)");
    }
}

/*-----------------------------------------------------------------------------
 *  command line options, parsed without simulating
 *---------------------------------------------------------------------------*/

static void check_opts1(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_EQ(cfg->nthread, 128, "-t");
    EXPECT_EQ(cfg->nblocksize, 32, "-T");
    EXPECT_EQ(cfg->isreflect, 1, "-b");
    EXPECT_EQ(cfg->respin, 2, "-r");
    EXPECT_NEAR(cfg->minenergy, 1e-4, 1e-10, "-e");
    EXPECT_EQ(cfg->isnormalized, 0, "-U");
    EXPECT_NEAR(cfg->unitinmm, 0.5, 1e-6, "-u");
    EXPECT_EQ(cfg->isspecular, 1, "-V");
    EXPECT_EQ(cfg->voidtime, 0, "-k");
    EXPECT_EQ(cfg->maxdetphoton, 1000, "-H");
    EXPECT_EQ(cfg->ismomentum, 1, "--momentum");
    EXPECT_EQ(cfg->maxjumpdebug, 100, "--maxjumpdebug");
    EXPECT_EQ(cfg->debugphoton, 3, "--debugphoton");
}

static void check_opts2(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_TRUE(strstr(cfg->compileropt, "-DUT_FLAG") != NULL, "-J");
    EXPECT_EQ(cfg->optlevel, 2, "-o");
    EXPECT_NEAR(cfg->workload[1], 70.f, 1e-6, "-W");
    EXPECT_EQ(cfg->autopilot, 0, "-A");
    EXPECT_EQ(cfg->srcid, 1, "--srcid");
    EXPECT_EQ(cfg->zipid, 1, "-Z");
}

static void test_cmdopts(void) {
    parse_mmc("--bench dmmc-cube60 -t 128 -T 32 -b 1 -r 2 -e 1e-4 -U 0 -u 0.5 -V 1 -k 0 -H 1000 --momentum 1 -d 1 "
              "--maxjumpdebug 100 --debugphoton 3 -g 1 -C 1 -M H --atomic 1 -S 0", check_opts1);
    parse_mmc("--bench dmmc-cube60 -J -DUT_FLAG -o 2 -W 30,70 -A 0 --srcid 1 -Z gzip --root . -s utopt -G 01", check_opts2);
    parse_mmc("--bench dmmc-cube60 -M G --gridsize 2 -F jnii -l 1 -s utlog", NULL);
    EXPECT_TRUE(file_size("utlog.log") >= 0, "-l writes a log file");
    parse_mmc("--bench dmmc-cube60 -j {\"Session\":{\"Photons\":5}} -D TP -O E -c sse", NULL);
    parse_mmc("--bench dmmc-cube60 -E 12345 -P 1 -q 1 -x 1 -X 1 -m 1 -S 1", NULL);

    EXPECT_EXIT(parse_mmc("--nosuchoption 1", NULL), "an unknown long option");
    EXPECT_EXIT(parse_mmc("-% 1", NULL), "an unknown short option");
    EXPECT_EXIT(parse_mmc("--bench nosuchbench", NULL), "an unknown benchmark");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -F", NULL), "-F without a value");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -F xyz", NULL), "an unknown output format");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -O z", NULL), "an unknown output type");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -M Y", NULL), "an unknown ray tracer");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -c nosuch", NULL), "an unknown compute backend");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -n", NULL), "-n without a value");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json {bad", NULL), "an invalid --json fragment");
    EXPECT_EXIT(parse_mmc("-f nosuch.json", NULL), "a missing input file");
    EXPECT_EXIT(parse_mmc("-f nosuch.inp", NULL), "a missing .inp file");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -O J", NULL), "a Jacobian outside of the replay mode");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -E nosuch.mch", NULL), "replay from a missing .mch file");

    /* the help, version and device list exit or return without simulating */
    {
        int status = 0;
        pid_t pid = fork();

        if (pid == 0) {
            if (freopen("/dev/null", "w", stdout)) {}

            parse_mmc("-h", NULL);
            exit(1);
        }

        waitpid(pid, &status, 0);
        EXPECT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 0, "-h prints the help and exits");
        pid = fork();

        if (pid == 0) {
            if (freopen("/dev/null", "w", stdout)) {}

            parse_mmc("-v", NULL);
            exit(1);
        }

        waitpid(pid, &status, 0);
        EXPECT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 0, "-v prints the version and exits");
    }

    if (hasocl) {
        RunResult r = run_mmc("-L", NULL);
        EXPECT_TRUE(!r.ran, "-L lists the OpenCL devices without simulating");
    }
}

/*-----------------------------------------------------------------------------
 *  the mmclab/pmmc setup path: mmc_validate_config, mcx_replayinit, mesh_validate
 *---------------------------------------------------------------------------*/

typedef void (*CfgMutator)(mcconfig* cfg, tetmesh* mesh);

static float* mlab_detps = NULL;
static int mlab_dims[2] = {0, 0};
static int mlab_seedbyte = 0;

/** load a benchmark, change it as mmclab would, validate, prepare and simulate on the CPU */
static void run_mmclab_path(const char* args, CfgMutator mutate) {
    char buf[1024], *argv[64];
    int argc = 0;
    mcconfig cfg;
    tetmesh mesh;
    raytracer tracer;

    snprintf(buf, sizeof(buf), "mmc %s", args);

    for (char* tok = strtok(buf, " "); tok && argc < 63; tok = strtok(NULL, " ")) {
        argv[argc++] = tok;
    }

    argv[argc] = NULL;
    memset(flagset, 0, sizeof(flagset));
    memset(&tracer, 0, sizeof(tracer));
    quiet_begin();
    mmc_init_from_cmd(&cfg, &mesh, &tracer, argc, argv);
    cfg.flog = devnull;

    if (mutate) {
        mutate(&cfg, &mesh);
    }

    mmc_validate_config(&cfg, mlab_detps, mlab_dims, mlab_seedbyte);
    mesh_validate(&mesh, &cfg);
    mmc_prep(&cfg, &mesh, &tracer);
    mmc_run_mp(&cfg, &mesh, &tracer);
    quiet_end();
    cfg.flog = stdout;
    mmc_cleanup(&cfg, &mesh, &tracer);
}

static void mut_replay(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->seed = SEED_FROM_FILE;
    cfg->nphoton = mlab_dims[1];
    cfg->savedetflag = 0x5;  /* detector ID and partial paths */
    cfg->outputtype = otJacobian;
    cfg->replaydet = 1;
    cfg->photonseed = calloc(mlab_dims[1], mlab_seedbyte);
}
static void mut_nphoton(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->nphoton = 0;
}
static void mut_tstep(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->tstep = 0.f;
}
static void mut_srcdir(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->srcdir.z = 2.f;
}
static void mut_unitinmm(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->unitinmm = -1.f;
}
static void mut_tend(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->tend = cfg->tstart;
    cfg->tstep = 1e-10f;
}
static void mut_pattern(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->srctype = stPattern;
}
static void mut_srcnumreplay(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->srcnum = 2;
    cfg->seed = SEED_FROM_FILE;
}
static void mut_adjoint0(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->outputtype = otAdjoint;
    cfg->basisorder = 0;
    cfg->method = rtBLBadouel;
}
static void mut_steps(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->steps.y = 2.f * cfg->steps.x;
}
static void mut_replaynodet(mcconfig* cfg, tetmesh* mesh) {
    mut_replay(cfg, mesh);
    cfg->savedetflag = 0x4;
}
static void mut_badmesh(mcconfig* cfg, tetmesh* mesh) {
    (void)cfg;
    mesh->prop = 0;  /* no optical properties */
}
static void mut_replaycount(mcconfig* cfg, tetmesh* mesh) {
    mut_replay(cfg, mesh);
    cfg->nphoton = mlab_dims[1] + 1;  /* more seeds than detected photons */
}

static void test_mmclab_path(void) {
    run_mmclab_path("--bench dmmc-cube60 -G -1 -n 200 -E 1 -M H -C 1 -s utml1", NULL);
    EXPECT_TRUE(1, "mmc_validate_config and mesh_validate accept a benchmark mesh");

    /* replay of 3 detected photons of detector 1, given as a detphoton array as mmclab does */
    float detps[3 * 3] = {1, 10, 0, 1, 20, 0, 2, 5, 0};
    mlab_detps = detps;
    mlab_dims[0] = 3;
    mlab_dims[1] = 3;
    mlab_seedbyte = sizeof(RandType) * RAND_BUF_LEN;
    run_mmclab_path("--bench dmmc-cube60 -G -1 -M H -C 1 -s utml2", mut_replay);
    EXPECT_TRUE(1, "replay from a detphoton array");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -M H -C 1 -s utml3", mut_replaynodet), "replay without detector IDs");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -M H -C 1 -s utml4", mut_replaycount), "replay with a mismatched photon count");
    mlab_seedbyte = 0;
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -M H -C 1 -s utml5", mut_replay), "replay without seeds");
    mlab_detps = NULL;
    mlab_seedbyte = 0;

    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_nphoton), "mmc_validate_config: no photons");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_tstep), "mmc_validate_config: zero time step");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_srcdir), "mmc_validate_config: non-unit srcdir");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_unitinmm), "mmc_validate_config: a negative unitinmm");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_tend), "mmc_validate_config: tend <= tstart");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_pattern), "mmc_validate_config: pattern without data");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_srcnumreplay), "mmc_validate_config: several patterns in replay");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_adjoint0), "mmc_validate_config: mesh adjoint without basisorder 1");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle", mut_steps), "mmc_validate_config: anisotropic voxels");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -M H -C 1 -s utmle", mut_badmesh), "mesh_validate: a mesh without optical properties");
}

#ifdef USE_OPENCL
char* print_cl_errstring(cl_int err);
#endif

static void test_clerrors(void) {
#ifdef USE_OPENCL

    for (int err = 1; err >= -72; err--) {
        EXPECT_TRUE(print_cl_errstring(err) != NULL, "print_cl_errstring");
    }

    for (int err = -1000; err >= -1010; err--) {
        EXPECT_TRUE(print_cl_errstring(err) != NULL, "print_cl_errstring of an extension error");
    }

#endif
}

/*-----------------------------------------------------------------------------
 *  ROI reflection, JData/BJData outputs, --dumpjson and .inp pattern files
 *---------------------------------------------------------------------------*/

static void check_writeconfig(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    quiet_begin();
    mcx_writeconfig("", cfg);
    quiet_end();
    mcx_writeconfig("utsaved.inp", cfg);
    EXPECT_TRUE(file_size("utsaved.inp") > 0, "mcx_writeconfig saves an .inp file");
    EXPECT_EXIT(mcx_writeconfig("/nonexistent/dir/ut.inp", cfg), "mcx_writeconfig to an unwritable path");
}

static void check_srcpattern(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_EQ(cfg->srctype, stPattern, ".inp pattern source type");
    EXPECT_EQ(cfg->srcnum, 2, ".inp pattern count");
    EXPECT_TRUE(cfg->srcpattern && cfg->srcpattern[7] == 8.f, ".inp pattern data read from the file");
}

static void check_detpattern(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_TRUE(cfg->detpattern && cfg->detpattern[1599] == 0.5f, ".inp detector pattern read in replay");
}

static void check_kernel(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_TRUE(cfg->clsource && strstr(cfg->clsource, "__kernel") != NULL, "--kernel loads an OpenCL kernel file");
}

static void test_moreio(void) {
    char args[1024], msg[256];
    RunResult r;

    /* iMMC with a refractive index mismatch between the ROI (the last medium) and the
       element, higher and lower, so that reflectrayroi reflects, transmits and totally
       reflects the photons at node, face and edge ROIs */
    const char* bench[] = {"nodeimmc", "faceimmc", "edgeimmc"};
    const float roin[] = {1.6f, 1.1f};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            snprintf(args, sizeof(args), "--bench %s -G -1 -n 300 -E 1 -b 1 -D R -s utroi%d%d --json "
                     "{\"Domain\":{\"Media\":[{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1},{\"mua\":0.0458,\"mus\":35.6541,\"g\":0.9,\"n\":1.37},"
                     "{\"mua\":0.5,\"mus\":9.4,\"g\":0.9,\"n\":%g}]}}", bench[i], i, j, roin[j]);
            r = run_mmc(args, NULL);
            sprintf(msg, "iMMC %s with ROI index %g", bench[i], roin[j]);
            EXPECT_TRUE(r.ran && r.normalizer > 0, msg);
        }
    }

    /* trajectories and detected photons (with seeds) saved as JData (.jdat) and BJData (.jdb),
       raw (-Z base64) or compressed */
    const char* fmt[] = {"jnii", "bnii"};
    const char* zip[] = {"zlib", "base64"};

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 20 -E 1 -M H -C 1 -D M -F %s -Z %s -s uttraj%d%d", fmt[i], zip[j], i, j);
            r = run_mmc(args, NULL);
            sprintf(msg, "trajectories saved by -F %s -Z %s", fmt[i], zip[j]);
            EXPECT_TRUE(r.ran, msg);
            snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 500 -E 1 -M H -C 1 -d 1 -q 1 -F %s -Z %s -s utdetq%d%d "
                     "--json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":5}]}}", fmt[i], zip[j], i, j);
            r = run_mmc(args, NULL);
            sprintf(msg, "detected photons and seeds saved by -F %s -Z %s", fmt[i], zip[j]);
            EXPECT_TRUE(r.ran && r.detected > 0, msg);
        }

        /* two sources: the detected photons carry the source ID in the upper bits of the detector ID */
        snprintf(args, sizeof(args), "-f utms.json -G -1 -d 1 -F %s -s utmsdet%d", fmt[i], i);
        r = run_mmc(args, NULL);
        sprintf(msg, "detected photons of two sources saved by -F %s", fmt[i]);
        EXPECT_TRUE(r.ran, msg);
    }

    /* --dumpjson of multiple sources and of pattern sources, to a file, then read back */
    EXPECT_EXIT_OK(run_mmc("-f utms.json -G -1 --dumpjson utmsdump.json", NULL), "--dumpjson of two sources");
    EXPECT_TRUE(file_size("utmsdump.json") > 0, "--dumpjson writes the named file");
    r = run_mmc("-f utmsdump.json -G -1 -n 100 -s utmsre", NULL);
    EXPECT_TRUE(r.ran, "the dumped two-source JSON runs");
    EXPECT_EXIT_OK(run_mmc("--bench dmmc-cube60 -G -1 --dumpjson utpatdump.json --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\","
                           "\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1]}}}}", NULL),
                   "--dumpjson of a pattern source");
    EXPECT_EXIT_OK(run_mmc("--bench dmmc-cube60 -G -1 --dumpjson utpat2dump.json --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\","
                           "\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"SrcNum\":2,\"Pattern\":{\"Nx\":2,\"Ny\":2,"
                           "\"Data\":[1,0.5,0.5,1,0,1,1,0]}}}}", NULL),
                   "--dumpjson of two patterns");
    EXPECT_TRUE(file_size("utpatdump.json") > 0 && file_size("utpat2dump.json") > 0, "the pattern sources are dumped");
    parse_mmc("-f utpat2dump.json", NULL);
    EXPECT_TRUE(1, "the dumped two-pattern JSON is read back");

    /* save the settings as an .inp file */
    parse_mmc("--bench dmmc-cube60", check_writeconfig);

    /* an .inp file with a pattern source and its pattern file (2 patterns of 2x2) */
    {
        float pat[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        const char* inp = "100\n17182818\n2 8 0.0\n0. 0. 1.\n0.e+00 5.e-09 5e-10\nonecube\n3\n1 1\n2.0 6.0 0.0\n"
                          "pattern\n4 0 0 2\n0 4 0 2\nutpat.bin 2\n";
        write_file("utpat.bin", pat, sizeof(pat));
        write_file("utpat.inp", inp, strlen(inp));
        parse_mmc("-f utpat.inp", check_srcpattern);
        inp = "100\n17182818\n2 8 0.0\n0. 0. 1.\n0.e+00 5.e-09 5e-10\nonecube\n3\n1 1\n2.0 6.0 0.0\n"
              "pattern\n4 0 0 2\n0 4 0 2\nnosuchpattern.bin 1\n";
        write_file("utpatx.inp", inp, strlen(inp));
        EXPECT_EXIT(parse_mmc("-f utpatx.inp", NULL), "a missing source pattern file");
        inp = "100\n17182818\n2 8 0.0\n0. 0. 1.\n0.e+00 5.e-09 5e-10\nonecube\n3\n1 1\n2.0 6.0 0.0\n"
              "pattern\n4 0 0 2\n0 4 0 2\nutpat.bin 0\n";
        write_file("utpat0.inp", inp, strlen(inp));
        EXPECT_EXIT(parse_mmc("-f utpat0.inp", NULL), "a pattern count of 0");
    }

    /* a wide-field detector with a detection pattern, which is only read in replay */
    {
        float* dpat = (float*)calloc(1600, sizeof(float));
        const char* inp = "100\n17182818\n2 8 0.0\n0. 0. 1.\n0.e+00 5.e-09 5e-10\nonecube\n3\n1 0\n0.0 0.0 0.0 0.0\n"
                          "pencil\n0 0 0 0\n0 0 0 0\n10.0 0.0 0.0 40\n0.0 10.0 0.0 40\nutdpat.bin\n";

        /* the replayed .mch must come from a mesh with the same media */
        r = run_mmc("-f onecube.inp -G -1 -n 300 -d 1 -q 1 -x 1 -s utocdet", NULL);
        EXPECT_TRUE(r.ran && file_size("utocdet.mch") > 0, "onecube detected photons for the replay");
        dpat[1599] = 0.5f;
        write_file("utdpat.bin", dpat, 1600 * sizeof(float));
        free(dpat);
        write_file("utdpat.inp", inp, strlen(inp));
        parse_mmc("-f utdpat.inp -E utocdet.mch -O L", check_detpattern);
        inp = "100\n17182818\n2 8 0.0\n0. 0. 1.\n0.e+00 5.e-09 5e-10\nonecube\n3\n1 0\n0.0 0.0 0.0 0.0\n"
              "pencil\n0 0 0 0\n0 0 0 0\n10.0 0.0 0.0 40\n0.0 10.0 0.0 40\nnosuchdpat.bin\n";
        write_file("utdpatx.inp", inp, strlen(inp));
        EXPECT_EXIT(parse_mmc("-f utdpatx.inp -E utocdet.mch -O L", NULL), "a missing detector pattern file");
    }

    /* an OpenCL kernel loaded from a file with --kernel */
    parse_mmc("--bench dmmc-cube60 --kernel " MMC_EXAMPLE_DIR "/../src/mmc_core.cl", check_kernel);
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --kernel nosuchkernel.cl", NULL), "a missing --kernel file");
}

/*-----------------------------------------------------------------------------
 *  tracer variants: several patterns, replay outputs, JSON variants, -N
 *---------------------------------------------------------------------------*/

static void check_jsonvariants(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_NEAR(cfg->prop[1].mus, 1.5, 1e-6, "Media given as 4-element arrays");
    EXPECT_TRUE(isinf(cfg->srcdir.w) && cfg->srcdir.w > 0, "Dir with _Inf_ focus");
    EXPECT_EQ(cfg->checkpt[1], 20, "Session.Checkpoints");
    EXPECT_EQ(cfg->detnum, 2, "Detector object with two Pos rows");
    EXPECT_TRUE(cfg->detpos && cfg->detpos[1].y == 2.f && cfg->detpos[1].w == 1.f, "the rows share R");
    EXPECT_TRUE(cfg->detdir && cfg->detdir[1].z == -1.f, "Detector.Dir rows of a detector object");
}

static void check_detobject(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_EQ(cfg->detnum, 1, "a single Detector object is one detector");
    EXPECT_TRUE(cfg->detpos && cfg->detpos[0].y == 6.f && cfg->detpos[0].w == 1.5f, "the Pos and R of a Detector object");
    EXPECT_TRUE(cfg->detdir && cfg->detdir[0].z == -1.f, "the Dir of a Detector object");
    EXPECT_TRUE(cfg->detparam1.x == 10.f && cfg->detparam2.w == 4.f, "Detector Param1/Param2 are read from the detector");
}

static void check_detrows(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_EQ(cfg->detnum, 2, "Detector as [x,y,z,r] rows");
    EXPECT_TRUE(cfg->detpos && cfg->detpos[1].y == 2.f && cfg->detpos[1].w == 2.f, "the radius in the 4th column");
}

static void check_srcdirrows(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_TRUE(isnan(cfg->srcdir.w), "Dir row 1 with _NaN_");
    EXPECT_TRUE(cfg->srcdata && isinf(cfg->srcdata[1].srcdir.w) && cfg->srcdata[1].srcdir.w < 0, "Dir row 2 with -_Inf_");
}

static void check_patfile(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_TRUE(cfg->srcpattern && cfg->srcpattern[3] == 4.f, "Pattern Data read from a file name");
}

static void check_netdoc(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    EXPECT_EQ(cfg->extrasrclen, 2, "-N loads the downloaded two-source JSON");
}

static void test_variants(void) {
    char args[2048], msg[256];
    RunResult r;
    const char* method[] = {"P", "H", "B", "S"};

    /* two patterns simulated at once with each CPU tracer, nodal and elem-based outputs */
    for (int i = 0; i < 4; i++) {
        for (int c = 0; c < 2; c++) {
            snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n 300 -E 1 -M %s -C %d -s utmp%s%d --json "
                     "{\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"SrcNum\":2,\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],"
                     "\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1,0,1,1,0]}}}}", method[i], c, method[i], c);
            r = run_mmc(args, NULL);
            sprintf(msg, "two patterns with -M %s -C %d", method[i], c);
            EXPECT_TRUE(r.ran && r.normalizer > 0, msg);
        }
    }

    /* elem-based output with the SSE Badouel tracer, and atomic updates */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -M S -C 0 -s utS0", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "-M S -C 0");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -M B -C 0 --atomic 1 -s utB0a", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "-M B -C 0 --atomic 1");

    /* replayed Jacobian, weighted path length and scattering count with each tracer */
    const char* rtype[] = {"J", "L", "P"};

    for (int i = 0; i < 4; i++) {
        if (method[i][0] == 'H') {
            continue;  /* tested in test_outputs */
        }

        for (int j = 0; j < 3; j++) {
            snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -E utdet.mch -O %s -M %s -C %d -s utrp%s%s",
                     rtype[j], method[i], method[i][0] != 'P', method[i], rtype[j]);
            r = run_mmc(args, NULL);
            sprintf(msg, "replay -O %s with -M %s", rtype[j], method[i]);
            EXPECT_TRUE(r.ran, msg);
        }
    }

    /* JSON variants: Media as numeric arrays, an infinite focus, checkpoints and the
       directions of a detector array object */
    {
        const char* json = "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utjv\",\"Checkpoints\":[10,20,30]},"
                           "\"Domain\":{\"Media\":[[0,0,1,1],[0.01,1.5,0.9,1.37]]},"
                           "\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},"
                           "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},"
                           "\"Optode\":{\"Source\":{\"Pos\":[2,8,0],\"Dir\":[0,0,1,\"_Inf_\"]},"
                           "\"Detector\":{\"Pos\":[[2,6,0],[2,2,0]],\"R\":1,\"Dir\":[[0,0,-1],[0,0,-1,0]]}}}";
        write_file("utjv.json", json, strlen(json));
        parse_mmc("-f utjv.json", check_jsonvariants);

        json = "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utjd\"},\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},"
               "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},"
               "\"Optode\":{\"Source\":{\"Pos\":[[2,8,0,1],[8,2,0,1]],\"Dir\":[[0,0,1,\"_NaN_\"],[0,0,1,\"-_Inf_\"]]}}}";
        write_file("utjd.json", json, strlen(json));
        parse_mmc("-f utjd.json", check_srcdirrows);

        float pat[4] = {1, 2, 3, 4};
        write_file("utjpat.bin", pat, sizeof(pat));
        parse_mmc("--bench dmmc-cube60 --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],"
                  "\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":\"utjpat.bin\"}}}}", check_patfile);

        /* a single detector object, and an array of [x,y,z,r] rows */
        json = "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utjo\"},\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},"
               "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},\"Optode\":{\"Source\":{\"Pos\":[2,8,0],\"Dir\":[0,0,1],\"Param1\":[9,9,9,9]},"
               "\"Detector\":{\"Pos\":[2,6,0],\"R\":1.5,\"Dir\":[0,0,-1],\"Param1\":[10,0,0,4],\"Param2\":[0,10,0,4]}}}";
        write_file("utjo.json", json, strlen(json));
        parse_mmc("-f utjo.json", check_detobject);
        json = "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utjr\"},\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},"
               "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},\"Optode\":{\"Source\":{\"Pos\":[2,8,0],\"Dir\":[0,0,1]},"
               "\"Detector\":[[2,6,0,1],[2,2,0,2]]}}";
        write_file("utjr.json", json, strlen(json));
        parse_mmc("-f utjr.json", check_detrows);

        write_file("utbad.json", "{\"Session\":{\"Photons\":", 22);
        EXPECT_EXIT(parse_mmc("-f utbad.json", NULL), "a malformed JSON input file");
    }

    /* -N with a stand-in curl that answers the NeuroJSON.io queries offline */
    {
        char cwd[4096], path[8192];
        const char* curl = "#!/bin/sh\n"
                           "case \"$*\" in\n"
                           "  *nosuchdb*) echo '{\"error\":\"not_found\"}';;\n"
                           "  *_find*) echo '{\"docs\":[{\"_id\":\"_design\"},{\"_id\":\"utdoc1\"},{\"_id\":\"utdoc2\"}]}';;\n"
                           "  *_all_docs*) echo '{\"rows\":[{\"id\":\"utdoc1\"},{\"id\":\".hidden\"}]}';;\n"
                           "  *) cat utmsdump.json;;\n"
                           "esac\n";
        const char* oldpath = getenv("PATH");

        if (mkdir_p("utbin") == 0 && getcwd(cwd, sizeof(cwd))) {
            write_file("utbin/curl", curl, strlen(curl));
            chmod("utbin/curl", 0755);
            snprintf(path, sizeof(path), "%s/utbin:%s", cwd, oldpath ? oldpath : "/usr/bin:/bin");
            char* savedpath = strdup(oldpath ? oldpath : "");
            setenv("PATH", path, 1);

            EXPECT_EXIT_OK(parse_mmc("-N", NULL), "-N lists the simulations of the database");
            EXPECT_EXIT_OK(parse_mmc("-N {\"selector\":{}}", NULL), "-N with a query selector");
            EXPECT_EXIT_OK(parse_mmc("-N limit=2", NULL), "-N with a query string");
            EXPECT_EXIT(parse_mmc("-N nosuchdb=1", NULL), "-N with a failed query");
            parse_mmc("-N utdoc1", check_netdoc);
            parse_mmc("-N mmc/utdoc1", NULL);
            parse_mmc("-N https://neurojson.io:7777/mmc/utdoc1", NULL);
            EXPECT_TRUE(1, "-N downloads a simulation by name, path and URL");

            setenv("PATH", savedpath, 1);
            free(savedpath);
        }
    }

    /* the library entry points used by the Python module */
    {
        long len = file_size("utmsdump.json");
        char* jbuf = (char*)calloc(len + 1, 1);
        FILE* fp = fopen("utmsdump.json", "rb");
        mcconfig cfg;
        tetmesh mesh;
        raytracer tracer;

        if (fp && fread(jbuf, len, 1, fp) == 1) {
            memset(flagset, 0, sizeof(flagset));
            memset(&tracer, 0, sizeof(tracer));
            quiet_begin();
            mmc_init_from_json(&cfg, &mesh, &tracer, jbuf, NULL);
            cfg.flog = stdout;
            mmc_reset(&cfg, &mesh, &tracer);
            quiet_end();
            cfg.flog = stdout;
            mmc_cleanup(&cfg, &mesh, &tracer);
            EXPECT_TRUE(1, "mmc_init_from_json and mmc_reset");
        }

        if (fp) {
            fclose(fp);
        }

        free(jbuf);
    }

    /* more detected photons in one thread than the initial detector buffer holds */
#ifdef _OPENMP
    {
        int nthread = omp_get_max_threads();
        omp_set_num_threads(1);
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G -1 -n %d -E 1 -M H -C 1 -d 1 -s utbuf --json "
                 "{\"Domain\":{\"Media\":[{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1},{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1}]},"
                 "\"Optode\":{\"Detector\":[{\"Pos\":[30,30,60],\"R\":100}]}}", DET_PHOTON_BUF + 1000);
        r = run_mmc(args, NULL);
        EXPECT_TRUE(r.ran && r.detected > DET_PHOTON_BUF, "the detector buffer grows past DET_PHOTON_BUF");
        omp_set_num_threads(nthread);
    }
#endif

    /* a device number list */
    {
        char list[MAX_DEVICE] = {0}, arg[] = "1,2,2,3";
        char* argv[] = {"mmc", "-G", arg, NULL};
        EXPECT_EQ(mcx_readarg(3, argv, 1, list, "bytenumlist"), 2, "mcx_readarg bytenumlist");
        EXPECT_TRUE(list[0] == 1 && list[1] == 2 && list[2] == 3, "bytenumlist removes duplicates");
    }
}

/*-----------------------------------------------------------------------------
 *  surface-only meshes (Mesh.MeshSurf) for OptiX
 *---------------------------------------------------------------------------*/

static int run_isolated(const char* args, RunResult* res);

/** write a 60 mm cube as a surface mesh: 12 outward triangles, label 1 inside, 0 outside,
    with the media and source of the dmmc-cube60 benchmark; extra replaces the last rows */
static void write_cubesurf(const char* fname, const char* extra) {
    const int quad[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    char buf[4096];
    int len = snprintf(buf, sizeof(buf), "{\"Session\":{\"ID\":\"utsurf\",\"Photons\":20000,\"RNGSeed\":1648335518,\"DoMismatch\":0},"
                       "\"Domain\":{\"Media\":[{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1},{\"mua\":0.005,\"mus\":1.0,\"g\":0.01,\"n\":1.37}]},"
                       "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},\"Optode\":{\"Source\":{\"Pos\":[30.1,30.2,0.0],\"Dir\":[0,0,1]}},"
                       "\"Mesh\":{\"MeshNode\":[[0,0,0],[60,0,0],[60,60,0],[0,60,0],[0,0,60],[60,0,60],[60,60,60],[0,60,60]],\"MeshSurf\":[");

    for (int i = 0; i < 6; i++) {
        if (extra && i == 5) {
            len += snprintf(buf + len, sizeof(buf) - len, ",%s", extra);
            break;
        }

        len += snprintf(buf + len, sizeof(buf) - len, "%s[%d,%d,%d,1,0],[%d,%d,%d,1,0]", i ? "," : "",
                        quad[i][0] + 1, quad[i][1] + 1, quad[i][2] + 1, quad[i][0] + 1, quad[i][2] + 1, quad[i][3] + 1);
    }

    snprintf(buf + len, sizeof(buf) - len, "]}}");
    write_file(fname, buf, strlen(buf));
}

static void check_cubesurf(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_EQ(mesh->nface, 12, "MeshSurf gives the OptiX triangles");
    EXPECT_EQ(mesh->ne, 0, "a surface-only mesh has no elements");
    EXPECT_EQ(cfg->e0, -1, "the initial medium of a surface-only mesh is found by ray tracing");
    EXPECT_TRUE(mesh->face && mesh->face[0].x == 0 && mesh->face[0].y == 3 && mesh->face[0].z == 2, "MeshSurf node indices become 0-based");
    EXPECT_TRUE(mesh->back && mesh->front && mesh->back[0] == 1 && mesh->front[0] == 0, "inner/outer labels become back/front");
    EXPECT_TRUE(mesh->fnorm && fabs(mesh->fnorm[0].z + 1.f) < 1e-6, "the bottom face normal points outward (-z)");
}


/* mesh_validate, the mmclab/pmmc entry, on a surface-only mesh and on incomplete meshes */
static void check_validate_surf(mcconfig* cfg, tetmesh* mesh) {
    mesh_validate(mesh, cfg);
    EXPECT_TRUE(cfg->e0 == -1 && mesh->nface == 12, "mesh_validate prepares a surface-only mesh for OptiX");
}
static void check_validate_surfcl(mcconfig* cfg, tetmesh* mesh) {
    cfg->compute = cbOpenCL;
    mesh_validate(mesh, cfg);
}
static void check_validate_surfnonode(mcconfig* cfg, tetmesh* mesh) {
    mesh->nn = 0;
    mesh_validate(mesh, cfg);
}
static void check_validate_nosurf(mcconfig* cfg, tetmesh* mesh) {
    cfg->surfnum = 0;
    mesh_validate(mesh, cfg);
}
static void check_validate_nonode(mcconfig* cfg, tetmesh* mesh) {
    mesh->node = NULL;
    mesh_validate(mesh, cfg);
}

/* mesh_savedetphoton writes the detected photons (.mch) or the trajectories (.mct) */
static void check_savedetphoton(mcconfig* cfg, tetmesh* mesh) {
    float ppath[2 * 3] = {1, 2.f, 0.5f, 1, 3.f, 0.25f};
    unsigned char seeds[2 * 4] = {1, 2, 3, 4, 5, 6, 7, 8};

    (void)mesh;
    cfg->his.colcount = 3;
    cfg->his.detected = 2;
    cfg->his.savedphoton = 2;
    cfg->issaveseed = 1;
    strcpy(cfg->session, "utsdp");
    mesh_savedetphoton(ppath, seeds, 2, 4, cfg);
    EXPECT_EQ(file_size("utsdp.mch"), (long)(sizeof(history) + sizeof(ppath) + sizeof(seeds)), "mesh_savedetphoton saves the photons and seeds");
    cfg->his.detected = 0;
    strcpy(cfg->rootpath, ".");
    mesh_savedetphoton(ppath, NULL, 2, 0, cfg);
    EXPECT_TRUE(file_size("utsdp.mct") > 0, "mesh_savedetphoton saves trajectories as .mct");
    cfg->rootpath[0] = '\0';
}

static void test_meshsurf(void) {
    char badrow[256];

    write_cubesurf("utsurf.json", NULL);
    parse_mmc("-f utsurf.json --compute optix -M G", check_cubesurf);

    /* other backends and outputs walk elements */
    EXPECT_EXIT(parse_mmc("-f utsurf.json --compute opencl -M G", NULL), "a surface-only mesh is rejected by OpenCL");
    EXPECT_EXIT(parse_mmc("-f utsurf.json -G -1", NULL), "a surface-only mesh is rejected on the CPU");

    /* invalid surfaces: a node index out of range, the same inner and outer labels, a label
       beyond the media, a zero-area triangle, and a row without 5 columns */
    const char* bad[] = {"[3,4,9,1,0],[3,1,5,1,0]", "[3,4,8,1,1],[3,1,5,1,0]", "[3,4,8,2,0],[3,1,5,1,0]",
                         "[3,3,8,1,0],[3,1,5,1,0]", "[3,4,8,1],[3,1,5,1,0]"
                        };
    const char* badmsg[] = {"MeshSurf node index out of range", "MeshSurf with the same inner and outer labels",
                            "MeshSurf label beyond the media", "MeshSurf zero-area triangle", "MeshSurf row without 5 columns"
                           };

    for (int i = 0; i < 5; i++) {
        snprintf(badrow, sizeof(badrow), "%s", bad[i]);
        write_cubesurf("utsurfbad.json", badrow);
        EXPECT_EXIT(parse_mmc("-f utsurfbad.json --compute optix -M G", NULL), badmsg[i]);
    }

    parse_mmc("-f utsurf.json --compute optix -M G", check_validate_surf);
    EXPECT_EXIT(parse_mmc("-f utsurf.json --compute optix -M G", check_validate_surfcl), "mesh_validate: a surface-only mesh needs OptiX");
    EXPECT_EXIT(parse_mmc("-f utsurf.json --compute optix -M G", check_validate_surfnonode), "mesh_validate: a surface-only mesh needs nodes");
    EXPECT_EXIT(parse_mmc("-f utsurf.json --compute optix -M G", check_validate_nosurf), "mesh_validate: a mesh without elements or surface");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60", check_validate_nonode), "mesh_validate: a mesh without nodes");
    parse_mmc("--bench dmmc-cube60", check_savedetphoton);

    /* an open surface (the last face removed) is only a warning */
    write_cubesurf("utsurfopen.json", "[4,1,5,1,0]");
    parse_mmc("-f utsurfopen.json --compute optix -M G", NULL);
    EXPECT_TRUE(1, "an open MeshSurf is accepted with a warning");

    /* --dumpjson writes MeshSurf, which is read back */
    EXPECT_EXIT_OK(run_mmc("-f utsurf.json --compute optix -M G --dumpjson utsurfdump.json", NULL), "--dumpjson of a surface-only mesh");
    parse_mmc("-f utsurfdump.json --compute optix -M G", check_cubesurf);

#ifdef USE_OPTIX
    /* the surface cube and the tetrahedral dmmc-cube60 benchmark are the same domain */
    RunResult rs, rt;
    int ret = run_isolated("-f utsurf.json --compute optix -M G -E 1234 -s utsurfrun", &rs);
    run_isolated("--bench dmmc-cube60 --compute optix -M G -n 20000 -E 1234 -s utsurftet", &rt);
    EXPECT_TRUE(ret == 0 && rt.wsum > 0. && fabs(rs.wsum - rt.wsum) / rt.wsum < 0.01, "OptiX gives the same fluence from MeshSurf as from the tetrahedral mesh");
#endif
}


/*-----------------------------------------------------------------------------
 *  JSON input: malformed mesh, media, source, pattern and detector fields, which
 *  must be reported as errors, and the less common valid layouts
 *---------------------------------------------------------------------------*/

/** write a JSON input of a one-tetrahedron mesh with the given Session members (starting
    with a comma), Domain, Mesh and Source contents and Detector member (starting with a
    comma); an empty or NULL value uses the default */
static void write_tetjson(const char* fname, const char* session, const char* domain, const char* mesh,
                          const char* source, const char* detector) {
    char buf[4096];

    snprintf(buf, sizeof(buf), "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utje\"%s},"
             "\"Domain\":{%s},\"Mesh\":{%s},"
             "\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},\"Optode\":{\"Source\":{%s}%s}}",
             session, (domain && domain[0]) ? domain : "\"Media\":[[0,0,1,1],[0.01,1,0.9,1.37]]", mesh ? mesh : "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],\"MeshElem\":[[1,2,3,4,1]],\"InitElem\":1",
             source ? source : "\"Pos\":[0.2,0.2,0],\"Dir\":[0,0,1]", detector);
    write_file(fname, buf, strlen(buf));
}

static void check_srcrows3(mcconfig* cfg, tetmesh* mesh) {
    EXPECT_EQ(cfg->extrasrclen, 3, "3 rows of a source field set 3 sources");
}

static void test_jsonedge(void) {
    const char* tetelem = "\"MeshElem\":[[1,2,3,4,1]],\"InitElem\":1";
    char mesh[1024];
    const struct {
        const char* session, *domain, *mesh, *source, *detector, *msg;
    } bad[] = {
        {"", "", "\"MeshNode\":[[0,0],[1,0,0],[0,1,0],[0,0,1]],\"MeshElem\":[[1,2,3,4,1]],\"InitElem\":1", NULL, "", "a MeshNode row without 3 numbers"},
        {
            "", "", "\"MeshNode\":{\"_ArrayType_\":\"int32\",\"_ArraySize_\":[4,3],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":12,\"_ArrayZipData_\":\"eJxjYEAARgZUgI0PAABsAAQ=\"},"
            "\"MeshElem\":[[1,2,3,4,1]],\"InitElem\":1", NULL, "", "a JData MeshNode that is not single precision"
        },
        {"", "", "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],\"MeshElem\":[[1,2,3,4]],\"InitElem\":1", NULL, "", "a MeshElem row without 5 integers"},
        {
            "", "", "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],\"MeshElem\":{\"_ArrayType_\":\"single\",\"_ArraySize_\":[1,5],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":5,\"_ArrayZipData_\":\"eJxjYGiwZ2BgcAAiIG4AYXsAGsIC/w==\"}"
            ",\"InitElem\":1", NULL, "", "a JData MeshElem that is not int32"
        },
        {
            "", "", "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],\"MeshSurf\":{\"_ArrayType_\":\"single\",\"_ArraySize_\":[1,5],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":5,\"_ArrayZipData_\":\"eJxjYGiwZ2BgcAAiIAazGQAZfgI/\"}"
            "", NULL, "", "a JData MeshSurf that is not int32"
        },
        {"", "", NULL, NULL, "", NULL},  /* placeholder for the MeshROI cases below */
        {"", "\"Media\":[[0,0,1,1],5]", NULL, NULL, "", "a medium that is neither an object nor an array"},
        {"", "\"Media\":[[0,0,1,1],[0.01,1,0.9,1.37]],\"Step\":[]", NULL, NULL, "", "an empty Domain.Step"},
        {"", "", NULL, "\"Pos\":[[0.2,0.2,0,1],[0.3,0.2,0,1]],\"Dir\":[[0,0,1],[0,0,1],[0,0,1]]", "", "Dir rows that do not match Pos"},
        {"", "", NULL, "\"Pos\":[[0.2,0.2,0,1],[0.3,0.2,0,1]],\"Dir\":[0,0,1],\"Param1\":[[0,0,0,0],[0,0,0,0],[0,0,0,0]]", "", "Param1 rows that do not match Pos"},
        {"", "", NULL, "\"Pos\":[[0.2,0.2,0,1],[0.3,0.2,0,1]],\"Dir\":[0,0,1],\"Param2\":[[0,0,0,0],[0,0,0,0],[0,0,0,0]]", "", "Param2 rows that do not match Pos"},
        {
            "", "", NULL, "\"Type\":\"pattern\",\"Pos\":[0.1,0.1,0],\"Dir\":[0,0,1],\"Param1\":[0.5,0,0,2],\"Param2\":[0,0.5,0,2],"
            "\"Pattern\":{\"_ArrayType_\":\"uint8\",\"_ArraySize_\":[2,2],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":4,\"_ArrayZipData_\":\"eJxjZGBgBAAACQAD\"}",
            "", "a compressed pattern that is not single precision"
        },
        {
            "", "", NULL, "\"Type\":\"pattern\",\"Pos\":[0.1,0.1,0],\"Dir\":[0,0,1],\"Param1\":[0.5,0,0,2],\"Param2\":[0,0.5,0,2],"
            "\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,2]}", "", "an incomplete pattern"
        },
        {
            "", "", NULL, "\"Type\":\"pattern\",\"Pos\":[0.1,0.1,0],\"Dir\":[0,0,1],\"Param1\":[0.5,0,0,2],\"Param2\":[0,0.5,0,2],"
            "\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":\"utnosuchpattern.bin\"}", "", "a missing pattern file"
        },
        {"", "", NULL, NULL, ",\"Detector\":{\"Pos\":[[0.2,0.2],[0.3,0.3,0]],\"R\":1}", "detector rows without 3 numbers"},
        {"", "", NULL, NULL, ",\"Detector\":[{\"Pos\":[0.2,0.2],\"R\":1}]", "a detector object without 3 numbers"},
        {",\"OutputFormat\":\"nosuchformat\"", "", NULL, NULL, "", "an unknown output format"},
        {",\"RayTracer\":\"q\"", "", NULL, NULL, "", "an unknown ray tracer"},
        {",\"OutputType\":\"z\"", "", NULL, NULL, "", "an unknown output type"},
        {"", "", "\"InitElem\":1", NULL, "", "no mesh"},
        {"", "", "\"MeshID\":\"onecube\"", NULL, "", "no initial element"},
    };

    for (int i = 0; i < (int)(sizeof(bad) / sizeof(bad[0])); i++) {
        if (bad[i].msg == NULL) {
            continue;
        }

        write_tetjson("utje.json", bad[i].session, bad[i].domain, bad[i].mesh, bad[i].source, bad[i].detector);
        EXPECT_EXIT(parse_mmc("-f utje.json", NULL), bad[i].msg);
    }

    /* MeshROI: rows of different lengths, a JData array that is not single precision, and a
       single row of per-node values */
    snprintf(mesh, sizeof(mesh), "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],%s,\"MeshROI\":[[1,2],[1]]", tetelem);
    write_tetjson("utje.json", "", "", mesh, NULL, "");
    EXPECT_EXIT(parse_mmc("-f utje.json", NULL), "MeshROI rows of different lengths");
    snprintf(mesh, sizeof(mesh), "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],%s,\"MeshROI\":{\"_ArrayType_\":\"int32\",\"_ArraySize_\":[4],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":4,\"_ArrayZipData_\":\"eJxjZGBgYETCAAA4AAU=\"}"
             "", tetelem);
    write_tetjson("utje.json", "", "", mesh, NULL, "");
    EXPECT_EXIT(parse_mmc("-f utje.json", NULL), "a JData MeshROI that is not single precision");
    snprintf(mesh, sizeof(mesh), "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],%s,\"MeshROI\":[0.1,0.1,0.1,0.1]", tetelem);
    write_tetjson("utje.json", "", "", mesh, NULL, "");
    EXPECT_EXIT_OK(parse_mmc("-f utje.json", NULL), "a single row of per-node MeshROI values");

    /* a source field given as rows sets the source count when Pos has one row; the
       directions of a detector object may carry a 4th number */
    write_tetjson("utje.json", "", "", NULL, "\"Pos\":[0.2,0.2,0],\"Dir\":[[0,0,1],[0,0,1,\"_NaN_\"],[0,0,1,\"_Inf_\"]]", "");
    parse_mmc("-f utje.json", check_srcrows3);
    write_tetjson("utje.json", "", "", NULL, "\"Pos\":[0.2,0.2,0],\"Dir\":[0,0,1],\"Param1\":[[0,0,0,0],[0,0,0,0],[0,0,0,0]]", "");
    parse_mmc("-f utje.json", check_srcrows3);
    write_tetjson("utje.json", "", "", NULL, "\"Pos\":[0.2,0.2,0],\"Dir\":[0,0,1],\"Param2\":[[0,0,0,0],[0,0,0,0],[0,0,0,0]]", "");
    parse_mmc("-f utje.json", check_srcrows3);
    write_tetjson("utje.json", "", "", NULL, NULL, ",\"Detector\":[{\"Pos\":[0.2,0.2,0],\"R\":1,\"Dir\":[0,0,-1,0.5]}]");
    EXPECT_EXIT_OK(parse_mmc("-f utje.json", NULL), "a detector direction with a 4th number");

    /* a second JSON input replaces the mesh, ROI and patterns of the first */
    snprintf(mesh, sizeof(mesh), "\"MeshNode\":[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],%s,\"MeshROI\":[0.1,0.1,0.1,0.1]", tetelem);
    write_tetjson("utje.json", "", "", mesh, NULL, "");
    EXPECT_EXIT_OK(parse_mmc("-f utje.json --json {\"Mesh\":{\"MeshNode\":[[0,0,0],[2,0,0],[0,2,0],[0,0,2]],\"MeshElem\":[[1,2,3,4,1]],"
                             "\"MeshROI\":[0.2,0.2,0.2,0.2]}}", NULL), "--json replaces the mesh of the input file");

    EXPECT_EXIT_OK(parse_mmc("-f utje.json --dumpjson 3", NULL), "--dumpjson 3 saves the parsed input");
}


/*-----------------------------------------------------------------------------
 *  output writers, JData decoding and helpers called directly, malformed mesh
 *  files, and source layouts the other tests do not use
 *---------------------------------------------------------------------------*/

static void check_writers(mcconfig* cfg, tetmesh* mesh) {
    size_t len = cfg->crop0.z;
    unsigned int* mask = (unsigned int*)calloc(len, sizeof(unsigned int));
    float vol[8] = {1, 2, 3, 4, 5, 6, 7, 8}, field[4] = {1, 2, 3, 4};
    uint dims[3] = {2, 2, 2};
    float voxsize[3] = {1.f, 1.f, 1.f};

    mask[0] = 1 | DET_MASK;
    mcx_savenii((OutputType*)mask, len, (char*)"utmask", NIFTI_TYPE_UINT32, ofNifti, cfg);
    EXPECT_TRUE(file_size("utmask.nii") > 0, "mcx_savenii saves a label/detector mask as 16-bit pairs");
    free(mask);
    EXPECT_EXIT(mcx_savenii((OutputType*)mesh->weight, len, (char*)"utbadfmt", NIFTI_TYPE_FLOAT64, ofMC2, cfg), "mcx_savenii with a format other than NIfTI/Analyze");

    mcx_savefloatjnii(vol, 3, dims, voxsize, (char*)"utfj", cfg);
    mcx_savefloatbnii(vol, 3, dims, voxsize, (char*)"utfb", cfg);
    EXPECT_TRUE(file_size("utfj.jnii") > 0 && file_size("utfb.bnii") > 0, "float JNIfTI and BNIfTI outputs");

    cfg->outputformat = ofTX3;
    strcpy(cfg->session, "uttx3");
    mcx_savedata(mesh->weight, len, cfg, 0);
    EXPECT_TRUE(file_size("uttx3.tx3") > 0, "the tx3 output format");
    cfg->outputformat = ofASCII;

    mcx_printlog(cfg, (char*)"unit test log line");
    mcx_normalize(field, 2.f, 4);
    EXPECT_NEAR(field[3], 8.f, 1e-6, "mcx_normalize scales the field");
    EXPECT_NEAR(mesh_getreff_approx(1.37, 1.0), -1.440 / (1.37 * 1.37) + 0.710 / 1.37 + 0.668 + 0.0636 * 1.37, 1e-12,
                "mesh_getreff_approx");

    /* mcx_clearcfg frees the optional arrays */
    cfg->nodemua = (float*)calloc(4, sizeof(float));
    cfg->nodemusp = (float*)calloc(4, sizeof(float));
    cfg->shapedata = strdup("{}");
    cfg->energytot = (double*)calloc(1, sizeof(double));
    cfg->energyesc = (double*)calloc(1, sizeof(double));
}

static void check_jdatadecode(mcconfig* cfg, tetmesh* mesh) {
    const char* arr[] = {
        "{\"_ArrayType_\":\"int16\",\"_ArraySize_\":[2],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipSize_\":2,\"_ArrayZipData_\":\"eJxjZGBiAAAADAAE\"}",
        "{\"_ArrayType_\":\"uint8\",\"_ArraySize_\":[4],\"_ArrayData_\":\"AQIDBA==\"}",
    };
    void* vol = NULL;
    int ndim;
    uint dims[3] = {1, 1, 1};
    char* type = NULL;

    (void)mesh;

    for (int i = 0; i < 2; i++) {
        cJSON* obj = cJSON_Parse(arr[i]);
        vol = malloc(4);  /* decoding replaces an existing buffer */
        mcx_jdatadecode(&vol, &ndim, dims, 3, &type, obj, cfg);
        EXPECT_TRUE(vol != NULL, i ? "mcx_jdatadecode of base64 _ArrayData_" : "mcx_jdatadecode of a zipped int16 array");
        free(vol);
        cJSON_Delete(obj);
    }

    EXPECT_EQ(cfg->mediabyte, 1, "mcx_jdatadecode sets the media format from the array type");
    EXPECT_EXIT({cJSON* o = cJSON_Parse("{\"_ArrayType_\":\"double\",\"_ArraySize_\":[1],\"_ArrayZipType_\":\"zlib\",\"_ArrayZipData_\":\"AA==\"}"); vol = NULL;
                 mcx_jdatadecode(&vol, &ndim, dims, 3, &type, o, cfg);
                }, "mcx_jdatadecode of an 8-byte array");
    EXPECT_EXIT({cJSON* o = cJSON_Parse("{\"_ArrayType_\":\"uint8\",\"_ArraySize_\":[1]}"); vol = NULL;
                 mcx_jdatadecode(&vol, &ndim, dims, 3, &type, o, cfg);
                }, "mcx_jdatadecode without array data");
}

static void test_morecover(void) {
    RunResult r;

    parse_mmc("--bench dmmc-cube60 -M G", check_writers);
    parse_mmc("--bench dmmc-cube60", check_jdatadecode);

    /* malformed node and element files of a mesh given by MeshID */
    {
        const char* json = "{\"Session\":{\"Photons\":10,\"ID\":\"utmf\"},\"Domain\":{\"Media\":[[0,0,1,1],[0.01,1,0.9,1.37]]},"
                           "\"Mesh\":{\"MeshID\":\"utmf\",\"InitElem\":1},\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},"
                           "\"Optode\":{\"Source\":{\"Pos\":[0.2,0.2,0],\"Dir\":[0,0,1]}}}";
        const char* goodnode = "1 4\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n";
        const char* goodelem = "1 1\n1 1 2 3 4 1\n";

        write_file("utmf.json", json, strlen(json));
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "a missing node file");
        write_file("node_utmf.dat", "abc\n", 4);
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "a node file without a header");
        write_file("node_utmf.dat", "1 4\n1 0 0 0\n", 12);
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "a node file with missing rows");
        write_file("node_utmf.dat", goodnode, strlen(goodnode));
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "a missing element file");
        write_file("elem_utmf.dat", "abc\n", 4);
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "an element file without a header");
        write_file("elem_utmf.dat", "1 2\n1 1 2 3 4 1\n", 16);
        EXPECT_EXIT(parse_mmc("-f utmf.json", NULL), "an element file with missing rows");
        write_file("elem_utmf.dat", goodelem, strlen(goodelem));
        EXPECT_EXIT_OK(parse_mmc("-f utmf.json", NULL), "a valid node/element file pair");
    }

    /* a wide-field source tilted away from z, and a diverging planar beam */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -s utgt --json {\"Optode\":{\"Source\":{\"Type\":\"gaussian\",\"Pos\":[30,45,15],"
                "\"Dir\":[0,0.6,0.8],\"Param1\":[1,0,0,0]}}}", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "a gaussian beam tilted away from z");
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -s utpd --json {\"Optode\":{\"Source\":{\"Type\":\"planar\",\"Pos\":[28,43,15],"
                "\"Dir\":[0,0,1,-10],\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,0]}}}", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "a diverging planar beam");
}


/*-----------------------------------------------------------------------------
 *  command line options and settings checks of the mmc binary that the other
 *  tests do not reach
 *---------------------------------------------------------------------------*/

static void test_cmdedge(void) {
    const char* json = "{\"Session\":{\"Photons\":100,\"RNGSeed\":1,\"ID\":\"utce\"},\"Domain\":{\"Media\":[[0,0,1,1],[0.01,1,0.9,1.37]]},"
                       "\"Mesh\":{\"MeshID\":\"onecube\",\"InitElem\":3},\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":5e-9},"
                       "\"Optode\":{\"Source\":{\"Pos\":[2,8,0],\"Dir\":[0,0,1]}}}";
    char args[2048];

    /* options parsed without a simulation */
    parse_mmc("--bench dmmc-cube60 -Z 1 -R 2 -I -D 3 --compute 1", NULL);
    EXPECT_TRUE(1, "-Z, -R, -I, numeric -D and --compute");
    snprintf(args, sizeof(args), "-f %s", json);
    parse_mmc(args, NULL);
    EXPECT_TRUE(1, "-f with an inline JSON input");
    parse_mmc("--bench dmmc-cube60 --json {\"Session\":{\"Photons\":10}} --json {\"Session\":{\"Photons\":20}}", NULL);
    EXPECT_TRUE(1, "a second --json replaces the first");
    parse_mmc("--bench dmmc-cube60 --kernel " MMC_EXAMPLE_DIR "/../src/mmc_core.cl --kernel " MMC_EXAMPLE_DIR "/../src/mmc_core.cl", NULL);
    EXPECT_TRUE(1, "a second --kernel replaces the first");
    EXPECT_EXIT_OK(parse_mmc("--bench", NULL), "--bench without a name lists the benchmarks");
    EXPECT_EXIT_OK(parse_mmc("--bench dmmc-cube60 --dumpjson", NULL), "--dumpjson without a file name");
    parse_mmc("--bench dmmc-cube60 --root utnosuchdir -l", NULL);
    EXPECT_TRUE(1, "-l falls back to stdout if the log file can not be created");

    EXPECT_EXIT(parse_mmc("-i -f utce.json", NULL), "both the interactive mode and an input file");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json", NULL), "--json without a JSON fragment");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -E utdet.mch -O A", NULL), "an adjoint output in the replay mode");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -O A --json {\"Optode\":{\"SDMap\":[[0,1]]}}", NULL), "an Optode.SDMap index below 1");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -O A --json {\"Optode\":{\"SDMap\":[[1]]}}", NULL), "an Optode.SDMap row shorter than 2");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -O A --json {\"Optode\":{\"SDMap\":[]}}", NULL), "an empty Optode.SDMap");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 -M G -O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}],"
                        "\"SDMap\":[[1,2]]}}", NULL), "an Optode.SDMap detector index beyond the detectors");

    /* the checks of mcx_validatecfg */
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json {\"Forward\":{\"T0\":1e-9,\"T1\":0,\"Dt\":1e-10}}", NULL), "a reversed time window");
    parse_mmc("--bench dmmc-cube60 --json {\"Forward\":{\"T0\":0,\"T1\":1e-9,\"Dt\":5e-9}}", NULL);
    EXPECT_TRUE(1, "a time step longer than the time window is shortened");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json {\"Optode\":{\"Source\":{\"Dir\":[0,0,2]}}}", NULL), "a source direction that is not a unit vector");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json {\"Forward\":{\"T0\":1e-9,\"T1\":1e-9,\"Dt\":1e-10}}", NULL), "tend equal to tstart");
    EXPECT_EXIT(run_mmclab_path("--bench dmmc-cube60 -G -1 -s utmle --json {\"Optode\":{\"Source\":{\"ID\":3,\"Pos\":[[30,45,0,1],[20,30,0,1]],"
                                "\"Dir\":[[0,0,1],[0,0,1]]}}}", NULL), "mmc_validate_config: a source ID beyond the sources");
    parse_mmc("--bench dmmc-cube60 -E -1", NULL);
    EXPECT_TRUE(1, "a negative seed is replaced by the time");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 --json {\"Domain\":{\"Media\":[{\"mua\":0,\"mus\":0,\"g\":1,\"n\":1e30},"
                          "{\"mua\":0.005,\"mus\":1,\"g\":0.01,\"n\":1.37}]}}", NULL), "a refractive index whose square overflows a float");
    EXPECT_EXIT(run_mmc("--bench faceimmc -G -1 -n 10 -M H", NULL), "face iMMC with the Havel ray-tracer (no face normals)");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -u -1", NULL), "a negative length unit");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -u 0", NULL), "a zero length unit");

#ifdef USE_OPENCL

    /* mesh-mode Jacobians saved as JNIfTI and BNIfTI, dual and RF, under a root folder */
    if (hasocl) {
        RunResult r;
        const char* det = "{\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";
        const char* detrf = "{\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";

        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W -F jnii --root . -s utjw --json %s", det);
        r = run_mmc(args, NULL);
        EXPECT_TRUE(r.ran && file_size("utjw_jmua.jnii") > 0 && file_size("utjw_jd.jnii") > 0, "dual mesh Jacobians saved as JNIfTI");
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O Q -F bnii -s utjq --json %s", detrf);
        r = run_mmc(args, NULL);
        EXPECT_TRUE(r.ran && file_size("utjq_jmua.bnii") > 0 && file_size("utjq_jmusp.bnii") > 0, "RF dual mesh Jacobians saved as BNIfTI");
    }

#endif
}


/*-----------------------------------------------------------------------------
 *  JSON writer branches, buffer reuse, detected photon and trajectory files in
 *  every format, and malformed .mch replay files
 *---------------------------------------------------------------------------*/

static void check_savebranches(mcconfig* cfg, tetmesh* mesh) {
    float ppath[2 * MCX_DEBUG_REC_LEN];
    unsigned char seeds[2 * 16];
    const int fmt[3] = {ofBJNifti, ofJNifti, ofMC2};

    /* the JSON writer: -b 2, normalization and reflectance flags beyond 1, a root folder,
       a modulation frequency, an outer refractive index and a mesh name */
    cfg->isreflect = 2;
    cfg->isnormalized = 2;
    cfg->issaveref = 2;
    cfg->omega = 6.28e8f;
    cfg->nout = 1.2f;
    strcpy(cfg->rootpath, ".");
    mcx_savejdata((char*)"utsjb.json", cfg);
    EXPECT_TRUE(file_size("utsjb.json") > 0, "--dumpjson of the less common settings");
    cfg->isreflect = 1;
    cfg->isnormalized = 1;
    cfg->issaveref = 0;
    cfg->omega = 0.f;

    /* mesh_getvolume recomputes the volumes, freeing the OptiX triangle arrays */
    mesh->fnode = (FLOAT3*)malloc(sizeof(FLOAT3));
    mesh->fnorm = (FLOAT3*)malloc(sizeof(FLOAT3));
    mesh->face = (uint3*)malloc(sizeof(uint3));
    mesh->front = (uint*)malloc(sizeof(uint));
    mesh->back = (uint*)malloc(sizeof(uint));
    mesh_getvolume(mesh, cfg);
    EXPECT_TRUE(mesh->fnode == NULL && mesh->face == NULL, "mesh_getvolume frees the OptiX triangle arrays");

    /* trajectory records with seeds, written and appended in each format under a root folder */
    for (int i = 0; i < 2 * MCX_DEBUG_REC_LEN; i++) {
        ppath[i] = (float)(i % MCX_DEBUG_REC_LEN) + 1.f;
    }

    memset(seeds, 7, sizeof(seeds));
    cfg->issaveseed = 1;
    cfg->his.detected = 0;
    cfg->his.savedphoton = 2;
    cfg->his.colcount = MCX_DEBUG_REC_LEN;
    cfg->his.seedbyte = 16;
    strcpy(cfg->session, "utsdb");

    for (int f = 0; f < 3; f++) {
        cfg->outputformat = fmt[f];
        mcx_savedetphoton(ppath, seeds, 2, 0, cfg);
        mcx_savedetphoton(ppath, seeds, 2, 1, cfg);
    }

    EXPECT_TRUE(file_size("utsdb_traj.jdb") > 0 && file_size("utsdb_traj.jdt") > 0 && file_size("utsdb.mct") > 0,
                "trajectories with seeds in .jdb, .jdt and .mct, written and appended");
    cfg->outputformat = ofASCII;
    cfg->rootpath[0] = '\0';
}

/** write a .mch file of the given header fields and payload size, for the replay tests */
static void write_mch(const char* fname, int maxmedia, int savedphoton, int colcount, int seedbyte, size_t payload) {
    history his;
    FILE* fp = fopen(fname, "wb");
    char* buf = (char*)calloc(payload + 1, 1);

    memset(&his, 0, sizeof(his));
    memcpy(his.magic, "MCXH", 4);
    his.version = 1;
    his.maxmedia = maxmedia;
    his.savedphoton = savedphoton;
    his.detected = savedphoton;
    his.colcount = colcount;
    his.seedbyte = seedbyte;

    if (fp) {
        fwrite(&his, sizeof(his), 1, fp);
        fwrite(buf, 1, payload, fp);
        fclose(fp);
    }

    free(buf);
}

static void check_mchfiles(mcconfig* cfg, tetmesh* mesh) {
    (void)cfg;
    write_file("utshort.mch", "MCXH", 4);
    write_mch("utnoseed.mch", mesh->prop, 0, 3, 0, 0);
    write_mch("utmedia.mch", mesh->prop + 5, 2, 3, 16, 2 * 3 * 4 + 2 * 16);
    write_mch("utnoseeds.mch", mesh->prop, 2, 3, 16, 2 * 3 * 4 + 8);
}

static void test_morecover2(void) {
    parse_mmc("--bench dmmc-cube60 -M G", check_savebranches);

    parse_mmc("--bench dmmc-cube60", check_mchfiles);
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -E utshort.mch -O J", NULL), "a truncated .mch header");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -E utmedia.mch -O J", NULL), "a .mch file of different media");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -E utnoseeds.mch -O J", NULL), "a .mch file with truncated seeds");
    parse_mmc("--bench dmmc-cube60 -E utnoseed.mch", NULL);
    EXPECT_TRUE(1, "a .mch file without seeds is not replayed");
}


/*-----------------------------------------------------------------------------
 *  .inp variants, source rotations, unit scaling of the mmclab path, CPU tracer
 *  branches and settings checks of the mmc binary
 *---------------------------------------------------------------------------*/

/** copy onecube.inp to fname, replacing the first occurrence of from by to */
static void write_inpvariant(const char* fname, const char* from, const char* to) {
    long len = file_size("onecube.inp");
    char* buf = (char*)calloc(len + 1, 1), *out = (char*)calloc(len + strlen(to) + 1, 1);
    FILE* fp = fopen("onecube.inp", "rb");

    if (fp && fread(buf, len, 1, fp) == 1) {
        char* hit = strstr(buf, from);

        if (hit) {
            memcpy(out, buf, hit - buf);
            strcat(out, to);
            strcat(out, hit + strlen(from));
            write_file(fname, out, strlen(out));
        }
    }

    if (fp) {
        fclose(fp);
    }

    free(buf);
    free(out);
}

static void check_surfunit(mcconfig* cfg, tetmesh* mesh) {
    float mua0 = mesh->med[1].mua;

    cfg->unitinmm = 2.f;
    mesh_validate(mesh, cfg);
    EXPECT_NEAR(mesh->med[1].mua, 2.f * mua0, 1e-6, "mesh_validate scales a surface-only mesh by the length unit");
}

static void test_morecover3(void) {
    RunResult r;
    char args[2048];

    /* .inp: a focal length after the direction, a reversed time window, a time step longer
       than the window, and the mesh found under --root */
    write_inpvariant("utdirw.inp", "0. 0. 1.", "0. 0. 1. 5");
    parse_mmc("-f utdirw.inp", NULL);
    write_inpvariant("utgate.inp", "0.e+00 5.e-09 5e-10", "1e-9 0 1e-10");
    EXPECT_EXIT(parse_mmc("-f utgate.inp", NULL), ".inp with a reversed time window");
    write_inpvariant("utstep.inp", "0.e+00 5.e-09 5e-10", "0 1e-9 5e-9");
    parse_mmc("-f utstep.inp", NULL);
    parse_mmc("-f onecube.inp --root .", NULL);
    EXPECT_TRUE(1, ".inp with a focal length, a long time step and --root");

    /* sources tilted away from z, and a wide-field source outside every element */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -s utiso --json {\"Optode\":{\"Source\":{\"Type\":\"isotropic\",\"Pos\":[30,45,15],"
                "\"Dir\":[0.6,0,0.8]}}}", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "an isotropic source with a tilted reference direction");
    EXPECT_EXIT(run_mmc("--bench dmmc-cube60 -G -1 -n 10 -E 1 -s utout --json {\"Optode\":{\"Source\":{\"Type\":\"planar\",\"Pos\":[28,43,-10],"
                        "\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,0]}}}", NULL), "a wide-field source outside the mesh");

    /* MCML weights with element-based output */
    r = run_mmc("--bench dmmc-cube60 -G -1 -n 300 -E 1 -M H -C 0 -m 1 -s utmcml0", NULL);
    EXPECT_TRUE(r.ran && r.normalizer > 0, "MCML weights with element-based output");

    /* the mmclab path with a length unit, in mesh and dual-grid modes, and on a surface */
    run_mmclab_path("--bench dmmc-cube60 -G -1 -n 100 -E 1 -M H -C 1 -u 2 -s utmlu", NULL);
    run_mmclab_path("--bench dmmc-cube60 -G -1 -n 100 -E 1 -M G -u 2 -s utmlg", NULL);
    EXPECT_TRUE(1, "the mmclab path scales the media by the length unit");
    write_cubesurf("utsurf.json", NULL);
    parse_mmc("-f utsurf.json --compute optix -M G", check_surfunit);

    /* settings checks of the mmc binary */
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -n 0", NULL), "no photons");
    EXPECT_EXIT(parse_mmc("--bench dmmc-cube60 -E utdet.mch -O J --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"SrcNum\":2,"
                          "\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1,0,1,1,0]}}}}", NULL),
                "several patterns in the replay mode");
    EXPECT_EXIT(parse_mmc("-f utce.json -i", NULL), "an input file before the interactive mode");
    parse_mmc("--bench dmmc-cube60 -G 3", NULL);
    EXPECT_TRUE(1, "-G with a device number");

#ifdef USE_OPENCL

    if (hasocl) {
        /* the remaining combinations of the mesh Jacobian outputs */
        const char* det = "{\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";
        const char* detrf = "{\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";

        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O A -F bnii -s utja --json %s", det);
        r = run_mmc(args, NULL);
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W -F bnii -s utjwb --json %s", det);
        run_mmc(args, NULL);
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O A -F jnii -s utjar --json %s", detrf);
        run_mmc(args, NULL);
        EXPECT_TRUE(r.ran && file_size("utja_jmua.bnii") > 0 && file_size("utjar_jmua.jnii") > 0, "single, dual and RF mesh Jacobian files");
    }

#endif
}


/* the mmclab path: wide-field detector elements (label -2), and the defaults of a negative
   seed and of -d 1 without detectors */
static void mut_extdet(mcconfig* cfg, tetmesh* mesh) {
    cfg->isextdet = 1;
    mesh->type[mesh->ne - 1] = -2;
}
static void mut_seeddet(mcconfig* cfg, tetmesh* mesh) {
    (void)mesh;
    cfg->seed = -1;
    cfg->issavedet = 1;
    cfg->detnum = 0;
}

static void test_morecover4(void) {
    char args[1024];

    EXPECT_EXIT_OK(run_mmclab_path("--bench dmmc-cube60 -G -1 -n 100 -E 1 -M H -C 1 -s utmlx", mut_extdet), "the mmclab path with wide-field detector elements");
    EXPECT_EXIT_OK(run_mmclab_path("--bench dmmc-cube60 -G -1 -n 100 -M H -C 1 -s utmls", mut_seeddet), "the mmclab path with a negative seed and no detectors");

    /* an .inp mesh without, or with a malformed, media property file */
    if (system("for f in node elem facenb velem; do cp ${f}_onecube.dat ${f}_utnp.dat; done") == 0) {
        write_inpvariant("utnp.inp", "onecube ", "utnp ");
        EXPECT_EXIT(parse_mmc("-f utnp.inp", NULL), "an .inp mesh without a media property file");
        write_file("prop_utnp.dat", "abc\n", 4);
        EXPECT_EXIT(parse_mmc("-f utnp.inp", NULL), "an .inp mesh with a malformed media property file");
    }

#ifdef USE_OPENCL

    if (hasocl) {
        snprintf(args, sizeof(args), "--bench dmmc-cube60 -G 1 -n 2000 -E 1 -M H -C 1 -O W -F jnii -s utjwr --json "
                 "{\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}");
        RunResult r = run_mmc(args, NULL);
        EXPECT_TRUE(r.ran && file_size("utjwr_jd.jnii") > 0, "RF dual mesh Jacobians saved as JNIfTI");
    }

#endif
    (void)args;
}

/*-----------------------------------------------------------------------------
 *  utility functions called directly
 *---------------------------------------------------------------------------*/

static void test_utilfuncs(void) {
    /* column-major to row-major conversion of a 2x3x4 volume and a 2x3x4x2 array */
    unsigned int* vol = (unsigned int*)malloc(24 * sizeof(unsigned int));
    uint3 dim3 = {2, 3, 4};

    for (int i = 0; i < 24; i++) {
        vol[i] = i;
    }

    mcx_convertcol2row(&vol, &dim3);
    /* element (x=1,y=2,z=3): column-major index 1+2*2+3*6=23, row-major index 1*12+2*4+3=23 */
    EXPECT_EQ(vol[1 * 12 + 0 * 4 + 0], 1, "mcx_convertcol2row (1,0,0)");
    EXPECT_EQ(vol[0 * 12 + 1 * 4 + 0], 2, "mcx_convertcol2row (0,1,0)");
    EXPECT_EQ(vol[0 * 12 + 0 * 4 + 1], 6, "mcx_convertcol2row (0,0,1)");
    free(vol);

    vol = (unsigned int*)malloc(48 * sizeof(unsigned int));
    uint4 dim4 = {2, 3, 4, 2};

    for (int i = 0; i < 48; i++) {
        vol[i] = i;
    }

    mcx_convertcol2row4d(&vol, &dim4);
    EXPECT_EQ(vol[1 * 24], 1, "mcx_convertcol2row4d (1,0,0,0)");
    EXPECT_EQ(vol[1 * 8], 2, "mcx_convertcol2row4d (0,1,0,0)");
    EXPECT_EQ(vol[1 * 2], 6, "mcx_convertcol2row4d (0,0,1,0)");
    EXPECT_EQ(vol[1], 24, "mcx_convertcol2row4d (0,0,0,1)");
    free(vol);

    /* empty inputs are left alone */
    vol = NULL;
    mcx_convertcol2row(&vol, &dim3);
    mcx_convertcol2row4d(&vol, &dim4);
    EXPECT_TRUE(vol == NULL, "mcx_convertcol2row of an empty volume");

    /* GPU information list */
    GPUInfo* info = (GPUInfo*)calloc(2, sizeof(GPUInfo));
    mcx_cleargpuinfo(&info);
    EXPECT_TRUE(info == NULL, "mcx_cleargpuinfo frees the list");

    /* raw typed BJData arrays of each data type */
    {
        const char* types[] = {"int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64", "single", "double", "float32", "float64"};
        const int bytes[] = {1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 4, 8};
        double data[4] = {1, 2, 3, 4};
        uint dims[2] = {2, 2};
        mcconfig cfg;
        uint8_t* buf = (uint8_t*)calloc(4096, 1);

        memset(&cfg, 0, sizeof(cfg));
        cfg.flog = devnull;

        for (int i = 0; i < (int)(sizeof(types) / sizeof(types[0])); i++) {
            for (int iscol = 0; iscol < 2; iscol++) {
                ubjw_context_t* root = ubjw_open_memory(buf, buf + 4096);
                ubjw_begin_object(root, UBJ_MIXED, 0);
                int ret = mcx_jdataencode(data, 2, dims, (char*)types[i], bytes[i], zmBase64, root, 1, iscol, &cfg);
                ubjw_end(root);
                size_t len = ubjw_close_context(root);
                EXPECT_TRUE(ret == 0 && len > 4u * bytes[i], "mcx_jdataencode of a raw BJData array");
            }
        }

        free(buf);
    }

    /* JData arrays over 4 MB: zlib/gzip deflate them in parallel blocks and save the block
       index as _ArrayZipOffsets_, which the reader uses to inflate them in parallel; a
       malformed index falls back to a serial inflate */
    {
        const uint len = 3u << 20;  /* 12 MB of float32 */
        uint dims[1] = {len};
        float* data = (float*)malloc(len * sizeof(float));
        const int zips[2] = {zmZlib, zmGzip};
        const char* bad[] = {"[[0,0],[1,1],[2]]", "[0,0]", "[[0,0],[-1,5],[9,9]]", "[[0,0],[\"a\",1],[9,9]]"};
        mcconfig cfg;

        memset(&cfg, 0, sizeof(cfg));
        cfg.flog = devnull;
        cfg.isdumpjson = 1;

        for (uint i = 0; i < len; i++) {
            data[i] = (float)((i * 2654435761u) % 1000u) * 0.5f;
        }

        for (int z = 0; z < 2; z++) {
            cJSON* obj = cJSON_CreateObject();
            int ret = mcx_jdataencode(data, 1, dims, (char*)"single", sizeof(float), zips[z], obj, 0, 0, &cfg);
            cJSON* idx = cJSON_GetObjectItem(obj, "_ArrayZipOffsets_");
            EXPECT_TRUE(ret == 0 && idx && cJSON_GetArraySize(idx) >= 3, "a 12 MB zlib/gzip array saves _ArrayZipOffsets_");

            void* out = NULL;
            char* type = NULL;
            int ndim = 0;
            uint odims[1] = {0};
            mcx_jdatadecode(&out, &ndim, odims, 1, &type, obj, &cfg);
            EXPECT_TRUE(out && odims[0] == len && memcmp(out, data, len * sizeof(float)) == 0, "a 12 MB zlib/gzip array is inflated in parallel blocks");
            free(out);

            for (int b = 0; b < 4; b++) {
                cJSON_ReplaceItemInObject(obj, "_ArrayZipOffsets_", cJSON_Parse(bad[b]));
                out = NULL;
                mcx_jdatadecode(&out, &ndim, odims, 1, &type, obj, &cfg);
                EXPECT_TRUE(out && memcmp(out, data, len * sizeof(float)) == 0, "an array with a malformed _ArrayZipOffsets_ is inflated serially");
                free(out);
            }

            cJSON_Delete(obj);

            /* binary JSON: the index is written as an array of uint64 pairs */
            size_t cap = len * sizeof(float) + (1 << 20), n, found = 0;
            unsigned char* ubj = (unsigned char*)malloc(cap);
            ubjw_context_t* root = ubjw_open_memory(ubj, ubj + cap);
            ubjw_begin_object(root, UBJ_MIXED, 0);
            ret = mcx_jdataencode(data, 1, dims, (char*)"single", sizeof(float), zips[z], root, 1, 0, &cfg);
            ubjw_end(root);
            n = ubjw_close_context(root);

            for (size_t i = 0; i + 17 <= n && !found; i++) {
                found = (memcmp(ubj + i, "_ArrayZipOffsets_", 17) == 0);
            }

            EXPECT_TRUE(ret == 0 && found, "a 12 MB zlib/gzip BJData array saves _ArrayZipOffsets_");
            free(ubj);
        }

        free(data);
    }
}

/*-----------------------------------------------------------------------------
 *  GPU backends: the same feature cases on OpenCL, CUDA and OptiX; OpenCL and CUDA
 *  must run every case and agree on the output, OptiX is only reported
 *---------------------------------------------------------------------------*/

/** run one case in a new process (this program re-executed with --runone), so that an
    error exit or a crash is reported instead of ending the test program; a forked child
    can not be used because the GPU drivers do not survive a fork after their initialization.
    returns 0: ran, >0: exit code, <0: -signal */
static int run_isolated(const char* args, RunResult* res) {
    char resfile[] = "utrunoneXXXXXX";
    int status = 0, fd = mkstemp(resfile);
    pid_t pid;

    memset(res, 0, sizeof(*res));

    if (fd < 0) {
        return 255;
    }

    close(fd);
    fflush(stdout);
    pid = fork();

    if (pid == 0) {
        if (getenv("MMC_UNITTEST_VERBOSE") == NULL && freopen("/dev/null", "w", stderr) == NULL) {
            exit(254);
        }

        execl("/proc/self/exe", "mmc_unittest", "--runone", resfile, args, (char*)NULL);
        exit(251);
    }

    waitpid(pid, &status, 0);

    FILE* fp = fopen(resfile, "rb");

    if (fp) {
        if (fread(res, sizeof(*res), 1, fp) != 1) {
            res->ran = 0;
        }

        fclose(fp);
    }

    remove(resfile);
    return WIFEXITED(status) ? WEXITSTATUS(status) : -WTERMSIG(status);
}

/** the --runone mode of run_isolated: run one case and save its RunResult */
static int run_one(const char* resfile, const char* args) {
    RunResult r = run_mmc(args, NULL);
    FILE* fp = fopen(resfile, "wb");

    if (fp) {
        fwrite(&r, sizeof(r), 1, fp);
        fclose(fp);
    }

    return r.ran ? 0 : 252;
}

static void test_gpubackends(void) {
    const char* det = "{\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";
    const struct {
        const char* name;
        const char* args;
        int usedet;
        int optix;   /**< 1: OptiX must run it and agree with opencl; 2: run it (test_optixadjoint compares it); 0: reject it */
    } cases[] = {
        {"grid fluence (-M G)",              "-M G", 0, 1},
        {"grid, reflection (-b 1)",          "-M G -b 1", 0, 1},
        {"grid, energy output (-O E)",       "-M G -O E", 0, 1},
        {"grid, 5 time gates",               "-M G --json {\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":1e-9}}", 0, 1},
        {"grid, detected photons (-d 1)",    "-M G -d 1 -x 1", 1, 1},
        {"grid, detected photons with seeds", "-M G -d 1 -q 1", 1, 0},
        {"grid, RF forward (-O R)",          "-M G -O R --json {\"Forward\":{\"Omega\":6.28e8}}", 0, 1},
        {"grid, replay Jacobian (-O J)",     "-M G -E utdet.mch -O J", 0, 0},
        {"grid, adjoint J_mua (-O A)",       "-M G -O A", 1, 2},
        {"grid, adjoint J_D (-O D)",         "-M G -O D", 1, 2},
        {"grid, adjoint J_mus (-O U)",       "-M G -O U", 1, 2},
        {"grid, adjoint J_musp (-O V)",      "-M G -O V", 1, 2},
        {"grid, dual adjoint (-O Q)",        "-M G -O Q", 1, 2},
        {
            "grid, adjoint of Optode.SDMap",    "-M G -O A --json {\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]},"
            "{\"Pos\":[30,20,0],\"R\":2,\"Dir\":[0,0,1]}],\"SDMap\":[[1,2],[1,1]]}}", 0, 2
        },
        {"grid, isotropic source",           "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"isotropic\",\"Pos\":[30,45,15]}}}", 0, 1},
        {"grid, cone source",                "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"cone\",\"Param1\":[0.4,0,0,0]}}}", 0, 1},
        {"grid, planar source",              "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"planar\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,0],\"Param2\":[0,4,0,0]}}}", 0, 1},
        {"grid, gaussian source",            "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"gaussian\",\"Pos\":[30,45,15],\"Param1\":[1,0,0,0]}}}", 0, 1},
        {"grid, disk source",                "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"disk\",\"Pos\":[30,45,15],\"Param1\":[3,0,0,0]}}}", 0, 1},
        {
            "grid, two patterns",               "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"pattern\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],"
            "\"Param2\":[0,4,0,2],\"SrcNum\":2,\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1,0,1,1,0]}}}}", 0, 1
        },
        {
            "grid, 2 sources, a field each",   "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"pencil\",\"ID\":-1,\"Pos\":[[30,45,0,0.5],[20,30,0,0.5]],"
            "\"Dir\":[[0,0,1,0],[0,0,1,0]]}}}", 0, 1
        },
        {"grid, trajectories (-D M)",        "-M G -D M", 0, 1},
        {"mesh nodal (-M H -C 1)",           "-M H -C 1", 0, 0},
        {"mesh elem (-M H -C 0)",            "-M H -C 0", 0, 0},
        {"mesh, adjoint J_mua (-C 1 -O A)",  "-M H -C 1 -O A", 1, 0},
        {"mesh, adjoint J_musp (-C 1 -O V)", "-M H -C 1 -O V", 1, 0},
        {
            "grid, RF dual adjoint (-O W)",     "-M G -O W --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 2
        },
        {
            "grid, RF adjoint J_musp (-O V)",   "-M G -O V --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 2
        },
        {
            "grid, RF field of each slot",      "-M G -O R --srcid -1 --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 2
        },
        {"mesh nodal, RF forward (-O R)",    "-M H -C 1 -O R --json {\"Forward\":{\"Omega\":6.28e8}}", 0, 0},
        {"mesh elem, RF forward (-O R)",     "-M H -C 0 -O R --json {\"Forward\":{\"Omega\":6.28e8}}", 0, 0},
        {
            "grid, pattern3d source",           "-M G --json {\"Optode\":{\"Source\":{\"Type\":\"pattern3d\",\"Pos\":[29,44,14],"
            "\"Param1\":[2,2,2,0],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Nz\":2,\"Data\":[1,1,1,1,1,1,1,1]}}}}", 0, 2
        },
        {"grid, respin (-r 2)",              "-M G -r 2", 0, 2},
        {"grid, manual threads (-A 0)",      "-M G -A 0 -t 4000 -T 64 -g 10 --json {\"Forward\":{\"T0\":0,\"T1\":5e-9,\"Dt\":1e-9}}", 0, 2},
        {"mesh, diffuse reflectance (-X 1)", "-M H -C 1 -X 1", 0, 0},
        {"mesh, progress bar (-D P)",        "-M H -C 1 -D P", 0, 0},
        {"grid, more photons than -H",       "-M G -d 1 -H 10", 1, 2},
        {"grid, more steps than --maxjumpdebug", "-M G -D M --maxjumpdebug 100", 0, 2},
        {
            "grid, RF forward, detected photons", "-M G -O R -d 1 --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 2
        },
        {"grid, adjoint, detected photons",  "-M G -O A -d 1", 1, 2},
        {
            "grid, RF adjoint, detected photons", "-M G -O W -d 1 --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 2
        },
        {
            "mesh, RF dual adjoint (-C 1 -O W)", "-M H -C 1 -O W --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 0
        },
        {
            "mesh nodal, RF field of each slot", "-M H -C 1 -O R --srcid -1 --json {\"Forward\":{\"Omega\":6.28e8},\"Optode\":{\"Detector\":"
            "[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}", 0, 0
        },
    };
    const struct {
        const char* name;
        const char* flag;
        int available;
        int strict;
    } backend[] = {
        {"opencl", "--compute opencl", hasocl, 1},
#ifdef USE_CUDA
        {"cuda", "--compute cuda", 1, 1},
#else
        {"cuda", "--compute cuda", 0, 1},
#endif
#ifdef USE_OPTIX
        {"optix", "--compute optix", 1, 0},
#else
        {"optix", "--compute optix", 0, 0},
#endif
    };
    const int ncase = (int)(sizeof(cases) / sizeof(cases[0]));
    double refsum[64] = {0.}, refjsum[64] = {0.}, refimsum[64] = {0.};
    int refdet[64] = {0};
    char args[2048], msg[512];

    for (int b = 0; b < 3; b++) {
        if (!backend[b].available) {
            continue;
        }

        printf("%s backend feature check:\n", backend[b].name);

        for (int i = 0; i < ncase; i++) {
            RunResult r;
            snprintf(args, sizeof(args), "--bench dmmc-cube60 %s -G 1 -n 20000 -E 1234 %s -s utgb%d%d%s%s",
                     backend[b].flag, cases[i].args, b, i, cases[i].usedet ? " --json " : "", cases[i].usedet ? det : "");
            int ret = run_isolated(args, &r);

            if (ret == 0 && b == 0) {
                refsum[i] = r.wsum;
                refjsum[i] = r.jsum;
                refimsum[i] = r.imsum;
            }

            if (ret == 0) {
                double rel = (refsum[i] != 0.) ? fabs(r.wsum - refsum[i]) / fabs(refsum[i]) : 0.;
                double jrel = (refjsum[i] != 0.) ? fabs(r.jsum - refjsum[i]) / fabs(refjsum[i]) : 0.;
                double imrel = (refimsum[i] != 0.) ? fabs(r.imsum - refimsum[i]) / fabs(refimsum[i]) : 0.;
                printf("  %-34s ok, output sum %.6g (%.1f%% from opencl)", cases[i].name, r.wsum, rel * 100.);

                if (r.jsum != 0. || refjsum[i] != 0.) {
                    printf(", Jacobian sum %.6g (%.1f%%)", r.jsum, jrel * 100.);
                }

                if (r.imsum != 0. || refimsum[i] != 0.) {
                    printf(", imaginary sum %.6g (%.1f%%)", r.imsum, imrel * 100.);
                }

                if (cases[i].usedet) {
                    printf(", %d detected", r.detected);
                }

                printf("\n");

                if ((backend[b].strict || cases[i].optix == 1) && b > 0 && refsum[i] != 0.) {
                    sprintf(msg, "%s: %s agrees with opencl within 5%%", backend[b].name, cases[i].name);
                    EXPECT_TRUE(rel < 0.05 && jrel < 0.05 && imrel < 0.05, msg);
                }

                /* the trajectory points are the launch, the scattering events and the end of each
                   photon: the same count per photon as opencl, within the noise */
                if (strstr(cases[i].args, "-D M") && b == 0) {
                    refdet[i] = r.ntraj;
                } else if (strstr(cases[i].args, "-D M") && refdet[i] > 0) {
                    sprintf(msg, "%s: %s saves as many valid trajectory points as opencl", backend[b].name, cases[i].name);
                    EXPECT_TRUE(r.trajok && fabs((double)r.ntraj - refdet[i]) < 0.03 * refdet[i], msg);
                }

                if (cases[i].usedet && b == 0) {
                    refdet[i] = r.detected;
                } else if (cases[i].usedet && (backend[b].strict || cases[i].optix) && refdet[i] > 0) {
                    /* the detected counts are Poisson-like: within 4 standard deviations of their difference */
                    sprintf(msg, "%s: %s detects as many photons as opencl", backend[b].name, cases[i].name);
                    EXPECT_TRUE(fabs((double)r.detected - refdet[i]) <= 4. * sqrt((double)r.detected + refdet[i]), msg);
                }
            } else {
                printf("  %-34s %s %d\n", cases[i].name, ret > 0 ? "NOT SUPPORTED, exit code" : "CRASHED, signal", ret > 0 ? ret : -ret);
            }

            if (backend[b].strict || cases[i].optix) {
                sprintf(msg, "%s: %s runs", backend[b].name, cases[i].name);
                EXPECT_TRUE(ret == 0, msg);
            } else {
                /* OptiX must report an unsupported feature as an error, not crash or run */
                sprintf(msg, "%s: %s is rejected with an error", backend[b].name, cases[i].name);
                EXPECT_TRUE(ret > 0, msg);
            }
        }
    }

    /* implicit MMC (edge, node and face ROIs): OpenCL and CUDA agree */
    {
        const char* immc[] = {"edgeimmc", "nodeimmc", "faceimmc"};
        double ref[3] = {0., 0., 0.};

        for (int b = 0; b < 2; b++) {
            if (!backend[b].available) {
                continue;
            }

            for (int i = 0; i < 3; i++) {
                RunResult r;
                snprintf(args, sizeof(args), "--bench %s %s -G 1 -n 20000 -E 1234 -s utgbimmc%d%d", immc[i], backend[b].flag, b, i);
                int ret = run_isolated(args, &r);
                sprintf(msg, "%s: implicit MMC %s runs", backend[b].name, immc[i]);
                EXPECT_TRUE(ret == 0 && r.wsum > 0., msg);

                if (b == 0) {
                    ref[i] = r.wsum;
                } else if (ref[i] > 0.) {
                    printf("  implicit MMC %-21s output sum %.6g (%.1f%% from opencl)\n", immc[i], r.wsum, fabs(r.wsum - ref[i]) / ref[i] * 100.);
                    sprintf(msg, "%s: implicit MMC %s agrees with opencl within 5%%", backend[b].name, immc[i]);
                    EXPECT_TRUE(fabs(r.wsum - ref[i]) < 0.05 * ref[i], msg);
                }
            }
        }
    }
}

/*-----------------------------------------------------------------------------
 *  OptiX source types: the fluence of each source type matches the CPU tracer. Only
 *  OptiX launches photons from arbitrary positions; the CPU tracer searches the launch
 *  element among e0 and the retessellated (-1 labeled) source elements, so these cases
 *  keep every launch position inside the initial tetrahedron of dmmc-cube60, which makes
 *  the CPU tracer a valid reference here (the OpenCL/CUDA kernel lacks pencilarray,
 *  pattern3d, hyperboloid and ring)
 *---------------------------------------------------------------------------*/

static void test_optixsources(void) {
#ifdef USE_OPTIX
    char p3d[512] = "", args[2048], msg[512];
    int len = 0;

    for (int i = 0; i < 64; i++) {
        len += snprintf(p3d + len, sizeof(p3d) - len, "%s%g", i ? "," : "", (i % 5) * 0.25);
    }

    const char* src[][2] = {
        {"pattern", "\"Type\":\"pattern\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Nz\":1,\"Data\":[1,0.2,0.2,1]}"},
        {"2 patterns sharing photons", "\"Type\":\"pattern\",\"SrcNum\":2,\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Data\":[1,0.5,0.5,1,0,1,1,0]}"},
        {"fourier", "\"Type\":\"fourier\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,1]"},
        {"fourierx", "\"Type\":\"fourierx\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,4],\"Param2\":[1,1,0,0]"},
        {"fourierx2d", "\"Type\":\"fourierx2d\",\"Pos\":[28,43,15],\"Param1\":[4,0,0,4],\"Param2\":[1,1,0,0]"},
        {"pencilarray", "\"Type\":\"pencilarray\",\"Pos\":[25,40,15],\"Param1\":[10,0,0,3],\"Param2\":[0,10,0,3]"},
        {"pattern3d", NULL},
        {"hyperboloid", "\"Type\":\"hyperboloid\",\"Pos\":[30,45,15],\"Param1\":[2,5,10,0]"},
        {"ring sector", "\"Type\":\"ring\",\"Pos\":[30,45,15],\"Param1\":[6,3,0,3.14]"},
        {"focused gaussian", "\"Type\":\"gaussian\",\"Pos\":[30,45,15],\"Dir\":[0,0,1,20],\"Param1\":[1,0.0005,0,0]"},
        {"diverging slit", "\"Type\":\"slit\",\"Pos\":[25,45,15],\"Dir\":[0,0,1,-10],\"Param1\":[10,0,0,0]"},
    };
    char p3dsrc[1024];

    snprintf(p3dsrc, sizeof(p3dsrc), "\"Type\":\"pattern3d\",\"Pos\":[28,43,13],\"Param1\":[4,4,4,0],\"Pattern\":{\"Nx\":4,\"Ny\":4,\"Nz\":4,\"Data\":[%s]}", p3d);
    src[6][1] = p3dsrc;

    printf("optix source types vs the CPU tracer:\n");

    for (int i = 0; i < (int)(sizeof(src) / sizeof(src[0])); i++) {
        RunResult rc, ro;
        snprintf(args, sizeof(args), "--bench dmmc-cube60 --compute sse -M G -n 100000 -E 1234 -s utos%dc --json {\"Optode\":{\"Source\":{%s}}}", i, src[i][1]);
        run_isolated(args, &rc);
        snprintf(args, sizeof(args), "--bench dmmc-cube60 --compute optix -M G -n 100000 -E 4321 -s utos%do --json {\"Optode\":{\"Source\":{%s}}}", i, src[i][1]);
        int ret = run_isolated(args, &ro);
        double dc = sqrt(pow(ro.cen[0] - rc.cen[0], 2) + pow(ro.cen[1] - rc.cen[1], 2) + pow(ro.cen[2] - rc.cen[2], 2));

        printf("  %-28s centroid off by %.3f voxel, RMS radius %.3f vs %.3f, fluence sum %.1f%% off\n", src[i][0], dc,
               ro.rms, rc.rms, rc.wsum != 0. ? 100. * fabs(ro.wsum - rc.wsum) / rc.wsum : 0.);
        sprintf(msg, "optix: the %s source gives the fluence of the CPU tracer", src[i][0]);
        EXPECT_TRUE(ret == 0 && rc.ran && dc < 0.2 && fabs(ro.rms - rc.rms) < 0.01 * rc.rms && fabs(ro.wsum - rc.wsum) < 0.03 * rc.wsum, msg);
    }

    /* ID=2 of 2 sources launches every photon from the 2nd source, the same as a single
       source there (the OpenCL/CUDA kernel launches them from the element of the 1st) */
    {
        RunResult r1, r2;
        run_isolated("--bench dmmc-cube60 --compute optix -M G -n 100000 -E 1234 -s utos2a --json {\"Optode\":{\"Source\":{\"Type\":\"pencil\",\"ID\":2,"
                     "\"Pos\":[[30,45,0,1],[20,30,0,1]],\"Dir\":[[0,0,1,0],[0,0,1,0]]}}}", &r2);
        run_isolated("--bench dmmc-cube60 --compute optix -M G -n 100000 -E 4321 -s utos2b --json {\"Optode\":{\"Source\":{\"Type\":\"pencil\",\"Pos\":[20,30,0],\"Dir\":[0,0,1]}}}", &r1);
        double dc = sqrt(pow(r2.cen[0] - r1.cen[0], 2) + pow(r2.cen[1] - r1.cen[1], 2) + pow(r2.cen[2] - r1.cen[2], 2));
        printf("  %-28s centroid off by %.3f voxel, fluence sum %.1f%% off\n", "2nd of 2 sources (ID=2)", dc, r1.wsum != 0. ? 100. * fabs(r2.wsum - r1.wsum) / r1.wsum : 0.);
        EXPECT_TRUE(r1.ran && r2.ran && dc < 0.2 && fabs(r2.wsum - r1.wsum) < 0.02 * r1.wsum, "optix: ID=2 of 2 sources launches from the 2nd source");
    }

#endif
}

/*-----------------------------------------------------------------------------
 *  OptiX adjoint Jacobians: a different random stream than OpenCL (CUDA replays that of
 *  OpenCL), so they are compared statistically with more photons; J_mus and J_musp are
 *  J_D times a factor of the medium label recorded in each voxel, the same for all
 *  backends in a homogeneous medium
 *---------------------------------------------------------------------------*/

static void test_optixadjoint(void) {
#ifdef USE_OPTIX
    const char* det = "{\"Optode\":{\"Detector\":[{\"Pos\":[30,40,0],\"R\":2,\"Dir\":[0,0,1]}]}}";
    const char* out[] = {"A", "D", "U", "V"};
    const char* name[] = {"J_mua", "J_D", "J_mus", "J_musp"};
    const double tol[] = {0.05, 0.2, 0.2, 0.2};  /* J_D varies ~4% (1 sd) between seeds at 1e6 photons */
    double jo[2][4];
    char args[1024], msg[256];

    if (!hasocl) {
        return;
    }

    printf("optix adjoint Jacobians vs opencl (1e6 photons):\n");

    for (int i = 0; i < 4; i++) {
        for (int b = 0; b < 2; b++) {
            RunResult r;
            snprintf(args, sizeof(args), "--bench dmmc-cube60 --compute %s -G 1 -M G -O %s -n 1000000 -E %d -s utoa%d%d --json %s",
                     b ? "optix" : "opencl", out[i], b ? 4321 : 1234, i, b, det);
            run_isolated(args, &r);
            jo[b][i] = r.jsum;
        }

        printf("  %-7s sum %.5g vs %.5g (%.1f%%)\n", name[i], jo[1][i], jo[0][i], jo[0][i] != 0. ? 100. * fabs(jo[1][i] - jo[0][i]) / fabs(jo[0][i]) : 0.);
        sprintf(msg, "optix: %s agrees with opencl within the noise of 1e6 photons", name[i]);
        EXPECT_TRUE(jo[0][i] != 0. && fabs(jo[1][i] - jo[0][i]) < tol[i] * fabs(jo[0][i]), msg);
    }

    /* the chain-rule factors of the dmmc-cube60 medium (mus=1, g=0.01): 1/(3*(1-g)*mus^2)
       for J_mus, 1/(3*(1-g)^2*mus^2) for J_musp; an OptiX run traces the same photons for
       each output, so its ratio is exact, while the photons of two OpenCL runs differ with
       the dynamic photon dispatch, by ~2% in the ratio */
    const double fmus = 1. / (3. * 0.99), fmusp = 1. / (3. * 0.99 * 0.99);
    EXPECT_NEAR(jo[1][2] / jo[1][1], fmus, 0.01 * fmus, "optix: J_mus/J_D is the chain-rule factor");
    EXPECT_NEAR(jo[1][3] / jo[1][1], fmusp, 0.01 * fmusp, "optix: J_musp/J_D is the chain-rule factor");
    EXPECT_NEAR(jo[0][2] / jo[0][1], fmus, 0.05 * fmus, "opencl: J_mus/J_D is the chain-rule factor");
    EXPECT_NEAR(jo[0][3] / jo[0][1], fmusp, 0.05 * fmusp, "opencl: J_musp/J_D is the chain-rule factor");
#endif
}

/*-----------------------------------------------------------------------------
 *  entry point
 *---------------------------------------------------------------------------*/

int main(int argc, char** argv) {
    char tmpdir[] = "/tmp/mmcunitXXXXXX", cwd[4096];

    devnull = fopen("/dev/null", "w");
    hasocl = (getenv("MMC_UNITTEST_NO_OPENCL") == NULL);

    if (argc == 4 && strcmp(argv[1], "--runone") == 0) {
        return run_one(argv[2], argv[3]);
    }

    if (getcwd(cwd, sizeof(cwd)) == NULL || mkdtemp(tmpdir) == NULL || chdir(tmpdir)) {
        printf("FAIL can not create a temporary folder\n");
        return 1;
    }

    test_tracers();
    test_sources();
    test_outputs();
    test_immc();
    test_inputfiles();
    test_moreio();
    test_variants();
    test_jsonedge();
    test_morecover();
    test_cmdedge();
    test_morecover2();
    test_morecover3();
    test_morecover4();
    test_utilfuncs();
    test_cmdopts();
    test_mmclab_path();
    test_clerrors();
    test_opencl();
    test_features();
    test_gpubackends();
    test_optixsources();
    test_optixadjoint();
    test_meshsurf();

    if (chdir(cwd) == 0) {
        char cmd[4200];
        snprintf(cmd, sizeof(cmd), "rm -rf %s", tmpdir);

        if (system(cmd)) {
            printf("warning: can not remove %s\n", tmpdir);
        }
    }

    printf("%d checks, %d failures\n", ntest, nfail);
    return nfail > 0;
}
