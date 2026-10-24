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

#include "mmc_host.h"
#include "mmc_utils.h"
#include "mmc_mesh.h"
#include "mmc_raytrace.h"
#include "zmat/zmatlib.h"
#include "ubj/ubj.h"

#ifdef USE_OPENCL
    #include "mmc_cl_host.h"
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

static void write_file(const char* fname, const void* buf, size_t len) {
    FILE* fp = fopen(fname, "wb");

    if (fp) {
        fwrite(buf, 1, len, fp);
        fclose(fp);
    }
}

/** results of one short MMC call */
typedef struct {
    double normalizer;
    int detected;
    int ran;
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

    if (cfg.isgpuinfo == 0) {
        mmc_prep(&cfg, &mesh, &tracer);

        if (cfg.compute == cbSSE || cfg.gpuid > MAX_DEVICE) {
            mmc_run_mp(&cfg, &mesh, &tracer);
        }

#ifdef USE_OPENCL
        else {
            mmc_run_cl(&cfg, &mesh, &tracer);
        }

#endif
        res.ran = 1;
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

static void test_immc(void) {
    const char* bench[] = {"edgeimmc", "nodeimmc", "faceimmc"};
    char args[512], msg[128];

    for (int i = 0; i < 3; i++) {
        snprintf(args, sizeof(args), "--bench %s -G -1 -n 500 -E 1 -b 1 -s utimmc%d", bench[i], i);
        RunResult r = run_mmc(args, NULL);
        sprintf(msg, "iMMC %s with reflection runs", bench[i]);
        EXPECT_TRUE(r.ran, msg);
        snprintf(args, sizeof(args), "--bench %s -G -1 -n 500 -E 1 -b 0 -m 1 -s utimmcb%d", bench[i], i);
        r = run_mmc(args, NULL);
        sprintf(msg, "iMMC %s without reflection, MCML-style", bench[i]);
        EXPECT_TRUE(r.ran, msg);
    }
}

/*-----------------------------------------------------------------------------
 *  OpenCL runs: tracer variants, RF and the adjoint Jacobians
 *---------------------------------------------------------------------------*/

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
                           "\"Pos\":[28,43,15],\"Param1\":[4,0,0,2],\"Param2\":[0,4,0,2],\"Pattern\":{\"Nx\":2,\"Ny\":2,\"Nz\":2,"
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
}

/*-----------------------------------------------------------------------------
 *  entry point
 *---------------------------------------------------------------------------*/

int main(int argc, char** argv) {
    char tmpdir[] = "/tmp/mmcunitXXXXXX", cwd[4096];

    devnull = fopen("/dev/null", "w");
    hasocl = (getenv("MMC_UNITTEST_NO_OPENCL") == NULL);

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
    test_utilfuncs();
    test_cmdopts();
    test_mmclab_path();
    test_clerrors();
    test_opencl();

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
