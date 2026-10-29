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
\file    pmmc.cpp

@brief   Python interface using Pybind11 for MMC
*******************************************************************************/
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <vector>
#include <cstdarg>
#include <streambuf>
#include <ostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <pybind11/iostream.h>

#include "mmc_const.h"
#include "mmc_mesh.h"
#include "mmc_host.h"
#ifdef USE_OPENCL
    #include "mmc_cl_host.h"
#endif
#ifdef USE_CUDA
    #include "mmc_cu_host.h"
    int mcx_list_cu_gpu(mcconfig* cfg, GPUInfo** info); /* defined in mmc_cu_host.cu; mmc_cuda_query_gpu.h clashes with the vector types */
#endif
#include "mmc_tictoc.h"
#include "mmc_raytrace.h"
#include "mmc_highorder.h"

// Python binding for runtime_error exception in Python.
namespace pybind11 {
PYBIND11_RUNTIME_EXCEPTION(runtime_error, PyExc_RuntimeError);
}

namespace py = pybind11;

float* det_ps = nullptr;     //! buffer to receive data from cfg.detphotons field
int dim_det_ps[2] = {0, 0};  //! dimensions of the cfg.detphotons array
int seed_byte = 0;

/**
 * Macro to find and extract a scalar property from a source Python dictionary configuration and assign it in a destination
 * MCX mcconfig. The scalar is cast to the python type before assignment.
 */
#define GET_SCALAR_FIELD(src_pydict, dst_mcx_config, property, py_type) if ((src_pydict).contains(#property))\
    {try {(dst_mcx_config).property = py_type((src_pydict)[#property]);\
            pyecho() << #property << ": " << (float) (dst_mcx_config).property << std::endl;} \
        catch (const std::runtime_error &err)\
        {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#property) + ". Reason: " + err.what()));}\
    }

#define GET_VEC3_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>()};\
            pyecho() << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << "]\n";} \
        catch (const std::runtime_error &err ) {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#prop) + ". Reason: " + err.what()));}}

#define GET_VEC4_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>(), list[3].cast<type>()}; \
            pyecho() << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << ", " << dst.prop.w << "]\n";} \
        catch (const std::runtime_error &err ) {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#prop) + ". Reason: " + err.what()));}}

#define GET_VEC34_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>(), list.size() == 4 ? list[3].cast<type>() : dst.prop.w}; \
            pyecho() << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << ", " << dst.prop.w << "]\n";} \
        catch (const std::runtime_error &err ) {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#prop) + ". Reason: " + err.what()));}                                                                 \
    }

#define MEXERROR(a)  mcx_error(999,a,__FILE__,__LINE__)   //! Macro to add unit name and line number in error printing

extern const char debugflag[];

/**
 * @brief Force matlab refresh the command window to print all buffered messages
 */

extern "C" void mcx_python_flush() {
    std::cout.flush();
}

/** prints the log of cfg.flog=2 (stderr, the default) to Python's sys.stdout, as mexPrintf does in
    MATLAB: the log then shows in Jupyter, in order with Python's own output. Called by the thread
    that holds the GIL (OpenMP thread 0, see MCX_FPRINTF); without a Python stdout it uses stderr */
extern "C" int mcx_python_printf(const char* format, ...) {
    va_list args, args2;
    va_start(args, format);
    va_copy(args2, args);
    int len = vsnprintf(NULL, 0, format, args);
    va_end(args);

    if (len < 0) {
        va_end(args2);
        return len;
    }

    std::vector<char> buf(len + 1);
    vsnprintf(buf.data(), buf.size(), format, args2);
    va_end(args2);

    try {
        py::gil_scoped_acquire gil;
        py::object out = py::module_::import("sys").attr("stdout");

        if (out.is_none()) {
            throw std::runtime_error("no sys.stdout");
        }

        out.attr("write")(py::str(buf.data(), len));
        out.attr("flush")();
    } catch (...) {
        fputs(buf.data(), stderr);
    }

    return len;
}


/**
 * The echo of the cfg fields (and the warnings printed while reading them) goes where cfg.flog sends
 * the simulation log: Python's sys.stdout (std::cout, redirected by pybind11) for 2 (default), or
 * the stdout of the process for 1 and 0, the log file, or nowhere for os.devnull
 */
static FILE* pyecho_logfp = NULL;   /**< NULL or stderr: std::cout, i.e. sys.stdout */

class PyEchoFileBuf : public std::streambuf {
  public:
    FILE* fp = NULL;
  protected:
    int overflow(int c) override {
        return (c == EOF || fputc(c, fp) == EOF) ? EOF : c;
    }
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        return (std::streamsize)fwrite(s, 1, (size_t)n, fp);
    }
};

static std::ostream& pyecho() {
    static PyEchoFileBuf filebuf;
    static std::ostream fileout(&filebuf);

    if (pyecho_logfp == NULL || pyecho_logfp == stderr) {
        return std::cout;
    }

    filebuf.fp = pyecho_logfp;
    return fileout;
}

/**
 * Parse user input cfg object and convert to to MMC's mcconfig object.
 * @param user_cfg
 * @param mcx_config
 */


void parse_config(const py::dict& user_cfg, mcconfig& mcx_config, tetmesh& mesh) {
    int patndim = 0, patdim0 = 1;  /* the shape of srcpattern, read before srctype may be known */

    mcx_initcfg(&mcx_config);
    mesh_init(&mesh);

    mcx_config.flog = stderr;

    /* cfg.flog first: the echo of the other fields goes where it sends the log */
    if (user_cfg.contains("flog")) {
        auto logfile_id_value = user_cfg["flog"];

        if (py::int_::check_(logfile_id_value)) {
            int logid = py::int_(logfile_id_value);
            mcx_config.flog = (logid >= 2 ? stderr : (logid == 1 ? stdout : (mcx_config.printnum = -1, stdout)));
        } else if (py::str::check_(logfile_id_value)) {
            std::string logfile_id_string_value = py::str(logfile_id_value);

            if (logfile_id_string_value.empty()) {
                throw py::value_error("the 'flog' field must be an integer or non-empty string");
            }

            mcx_config.flog = fopen(logfile_id_string_value.c_str(), "a+");

            if (mcx_config.flog == NULL) {
                throw py::value_error("Log output file can not be written");
            }
        }
    }

    pyecho_logfp = mcx_config.flog;

    GET_SCALAR_FIELD(user_cfg, mcx_config, nphoton, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, nblocksize, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, nthread, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tstart, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tstep, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tend, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isreflect, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isspecular, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, ismomentum, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issaveexit, py::int_); /* 2 saves a detector image */
    GET_SCALAR_FIELD(user_cfg, mcx_config, isextdet, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issave2pt, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issavedet, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issaveseed, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, optlevel, py::int_);

    if (user_cfg.contains("isatomic") && !py::bool_(user_cfg["isatomic"])) {
        pyecho() << "WARNING: isatomic=0 is deprecated and ignored; non-atomic accumulation raced and lost part of the deposited energy, at no speed benefit" << std::endl;
    }

    GET_SCALAR_FIELD(user_cfg, mcx_config, basisorder, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, roulettesize, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, nout, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isref3, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isnormalized, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issaveref, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, debugphoton, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, minenergy, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, replaydet, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, unitinmm, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, voidtime, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, mcmethod, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, maxdetphoton, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, maxjumpdebug, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, e0, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, srcid, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, adjointmode, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isnodalmua, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isnodalmusp, py::int_);
    GET_VEC3_FIELD(user_cfg, mcx_config, steps, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, srcparam1, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, srcparam2, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, detparam1, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, detparam2, float);

    /* srcpos / srcdir can be either a single 3- or 4-vector (single source)
     * or an (Nsrc, 3) / (Nsrc, 4) matrix (multi-source). The matrix form
     * mirrors mmclab.cpp:783-821: row 0 fills cfg.srcpos/srcdir (the
     * "main" slot the kernel falls back to in single-source mode); all
     * rows are then copied into cfg.srcdata[0..Nsrc-1] and cfg.extrasrclen
     * is set to Nsrc, so the multi-source kernel dispatch picks them up.
     *
     * srcdir's optional 4th column (focal length / RF NaN/Inf markers) is
     * read whenever the row has 4 entries; the floats parse cleanly even
     * for NaN/Inf if the caller passes math.nan / math.inf in the numpy
     * array (no string-literal handling like mmclab.cpp because pybind11
     * arrays are typed). */
    {
        auto parse_multisrc_matrix = [&](const char* name, float4 & dst_single,
        float4 ExtraSrc::* dst_field, bool is_srcdir) {
            if (!user_cfg.contains(name)) {
                return;
            }

            auto arr = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg[name]);

            if (!arr) {
                throw py::value_error(std::string("Invalid ") + name + " field value");
            }

            auto info = arr.request();

            int nrows = 1, ncols = 0;

            if (info.ndim == 1) {
                ncols = (int)info.shape[0];
            } else if (info.ndim == 2) {
                nrows = (int)info.shape[0];
                ncols = (int)info.shape[1];
            } else {
                throw py::value_error(std::string(name) + " must be a 1D or 2D array");
            }

            if (ncols < 3 || ncols > 4) {
                throw py::value_error(std::string(name) + " must have 3 or 4 columns");
            }

            float* val = (float*)info.ptr;

            /* numpy passes f_style: linear index = col * nrows + row */
            for (int c = 0; c < ncols; c++) {
                ((float*)(&dst_single))[c] = val[c * nrows + 0];
            }

            pyecho() << name << ": [" << dst_single.x << ", " << dst_single.y << ", "
                     << dst_single.z << ", " << dst_single.w << "]" << std::endl;

            if (nrows == 1 && mcx_config.extrasrclen == 0) {
                return;
            }

            /* multi-source: populate srcdata[].srcpos or srcdata[].srcdir */
            if (mcx_config.extrasrclen != 0 && mcx_config.extrasrclen != nrows) {
                throw py::value_error("Length of sub-elements of srcpos/srcdir/srcparam1/srcparam2 must match");
            }

            mcx_config.extrasrclen = nrows;

            if (mcx_config.srcdata == NULL) {
                mcx_config.srcdata = (ExtraSrc*)calloc(sizeof(ExtraSrc), mcx_config.extrasrclen);
            }

            for (int c = 0; c < ncols; c++)
                for (int r = 0; r < nrows; r++) {
                    ((float*)(&(mcx_config.srcdata[r].*dst_field)))[c] = val[c * nrows + r];
                }

            pyecho() << name << " multi-source rows=" << nrows << std::endl;
        };

        parse_multisrc_matrix("srcpos", mcx_config.srcpos, &ExtraSrc::srcpos, false);
        parse_multisrc_matrix("srcdir", mcx_config.srcdir, &ExtraSrc::srcdir, true);
    }

    /* without node or elem (or the face of an OptiX surface mesh), the mesh setup dereferenced NULL */
    if (!user_cfg.contains("node") || user_cfg["node"].is_none()) {
        throw py::value_error("the 'node' field (the mesh nodes) is missing");
    }

    if ((!user_cfg.contains("elem") || user_cfg["elem"].is_none()) && !user_cfg.contains("face")) {
        throw py::value_error("the 'elem' field (the tetrahedra) is missing; an OptiX surface-only mesh gives 'face' instead");
    }

    if (user_cfg.contains("node")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["node"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid node field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && (buffer_info.shape.at(0) < 4 || buffer_info.shape.at(1) != 3)) || (buffer_info.shape.size() == 1)) {
            throw py::value_error("the 'node' field must have 3 columns (x,y,z) and minimum 4 nodes");
        }

        mesh.nn = buffer_info.shape.at(0);

        if (mesh.node) {
            free(mesh.node);
        }

        mesh.node = (FLOAT3*) malloc(mesh.nn * sizeof(FLOAT3));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 3; j++)
            for (int i = 0; i < mesh.nn; i++) {
                ((float*) (&mesh.node[i]))[j] = val[j * mesh.nn + i];
            }
    }


    if (user_cfg.contains("elem")) {
        auto f_style_volume = py::array_t < int, py::array::f_style | py::array::forcecast >::ensure(user_cfg["elem"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid elem field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && (buffer_info.shape.at(0) == 0 || (buffer_info.shape.at(1) != 4 && buffer_info.shape.at(1) != 10))) || (buffer_info.shape.size() == 1)) {
            throw py::value_error("the 'elem' field must have 4 or 10 columns");
        }

        mesh.ne = buffer_info.shape.at(0);
        mesh.elemlen = buffer_info.shape.at(1);

        if (mesh.elem) {
            free(mesh.elem);
        }

        mesh.elem = (int*) malloc(mesh.ne * mesh.elemlen * sizeof(int));
        auto val = static_cast<int*>(buffer_info.ptr);

        for (int j = 0; j < mesh.elemlen; j++)
            for (int i = 0; i < mesh.ne; i++) {
                mesh.elem[i * mesh.elemlen + j] = val[j * mesh.ne + i];
            }
    }


    if (user_cfg.contains("noderoi")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["noderoi"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid noderoi field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && (buffer_info.shape.at(0) != 1 && buffer_info.shape.at(1) != 1))) {
            throw py::value_error("the 'noderoi' field must have 1 row or 1 column");
        }

        mesh.nn = (buffer_info.shape.size() == 1) ? buffer_info.shape.at(0) : buffer_info.shape.at(0) * buffer_info.shape.at(1);

        if (mesh.noderoi) {
            free(mesh.noderoi);
        }

        mesh.noderoi = (float*) malloc(mesh.nn * sizeof(float));
        auto val = static_cast<float*>(buffer_info.ptr);
        memcpy(mesh.noderoi, val, mesh.nn * sizeof(float));
        mcx_config.implicit = 1; /* enable implicit MMC (iMMC), as mmclab */
    }


    if (user_cfg.contains("edgeroi")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["edgeroi"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid edgeroi field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(1) != 6) || (buffer_info.shape.size() == 1)) {
            throw py::value_error("the 'edgeroi' field must have 6 columns");
        }

        mesh.ne = buffer_info.shape.at(0);

        if (mesh.edgeroi) {
            free(mesh.edgeroi);
        }

        mesh.edgeroi = (float*) malloc(mesh.ne * 6 * sizeof(float));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 6; j++)
            for (int i = 0; i < mesh.ne; i++) {
                mesh.edgeroi[i * 6 + j] = val[j * mesh.ne + i];
            }

        mcx_config.implicit = 1; /* edge-based iMMC */
    }


    if (user_cfg.contains("faceroi")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["faceroi"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid faceroi field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(1) != 4) || (buffer_info.shape.size() == 1)) {
            throw py::value_error("the 'faceroi' field must have 6 columns");
        }

        mesh.ne = buffer_info.shape.at(0);

        if (mesh.faceroi) {
            free(mesh.faceroi);
        }

        mesh.faceroi = (float*) malloc(mesh.ne * 4 * sizeof(float));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 4; j++)
            for (int i = 0; i < mesh.ne; i++) {
                mesh.faceroi[i * 4 + j] = val[j * mesh.ne + i];
            }

        mcx_config.implicit = 2; /* face-based iMMC */
    }


    if (user_cfg.contains("face") && !(user_cfg.contains("elem") && !user_cfg["elem"].is_none())) {
        /* a surface-only mesh for OptiX: 3 node indices (1-based), the inner and the outer labels,
           given instead of cfg["elem"]; with "elem", "face" is not mmc's (it is often the boundary
           of the tetrahedral mesh) and is ignored, as it was before surface-only meshes */
        auto f_style_volume = py::array_t < int, py::array::f_style | py::array::forcecast >::ensure(user_cfg["face"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid face field value");
        }

        auto buffer_info = f_style_volume.request();

        if (buffer_info.shape.size() != 2 || buffer_info.shape.at(0) == 0 || buffer_info.shape.at(1) != 5) {
            throw py::value_error("the 'face' field (a surface-only mesh for OptiX) must have 5 columns: 3 node indices, the inner and the outer labels");
        }

        mcx_config.surfnum = buffer_info.shape.at(0);

        if (mcx_config.surf) {
            free(mcx_config.surf);
        }

        mcx_config.surf = (int*) malloc(mcx_config.surfnum * 5 * sizeof(int));
        auto val = static_cast<int*>(buffer_info.ptr);

        for (int j = 0; j < 5; j++)
            for (int i = 0; i < (int)mcx_config.surfnum; i++) {
                mcx_config.surf[i * 5 + j] = val[j * mcx_config.surfnum + i];
            }
    }

    if (user_cfg.contains("facenb")) {
        auto f_style_volume = py::array_t < int, py::array::f_style | py::array::forcecast >::ensure(user_cfg["facenb"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid facenb field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && (buffer_info.shape.at(0) == 0 || (buffer_info.shape.at(1) != 4 && buffer_info.shape.at(1) != 10))) || (buffer_info.shape.size() == 1)) {
            throw py::value_error("the 'facenb' field must have 4 or 10 columns");
        }

        mesh.ne = buffer_info.shape.at(0);
        mesh.elemlen = buffer_info.shape.at(1);

        if (mesh.facenb) {
            free(mesh.facenb);
        }

        mesh.facenb = (int*) malloc(mesh.ne * mesh.elemlen * sizeof(int));
        auto val = static_cast<int*>(buffer_info.ptr);

        for (int j = 0; j < mesh.elemlen; j++)
            for (int i = 0; i < mesh.ne; i++) {
                mesh.facenb[i * mesh.elemlen + j] = val[j * mesh.ne + i];
            }
    }


    if (user_cfg.contains("elemprop")) {
        auto f_style_volume = py::array_t < int, py::array::f_style | py::array::forcecast >::ensure(user_cfg["elemprop"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid elemprop field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && (buffer_info.shape.at(0) != 1 && buffer_info.shape.at(1) != 1))) {
            throw py::value_error("the 'elemprop' field must have 1 row or 1 column");
        }

        mesh.ne = (buffer_info.shape.size() == 1) ? buffer_info.shape.at(0) : buffer_info.shape.at(0) * buffer_info.shape.at(1);

        if (mesh.type) {
            free(mesh.type);
        }

        mesh.type = (int*) malloc(mesh.ne * sizeof(int));
        auto val = static_cast<int*>(buffer_info.ptr);
        memcpy(mesh.type, val, mesh.ne * sizeof(int));
    }


    if (user_cfg.contains("detpos")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["detpos"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid detpos field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(0) > 0 && buffer_info.shape.at(1) != 4) || (buffer_info.shape.size() == 1 && buffer_info.shape.at(0) != 4)) {
            throw py::value_error("the 'detpos' field must have 4 columns (x,y,z,radius)");
        }

        mcx_config.detnum = (buffer_info.shape.size() == 1) ? 1 : buffer_info.shape.at(0);

        if (mcx_config.detpos) {
            free(mcx_config.detpos);
        }

        mcx_config.detpos = (float4*) malloc(mcx_config.detnum * sizeof(float4));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 4; j++)
            for (int i = 0; i < mcx_config.detnum; i++) {
                ((float*) (&mcx_config.detpos[i]))[j] = val[j * mcx_config.detnum + i];
            }
    }


    if (user_cfg.contains("prop")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["prop"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid prop field format");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(0) > 0 && buffer_info.shape.at(1) != 4) || (buffer_info.shape.size() == 1 && buffer_info.shape.at(0) != 4)) {
            throw py::value_error("the 'prop' field must have 4 columns (mua,mus,g,n)");
        }

        mesh.prop = (buffer_info.shape.size() == 1) ? 0 : buffer_info.shape.at(0) - 1;

        if (mesh.med) {
            free(mesh.med);
        }

        mesh.med = (medium*) malloc((mesh.prop + 1) * sizeof(medium));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 4; j++)
            for (int i = 0; i <= mesh.prop; i++) {
                ((float*) (&mesh.med[i]))[j] = val[j * (mesh.prop + 1) + i];
            }

        mcx_config.his.maxmedia = mesh.prop;
    }



    if (user_cfg.contains("compileropt")) {
        std::string compileropt = py::str(user_cfg["compileropt"]);

        if (compileropt.empty()) {
            throw py::value_error("the 'compileropt' field must be a non-empty string");
        }

        if (compileropt.size() > MAX_PATH_LENGTH) {
            throw py::value_error("the 'compileropt' field is too long");
        }

        strncpy(mcx_config.compileropt, compileropt.c_str(), MAX_PATH_LENGTH);
    }


    if (user_cfg.contains("kernelfile")) {
        std::string kernelfile = py::str(user_cfg["kernelfile"]);

        if (kernelfile.empty()) {
            throw py::value_error("the 'kernelfile' field must be a non-empty string");
        }

        if (kernelfile.size() >= MAX_PATH_LENGTH) {
            throw py::value_error("the 'kernelfile' field is too long");
        }

        strncpy(mcx_config.kernelfile, kernelfile.c_str(), MAX_PATH_LENGTH - 1);
    }


    if (user_cfg.contains("session")) {
        std::string session = py::str(user_cfg["session"]);

        if (session.empty()) {
            throw py::value_error("the 'session' field must be a non-empty string");
        }

        if (session.size() > MAX_SESSION_LENGTH) {
            throw py::value_error("the 'session' field is too long");
        }

        strncpy(mcx_config.session, session.c_str(), MAX_SESSION_LENGTH);
    }


    if (user_cfg.contains("srctype")) {
        std::string src_type = py::str(user_cfg["srctype"]);
        const char* srctypeid[] = {"pencil", "isotropic", "cone", "gaussian", "planar",
                                   "pattern", "fourier", "arcsine", "disk", "fourierx", "fourierx2d", "zgaussian",
                                   "line", "slit", "pencilarray", "pattern3d", "hyperboloid", "ring", ""
                                  };

        if (src_type.empty()) {
            throw py::value_error("the 'srctype' field must be a non-empty string");
        }

        mcx_config.srctype = mcx_keylookup((char*)(src_type.c_str()), srctypeid);

        if (mcx_config.srctype == -1) {
            throw py::value_error("the specified source type is not supported");
        }
    }


    if (user_cfg.contains("method")) {
        std::string method_str = py::str(user_cfg["method"]);
        const char* methods[] = {"plucker", "havel", "badouel", "elem", "grid", ""};

        if (method_str.empty()) {
            throw py::value_error("the 'method' field must be a non-empty string");
        }

        mcx_config.method = mcx_keylookup((char*)(method_str.c_str()), methods);

        if (mcx_config.method == -1) {
            throw py::value_error("the specified ray-tracing method is not supported");
        }
    }


    if (user_cfg.contains("outputtype")) {
        std::string output_type_str = py::str(user_cfg["outputtype"]);
        const char* outputtype[] = {"flux", "fluence", "energy", "jacobian", "wl", "wp",
                                    "rf", "rfmus", "adjoint", "adjoint_dcoeff", "adjoint_mus", "adjoint_musp", "adjoint_mua_d", "adjoint_mua_musp", ""
                                   };

        if (output_type_str.empty()) {
            throw py::value_error("the 'outputtype' field must be a non-empty string");
        }

        mcx_config.outputtype = mcx_keylookup((char*)(output_type_str.c_str()), outputtype);

        if (mcx_config.outputtype == -1) {
            throw py::value_error("the specified output type is not supported");
        }
    }

    if (user_cfg.contains("omega")) {
        mcx_config.omega = py::float_(user_cfg["omega"]);
        pyecho() << "omega: " << mcx_config.omega << std::endl;
    }

    if (user_cfg.contains("detdir")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["detdir"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid detdir field value");
        }

        auto buffer_info = f_style_volume.request();

        if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(0) > 0 && buffer_info.shape.at(1) != 4) || (buffer_info.shape.size() == 1 && buffer_info.shape.at(0) != 4)) {
            throw py::value_error("the 'detdir' field must have 4 columns (dx,dy,dz,focal)");
        }

        int nd = (buffer_info.shape.size() == 1) ? 1 : buffer_info.shape.at(0);

        if (mcx_config.detdir) {
            free(mcx_config.detdir);
        }

        mcx_config.detdir = (float4*) malloc(nd * sizeof(float4));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int j = 0; j < 4; j++)
            for (int i = 0; i < nd; i++) {
                ((float*) (&mcx_config.detdir[i]))[j] = val[j * nd + i];
            }

        pyecho() << "detdir: [" << nd << ",4]" << std::endl;
    }

    if (user_cfg.contains("sdmap")) {
        auto f_style_volume = py::array_t < double, py::array::f_style | py::array::forcecast >::ensure(user_cfg["sdmap"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid sdmap field value");
        }

        auto buffer_info = f_style_volume.request();

        if (mcx_config.sdmap) {
            free(mcx_config.sdmap);
            mcx_config.sdmap = nullptr;
        }

        mcx_config.sdmaplen = 0;

        if (buffer_info.size > 0) {
            if ((buffer_info.shape.size() > 1 && buffer_info.shape.at(1) < 2) || (buffer_info.shape.size() == 1 && buffer_info.shape.at(0) < 2) || buffer_info.shape.size() > 2) {
                throw py::value_error("the 'sdmap' field must have 2 or more columns (srcid,detid,...)");
            }

            unsigned int nch = (buffer_info.shape.size() == 1) ? 1 : buffer_info.shape.at(0);
            auto val = static_cast<double*>(buffer_info.ptr);

            mcx_config.sdmap = (unsigned int*) malloc(sizeof(unsigned int) * 2 * nch);
            mcx_config.sdmaplen = nch;

            for (unsigned int i = 0; i < nch; i++) {
                if (val[i] < 1 || val[nch + i] < 1) {
                    throw py::value_error("the source and detector indices in 'sdmap' start from 1");
                }

                mcx_config.sdmap[i << 1] = (unsigned int)val[i] - 1;
                mcx_config.sdmap[(i << 1) + 1] = (unsigned int)val[nch + i] - 1;
            }
        }
    }

    /* per-node optical properties: cfg.nodemua / cfg.nodemusp are length-nn float arrays.
     * Either array's presence auto-sets the corresponding isnodal* flag. */
    if (user_cfg.contains("nodemua")) {
        auto arr = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["nodemua"]);

        if (!arr) {
            throw py::value_error("Invalid nodemua field value");
        }

        auto buffer_info = arr.request();
        size_t nn = (size_t)buffer_info.size;

        if (mcx_config.nodemua) {
            free(mcx_config.nodemua);
        }

        mcx_config.nodemua = (float*) malloc(nn * sizeof(float));
        memcpy(mcx_config.nodemua, buffer_info.ptr, nn * sizeof(float));
        mcx_config.isnodalmua = 1;
        pyecho() << "nodemua: [" << nn << "]" << std::endl;
    }

    if (user_cfg.contains("nodemusp")) {
        auto arr = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["nodemusp"]);

        if (!arr) {
            throw py::value_error("Invalid nodemusp field value");
        }

        auto buffer_info = arr.request();
        size_t nn = (size_t)buffer_info.size;

        if (mcx_config.nodemusp) {
            free(mcx_config.nodemusp);
        }

        mcx_config.nodemusp = (float*) malloc(nn * sizeof(float));
        memcpy(mcx_config.nodemusp, buffer_info.ptr, nn * sizeof(float));
        mcx_config.isnodalmusp = 1;
        pyecho() << "nodemusp: [" << nn << "]" << std::endl;
    }


    if (user_cfg.contains("compute")) {
        std::string compute_str = py::str(user_cfg["compute"]);
        const char* computebackend[] = {"sse", "opencl", "cuda", "optix", ""};  /* as mmclab; without OptiX, "optix" runs OpenCL */

        if (compute_str.empty()) {
            throw py::value_error("the 'compute' field must be a non-empty string");
        }

        mcx_config.compute = mcx_keylookup((char*)(compute_str.c_str()), computebackend);

        if (mcx_config.compute == -1) {
            throw py::value_error("the specified compute backend is not supported");
        }
    }


    if (user_cfg.contains("debuglevel")) {
        std::string debug_level = py::str(user_cfg["debuglevel"]);
        char debuglevel[MAX_SESSION_LENGTH] = {'\0'};

        if (debug_level.empty()) {
            throw py::value_error("the 'debuglevel' field must be a non-empty string");
        }

        if (debug_level.size() > MAX_SESSION_LENGTH) {
            throw py::value_error("the 'debuglevel' field is too long");
        }

        strncpy(debuglevel, debug_level.c_str(), MAX_SESSION_LENGTH);
        /* the same flags as mmc -D and mmclab, e.g. 'T' timing, 'P' progress, 'M' trajectories (dlTraj) */
        mcx_config.debuglevel = mcx_parsedebugopt(debuglevel, debugflag);

        if (mcx_config.debuglevel == 0) {
            throw py::value_error("the specified debuglevel is not supported");
        }
    }


    if (user_cfg.contains("srcpattern")) {
        auto f_style_volume = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["srcpattern"]);

        if (!f_style_volume) {
            throw py::value_error("Invalid srcpattern field value");
        }

        auto buffer_info = f_style_volume.request();

        if (mcx_config.srcpattern) {
            free(mcx_config.srcpattern);
        }

        patndim = buffer_info.ndim;
        patdim0 = buffer_info.shape[0];

        mcx_config.srcpattern = (float*) malloc(buffer_info.size * sizeof(float));
        auto val = static_cast<float*>(buffer_info.ptr);

        for (int i = 0; i < buffer_info.size; i++) {
            mcx_config.srcpattern[i] = val[i];
        }
    }


    if (user_cfg.contains("shapes")) {
        std::string shapes_string = py::str(user_cfg["shapes"]);

        if (shapes_string.empty()) {
            throw py::value_error("the 'shapes' field must be a non-empty string");
        }

        mcx_config.shapedata = (char*) calloc(shapes_string.size() + 2, 1);
        strncpy(mcx_config.shapedata, shapes_string.c_str(), shapes_string.size() + 1);
    }

    // mcx uses detphotons raw-data, mmc uses replayweight/replaytime after post-processing, may adopt mcx's way
    if (user_cfg.contains("detphotons")) {
        auto detphotons = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["detphotons"]);

        if (!detphotons) {
            throw py::value_error("Invalid detphotons field value");
        }

        auto buffer_info = detphotons.request();

        det_ps = static_cast<float*>(buffer_info.ptr);
        dim_det_ps[0] = buffer_info.shape.at(0);
        dim_det_ps[1] = buffer_info.shape.at(1);
    }


    if (user_cfg.contains("seed")) {
        auto seed_value = user_cfg["seed"];

        // If the seed value is scalar (int or float), then assign it directly
        if (py::int_::check_(seed_value)) {
            mcx_config.seed = py::int_(seed_value);
        } else if (py::float_::check_(seed_value)) {
            mcx_config.seed = py::float_(seed_value).cast<int>();
        }
        // Set seed from array
        else {
            auto f_style_array = py::array_t < uint8_t, py::array::f_style | py::array::forcecast >::ensure(seed_value);

            if (!f_style_array) {
                throw py::value_error("Invalid seed field value");
            }

            auto buffer_info = f_style_array.request();
            seed_byte = buffer_info.shape.at(0);

            if (buffer_info.shape.at(0) != sizeof(RandType) * RAND_BUF_LEN) {
                throw py::value_error("the row number of cfg.seed does not match RNG seed byte-length");
            }

            mcx_config.photonseed = malloc(buffer_info.size);
            memcpy(mcx_config.photonseed, buffer_info.ptr, buffer_info.size);
            mcx_config.seed = SEED_FROM_FILE;
            mcx_config.nphoton = buffer_info.shape.at(1);
        }
    }


    if (user_cfg.contains("gpuid")) {
        auto gpu_id_value = user_cfg["gpuid"];

        if (py::int_::check_(gpu_id_value)) {
            int gpuid = py::int_(gpu_id_value);
            mcx_config.gpuid = gpuid;
            memset(mcx_config.deviceid, 0, MAX_DEVICE);

            if (gpuid >= MAX_DEVICE) {
                throw py::value_error("GPU id can not be more than 256");
            } else if (gpuid > 0) {
                memset(mcx_config.deviceid, '0', mcx_config.gpuid - 1);
                mcx_config.deviceid[mcx_config.gpuid - 1] = '1';
            } /* gpuid<=0, e.g. -1 to run on the CPU, is passed on as in mmclab */
        } else if (py::str::check_(gpu_id_value)) {
            std::string gpu_id_string_value = py::str(gpu_id_value);

            if (gpu_id_string_value.empty()) {
                throw py::value_error("the 'gpuid' field must be an integer or non-empty string");
            }

            if (gpu_id_string_value.size() > MAX_DEVICE) {
                throw py::value_error("the 'gpuid' field is too long");
            }

            strncpy(mcx_config.deviceid, gpu_id_string_value.c_str(), MAX_DEVICE);
        }

        for (int i = 0; i < MAX_DEVICE; i++)
            if (mcx_config.deviceid[i] == '0') {
                mcx_config.deviceid[i] = '\0';
            }
    }


    if (user_cfg.contains("workload")) {
        auto workload_value = py::array_t < float, py::array::f_style | py::array::forcecast >::ensure(user_cfg["workload"]);

        if (!workload_value) {
            throw py::value_error("Invalid workload field value");
        }

        auto buffer_info = workload_value.request();

        if (buffer_info.shape.size() < 2 && buffer_info.size > MAX_DEVICE) {
            throw py::value_error("the workload list can not be longer than 256");
        }

        for (int i = 0; i < buffer_info.size; i++) {
            mcx_config.workload[i] = static_cast<float*>(buffer_info.ptr)[i];
        }
    }


    /* a stack of patterns for photon sharing is [srcnum, Nx, Ny] for 'pattern' and
       [srcnum, Nx, Ny, Nz] for 'pattern3d', whose 3-D srcpattern is one volume; srcnum was
       taken from any 3-D srcpattern, which broke pattern3d */
    if (patndim == (mcx_config.srctype == stPattern3D ? 4 : 3)) {
        mcx_config.srcnum = patdim0;
    }

    /* save the detected photons when detectors are given, as the mmc command line and pmcx do;
       mmc's default issavedet=0 made cfg['detpos'] alone return no 'detp'. The adjoint outputs
       launch the detectors as sources instead */
    if (!user_cfg.contains("issavedet") && mcx_config.detnum > 0 && !MCX_IS_ADJOINT_TYPE(mcx_config.outputtype)) {
        mcx_config.issavedet = 1;
    }

    if (mesh.ne == 0 && mcx_config.surfnum > 0) {
        return;  /* a surface-only mesh for OptiX: no element data to derive */
    }

    if (mesh.type == NULL) {
        throw py::value_error("the 'elemprop' field is missing");  /* mesh_getvolume reads it */
    }

    /* an out-of-range (e.g. 0-based) elem aborted the process, and labels beyond the rows of
       prop silently gave a zero fluence */
    for (int i = 0; i < mesh.ne * mesh.elemlen; i++) {
        if (mesh.elem[i] < 1 || mesh.elem[i] > mesh.nn) {
            throw py::value_error("the 'elem' field must hold 1-based node indices, between 1 and the number of nodes");
        }
    }

    if (mesh.med != NULL) {
        for (int i = 0; i < mesh.ne; i++) {
            if (mesh.type[i] > mesh.prop) {
                throw py::value_error("the 'elemprop' field has labels beyond the media of 'prop' (its rows minus 1)");
            }
        }
    }

    if (mesh.evol == NULL) {
        mesh_getvolume(&mesh, &mcx_config);
    }

    /* a zero-volume element made the ray-tracer loop forever; mmclab rejects them too */
    if (mcx_config.compute != cbOptiX) {
        for (int i = 0; i < mesh.ne; i++) {
            if (!(mesh.evol[i] > 0.f)) {
                throw py::value_error("the mesh contains degenerated (zero-volume) elements, e.g. element " + std::to_string(i + 1));
            }
        }
    }

    if (mesh.facenb == NULL) {
        mesh_getfacenb(&mesh, &mcx_config);
    }

    if (!user_cfg.contains("e0")) {
        mesh_initelem(&mesh, &mcx_config);
    }

    // Flush the std::cout and std::cerr
    std::cout.flush();
    std::cerr.flush();
}

/**
 * Function that's called to cleanup any memory/configs allocated by PMMC. It is used in both normal and exceptional
 * termination of the application
 * @param gpu_info reference to an array of MCXGPUInfo data structure
 * @param mcx_config reference to mcconfig data structure
 */
inline void cleanup_configs(MCXGPUInfo*& gpu_info, mcconfig& mcx_config) {
    mcx_cleargpuinfo(&gpu_info);
    mcx_clearcfg(&mcx_config);
    pyecho_logfp = NULL;  /* mcx_clearcfg closed the log file */
}

inline void cleanup_mesh(mcconfig& mcx_config, tetmesh& mesh) {
    mesh_clear(&mesh, &mcx_config);
    mcx_clearcfg(&mcx_config);
    pyecho_logfp = NULL;  /* mcx_clearcfg closed the log file */
}


static int pmmc_task_guard(void (*task)(void*), void* arg);

py::dict pmmc_interface(const py::dict& user_cfg) {
    unsigned int hostdetreclen;
    mcconfig mcx_config;  /* mcx_config: structure to store all simulation parameters */
    tetmesh mesh;
    raytracer tracer = {NULL, 0, NULL, NULL, NULL};
    GPUInfo* gpu_info = nullptr;        /** gpuInfo: structure to store GPU information */
    unsigned int active_dev = 0;     /** activeDev: count of total active GPUs to be used */
    std::vector<std::string> exception_msgs;
    int thread_id = 0;
    size_t field_dim[5] = {0};
    py::dict output;

    try {
        /*
         * To start an MCX simulation, we first create a simulation configuration and set all elements to its default settings.
         */
        det_ps = nullptr;

        parse_config(user_cfg, mcx_config, mesh);

        if (mcx_config.compute == cbCUDA) {
#ifdef USE_CUDA
            char deviceid0[MAX_DEVICE];  /* mmc_run_cu lists the devices again, keep the unconverted mask */
            memcpy(deviceid0, mcx_config.deviceid, MAX_DEVICE);
            active_dev = mcx_list_cu_gpu(&mcx_config, &gpu_info);
            memcpy(mcx_config.deviceid, deviceid0, MAX_DEVICE);
#else
            throw py::value_error("compute='cuda' needs a pmmc built with CUDA (CMAKE_ARGS=\"-DBUILD_CUDA=ON\"); use 'opencl' or 'sse'");
#endif
        } else {
#ifdef USE_OPENCL
            mcx_list_cl_gpu(&mcx_config, &active_dev, NULL, &gpu_info);
#endif
        }

        /** The next step, we identify gpu number and query all GPU info */
        if (!active_dev) {
            mcx_error(-1, "No GPU device found\n", __FILE__, __LINE__);
        }

        mcx_python_flush();

        /** Find the source (-1) and wide-field detector (-2) elements first, as mmclab does: it sets
         *  cfg.isextdet, which mmc_validate_config (keeping issavedet) and mesh_validate (the medium
         *  of the detector elements) need */
        mesh_srcdetelem(&mesh, &mcx_config);

        /** Validate all input fields, and warn incompatible inputs */
        mmc_validate_config(&mcx_config, det_ps, dim_det_ps, seed_byte);

        /** Build srcdata from detectors for adjoint / srcid=-2 mode BEFORE
         *  mesh_validate runs.  mesh_validate allocates mesh.weight sized by
         *  the final nsrcslots = max(srcnum, extrasrclen); doing the
         *  detector-append later would leave mesh.weight under-sized for the
         *  Ns+Nd slots the kernel ends up writing to, producing a heap
         *  overflow at simulation end.  This matches mmclab.cpp's order:
         *  build-srcdata first, validate second.
         *
         *  Convention (matches mmclab.cpp and mmc_validate_config):
         *    slots 0..Ns-1        : forward sources (copy of cfg.srcpos)
         *    slots Ns..Ns+Nd-1    : detector-as-adjoint sources
         *
         *  If the multi-source srcpos parser already populated cfg.srcdata
         *  with M forward slots (cfg.extrasrclen == M), keep them and only
         *  append Nd detectors.  Otherwise replicate cfg.srcpos/srcdir Ns
         *  times (single-source / photon-sharing convention).
         */
        if ((MCX_IS_ADJOINT_TYPE(mcx_config.outputtype) || mcx_config.srcid == -2)
                && mcx_config.detnum > 0 && mcx_config.detdir != nullptr) {
            int Nd_pre = mcx_config.detnum;
            int Ns_pre;
            int already_populated_pre = (mcx_config.srcdata != nullptr && mcx_config.extrasrclen > 0);

            if (already_populated_pre) {
                Ns_pre = mcx_config.extrasrclen;
                mcx_config.srcdata = (ExtraSrc*)realloc(mcx_config.srcdata, (Ns_pre + Nd_pre) * sizeof(ExtraSrc));
                memset(mcx_config.srcdata + Ns_pre, 0, Nd_pre * sizeof(ExtraSrc));

                /* Uniform per-slot launch weight (MCX parity, mcx_utils.c:1822, :1985-1986):
                 * every slot launches with srcpos.w = 1 so the total simulated energy
                 * equals N_photon. Matches mmclab.cpp adjoint setup. */
                for (int is = 0; is < Ns_pre; is++) {
                    if (mcx_config.srcdata[is].srcpos.w == 0.f) {
                        mcx_config.srcdata[is].srcpos.w = 1.f;
                    }
                }
            } else {
                if (mcx_config.srcdata) {
                    free(mcx_config.srcdata);
                }

                Ns_pre = (mcx_config.srcnum > 0) ? mcx_config.srcnum : 1;
                mcx_config.srcdata = (ExtraSrc*)calloc(Ns_pre + Nd_pre, sizeof(ExtraSrc));

                /* Slots 0..Ns-1: forward sources (replicate the single main source).
                 * Uniform unit weight per slot (MCX parity); preserve cfg.srcdir.w
                 * (focal length / lens parameter). Detector slots below honor
                 * cfg.detdir.w. */
                for (int is = 0; is < Ns_pre; is++) {
                    mcx_config.srcdata[is].srcpos    = {mcx_config.srcpos.x, mcx_config.srcpos.y,
                                                        mcx_config.srcpos.z, 1.f
                                                       };
                    mcx_config.srcdata[is].srcdir    = {mcx_config.srcdir.x, mcx_config.srcdir.y,
                                                        mcx_config.srcdir.z, mcx_config.srcdir.w
                                                       };
                    mcx_config.srcdata[is].srcparam1 = mcx_config.srcparam1;
                    mcx_config.srcdata[is].srcparam2 = mcx_config.srcparam2;
                }
            }

            mcx_config.extrasrclen = Ns_pre + Nd_pre;

            /* Slots Ns..Ns+Nd-1: detector-as-reversed-source, unit weight per slot. */
            for (int id = 0; id < Nd_pre; id++) {
                mcx_config.srcdata[Ns_pre + id].srcpos  = {mcx_config.detpos[id].x, mcx_config.detpos[id].y,
                                                           mcx_config.detpos[id].z, 1.f
                                                          };
                mcx_config.srcdata[Ns_pre + id].srcdir  = {mcx_config.detdir[id].x, mcx_config.detdir[id].y,
                                                           mcx_config.detdir[id].z, mcx_config.detdir[id].w
                                                          };
                mcx_config.srcdata[Ns_pre + id].srcparam1 = {mcx_config.detpos[id].w, 0.f, 0.f, 0.f};
                mcx_config.srcdata[Ns_pre + id].srcparam2 = {0.f, 0.f, 0.f, 0.f};
            }

            if (MCX_IS_ADJOINT_TYPE(mcx_config.outputtype)) {
                mcx_config.srcid = -1;
            }
        }

        mesh_validate(&mesh, &mcx_config);

        hostdetreclen = (2 + ((mcx_config.ismomentum) > 0)) * mesh.prop + (mcx_config.issaveexit > 0) * 6 + 2;

        /** One must define the domain and properties */
        if (mesh.node == nullptr || mesh.prop == 0) {
            throw py::value_error("You must define 'node' and 'prop' field.");
        }

#if defined(MMC_LOGISTIC) || defined(MMC_SFMT)
        mcx_config.issaveseed = 0;
#endif

        if (mcx_config.issavedet >= 1) {
            mcx_config.exportdetected = NULL;
        }

        if (mcx_config.issaveseed == 1 && mcx_config.seed != SEED_FROM_FILE) {
            mcx_config.photonseed = malloc(mcx_config.maxdetphoton * sizeof(RandType) * RAND_BUF_LEN);
        }

        if (mcx_config.debuglevel & dlTraj) {
            mcx_config.exportdebugdata = (float*)malloc(mcx_config.maxjumpdebug * sizeof(float) * MCX_DEBUG_REC_LEN);
            mcx_config.debuglevel |= dlTraj;
        }

        /* (mesh_srcdetelem and the detector-as-adjoint-source append run above, before
         * mesh_validate, so mesh.weight is sized for the correct nsrcslots.) */

        if (mcx_config.isgpuinfo == 0) {
            mmc_prep(&mcx_config, &mesh, &tracer);
        }

        /** Enclose all simulation calls inside a try/catch construct for exception handling */
        mmc_task_guard = pmmc_task_guard;  /* catches the errors of the photons of mmc_run_mp */

        try {

            if (mcx_config.compute == cbSSE || mcx_config.gpuid > MAX_DEVICE) {
                mmc_run_mp(&mcx_config, &mesh, &tracer);
            }

#ifdef USE_CUDA
            else if (mcx_config.compute == cbCUDA) {
                mmc_run_cu(&mcx_config, &mesh, &tracer);
            }

#endif
#ifdef USE_OPENCL
            else {
                mmc_run_cl(&mcx_config, &mesh, &tracer);
            }

#endif
        } catch (const char* err) {
            exception_msgs.push_back("Error from thread (" + std::to_string(thread_id) + "): " + err);
        } catch (const std::exception& err) {
            exception_msgs.push_back("C++ Error from thread (" + std::to_string(thread_id) + "): " + err.what());
        } catch (...) {
            exception_msgs.push_back("Unknown Exception from thread (" + std::to_string(thread_id) + ")");
        }


        tracer_clear(&tracer);

        /** If error is detected, gracefully terminate the mex and return back to Python */
        if (!exception_msgs.empty()) {
            throw py::runtime_error("PMMC terminated due to an exception!");
        }

        if (mcx_config.debuglevel & dlTraj) {
            field_dim[0] = MCX_DEBUG_REC_LEN;
            field_dim[1] = mcx_config.debugdatalen; // his.savedphoton is for one repetition, should correct
            field_dim[2] = 0;
            field_dim[3] = 0;
            auto photon_traj_data = py::array_t<float, py::array::f_style>({field_dim[0], field_dim[1]});

            if (mcx_config.debuglevel & dlTraj) {
                memcpy(photon_traj_data.mutable_data(), mcx_config.exportdebugdata, field_dim[0] * field_dim[1] * sizeof(float));
            }

            if (mcx_config.exportdebugdata) {
                free(mcx_config.exportdebugdata);
            }

            mcx_config.exportdebugdata = nullptr;
            output["traj"] = photon_traj_data;
        }

        if (mcx_config.issaveseed == 1) {
            field_dim[0] = (mcx_config.issaveseed > 0) * sizeof(RandType) * RAND_BUF_LEN;
            field_dim[1] = mcx_config.detectedcount; // his.savedphoton is for one repetition, should correct
            field_dim[2] = 0;
            field_dim[3] = 0;
            auto detected_seeds = py::array_t<uint8_t, py::array::f_style>({field_dim[0], field_dim[1]});
            memcpy(detected_seeds.mutable_data(), mcx_config.exportseed, field_dim[0] * field_dim[1]);
            free(mcx_config.exportseed);
            mcx_config.exportseed = nullptr;
            output["seeds"] = detected_seeds;
        }

        if (mcx_config.issavedet >= 1) {
            if (mcx_config.issaveexit != 2) {
                field_dim[0] = hostdetreclen;
                field_dim[1] = mcx_config.detectedcount;
                field_dim[2] = 0;
                field_dim[3] = 0;

                if (mcx_config.detectedcount > 0) {
                    auto partial_path = py::array_t<float, py::array::f_style>(std::initializer_list<size_t>({field_dim[0], mcx_config.detectedcount}));
                    memcpy(partial_path.mutable_data(), mcx_config.exportdetected,
                           field_dim[0] * field_dim[1] * sizeof(float));
                    output["detp"] = partial_path;
                }
            } else {
                field_dim[0] = mcx_config.detparam1.w;
                field_dim[1] = mcx_config.detparam2.w;
                field_dim[2] = mcx_config.maxgate;
                field_dim[3] = 0;

                if (field_dim[0] * field_dim[1] > 0) {
                    /* bin the detected photons into a detparam1.w x detparam2.w x maxgate image, as mmclab */
                    auto detimage = py::array_t<float, py::array::f_style>(std::initializer_list<size_t>({field_dim[0], field_dim[1], field_dim[2]}));
                    float* detmap = static_cast<float*>(detimage.mutable_data());
                    memset(detmap, 0, field_dim[0] * field_dim[1] * field_dim[2] * sizeof(float));
                    mesh_getdetimage(detmap, mcx_config.exportdetected, mcx_config.detectedcount, &mcx_config, &mesh);
                    output["detp"] = detimage;
                }
            }

            free(mcx_config.exportdetected);
            mcx_config.exportdetected = NULL;
        }

        if (mcx_config.issave2pt) {
            int isrfforward = (mcx_config.omega > 0.f && mcx_config.seed != SEED_FROM_FILE);
            int isadjoint = MCX_IS_ADJOINT_TYPE(mcx_config.outputtype);

            /* When extrasrclen exceeds srcnum, the kernel runs nsrcslots = extrasrclen
             * slots in the weight buffer regardless of whether outputtype is adjoint
             * (e.g. cfg.srcid = -2 forward-only multi-source mode). */
            int nsrcslots = mesh_outputslots(&mcx_config);  /* one per source only for srcid<0 (also set by the adjoint outputs) */
            int ismultislot = (nsrcslots > mcx_config.srcnum);

            size_t datalen = (mcx_config.method == rtBLBadouelGrid) ? mcx_config.crop0.z : ( (mcx_config.basisorder) ? mesh.nn : mesh.ne);
            field_dim[0] = nsrcslots;
            field_dim[1] = datalen;
            field_dim[2] = mcx_config.maxgate;
            field_dim[3] = 0;
            field_dim[4] = 0;

            std::vector<size_t> array_dims;

            /* Always output forward fluence in output["flux"].
             * In adjoint mode, flux["flux"] has nsrcslots = srcnum+detnum slices (sources then
             * detectors-as-sources), and Jacobians go to output["jmua"], output["jd"], etc. */
            if (mcx_config.method == rtBLBadouelGrid) {
                /* multi-slot layouts use either:
                 *   detector-adjoint slots (nsrcslots > srcnum; adjoint or srcid==-2 forward):
                 *       voxel-fastest, gate-middle, slot-slowest -> [..., maxgate, nsrcslots]
                 *   pattern source (nsrcslots == srcnum > 1):
                 *       pidx-fastest                             -> [srcnum, ..., maxgate] */
                if (ismultislot) {
                    array_dims = {(size_t)mcx_config.dim.x, (size_t)mcx_config.dim.y,
                                  (size_t)mcx_config.dim.z, (size_t)mcx_config.maxgate, (size_t)nsrcslots
                                 };
                } else if (mcx_config.srcnum > 1) {
                    array_dims = {(size_t)mcx_config.srcnum, (size_t)mcx_config.dim.x,
                                  (size_t)mcx_config.dim.y, (size_t)mcx_config.dim.z, (size_t)mcx_config.maxgate
                                 };
                } else {
                    array_dims = {(size_t)mcx_config.dim.x, (size_t)mcx_config.dim.y,
                                  (size_t)mcx_config.dim.z, (size_t)mcx_config.maxgate
                                 };
                }
            } else {
                /* mesh-mode layout. Kernel write order:
                 *   detector-adjoint slots (nsrcslots > srcnum, srcnum==1):
                 *       field[node + gate*nn + slot*nn*maxgate]  -> shape [datalen, maxgate, nsrcslots]
                 *   pattern source (nsrcslots == srcnum > 1):
                 *       field[(gate*ne + eid)*srcnum + pidx]     -> shape [nsrcslots, datalen, maxgate]
                 *   single source:
                 *       field[node + gate*nn]                    -> shape [datalen, maxgate] */
                if (nsrcslots > 1) {
                    if (ismultislot) {
                        array_dims = {datalen, (size_t)mcx_config.maxgate, (size_t)nsrcslots};
                    } else {
                        array_dims = {(size_t)nsrcslots, datalen, (size_t)mcx_config.maxgate};
                    }
                } else {
                    array_dims = {datalen, (size_t)mcx_config.maxgate};
                }
            }

            auto data = py::array_t<double, py::array::f_style>(array_dims);
            memcpy(data.mutable_data(), mesh.weight, data.size() * sizeof(double));
            output["flux"] = data;

            /* RF forward: also return imaginary part in a separate field */
            if (isrfforward && mcx_config.exportadjoint) {
                auto im_data = py::array_t<float, py::array::f_style>(array_dims);
                auto* im_ptr = static_cast<float*>(im_data.mutable_data());
                memcpy(im_ptr, mcx_config.exportadjoint, data.size() * sizeof(float));
                output["fluximag"] = im_data;
            }

            if (mcx_config.issaveref) {
                field_dim[1] = mesh.nf;
                field_dim[2] = mcx_config.maxgate;
                array_dims = {field_dim[1], field_dim[2]};
                auto dref_array = py::array_t<double, py::array::f_style>(array_dims);
                auto* dref = static_cast<double*>(dref_array.mutable_data());
                memcpy(dref, mesh.dref, dref_array.size() * sizeof(double));

                output["dref"] = dref_array;
            }

            /* Output adjoint Jacobian in separate dict fields (output["jmua"], output["jd"], etc.).
             * Grid mode: shape [Nx, Ny, Nz, maxgate, Nch].
             * Mesh mode (basisorder=1): shape [nn, Nch] (CW, no maxgate dim).
             * Nch: Ns*Nd, or the rows of cfg['sdmap'] */
            if (isadjoint && mcx_config.exportjacob) {
                int nsrcpairs = (int)mcx_jacchannels(&mcx_config);
                int isdual = MCX_IS_DUAL_ADJOINT_TYPE(mcx_config.outputtype);
                size_t adjlen;
                std::vector<size_t> jdims;

                if (mcx_config.method == rtBLBadouelGrid) {
                    adjlen = (size_t)mcx_config.dim.x * mcx_config.dim.y * mcx_config.dim.z *
                             mcx_config.maxgate * nsrcpairs;
                    jdims = {(size_t)mcx_config.dim.x, (size_t)mcx_config.dim.y,
                             (size_t)mcx_config.dim.z, (size_t)mcx_config.maxgate,
                             (size_t)nsrcpairs
                            };
                } else {
                    /* mesh mode: nodal output, CW only */
                    adjlen = (size_t)mesh.nn * nsrcpairs;
                    jdims = {(size_t)mesh.nn, (size_t)nsrcpairs};
                }

                const char* jname1 = "jmua";
                const char* jname2 = "jd";

                switch (mcx_config.outputtype) {
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

                /* For RF: output real and imaginary parts in separate _re/_im arrays */
                float* re1 = mcx_config.exportjacob;
                float* re2 = isdual      ? mcx_config.exportjacob + adjlen                    : nullptr;
                float* im1 = isrfforward ? mcx_config.exportjacob + (isdual ? 2 : 1) * adjlen : nullptr;
                float* im2 = (isrfforward && isdual) ? mcx_config.exportjacob + 3 * adjlen    : nullptr;

                auto add_jac = [&](const char* fname, float * re_data, float * im_data) {
                    auto jre = py::array_t<float, py::array::f_style>(jdims);
                    memcpy(jre.mutable_data(), re_data, adjlen * sizeof(float));
                    output[fname] = jre;

                    if (isrfforward && im_data) {
                        std::string imname = std::string(fname) + "_im";
                        auto jim = py::array_t<float, py::array::f_style>(jdims);
                        memcpy(jim.mutable_data(), im_data, adjlen * sizeof(float));
                        output[imname.c_str()] = jim;
                    }
                };

                add_jac(jname1, re1, im1);

                if (isdual) {
                    add_jac(jname2, re2, im2);
                }
            }
        }
    } catch (const char* err) {
        cleanup_configs(gpu_info, mcx_config);
        throw py::runtime_error(err);
    } catch (const py::type_error& err) {
        cleanup_configs(gpu_info, mcx_config);
        throw err;
    } catch (const py::value_error& err) {
        cleanup_configs(gpu_info, mcx_config);
        throw err;
    } catch (const py::runtime_error& err) {
        cleanup_configs(gpu_info, mcx_config);
        std::string error_msg = err.what();

        for (const auto& m : exception_msgs) {
            error_msg += (m + "\n");
        }

        throw py::runtime_error(error_msg);
    } catch (const std::exception& err) {
        cleanup_configs(gpu_info, mcx_config);
        throw py::runtime_error(std::string("C++ Error: ") + err.what());
    } catch (...) {
        cleanup_configs(gpu_info, mcx_config);
        throw py::runtime_error("Unknown exception occurred");
    }

    /** Clear up simulation data structures by calling the destructors */
    cleanup_mesh(mcx_config, mesh);
    // return the MCX output dictionary
    return output;
}


/**
 * @brief Run one photon of mmc_run_mp in a try block, see mmc_task_guard in mmc_utils.h
 *
 * An exception can not leave the OpenMP parallel region of mmc_run_mp: the error of a photon
 * is recorded by mcx_error/mesh_error and raised again after the region; the exception is
 * caught here, in the thread that threw it.
 */

static int pmmc_task_guard(void (*task)(void*), void* arg) {
    try {
        task(arg);
    } catch (...) {
        mmc_defer_error(-1, "unknown exception in a simulation thread", __FILE__, __LINE__);
        return 1;
    }

    return 0;
}

/**
 * @brief Error reporting function in PMMC, equivalent to mcx_error in binary mode
 *
 * @param[in] id: a single integer for the types of the error
 * @param[in] msg: the error message string
 * @param[in] filename: the unit file name where this error is raised
 * @param[in] linenum: the line number in the file where this error is raised
 */

int mmc_throw_exception(const int id, const char* msg, const char* filename, const int linenum) {
    throw msg;
    return id;
}

void print_mcx_usage() {
    std::cout
            << "PMMC (" MMC_VERSION ")\nUsage:\n    output = pmmc.run(cfg);\n\nRun 'help(pmmc.run)' for more details.\n";
}

py::dict pmmc_interface_wargs(py::args args, const py::kwargs& kwargs) {
    if (py::len(kwargs) == 0) {
        print_mcx_usage();
        return {};
    }

    return pmmc_interface(kwargs);
}

py::str print_version() {
    mcconfig mcx_config;            /** mcxconfig: structure to store all simulation parameters */
    mcx_initcfg(&mcx_config);
    mcx_printheader(&mcx_config);
    mcx_clearcfg(&mcx_config);
    return py::str(MMC_VERSION);
}

py::list get_GPU_info() {
    mcconfig mcx_config;            /** mcxconfig: structure to store all simulation parameters */
    GPUInfo* gpu_info = nullptr;        /** gpuinfo: structure to store GPU information */
    mcx_initcfg(&mcx_config);
    mcx_config.isgpuinfo = 3;
    py::list output;
    unsigned int  workdev;

    try {

        if (mcx_config.compute == cbCUDA) {
#ifdef USE_CUDA
            workdev = mcx_list_cu_gpu(&mcx_config, &gpu_info);
#endif
        } else {
#ifdef USE_OPENCL
            mcx_list_cl_gpu(&mcx_config, &workdev, NULL, &gpu_info);
#endif
        }

    } catch (...) {
        std::cerr << "No CUDA-capable device was found." << std::endl;
        return output;
    }

    if (!workdev) {
        std::cerr << "no active GPU device found" << std::endl;
    }

    if (workdev > MAX_DEVICE) {
        workdev = MAX_DEVICE;
    }

    for (int i = 0; i < gpu_info[0].devcount; i++) {
        py::dict current_device_info;
        current_device_info["name"] = gpu_info[i].name;
        current_device_info["id"] = gpu_info[i].id;
        current_device_info["devcount"] = gpu_info[i].devcount;
        current_device_info["major"] = gpu_info[i].major;
        current_device_info["minor"] = gpu_info[i].minor;
        current_device_info["globalmem"] = gpu_info[i].globalmem;
        current_device_info["constmem"] = gpu_info[i].constmem;
        current_device_info["sharedmem"] = gpu_info[i].sharedmem;
        current_device_info["regcount"] = gpu_info[i].regcount;
        current_device_info["clock"] = gpu_info[i].clock;
        current_device_info["sm"] = gpu_info[i].sm;
        current_device_info["core"] = gpu_info[i].core;
        current_device_info["autoblock"] = gpu_info[i].autoblock;
        current_device_info["autothread"] = gpu_info[i].autothread;
        current_device_info["maxgate"] = gpu_info[i].maxgate;
        output.append(current_device_info);
    }

    cleanup_configs(gpu_info, mcx_config);
    return output;
}

/**
 * The docstring of run(), one literal per line and joined on first use: as one literal, it would
 * exceed the 4095 characters -pedantic accepts and the 16 KB limit of MSVC
 */
static const char* pmmc_run_doc_lines[] = {
    "Run an MMC simulation on a tetrahedral mesh.\n",
    "\n",
    "    res = pmmc.run(cfg)\n",
    "    res = pmmc.run(nphoton=1000000, node=node, elem=elem, elemprop=elemprop, ...)\n",
    "\n",
    "The settings are given as a dict, or as keyword arguments. Array values may be\n",
    "NumPy arrays (C or Fortran order) or nested lists, with the shapes of mmclab;\n",
    "element and node indices are 1-based. pmmc.run computes facenb, evol and e0 when\n",
    "they are absent, but does not reorient the elements or add the source domain of\n",
    "wide-field sources: pmmc.mmclab(cfg) does both, and post-processes the outputs as\n",
    "mmclab does.\n",
    "\n",
    "Parameters\n",
    "----------\n",
    "cfg : dict\n",
    "    The simulation settings. Keys marked with * are required; values in [] are the\n",
    "    defaults. Lengths are in the unit of the mesh nodes; cfg['unitinmm'] gives it in\n",
    "    mm.\n",
    "\n",
    "    Mesh and optical properties\n",
    "\n",
    "    * nphoton : int\n",
    "        Total number of photons to simulate.\n",
    "    * node : ndarray, shape (Nn, 3)\n",
    "        The mesh nodes.\n",
    "    * elem : ndarray, shape (Ne, 4)\n",
    "        The tetrahedra, 1-based node indices, with positive volumes (see\n",
    "        iso2mesh.meshreorient).\n",
    "    * elemprop : ndarray, shape (Ne,)\n",
    "        The medium label (from 1) of each element; elements labeled -1 form the\n",
    "        source domain of wide-field sources.\n",
    "    face : ndarray, shape (Nf, 5)\n",
    "        OptiX only: a surface mesh used instead of elem, [n1, n2, n3, inner label,\n",
    "        outer label].\n",
    "    * prop : ndarray, shape (N, 4)\n",
    "        One row per medium, [mua (1/mm), mus (1/mm), g, n]; row 0 is label 0\n",
    "        (outside of the mesh).\n",
    "    facenb : ndarray, shape (Ne, 4)\n",
    "        Element face neighbors [computed].\n",
    "    evol : ndarray, shape (Ne,)\n",
    "        Element volumes [computed].\n",
    "    e0 : int\n",
    "        The 1-based element enclosing the source [computed].\n",
    "    unitinmm : float\n",
    "        The length unit of the mesh in mm [1.0].\n",
    "    nout : float\n",
    "        Refractive index outside of the mesh [1.0].\n",
    "    nodemua, nodemusp : ndarray, shape (Nn,)\n",
    "        Per-node mua and mus', replacing those of prop.\n",
    "    noderoi, edgeroi, faceroi : ndarray\n",
    "        Implicit MMC (iMMC): sphere radii per node, (Nn,); cylinder radii per\n",
    "        element edge, (Ne, 6); slab thicknesses per element face, (Ne, 4).\n",
    "\n",
    "    Time gates\n",
    "\n",
    "    * tstart, tend, tstep : float\n",
    "        Start and end of the simulated time window, and width of each gate (s).\n",
    "\n",
    "    Source\n",
    "\n",
    "    * srcpos : array_like\n",
    "        [x, y, z] position of the source; several sources: one row per source.\n",
    "    * srcdir : array_like\n",
    "        [vx, vy, vz] unit launch direction; an optional 4th value is the focal\n",
    "        length: 0 collimated, >0 converging, <0 diverging.\n",
    "    srctype : str\n",
    "        ['pencil'], 'isotropic', 'arcsine', 'cone', 'zgaussian', 'gaussian',\n",
    "        'disk', 'planar', 'pattern', 'pattern3d', 'fourier', 'fourierx',\n",
    "        'fourierx2d', 'line', 'slit'; on the CPU and OptiX also 'pencilarray',\n",
    "        'hyperboloid' and 'ring'. See the mmclab help (help mmclab in MATLAB) for\n",
    "        srcparam1/srcparam2 of each.\n",
    "    srcparam1, srcparam2 : array_like\n",
    "        Four parameters each, see srctype.\n",
    "    srcpattern : ndarray\n",
    "        The intensity pattern of 'pattern' sources, (Nx, Ny), or (srcnum, Nx, Ny)\n",
    "        for several patterns at once; or the 3-D volume of 'pattern3d'.\n",
    "    srcid : int\n",
    "        With several sources: [0] all in one output, -1 each in its own output, -2\n",
    "        also run each detector as a source, N>0 only the N-th source.\n",
    "    voidtime : int\n",
    "        For a source outside of the mesh, [1]: start the clock at launch, 0: when\n",
    "        the photon enters the mesh.\n",
    "    isspecular : int\n",
    "        1: apply the specular reflection at entry [0].\n",
    "    omega : float\n",
    "        Modulation angular frequency 2*pi*f (rad/s) [0]; >0 makes an RF run.\n",
    "\n",
    "    Detectors\n",
    "\n",
    "    detpos : ndarray, shape (N, 4)\n",
    "        One row per disk detector, [x, y, z, radius].\n",
    "    detdir : ndarray, shape (N, 4)\n",
    "        The inward direction and focal length of each detector, required by the\n",
    "        adjoint outputs.\n",
    "    sdmap : ndarray, shape (Nch, 2)\n",
    "        [srcid, detid] rows (1-based) selecting the channels of the adjoint\n",
    "        Jacobians [all pairs, detector index fastest].\n",
    "    detparam1, detparam2 : array_like\n",
    "        [ux, uy, uz, Nx] and [vx, vy, vz, Ny]: the edges and pixel counts of a\n",
    "        wide-field camera detector at detpos, used with issaveexit=2.\n",
    "    maxdetphoton : int\n",
    "        The most detected photons kept [1000000].\n",
    "    issavedet : int\n",
    "        1: save the detected photons; 0: do not. [1 when detpos is given, except with the\n",
    "        adjoint outputs; 0 otherwise]\n",
    "    issaveexit : int\n",
    "        1: also save the exit position and direction [0]; 2: instead, add up the\n",
    "        exiting weight on the wide-field detector.\n",
    "    ismomentum : int\n",
    "        1: also save the momentum transfer [0].\n",
    "    issaveseed : int\n",
    "        1: save the RNG seed of each detected photon, for replay [0].\n",
    "\n",
    "    Monte Carlo settings\n",
    "\n",
    "    seed : int or ndarray\n",
    "        RNG seed [1648335518]; a uint8 array (the 'seeds' output of an earlier run)\n",
    "        replays those photons.\n",
    "    isreflect : int\n",
    "        Boundary condition at the mesh surface: [1] Fresnel, 0 index-matched, 2\n",
    "        absorbing, 3 mirror.\n",
    "    isnormalized : int\n",
    "        [1]: normalize the output to a unit source; 0: raw output.\n",
    "    minenergy : float\n",
    "        Weight below which Russian roulette starts [1e-6].\n",
    "    roulettesize : float\n",
    "        The Russian roulette survival factor [10].\n",
    "    mcmethod : int\n",
    "        [0]: MCX-style (Beer-Lambert absorption); 1: MCML-style (albedo weight).\n",
    "    method : str\n",
    "        Ray tracer: ['grid'] dual-grid MMC with a voxel output; 'elem' a mesh\n",
    "        output; 'plucker', 'havel', 'badouel' on the CPU only.\n",
    "    steps : array_like\n",
    "        [dx, dy, dz], the voxel size of the 'grid' output [1, 1, 1].\n",
    "    basisorder : int\n",
    "        Mesh output: [1] values at the nodes, 0 values in the elements.\n",
    "\n",
    "    Output\n",
    "\n",
    "    outputtype : str\n",
    "        What res['flux'] holds: ['flux'] fluence rate (1/mm^2/s), 'fluence',\n",
    "        'energy'; replay: 'jacobian' (or 'wl'), 'wp'; forward with omega>0: 'rf';\n",
    "        adjoint Jacobians (need detpos and detdir): 'adjoint', 'adjoint_dcoeff',\n",
    "        'adjoint_mus', 'adjoint_musp', 'adjoint_mua_d', 'adjoint_mua_musp'.\n",
    "    adjointmode : int\n",
    "        Mesh-output adjoint mua Jacobian: [0] the full FEM form, 1 the nodal\n",
    "        approximation.\n",
    "    issave2pt : int\n",
    "        [1]: compute res['flux']; 0: skip it.\n",
    "    issaveref : int\n",
    "        1: save the diffuse reflectance on the mesh surface in res['dref'] [0].\n",
    "    session : str\n",
    "        A name for the output files.\n",
    "    flog : int or str\n",
    "        Where the messages go: [2] sys.stdout, 1 the stdout of the process, 0 the\n",
    "        same without the banner, or a file name (appended); '/dev/null' (Linux,\n",
    "        macOS) or 'NUL' (Windows) silences everything. Errors are always raised.\n",
    "\n",
    "    Replay\n",
    "\n",
    "    detphotons : ndarray\n",
    "        The 'detp' output of the baseline run.\n",
    "    replaydet : int\n",
    "        [0]: replay all detectors into one output; N>0: only detector N.\n",
    "\n",
    "    Compute backend\n",
    "\n",
    "    compute : str\n",
    "        'sse' (CPU, OpenMP), ['opencl'], 'cuda' or 'optix'; 'cuda' and 'optix' need\n",
    "        a module built with them.\n",
    "    gpuid : int or str\n",
    "        The device to use, from 1 [1]; -1 runs on the CPU; a string of 0s and 1s\n",
    "        selects several, e.g. '1101'. pmmc.gpuinfo() lists them.\n",
    "    workload : array_like\n",
    "        The relative share of photons for each selected device.\n",
    "    nthread, nblocksize : int\n",
    "        Total GPU threads [8192] and threads per block [64].\n",
    "    optlevel : int\n",
    "        OpenCL kernel optimization level, 0 (none) to 3 [3].\n",
    "    compileropt : str\n",
    "        Extra OpenCL JIT compiler options.\n",
    "    kernelfile : str\n",
    "        An OpenCL kernel file to run instead of the built-in one.\n",
    "\n",
    "    Debugging\n",
    "\n",
    "    debuglevel : str\n",
    "        Debug flags, e.g. 'P' progress bar, 'T' timing, 'M' record photon\n",
    "        trajectories; see the mmclab help for the others.\n",
    "    debugphoton : int\n",
    "        Print the debug information of this photon only [-1: all].\n",
    "    maxjumpdebug : int\n",
    "        The most trajectory positions recorded [10000000].\n",
    "\n",
    "Returns\n",
    "-------\n",
    "res : dict\n",
    "    flux : ndarray, float64\n",
    "        The quantity of outputtype, normalized as set by isnormalized: with method\n",
    "        'grid', (Nx, Ny, Nz, Ngate) voxels of size steps covering the mesh; with\n",
    "        the mesh tracers, (Nn, Ngate) (basisorder 1) or (Ne, Ngate) (basisorder 0).\n",
    "        With srcnum>1, a first dimension of srcnum; with several source slots, a\n",
    "        last dimension. For RF runs, the real part.\n",
    "    fluximag : ndarray\n",
    "        The imaginary part of flux, for RF runs.\n",
    "    dref : ndarray, shape (Nf, Ngate)\n",
    "        The diffuse reflectance on the surface triangles (issaveref=1).\n",
    "    jmua, jd, jmus, jmusp : ndarray\n",
    "        The adjoint Jacobians, (Nx, Ny, Nz, Ngate, Nch) with method 'grid' or\n",
    "        (Nn, Nch) on the mesh; for RF, their imaginary parts are in jmua_im, jd_im\n",
    "        and so on.\n",
    "    detp : ndarray, float32\n",
    "        The detected photons (issavedet=1 and detpos given), (ncol, ndetected), one\n",
    "        column per photon: [detid, nscat (N), ppath (N), mom (N, ismomentum=1),\n",
    "        p (3) and v (3) (issaveexit=1), w0]; with issaveexit=2, the wide-field\n",
    "        detector image, (Nx, Ny, Ngate).\n",
    "    seeds : ndarray, uint8, shape (seedbyte, ndetected)\n",
    "        The RNG seed of each detected photon (issaveseed=1).\n",
    "    traj : ndarray, float32, shape (6, npos)\n",
    "        The photon trajectories (debuglevel containing 'M'), one column per\n",
    "        position, [id, x, y, z, weight, eid], the id and eid stored as uint32 bits.\n",
    "\n",
    "Raises\n",
    "------\n",
    "ValueError\n",
    "    If a setting is missing or invalid.\n",
    "RuntimeError\n",
    "    If the simulation fails.\n",
    "\n",
    "Examples\n",
    "--------\n",
    ">>> import numpy as np, iso2mesh as i2m, pmmc\n",
    ">>> node, face, elem = i2m.meshabox([0, 0, 0], [60, 60, 60], 10, 100)\n",
    ">>> res = pmmc.run(nphoton=1000000, node=node, elem=elem,\n",
    "...                elemprop=np.ones(elem.shape[0]), tstart=0, tend=5e-9, tstep=5e-9,\n",
    "...                srcpos=[30, 30, 0], srcdir=[0, 0, 1],\n",
    "...                prop=[[0, 0, 1, 1], [0.005, 1, 0.01, 1.37]])\n",
    nullptr
};

static const char* pmmc_run_doc() {
    static std::string doc;

    if (doc.empty()) {
        for (const char** line = pmmc_run_doc_lines; *line; line++) {
            doc += *line;
        }
    }

    return doc.c_str();
}

PYBIND11_MODULE(_pmmc, m) {
    m.doc() = "PMMC (" MMC_VERSION "): Python bindings for the Mesh-based Monte Carlo (MMC) photon transport simulator, https://mcx.space/#mmc";
    m.def("run", &pmmc_interface, pmmc_run_doc(), py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("run", &pmmc_interface_wargs, "The same as run(cfg), with the settings given as keyword arguments.", py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("gpuinfo",
          &get_GPU_info,
          "List the OpenCL devices.\n\n"
          "Returns\n-------\nlist of dict\n"
          "    One dict per device: name, id, devcount, major, minor, globalmem, constmem, sharedmem,\n"
          "    regcount, clock, sm, core, autoblock, autothread and maxgate.",
          py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("version",
          &print_version,
          "Print the MMC banner and return the version string.\n\n"
          "Returns\n-------\nstr\n    The version, e.g. 'v2026'.",
          py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
}
