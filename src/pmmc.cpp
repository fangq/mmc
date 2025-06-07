/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2025
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
            std::cout << #property << ": " << (float) (dst_mcx_config).property << std::endl;} \
        catch (const std::runtime_error &err)\
        {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#property) + ". Reason: " + err.what()));}\
    }

#define GET_VEC3_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>()};\
            std::cout << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << "]\n";} \
        catch (const std::runtime_error &err ) {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#prop) + ". Reason: " + err.what()));}}

#define GET_VEC4_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>(), list[3].cast<type>()}; \
            std::cout << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << ", " << dst.prop.w << "]\n";} \
        catch (const std::runtime_error &err ) {throw py::type_error(std::string("Failed to assign MCX property " + std::string(#prop) + ". Reason: " + err.what()));}}

#define GET_VEC34_FIELD(src, dst, prop, type) if (src.contains(#prop)) {try {auto list = py::list(src[#prop]);\
            dst.prop = {list[0].cast<type>(), list[1].cast<type>(), list[2].cast<type>(), list.size() == 4 ? list[3].cast<type>() : dst.prop.w}; \
            std::cout << #prop << ": [" << dst.prop.x << ", " << dst.prop.y << ", " << dst.prop.z << ", " << dst.prop.w << "]\n";} \
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


/**
 * Parse user input cfg object and convert to to MMC's mcconfig object.
 * @param user_cfg
 * @param mcx_config
 */


void parse_config(const py::dict& user_cfg, mcconfig& mcx_config, tetmesh& mesh) {
    mcx_initcfg(&mcx_config);
    mesh_init(&mesh);

    mcx_config.flog = stderr;

    GET_SCALAR_FIELD(user_cfg, mcx_config, nphoton, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, nblocksize, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, nthread, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tstart, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tstep, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, tend, py::float_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isreflect, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isspecular, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, ismomentum, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issaveexit, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issave2pt, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issavedet, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, issaveseed, py::bool_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, optlevel, py::int_);
    GET_SCALAR_FIELD(user_cfg, mcx_config, isatomic, py::bool_);
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
    GET_VEC3_FIELD(user_cfg, mcx_config, srcpos, float);
    GET_VEC34_FIELD(user_cfg, mcx_config, srcdir, float);
    GET_VEC3_FIELD(user_cfg, mcx_config, steps, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, srcparam1, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, srcparam2, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, detparam1, float);
    GET_VEC4_FIELD(user_cfg, mcx_config, detparam2, float);

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
            throw py::value_error("the 'methods' field must be a non-empty string");
        }

        mcx_config.method = mcx_keylookup((char*)(method_str.c_str()), methods);

        if (mcx_config.method == -1) {
            throw py::value_error("the specified ray-tracing method is not supported");
        }
    }


    if (user_cfg.contains("outputtype")) {
        std::string output_type_str = py::str(user_cfg["outputtype"]);
        const char* outputtype[] = {"flux", "fluence", "energy", "jacobian", "nscat", "wl", "wp", "wm", "rf", "length", "rfmus", "wltof", "wptof"};

        if (output_type_str.empty()) {
            throw py::value_error("the 'srctype' field must be a non-empty string");
        }

        mcx_config.outputtype = mcx_keylookup((char*)(output_type_str.c_str()), outputtype);

        if (mcx_config.outputtype >= 5) { // map wl to jacobian, wp to nscat
            mcx_config.outputtype -= 2;
        }

        if (mcx_config.outputtype == -1) {
            throw py::value_error("the specified output type is not supported");
        }
    }


    if (user_cfg.contains("compute")) {
        std::string compute_str = py::str(user_cfg["compute"]);
        const char* computebackend[] = {"sse", "opencl", "cuda", ""};

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
        const char debugflag[] = {'R', 'M', 'P', '\0'};
        char debuglevel[MAX_SESSION_LENGTH] = {'\0'};

        if (debug_level.empty()) {
            throw py::value_error("the 'debuglevel' field must be a non-empty string");
        }

        if (debug_level.size() > MAX_SESSION_LENGTH) {
            throw py::value_error("the 'debuglevel' field is too long");
        }

        strncpy(debuglevel, debug_level.c_str(), MAX_SESSION_LENGTH);
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

            if (buffer_info.shape.at(0) != sizeof(float) * RAND_BUF_LEN) {
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
            mcx_config.gpuid = py::int_(gpu_id_value);
            memset(mcx_config.deviceid, 0, MAX_DEVICE);

            if (mcx_config.gpuid > 0 && mcx_config.gpuid < MAX_DEVICE) {
                memset(mcx_config.deviceid, '0', mcx_config.gpuid - 1);
                mcx_config.deviceid[mcx_config.gpuid - 1] = '1';
            } else {
                throw py::value_error("GPU id must be positive and can not be more than 256");
            }
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

    //
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

    if (mesh.evol == NULL) {
        mesh_getvolume(&mesh, &mcx_config);
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
}

inline void cleanup_mesh(mcconfig& mcx_config, tetmesh& mesh) {
    mesh_clear(&mesh, &mcx_config);
    mcx_clearcfg(&mcx_config);
}


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
            mcx_list_cu_gpu(&mcx_config, &active_dev, NULL, &gpu_info);
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

        /** Validate all input fields, and warn incompatible inputs */
        mmc_validate_config(&mcx_config, det_ps, dim_det_ps, seed_byte);
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
            mcx_config.exportdetected = (float*) malloc(hostdetreclen * mcx_config.maxdetphoton * sizeof(float));
        }

        if (mcx_config.issaveseed == 1) {
            mcx_config.photonseed = malloc(mcx_config.maxdetphoton * sizeof(float) * RAND_BUF_LEN);
        }

        if (mcx_config.debuglevel & MCX_DEBUG_MOVE) {
            mcx_config.exportdebugdata = (float*)malloc(mcx_config.maxjumpdebug * sizeof(float) * MCX_DEBUG_REC_LEN);
            mcx_config.debuglevel |= dlTraj;
        }

        mesh_srcdetelem(&mesh, &mcx_config);

        if (mcx_config.isgpuinfo == 0) {
            mmc_prep(&mcx_config, &mesh, &tracer);
        }

        /** Enclose all simulation calls inside a try/catch construct for exception handling */
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

        if (mcx_config.debuglevel & MCX_DEBUG_MOVE) {
            field_dim[0] = MCX_DEBUG_REC_LEN;
            field_dim[1] = mcx_config.debugdatalen; // his.savedphoton is for one repetition, should correct
            field_dim[2] = 0;
            field_dim[3] = 0;
            auto photon_traj_data = py::array_t<float, py::array::f_style>({field_dim[0], field_dim[1]});

            if (mcx_config.debuglevel & MCX_DEBUG_MOVE) {
                memcpy(photon_traj_data.mutable_data(), mcx_config.exportdebugdata, field_dim[0] * field_dim[1] * sizeof(float));
            }

            if (mcx_config.exportdebugdata) {
                free(mcx_config.exportdebugdata);
            }

            mcx_config.exportdebugdata = nullptr;
            output["traj"] = photon_traj_data;
        }

        if (mcx_config.issaveseed == 1) {
            field_dim[0] = (mcx_config.issaveseed > 0) * RAND_BUF_LEN * sizeof(float);
            field_dim[1] = mcx_config.detectedcount; // his.savedphoton is for one repetition, should correct
            field_dim[2] = 0;
            field_dim[3] = 0;
            auto detected_seeds = py::array_t<uint8_t, py::array::f_style>({field_dim[0], field_dim[1]});
            memcpy(detected_seeds.mutable_data(), mcx_config.photonseed, field_dim[0] * field_dim[1]);
            free(mcx_config.photonseed);
            mcx_config.photonseed = nullptr;
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
                    auto partial_path = py::array_t<float, py::array::f_style>(std::initializer_list<size_t>({field_dim[0], field_dim[1], field_dim[2]}));
                    memcpy(partial_path.mutable_data(), mcx_config.exportdetected,
                           field_dim[0] * field_dim[1] * field_dim[2] * sizeof(float));
                    output["detp"] = partial_path;
                }
            }

            free(mcx_config.exportdetected);
            mcx_config.exportdetected = NULL;
        }

        if (mcx_config.issave2pt) {
            size_t datalen = (mcx_config.method == rtBLBadouelGrid) ? mcx_config.crop0.z : ( (mcx_config.basisorder) ? mesh.nn : mesh.ne);
            field_dim[0] = mcx_config.srcnum;
            field_dim[1] = datalen;
            field_dim[2] = mcx_config.maxgate;
            field_dim[3] = 0;
            field_dim[4] = 0;

            std::vector<size_t> array_dims;

            if (mcx_config.method == rtBLBadouelGrid) {
                field_dim[0] = mcx_config.srcnum;
                field_dim[1] = mcx_config.dim.x;
                field_dim[2] = mcx_config.dim.y;
                field_dim[3] = mcx_config.dim.z;
                field_dim[4] = mcx_config.maxgate;

                if (mcx_config.srcnum > 1) {
                    array_dims = {field_dim[0], field_dim[1], field_dim[2], field_dim[3], field_dim[4]};
                } else {
                    array_dims = {field_dim[1], field_dim[2], field_dim[3], field_dim[4]};
                }
            } else {
                if (mcx_config.srcnum > 1) {
                    array_dims = {field_dim[0], field_dim[1], field_dim[2]};
                } else {
                    array_dims = {field_dim[1], field_dim[2]};
                }
            }

            auto data = py::array_t<double, py::array::f_style>(array_dims);
            memcpy(data.mutable_data(), mesh.weight, data.size() * sizeof(double));
            output["flux"] = data;

            if (mcx_config.issaveref) {
                field_dim[1] = mesh.nf;
                field_dim[2] = mcx_config.maxgate;
                array_dims = {field_dim[1], field_dim[2]};
                auto dref_array = py::array_t<double, py::array::f_style>(array_dims);
                auto* dref = static_cast<double*>(dref_array.mutable_data());
                memcpy(dref, mesh.dref, dref_array.size() * sizeof(double));

                output["dref"] = dref_array;
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
            mcx_list_cu_gpu(&mcx_config, &workdev, NULL, &gpu_info);
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

PYBIND11_MODULE(_pmmc, m) {
    m.doc() = "PMMC (" MMC_VERSION "): Python bindings for Monte Carlo eXtreme photon transport simulator, https://mcx.space";
    m.def("run", &pmmc_interface, "Runs MCX with the given config.", py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("run", &pmmc_interface_wargs, "Runs MCX with the given config.", py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("gpuinfo",
          &get_GPU_info,
          "Prints out the list of CUDA-capable devices attached to this system.",
          py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
    m.def("version",
          &print_version,
          "Prints mcx version information.",
          py::call_guard<py::scoped_ostream_redirect,
          py::scoped_estream_redirect>());
}
