#!/bin/bash
# Builds the mmc-trinity packages on Windows, the Windows twin of
# .github/linux_release_builder/build_mmc_trinity.sh.  On Windows nvcc takes no
# host compiler other than MSVC, so everything is built with MSVC, from a Git
# Bash where the MSVC environment is set up (e.g. ilammy/msvc-dev-cmd with
# toolset 14.29: the STL of newer toolsets refuses CUDA 11.7):
#
#   OPTIXHOME=/path/to/optix-sdk [MATLAB_ROOT=...] [OCTAVE_HOME=.../mingw64] \
#   [PYTHONS="/c/.../python.exe ..."] bash .github/windows_release_builder/build_mmc_trinity.sh
#
# As on Linux, CUDA 11.7 comes from NVIDIA's redistributable archives (nvcc,
# cudart, cccl; cudart also carries the CL headers and OpenCL.lib), and
# OPTIXHOME holds the OptiX 7.5 SDK headers (include/optix*.h), e.g. a clone of
# https://github.com/NVIDIA/optix-dev at v7.5.0.  The objects are those of
# "make optix" (sse, opencl, cuda and optix backends), with the same defines;
# the CUDA code is built for sm_52 to sm_86, with the sm_86 PTX for newer GPUs.
#
# It builds
#   bin/mmc.exe          mmc-trinity
#   mmclab/mmc.mexw64    mmclab-trinity for MATLAB (with MATLAB_ROOT), by MATLAB's mex
#   mmclab/mmc.mex       mmclab-trinity for Octave 10+ (with OCTAVE_HOME): Octave's
#                        C mex API is in liboctmex, which MSVC links by an import library
#   pmmc/dist/*.whl      pmmc-trinity: pmmc with _rtmmc, for each Python in PYTHONS;
#                        _rtmmc also runs sse and opencl, so it is the only module
# Each has the C and CUDA runtimes linked in statically (the mex and Python
# modules: the C++ runtime DLL, as their hosts use); nvcuda.dll (the driver, at
# least R516 for CUDA 11.7) and OpenCL.dll are delay-loaded when a GPU backend
# is used, and the OpenMP runtime of MSVC, vcomp140.dll, is copied next to each.

set -e
cd "$(dirname "$0")/../.."
ROOT=$PWD

: "${OPTIXHOME:?set OPTIXHOME to the OptiX SDK folder (with include/optix.h)}"
[ -f "$OPTIXHOME/include/optix.h" ] || { echo "no $OPTIXHOME/include/optix.h" >&2; exit 1; }
OPTIXW=$(cygpath -m "$OPTIXHOME")

CUDAVER=11.7.99
CCCLVER=11.7.91
REDIST=https://developer.download.nvidia.com/compute/cuda/redist
CUDAHOME=${CUDAHOME:-$ROOT/../cuda-$CUDAVER}

if [ ! -x "$CUDAHOME/bin/nvcc.exe" ]; then
    mkdir -p "$CUDAHOME"
    DL=$(mktemp -d)
    for pkg in "cuda_nvcc:$CUDAVER" "cuda_cudart:$CUDAVER" "cuda_cccl:$CCCLVER"; do
        name=${pkg%%:*}
        ver=${pkg##*:}
        curl -sSfL --retry 3 -o "$DL/$name.zip" "$REDIST/$name/windows-x86_64/$name-windows-x86_64-$ver-archive.zip"
        unzip -q -o "$DL/$name.zip" -d "$DL"
        cp -af "$DL/$name-windows-x86_64-$ver-archive/." "$CUDAHOME/"
    done
    rm -rf "${DL:?}"
fi
CUDAW=$(cygpath -m "$CUDAHOME")
export PATH="$CUDAHOME/bin:$PATH"

# the options below are MSVC's, written with "-" rather than "/" (which Git
# Bash would take for paths); every path handed to cl, nvcc or link is mixed
# (C:/...) or relative, so MSYS has nothing to convert
export MSYS2_ARG_CONV_EXCL='*'

command -v cl > /dev/null || { echo "cl.exe is not on the PATH: set up the MSVC environment first" >&2; exit 1; }
# nvcc --use-local-env: the MSVC environment is set up already; without it nvcc
# runs vcvars64.bat again, which fails once MATLAB and Python are set up as well
# Git's own link.exe (coreutils) may come first on the PATH
MSBIN="$(dirname "$(command -v cl)")"
MSLINK="$MSBIN/link.exe"

cd src

# Git for Windows checks symbolic links out as small text files holding the target
for l in mmc_core.cu zmat/miniz/zlib.h; do
    if [ ! -L "$l" ] && [ "$(wc -c < "$l")" -lt 64 ]; then
        cp "$(dirname "$l")/$(cat "$l")" "$l"
    fi
done

rm -rf built ../bin/mmc.exe ../mmclab/mmc.mexw64 ../mmclab/mmc.mex ../pmmc/dist
mkdir -p built ../bin

# the OpenCL kernel and the OptiX PTX, as 0-terminated char arrays
xxd -i mmc_core.cl | sed 's/\([0-9a-f]\)$/\0, 0x00/' > mmc_core.clh

USERDEFS="-DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT -DUSE_CUDA -DUSE_OPTIX -DUSE_CUDA_VECTOR_TYPES"
GPUINC="-I$CUDAW/include -I$OPTIXW/include"

# the OptiX programs stay PTX of nvcc's default target: OptiX compiles them for the GPU
nvcc --use-local-env -O3 -ptx --expt-relaxed-constexpr -use_fast_math -m64 $USERDEFS $GPUINC -o built/mmc_optix_core.ptx mmc_optix_core.cu
(cd built && xxd -i mmc_optix_core.ptx | sed 's/\([0-9a-f]\)$/\0, 0x00/' > mmc_optix_core_ptx.h)

# native code from Maxwell (sm_52) to Ampere (sm_86), and the sm_86 PTX for newer GPUs to JIT;
# PTX of a pre-Volta GPU (nvcc's default sm_52) runs the mmc kernel 2-3x slower on newer GPUs
GENCODE=""
for sm in 52 60 70 75 80 86; do
    GENCODE="$GENCODE -gencode=arch=compute_$sm,code=sm_$sm"
done
GENCODE="$GENCODE -gencode=arch=compute_86,code=compute_86"

CLOPT="-nologo -c -O2 -openmp -EHsc -W1 -D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES"
MMCDEFS="-DMCX_EMBED_CL $USERDEFS -DMMC_USE_SSE -DHAVE_SSE2 -DUSE_SSE2 -DMMC_USE_SSE_MATH"
MMCINC="-I. -Ibuilt -Izmat -Izmat/easylzma -Iubj $GPUINC"
LIBPATHS="-LIBPATH:$CUDAW/lib/x64"
GPULIBS="cudart_static.lib cuda.lib OpenCL.lib delayimp.lib -DELAYLOAD:nvcuda.dll -DELAYLOAD:OpenCL.dll"

# objects <folder> <runtime: -MT or -MD> <defines> <sources...>: the C/C++ sources
# and the CUDA backend (with the flags of "make cuda", CUCCOPT) into a folder
objects() {
    local out=$1 rt=$2 defs=$3 f std
    shift 3
    mkdir -p "$out"
    for f in "$@"; do
        std=""
        [[ "$f" == *.cpp ]] && std="-std:c++14"
        cl $CLOPT $rt $std $MMCDEFS $defs $MMCINC -Fo"$out/$(basename "${f%.*}").obj" "$f"
    done
    nvcc --use-local-env -c -m64 -use_fast_math -Xptxas -O3 -Xcompiler $rt,-openmp -D_USE_MATH_DEFINES $GENCODE \
         -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL $defs \
         -o "$out/mmc_cu_host.obj" mmc_cu_host.cu
}

MMCSRC="mmc_rand_xorshift128p.c mmc_mesh.c mmc_raytrace.c mmc_utils.c mmc_tictoc.c mmc_host.c mmc_bench.c
        mmc_cl_utils.c mmc_cl_host.c cjson/cJSON.c ubj/ubjw.c mmc_highorder.cpp mmc_optix_utils.cpp mmc_optix_host.cpp"

echo "== mmc-trinity"
# zmat, as its Makefile builds it for mmc: miniz for zlib, the LZMA SDK, lz4; no
# blosc2/zstd; with MSVC its threads (block-parallel zlib, LZMA) are Win32's
ZMATDEFS="-DNO_BLOSC2 -DNO_ZSTD -DNO_ZLIB -D_LARGEFILE64_SOURCE=1 -DZMAT_USE_LZMA_SDK"
ZMATINC="-Izmat -Izmat/miniz -Izmat/easylzma -Izmat/easylzma/lzma -Izmat/lz4"
mkdir -p built/exe
for f in zmatlib miniz/miniz easylzma/compress easylzma/decompress easylzma/lzma_header \
         easylzma/lzip_header easylzma/common_internal lz4/lz4 lz4/lz4hc \
         easylzma/lzma/{LzmaEnc,LzmaDec,LzmaLib,LzFind,Bra,BraIA64,Alloc,7zCrc,7zCrcOpt,CpuArch} \
         easylzma/lzma/{Lzma2Dec,Lzma2DecMt,Lzma2Enc,MtCoder,MtDec,Xz,XzCrc64,XzCrc64Opt,XzEnc,XzDec} \
         easylzma/lzma/{Sha256,Sha256Opt,Delta,Bra86,7zStream,LzFindMt,LzFindOpt,Threads}; do
    cl $CLOPT -MT $ZMATDEFS $ZMATINC -Fo"built/exe/zmat_$(echo $f | tr / _).obj" "zmat/$f.c"
done
# -MT: the static C runtime, as the CUDA runtime (cudart_static) is linked
objects built/exe -MT "" $MMCSRC mmc.c mmc_neurojson.cpp
"$MSLINK" -nologo -OUT:../bin/mmc.exe built/exe/*.obj $LIBPATHS $GPULIBS
# the import library and export file link writes as some objects export symbols
rm -f ../bin/mmc.lib ../bin/mmc.exp

# the OpenMP runtime and the C++ runtime DLL it loads, from the MSVC
# redistributables (the LLVM OpenMP runtime of -openmp:llvm, libomp140, is
# not redistributable, and MSVC accepts no OpenMP 3 code without it)
for dll in vcomp140.dll vcruntime140.dll; do
    path=$(find "$(cygpath -u "$VCToolsRedistDir")" -ipath '*/x64/*' -iname $dll ! -ipath '*debug*' ! -ipath '*onecore*' 2>/dev/null | head -1)
    [ -n "$path" ] || { echo "$dll was not found in $VCToolsRedistDir" >&2; exit 1; }
    cp "$path" ../bin/
done

# the mex files and the Python module: -MD, the C++ runtime DLL MATLAB, Octave and Python use
if [ -n "$MATLAB_ROOT" ] || [ -n "$OCTAVE_HOME" ]; then
    objects built/mex -MD "-DMCX_CONTAINER -DMATLAB_MEX_FILE" $MMCSRC
fi

if [ -n "$MATLAB_ROOT" ]; then
    echo "== mmclab-trinity for MATLAB ($MATLAB_ROOT)"
    # built with cl and link as mex does, not by mex.bat: Git Bash runs a .bat by cmd with
    # its full path, and cmd splits a quoted "C:/Program Files/..." when more arguments are
    # quoted; c_mexapi_version exports the API version (R2017b, separate complex) MATLAB checks
    MLW=$(cygpath -m "$MATLAB_ROOT")
    mkdir -p built/mat
    cl $CLOPT -MD -std:c++14 $MMCDEFS -DMCX_CONTAINER -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD \
       -I"$MLW/extern/include" $MMCINC -Fo"built/mat/mmclab.obj" mmclab.cpp
    cl -nologo -c -O2 -MD -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2017b -I"$MLW/extern/include" \
       -Fo"built/mat/c_mexapi_version.obj" "$MLW/extern/version/c_mexapi_version.c"
    "$MSLINK" -nologo -DLL -EXPORT:mexFunction -EXPORT:mexfilerequiredapiversion -OUT:../mmclab/mmc.mexw64 \
              built/mat/*.obj built/mex/*.obj -LIBPATH:"$MLW/extern/lib/win64/microsoft" libmx.lib libmex.lib libmat.lib \
              $LIBPATHS $GPULIBS
    rm -f ../mmclab/mmc.lib ../mmclab/mmc.exp
fi

if [ -n "$OCTAVE_HOME" ]; then
    echo "== mmclab-trinity for Octave ($OCTAVE_HOME)"
    OCTW=$(cygpath -m "$OCTAVE_HOME")
    OCTINC=$(ls -d "$OCTW"/include/octave-*/octave | head -1)
    OCTMEX=$(ls "$OCTW"/bin/liboctmex*.dll | head -1)
    [ -f "$OCTMEX" ] || { echo "no liboctmex DLL in $OCTW/bin (Octave 10 or newer is needed)" >&2; exit 1; }
    # an import library of Octave's mex API, which MinGW built, from the DLL's exports
    {
        echo "LIBRARY $(basename "$OCTMEX")"
        echo "EXPORTS"
        "$MSBIN/dumpbin.exe" -nologo -exports "$OCTMEX" | awk 'f && NF >= 4 {print $4} /ordinal hint RVA/ {f = 1}'
    } > built/liboctmex.def
    "$MSBIN/lib.exe" -nologo -def:built/liboctmex.def -machine:x64 -out:built/liboctmex.lib
    mkdir -p built/oct
    cl $CLOPT -MD -std:c++14 $MMCDEFS -DMCX_CONTAINER -DMATLAB_MEX_FILE -I"$OCTINC" -I"$OCTINC/.." $MMCINC \
       -Fo"built/oct/mmclab.obj" mmclab.cpp
    # Octave 10 refuses a mex file without __octave_mex_soversion__ ("No SOVERSION found in .mex
    # file"), the version of liboctmex it was built for: the const int mkoctfile adds, here taken
    # from the DLL's name (liboctmex-1.dll)
    SOVER=$(basename "$OCTMEX" .dll | sed -n 's/^liboctmex-\([0-9]*\)$/\1/p')
    echo "__declspec(dllexport) const int __octave_mex_soversion__ = ${SOVER:-1};" > built/oct/soversion.c
    cl $CLOPT -MD -Fo"built/oct/soversion.obj" built/oct/soversion.c
    "$MSLINK" -nologo -DLL -EXPORT:mexFunction -OUT:../mmclab/mmc.mex built/oct/mmclab.obj built/oct/soversion.obj \
              built/mex/*.obj built/liboctmex.lib $LIBPATHS $GPULIBS
    rm -f ../mmclab/mmc.lib ../mmclab/mmc.exp
fi

if [ -n "$PYTHONS" ]; then
    echo "== pmmc-trinity"
    objects built/pyd -MD "-DMCX_CONTAINER -DPYBIND11_VERSION_MAJOR" $MMCSRC
    mkdir -p built/pyobj ../pmmc/dist
    rtmmc() {   # rtmmc <python>: _rtmmc for that Python, and a pmmc wheel with it
        local PY=$1 PYINC PYLIB EXT
        PYINC=$("$PY" -c "import sysconfig; print(sysconfig.get_paths()['include'])")
        PYLIB=$("$PY" -c "import sys, os; print(os.path.join(sys.base_prefix, 'libs'))")
        EXT=$("$PY" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
        echo "-- $PY ($EXT)"
        rm -rf built/pyobj/* built/whl
        cl $CLOPT -MD -std:c++14 -bigobj $MMCDEFS -DMCX_CONTAINER -DPYBIND11_VERSION_MAJOR -DPMMC_MODULE_NAME=_rtmmc \
           -I"$PYINC" -Ipybind11/include $MMCINC -Fo"built/pyobj/pmmc.obj" pmmc.cpp || return 1
        mkdir -p built/whl
        "$MSLINK" -nologo -DLL -OUT:"built/whl/_rtmmc$EXT" built/pyobj/pmmc.obj built/pyd/*.obj \
                  -LIBPATH:"$PYLIB" $LIBPATHS $GPULIBS || return 1
        rm -f built/whl/_rtmmc*.lib built/whl/_rtmmc*.exp
        cp ../bin/vcomp140.dll built/whl/
        "$PY" ../.github/windows_release_builder/pack_wheel.py built/whl ../pmmc ../pmmc/dist
    }
    # a failed Python version (e.g. one newer than the pinned pybind11 knows) skips its wheel only
    for PY in $PYTHONS; do
        rtmmc "$PY" || echo "pmmc-trinity: the wheel for $PY failed"
    done
fi

cd ..
dumpbin -nologo -dependents bin/*.exe bin/*.dll mmclab/mmc.mex* 2>/dev/null | grep -iE '\.dll|Dump of'
ls -l bin mmclab/mmc.mex* pmmc/dist 2>/dev/null || true
