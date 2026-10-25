#!/bin/bash
# Builds mmc with the CUDA and OptiX backends on Windows, the Windows twin of
# .github/linux_release_builder/build_mmc_optix.sh.  On Windows nvcc takes no
# host compiler other than MSVC, so the whole binary is built with MSVC, from
# a Git Bash where the MSVC environment is set up (e.g. ilammy/msvc-dev-cmd
# with toolset 14.29: the STL of newer toolsets refuses CUDA 11.7):
#
#   OPTIXHOME=/path/to/optix-sdk bash .github/windows_release_builder/build_mmc_optix.sh
#
# As on Linux, CUDA 11.7 comes from NVIDIA's redistributable archives (nvcc,
# cudart, cccl; cudart also carries the CL headers and OpenCL.lib), and
# OPTIXHOME holds the OptiX 7.5 SDK headers (include/optix*.h), e.g. a clone of
# https://github.com/NVIDIA/optix-dev at v7.5.0.  The objects are those of
# "make optix" (sse, opencl, cuda and optix backends), with the same defines.
#
# bin/mmc.exe has the C and CUDA runtimes linked in statically; it loads
# nvcuda.dll (the driver, at least R516 for CUDA 11.7) and OpenCL.dll when a
# GPU backend is used, and the OpenMP runtime of MSVC, vcomp140.dll, copied
# next to it from the redistributables with the vcruntime140.dll it needs.

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
# Git's own link.exe (coreutils) may come first on the PATH
MSLINK="$(dirname "$(command -v cl)")/link.exe"

cd src

# Git for Windows checks symbolic links out as small text files holding the target
for l in mmc_core.cu zmat/miniz/zlib.h; do
    if [ ! -L "$l" ] && [ "$(wc -c < "$l")" -lt 64 ]; then
        cp "$(dirname "$l")/$(cat "$l")" "$l"
    fi
done

mkdir -p built ../bin
rm -f built/*.obj ../bin/mmc.exe

# the OpenCL kernel and the OptiX PTX, as 0-terminated char arrays
xxd -i mmc_core.cl | sed 's/\([0-9a-f]\)$/\0, 0x00/' > mmc_core.clh

USERDEFS="-DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT -DUSE_CUDA -DUSE_OPTIX -DUSE_CUDA_VECTOR_TYPES"
GPUINC="-I$CUDAW/include -I$OPTIXW/include"

nvcc -O3 -ptx --expt-relaxed-constexpr -use_fast_math -m64 $USERDEFS $GPUINC -o built/mmc_optix_core.ptx mmc_optix_core.cu
(cd built && xxd -i mmc_optix_core.ptx | sed 's/\([0-9a-f]\)$/\0, 0x00/' > mmc_optix_core_ptx.h)

# -MT: the static C runtime, as the CUDA runtime (cudart_static) is linked
CLOPT="-nologo -c -O2 -MT -openmp -EHsc -W1 -D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES"
MMCDEFS="-DMCX_EMBED_CL $USERDEFS -DMMC_USE_SSE -DHAVE_SSE2 -DUSE_SSE2 -DMMC_USE_SSE_MATH"

# zmat, as its Makefile builds it for mmc: miniz for zlib, the LZMA SDK, lz4; no
# blosc2/zstd; with MSVC its threads (block-parallel zlib, LZMA) are Win32's
ZMATDEFS="-DNO_BLOSC2 -DNO_ZSTD -DNO_ZLIB -D_LARGEFILE64_SOURCE=1 -DZMAT_USE_LZMA_SDK"
ZMATINC="-Izmat -Izmat/miniz -Izmat/easylzma -Izmat/easylzma/lzma -Izmat/lz4"
for f in zmatlib miniz/miniz easylzma/compress easylzma/decompress easylzma/lzma_header \
         easylzma/lzip_header easylzma/common_internal lz4/lz4 lz4/lz4hc \
         easylzma/lzma/{LzmaEnc,LzmaDec,LzmaLib,LzFind,Bra,BraIA64,Alloc,7zCrc,7zCrcOpt,CpuArch} \
         easylzma/lzma/{Lzma2Dec,Lzma2DecMt,Lzma2Enc,MtCoder,MtDec,Xz,XzCrc64,XzCrc64Opt,XzEnc,XzDec} \
         easylzma/lzma/{Sha256,Sha256Opt,Delta,Bra86,7zStream,LzFindMt,LzFindOpt,Threads}; do
    cl $CLOPT $ZMATDEFS $ZMATINC -Fo"built/zmat_$(echo $f | tr / _).obj" "zmat/$f.c"
done

MMCINC="-I. -Ibuilt -Izmat -Izmat/easylzma -Iubj $GPUINC"
for f in mmc_rand_xorshift128p mmc_mesh mmc_raytrace mmc_utils mmc_tictoc mmc_host \
         mmc_bench mmc_cl_utils mmc_cl_host mmc cjson/cJSON ubj/ubjw; do
    cl $CLOPT $MMCDEFS $MMCINC -Fo"built/$(basename $f).obj" "$f.c"
done
for f in mmc_highorder mmc_neurojson mmc_optix_utils mmc_optix_host; do
    cl $CLOPT -std:c++14 $MMCDEFS $MMCINC -Fo"built/$f.obj" "$f.cpp"
done

# the CUDA backend: the flags of "make cuda" (CUCCOPT)
nvcc -c -m64 -use_fast_math -Xptxas -O3 -Xcompiler -MT,-openmp -D_USE_MATH_DEFINES \
     -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL \
     -o built/mmc_cu_host.obj mmc_cu_host.cu

# the driver and OpenCL DLLs are delay-loaded: the sse backend then runs (and
# mmc prints its help) on a computer with neither
"$MSLINK" -nologo -OUT:../bin/mmc.exe built/*.obj -LIBPATH:"$CUDAW/lib/x64" \
          cudart_static.lib cuda.lib OpenCL.lib delayimp.lib -DELAYLOAD:nvcuda.dll -DELAYLOAD:OpenCL.dll
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

cd ..
dumpbin -nologo -dependents bin/*.exe bin/*.dll | grep -i '\.dll'
