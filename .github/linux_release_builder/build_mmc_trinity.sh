#!/bin/bash
# Builds the mmc-trinity packages, mmc with the CUDA and OptiX backends next to
# the sse and opencl ones, inside the manylinux2014 builder:
#
#   docker build --build-arg WITH_CUDA=0 -t mcx-linux-release-builder .github/linux_release_builder
#   docker run --rm -v $PWD:/src -v /path/to/optix-sdk:/optix:ro \
#       [-v $MATLAB_ROOT:$MATLAB_ROOT:ro -e MATLAB_ROOT=$MATLAB_ROOT] \
#       [-v /usr/include/octave-X.Y.Z:/octinc:ro -e OCTAVE_INCDIR=/octinc] \
#       mcx-linux-release-builder bash .github/linux_release_builder/build_mmc_trinity.sh
#
# OptiX 7.5 needs CUDA 11, newer than the CUDA 10.2 the image can carry for mcx;
# only nvcc, the runtime and cccl are needed, so NVIDIA's redistributable
# archives of CUDA 11.7 (about 40 MB) are fetched rather than the full
# installer.  /optix holds the OptiX SDK headers (include/optix*.h), e.g. a
# clone of https://github.com/NVIDIA/optix-dev at v7.5.0.
#
# It builds
#   bin/mmc                     mmc-trinity: the sse, opencl, cuda and optix backends
#   mmclab/mmc.mexa64, mmc.mex  mmclab-trinity for MATLAB and Octave (when given)
#   pmmc/dist/*.whl             pmmc-trinity: pmmc wheels with _rtmmc, the cuda and
#                               optix backends, next to the opencl-only _pmmc
# The CUDA code is built for every GPU generation from sm_52 to sm_86, with the
# sm_86 PTX for newer GPUs (CUGENCODE in commons/Makefile_common.mk).  As with
# the other release builds, everything needs at most GLIBC_2.17; libstdc++,
# libgcc, libgomp (but in the mex files) and the CUDA runtime are linked
# statically, and libcuda.so.1 (the driver, at least R515 for CUDA 11.7) and
# libOpenCL.so.1 are the only GPU libraries loaded.

set -e
cd /src

git config --global --add safe.directory '*' 2>/dev/null || true

CUDAVER=11.7.99
CCCLVER=11.7.91
REDIST=https://developer.download.nvidia.com/compute/cuda/redist
CUDAHOME=/tmp/cuda-$CUDAVER

if [ ! -x "$CUDAHOME/bin/nvcc" ]; then
    mkdir -p "$CUDAHOME" /tmp/cudadl
    for pkg in "cuda_nvcc:$CUDAVER" "cuda_cudart:$CUDAVER" "cuda_cccl:$CCCLVER"; do
        name=${pkg%%:*}
        ver=${pkg##*:}
        wget -q -O /tmp/cudadl/$name.tar.xz "$REDIST/$name/linux-x86_64/$name-linux-x86_64-$ver-archive.tar.xz"
        tar -xJf /tmp/cudadl/$name.tar.xz -C "$CUDAHOME" --strip-components=1
    done
    rm -rf /tmp/cudadl
fi
export PATH=$CUDAHOME/bin:$PATH

[ -f /optix/include/optix.h ] || { echo "the OptiX SDK headers must be mounted at /optix" >&2; exit 1; }

export PATH=/src/.github/linux_release_builder:$PATH   # the mkoctfile stand-in
chmod +x /src/.github/linux_release_builder/mkoctfile

# the check at the end must only ever see what this run built
rm -rf bin/mmc mmclab/mmc.mex mmclab/mmc.mexa64 pmmc/dist pmmc/wheels lib/librtmmc-cudahost.a

# the CUDA runtime goes in statically, the driver (libcuda, from the stubs at link time) stays shared
CUDALIBS="-L$CUDAHOME/lib -L$CUDAHOME/lib/stubs -lcudart_static -lcuda -ldl -lrt"

cd src

if [ -n "$OCTAVE_INCDIR" ]; then
    echo "== mmclab-trinity for Octave (headers: $OCTAVE_INCDIR)"
    make clean
    make optixoct CUDAHOME=$CUDAHOME OPTIXHOME=/optix LIBCUDART="$CUDALIBS" MEXLINKOPT="-static-libgcc -static-libstdc++ -fopenmp"
fi

if [ -n "$MATLAB_ROOT" ]; then
    echo "== mmclab-trinity for MATLAB ($MATLAB_ROOT)"
    make clean
    make optixmex CUDAHOME=$CUDAHOME OPTIXHOME=/optix LIBCUDART="$CUDALIBS" MKMEX="$MATLAB_ROOT/bin/mex" \
        MEXLINKOPT="-static-libgcc -static-libstdc++ -fopenmp"
    strip ../mmclab/mmc.mex*
fi

echo "== mmc-trinity"
make clean
# libstdc++, libgcc and libgomp go in statically too; OpenCL and glibc stay shared;
# --as-needed drops the shared libgomp that -fopenmp adds after the static one
make optix AR=g++ CUDAHOME=$CUDAHOME OPTIXHOME=/optix \
    EXTRALIB="-static-libstdc++ -static-libgcc $CUDALIBS -lOpenCL \
              -Wl,-Bstatic -lgomp -Wl,-Bdynamic -lm -lpthread -Wl,--as-needed"

echo "== pmmc-trinity"
cd /src/pmmc
PYS=$(ls -d /opt/python/cp3*/bin | grep -v 't/bin')   # CPython; not the free-threaded builds
RTARGS="-DBUILD_OPTIX=ON	-DOPTIX_ROOT=/optix	-DCUDA_TOOLKIT_ROOT_DIR=$CUDAHOME"
# the CUDA code of _rtmmc, built once for all Python versions
PY0=$(echo "$PYS" | tail -1)
cmake -S ../src -B /tmp/rtmmc-cuda -DBUILD_PYTHON=ON -DBUILD_MEX=OFF -DBUILD_OPTIX=ON -DOPTIX_ROOT=/optix \
      -DCUDA_TOOLKIT_ROOT_DIR=$CUDAHOME -DPYTHON_EXECUTABLE=$PY0/python -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build /tmp/rtmmc-cuda --target rtmmc-cudahost -j "$(nproc)"
RTLIB=$(find /src/lib /tmp/rtmmc-cuda -name 'librtmmc-cudahost.a' | head -1)
[ -f "$RTLIB" ] || { echo "librtmmc-cudahost.a was not built" >&2; exit 1; }
cp "$RTLIB" /tmp/librtmmc-cudahost.a
for PYBIN in $PYS; do
    echo "-- $PYBIN"
    CMAKE_ARGS="$RTARGS	-DRTMMC_CUDAHOST_LIB=/tmp/librtmmc-cudahost.a" "$PYBIN/pip" wheel . --no-deps -w wheels/ \
        || echo "pmmc-trinity: the wheel for $PYBIN failed"
done
# libOpenCL is bundled into the wheels, as in the PyPI ones; libcuda is the driver's
for WHEEL in wheels/*.whl; do
    auditwheel repair --exclude libcuda.so.1 "$WHEEL" -w dist/
done
rm -rf wheels build

cd /src
bash .github/linux_release_builder/check_glibc.sh bin/mmc mmclab/mmc.mex mmclab/mmc.mexa64
ls -l pmmc/dist
