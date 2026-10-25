#!/bin/bash
# Builds mmc with the CUDA and OptiX backends inside the manylinux2014 builder:
#
#   docker build --build-arg WITH_CUDA=0 -t mcx-linux-release-builder .github/linux_release_builder
#   docker run --rm -v $PWD:/src -v /path/to/optix-sdk:/optix:ro \
#       mcx-linux-release-builder bash .github/linux_release_builder/build_mmc_optix.sh
#
# OptiX 7.5 needs CUDA 11, newer than the CUDA 10.2 the image can carry for mcx;
# only nvcc, the runtime and cccl are needed, so NVIDIA's redistributable
# archives of CUDA 11.7 (about 40 MB) are fetched rather than the full
# installer.  /optix holds the OptiX SDK headers (include/optix*.h), e.g. a
# clone of https://github.com/NVIDIA/optix-dev at v7.5.0.
#
# The binary, bin/mmc, has the sse, opencl, cuda and optix backends.  As with
# the other release builds it needs at most GLIBC_2.17; libstdc++, libgcc,
# libgomp and the CUDA runtime are linked statically, and libcuda.so.1 (the
# driver, at least R515 for CUDA 11.7) and libOpenCL.so.1 are the only GPU
# libraries it loads.

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

# the check at the end must only ever see what this run built
rm -f bin/mmc

cd src
make clean
# libstdc++, libgcc, libgomp and the CUDA runtime go in statically; the CUDA
# driver (libcuda, from the stubs at link time), OpenCL and glibc stay shared;
# --as-needed drops the shared libgomp that -fopenmp adds after the static one
make optix AR=g++ CUDAHOME=$CUDAHOME OPTIXHOME=/optix \
    EXTRALIB="-static-libstdc++ -static-libgcc -L$CUDAHOME/lib -L$CUDAHOME/lib/stubs -lcudart_static -lcuda -lOpenCL \
              -Wl,-Bstatic -lgomp -Wl,-Bdynamic -lm -lpthread -ldl -lrt -Wl,--as-needed"

cd /src
bash .github/linux_release_builder/check_glibc.sh bin/mmc
