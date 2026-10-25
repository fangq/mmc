#!/bin/bash
# Builds the Linux release binaries of mmc inside the manylinux2014 builder:
#
#   docker build --build-arg WITH_CUDA=0 -t mcx-linux-release-builder .github/linux_release_builder
#   docker run --rm -v $PWD:/src \
#       [-v $MATLAB_ROOT:$MATLAB_ROOT:ro -e MATLAB_ROOT=$MATLAB_ROOT] \
#       [-v /usr/include/octave-X.Y.Z:/octinc:ro -e OCTAVE_INCDIR=/octinc] \
#       mcx-linux-release-builder bash .github/linux_release_builder/build_mmc.sh
#
# MATLAB is the runner's (matlab-actions/setup-matlab), mounted read-only at
# the same path; the Octave headers are the runner's too.  Either may be left
# out to skip that mex file.  Every artifact is checked at the end: none may
# need a glibc symbol newer than 2.17, the oldest the image has, or link the
# GCC runtime libraries dynamically, or the build fails.

set -e
cd /src

# git refuses a checkout owned by another user, which a mounted one is
git config --global --add safe.directory '*' 2>/dev/null || true

export PATH=/src/.github/linux_release_builder:$PATH   # the mkoctfile stand-in
chmod +x /src/.github/linux_release_builder/mkoctfile

# the check at the end must only ever see what this run built
rm -f bin/mmc mmclab/mmc.mex mmclab/mmc.mexa64

cd src

if [ -n "$OCTAVE_INCDIR" ]; then
    echo "== mmclab for Octave (headers: $OCTAVE_INCDIR)"
    make clean
    make oct MEXLINKOPT="-static-libgcc -static-libstdc++ -fopenmp"
fi

if [ -n "$MATLAB_ROOT" ]; then
    echo "== mmclab for MATLAB ($MATLAB_ROOT)"
    make clean
    make mex MKMEX="$MATLAB_ROOT/bin/mex" MEXLINKOPT="-static-libgcc -static-libstdc++ -fopenmp"
    strip ../mmclab/mmc.mex*
fi

echo "== mmc"
make clean
# libstdc++, libgcc and libgomp go in statically, OpenCL and glibc stay shared;
# --as-needed drops the shared libgomp that -fopenmp adds after the static one
make AR=g++ EXTRALIB="-static-libstdc++ -static-libgcc -lOpenCL -Wl,-Bstatic -lgomp -Wl,-Bdynamic -lm -lpthread -ldl -Wl,--as-needed"

cd /src
bash .github/linux_release_builder/check_glibc.sh bin/mmc mmclab/mmc.mex mmclab/mmc.mexa64
