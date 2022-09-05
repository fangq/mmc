#!/bin/sh

nvcc -c -cubin -O3 -Xptxas -allow-expensive-optimizations -ptx --expt-relaxed-constexpr -arch=sm_60 -use_fast_math -maxrregcount=0 --use-local-env -Xcompiler  -m64 -Izmat -Izmat/easylzma -Iubj -I/usr/local/cuda/include -I/pub/optix-7.5/include -I/pub/optix-7.5/SDK -I/pub/optix-7.5/SDK/support  -o mmc_optix_core.ptx  mmc_optix_core.cu

g++ -c -Wall -g -DMCX_EMBED_CL -fno-strict-aliasing -m64 -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1 -O3 -fopenmp  -fPIC -DMCX_CONTAINER -DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT -DUSE_OPTIX -DNDEBUG -I../src -I../src/zmat/easylzma -I../src/ubj -I/usr/local/cuda/include -I/pub/optix-7.5/include -I/pub/optix-7.5/SDK -I/pub/optix-7.5/SDK/support  -o built/mmc_optix_host.o  mmc_optix_host.cpp

make mex
