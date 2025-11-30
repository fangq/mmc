ROOTDIR = ..
BINARY=mmc

FILES=mmc_rand_xorshift128p mmc_mesh mmc_raytrace mmc_utils mmc_tictoc mmc_host mmc_highorder mmc_bench mmc_cl_utils mmc_cl_host
CLPROGRAM=mmc_core

DOXYCFG=mmc_doxygen.cfg

USERCCFLAGS=-DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT

DUMMY:=$(shell mkdir -p built/cjson)

ifneq (,$(filter $(MAKECMDGOALS),cuda cudamex cudaoct optix optixmex optixoct))
    FILES+=mmc_cu_host
    USERCCFLAGS+=-DUSE_CUDA
    CUCCOPT= -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL -Xcompiler -fPIC
    EXTRALIB+=-lcudart
    LIBCUDART=-L$(LIBOPENCLDIR) -lcudart
endif

ifneq (,$(filter $(MAKECMDGOALS),optix optixmex optixoct))
    PTXSOURCE=mmc_optix_utils
    FILES+=$(PTXSOURCE) mmc_optix_host
    OPTIXPROGRAM=mmc_optix_core
    USERCCFLAGS+=-DUSE_OPTIX -DUSE_CUDA_VECTOR_TYPES -I$(CUDAHOME)/include -I$(OPTIXHOME)/include -I$(OPTIXHOME)/SDK
    EXTRALIB+=-lcuda -ldl
endif

include $(ROOTDIR)/commons/Makefile_common.mk
