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
    CUCCOPT= -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL -Xcompiler -fPIC $(USERCUCCOPT)
    EXTRALIB+=-lcudart
    LIBCUDART=-L$(LIBOPENCLDIR) -lcudart
endif

ifneq (,$(filter $(MAKECMDGOALS),optix optixmex optixoct))
    PTXSOURCE=mmc_optix_utils
    FILES+=$(PTXSOURCE) mmc_optix_host
    OPTIXPROGRAM=mmc_optix_core
    USERCCFLAGS+=-DUSE_OPTIX -DUSE_CUDA_VECTOR_TYPES -I$(CUDAHOME)/include -I$(OPTIXHOME)/include
    LIBCUDART+=-lcuda -ldl
    EXTRALIB+=-lcuda -ldl
endif

include $(ROOTDIR)/commons/Makefile_common.mk

## mmc_utils.c embeds the OpenCL kernel; recompile whenever mmc_core.cl changes ##
$(OBJDIR)/mmc_utils$(OBJSUFFIX): $(CLSOURCE)

##  Host unit tests (test/unit), linked with the same objects as mmc except its main() ##
##  and built with the flags of the default goal (ssemath, or omp on arm64 macOS)      ##
unittest: CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
unittest: ARFLAGS+= $(OPENMPLIB) $(FASTMATH)
ifeq ($(.DEFAULT_GOAL),ssemath)
unittest: CCFLAGS+=$(SSEFLAGS) -DUSE_SSE2 -DMMC_USE_SSE_MATH
endif
unittest: makedirs $(CLSOURCE) $(ZMATLIB) $(OBJS)
	$(CC) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -DMMC_EXAMPLE_DIR=\"$(abspath ../examples)\" -o $(OBJDIR)/test_mmc$(OBJSUFFIX) ../test/unit/test_mmc.c
	$(AR) $(ARFLAGS) $(AROUTPUT) $(BINDIR)/mmc_unittest $(OBJDIR)/test_mmc$(OBJSUFFIX) $(filter-out $(OBJDIR)/mmc$(OBJSUFFIX),$(OBJS)) $(USERARFLAGS) $(EXTRALIB)
	$(BINDIR)/mmc_unittest
