ROOTDIR = ..
BINARY=mmc

FILES=mmc_rand_xorshift128p mmc_mesh mmc_raytrace mmc_utils mmc_tictoc mmc_host mmc_highorder mmc_bench mmc_cl_utils mmc_cl_host
CLPROGRAM=mmc_core

#PLATFORM = $(shell uname -s)

#ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
#  EXTRALIB:=-static-libgcc -lgcc_eh
#  MEXLINKOPT:=$(EXTRALIB)
#endif

DOXYCFG=mmc_doxygen.cfg

USERCCFLAGS=-DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT

DUMMY:=$(shell mkdir -p built/cjson)

ifneq (,$(filter $(MAKECMDGOALS),cuda cudamex cudaoct))
    FILES+=mmc_cu_host
    USERCCFLAGS+=-DUSE_CUDA
    CUCCOPT= -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL -Xcompiler -fPIC
    EXTRALIB+=-lcudart
    LIBCUDART=-L$(LIBOPENCLDIR) -lcudart
endif

include $(ROOTDIR)/commons/Makefile_common.mk
