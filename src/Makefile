ROOTDIR = ..
BINARY=mmc

FILES=xorshift128p_rand simpmesh tettracing mcx_utils tictoc mmc cjson/cJSON mmc_host highordermesh mmc_cl_utils mmc_cl_host ubj/ubjw
CLPROGRAM=mmc_core

ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
  EXTRALIB:=-static-libgcc -lgcc_eh
  MEXLINKOPT:=$(EXTRALIB)
endif

DOXYCFG=mmcdoxy.cfg

USERCCFLAGS=-DUSE_OS_TIMER -DUSE_OPENCL -DMMC_XORSHIFT

DUMMY:=$(shell mkdir -p built/cjson)

ifeq ($(MAKECMDGOALS),cuda)
    FILES+=mmc_cu_host
    USERCCFLAGS+=-DUSE_CUDA
    CUCCOPT= -DUSE_ATOMIC -DMCX_SAVE_DETECTORS -DMCX_DO_REFLECTION -DUSE_DMMC -DUSE_BLBADOUEL
    EXTRALIB+=-lcudart
endif

include $(ROOTDIR)/commons/Makefile_common.mk

