ROOTDIR = ..
BINARY=mmc

FILES=posix_randr simpmesh tettracing mcx_utils tictoc mmc mmc_host cjson/cJSON

PLATFORM = $(shell uname -s)        
ifeq ($(findstring MINGW32,$(PLATFORM)), MINGW32)
  FILES+=drand48_r_libgw32c
endif
ifeq ($(findstring CYGWIN,$(PLATFORM)), CYGWIN)
  FILES+=drand48_r_libgw32c
endif

ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
  EXTRALIB:=-static-libgcc -lgcc_eh
  MEXLINKOPT:=$(EXTRALIB)
endif

DOXYCFG=mmcdoxy.cfg

USERCCFLAGS=-DUSE_OS_TIMER

DUMMY:=$(shell mkdir -p built/cjson)

include $(ROOTDIR)/commons/Makefile_common.mk

