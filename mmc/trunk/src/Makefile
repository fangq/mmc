ROOTDIR = ..
BINARY=mmc

PLATFORM = $(shell uname -s)
ifeq ($(findstring MINGW32,$(PLATFORM)), MINGW32)
  FILES=posix_randr simpmesh tettracing mcx_utils tictoc drand48_r_libgw32c mmc cjson/cJSON
else
  FILES=posix_randr simpmesh tettracing mcx_utils tictoc mmc cjson/cJSON
endif

ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
  EXTRALIB:=-static-gcc
endif

DOXYCFG=mmcdoxy.cfg

USERCCFLAGS=-DUSE_OS_TIMER

DUMMY:=$(shell mkdir -p built/cjson)

include $(ROOTDIR)/commons/Makefile_common.mk

