ROOTDIR = ..
BINARY=mmc

PLATFORM = $(shell uname -s)
ifeq ($(findstring MINGW32,$(PLATFORM)), MINGW32)
  FILES=posix_randr simpmesh tettracing mcx_utils tictoc drand48_r_libgw32c tetray
else
  FILES=posix_randr simpmesh tettracing mcx_utils tictoc tetray
endif

ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
  EXTRALIB:=-static-gcc /opt/local/lib/gcc44/libgomp.a
endif

DOXYCFG=mmcdoxy.cfg

USERCCFLAGS=-DUSE_OS_TIMER

include $(ROOTDIR)/commons/Makefile_common.mk

