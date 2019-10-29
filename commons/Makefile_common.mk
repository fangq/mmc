########################################################
#  MMC: Mesh-based Monte Carlo
#  Copyright (C) 2009 Qianqian Fang 
#                    <q.fang at neu.edu>
#
#  $Id$
########################################################

########################################################
# Base Makefile for all example/tests and main program
#
# This file specifies the compiler options for compiling
# and linking
########################################################

ROOTDIR ?= .
MMCDIR  ?= $(ROOTDIR)

MMCSRC :=$(MMCDIR)/src

CXX        := g++
AR         := $(CC)
BIN        := bin
BUILT      := built
BINDIR     := $(BIN)
OBJDIR 	   := $(BUILT)
CCFLAGS    += -c -Wall -g -DMCX_EMBED_CL -fno-strict-aliasing#-pedantic -std=c99 -mfpmath=sse -ffast-math -mtune=core2
INCLUDEDIR := $(MMCDIR)/src
AROUTPUT   += -o
MAKE       := make

LIBOPENCLDIR ?= /usr/local/cuda/lib64
LIBOPENCL=-lOpenCL
EXTRALIB   += -lm -lstdc++ -L$(LIBOPENCLDIR)

OPENMP     := -fopenmp
OPENMPLIB  := -fopenmp
FASTMATH   := #-ffast-math

ECHO	   := echo
MKDIR      := mkdir

ARCH = $(shell uname -m)
ifeq ($(findstring x86_64,$(ARCH)), x86_64)
     CCFLAGS+=-m64
endif

MEXLINKOPT +=$(OPENMPLIB)
MKMEX      :=mex
MKMEXOPT    =CC='$(CC)' CXX='$(CXX)' CXXLIBS='$$CXXLIBS $(LIBOPENCL)' CXXFLAGS='$(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='-L$$TMW_ROOT$$MATLABROOT/sys/os/$$ARCH $$LDFLAGS $(MEXLINKOPT)' $(FASTMATH) -cxx -outdir $(BINDIR)
MKOCT      :=mkoctfile

DLLFLAG=-fPIC

PLATFORM = $(shell uname -s)
ifeq ($(findstring MINGW32,$(PLATFORM)), MINGW32)
    MKMEX      :=cmd //c mex.bat
    MKMEXOPT    =COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LINKFLAGS='$$LINKFLAGS $(OPENMPLIB) $(MEXLINKOPT)' $(FASTMATH)
    DLLFLAG     =
endif
ifeq ($(findstring CYGWIN,$(PLATFORM)), CYGWIN)
    MKMEX      :=cmd /c mex.bat
    MKMEXOPT    =-f mexopts_cygwin64_gcc.bat COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LINKFLAGS='$$LINKFLAGS $(OPENMPLIB) $(MEXLINKOPT)' $(FASTMATH) -outdir ../mmclab
    DLLFLAG     =
else ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
  INCLUDEDIRS=-I/System/Library/Frameworks/OpenCL.framework/Headers
  LIBOPENCL=-framework OpenCL
  LIBOPENCLDIR=/System/Library/Frameworks/OpenCL.framework/Versions/A
  OPENMPLIB=-static-libgcc /usr/local/lib/libgomp.a
endif

INCLUDEDIR+=$(INCLUDEDIRS)
EXTRALIB+=$(LIBOPENCL)

NACL_SDK_ROOT ?= ../../../nacl
OSNAME := $(shell echo $(PLATFORM) | tr A-Z a-z)
PNACL_TC_PATH := $(abspath $(NACL_SDK_ROOT)/toolchain/$(OSNAME)_pnacl)

DOXY       := doxygen
DOCDIR     := $(MMCDIR)/doc

ifeq ($(CC),icc)
	OPENMP   := -qopenmp
	OPENMPLIB:= $(OPENMP)
	FASTMATH :=
	EXTRALIB :=
endif

ifeq ($(CC),clang)
	OPENMP   := -fopenmp
        OPENMPLIB:= -fopenmp=libiomp5
endif

ARFLAGS    := 

OBJSUFFIX  := .o
BINSUFFIX  := 
CLHEADER=.clh

OBJS       := $(addprefix $(OBJDIR)/, $(FILES))
OBJS       := $(subst $(OBJDIR)/$(MMCSRC)/,$(MMCSRC)/,$(OBJS))
OBJS       := $(addsuffix $(OBJSUFFIX), $(OBJS))
CLSOURCE  := $(addsuffix $(CLHEADER), $(CLPROGRAM))

release:   CCFLAGS+= -O3
sse ssemath mex oct mexsse octsse: CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1
sse ssemath omp mex oct mexsse octsse:   CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
sse ssemath omp:   ARFLAGS+= $(OPENMPLIB) $(FASTMATH)
ssemath mexsse octsse:   CCFLAGS+=-DUSE_SSE2 -DMMC_USE_SSE_MATH
mex mexsse:        ARFLAGS+=$(MKMEXOPT)
prof:      CCFLAGS+= -O3 -pg
prof:      ARFLAGS+= -O3 -g -pg

pnacl:     CC=$(PNACL_TC_PATH)/bin/pnacl-clang++
pnacl:     AR=$(PNACL_TC_PATH)/bin/pnacl-ar
pnacl:	   ARFLAGS= cr
pnacl:	   EXTRALIB   :=
pnacl:     INCLUDEDIR+= -I$(NACL_SDK_ROOT)/include/pnacl
pnacl:     BINARY=libmmc-pnacl.a

web: CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3
web: CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
web: ARFLAGS+= $(OPENMPLIB) $(FASTMATH) -DUSE_SSE2 -DMMC_USE_SSE_MATH
web: CFLAGS+=-D__SSE__ -D__SSE2__
web: CC=emcc
web: BINDIR:=webmmc
web: AR=emcc
web: EXTRALIB=-s SIMD=1 -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='["cwrap"]' -s FORCE_FILESYSTEM=1 -o $(BINDIR)/webmmc.html

mex oct mexsse octsse:   EXTRALIB=
mex oct mexsse octsse:   CCFLAGS+=$(DLLFLAG) -DMCX_CONTAINER
mex oct mexsse octsse:   CPPFLAGS+=-g $(DLLFLAG) -DMCX_CONTAINER
mex oct mexsse octsse:   BINDIR=../mmclab
mex mexsse:     AR=$(MKMEX)
mex mexsse:     AROUTPUT=-output
mex mexsse:     ARFLAGS+=mmclab.cpp -I$(INCLUDEDIR)
mexsse:         BINARY=mmc_sse

oct:            BINARY=mmc.mex
octsse:         BINARY=mmc_sse.mex
oct octsse:     ARFLAGS+=--mex mmclab.cpp -I$(INCLUDEDIR)
oct octsse:     AR=CC=$(CC) CXX=$(CXX) LDFLAGS='$(LFLAGS)' CPPFLAGS='$(CCFLAGS) $(USERCCFLAGS) -std=c++11' $(USEROCTOPT) $(MKOCT)
oct octsse:     USERARFLAGS=-o $(BINDIR)/mmc

TARGETSUFFIX:=$(suffix $(BINARY))

ifeq ($(TARGETSUFFIX),.so)
	CCFLAGS+= $(DLLFLAG) 
	ARFLAGS+= -shared -Wl,-soname,$(BINARY).1 
endif

ifeq ($(TARGETSUFFIX),.a)
        CCFLAGS+=
	AR         := ar
	ARFLAGS    := cr
	AROUTPUT   :=
	EXTRALIB   :=
	OPENMPLIB  :=
endif

all release sse ssemath prof omp mex oct mexsse octsse pnacl web: $(SUBDIRS) $(BINDIR)/$(BINARY)

$(SUBDIRS):
	$(MAKE) -C $@ --no-print-directory

makedirs:
	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); fi
	@if test ! -d $(BINDIR); then $(MKDIR) $(BINDIR); fi

makedocdir:
	@if test ! -d $(DOCDIR); then $(MKDIR) $(DOCDIR); fi

.SUFFIXES : $(OBJSUFFIX) .cpp

##  Compile .cpp files ##
$(OBJDIR)/%$(OBJSUFFIX): %.cpp
	@$(ECHO) Building $@
	$(CXX) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

##  Compile .cpp files ##
%$(OBJSUFFIX): %.cpp
	@$(ECHO) Building $@
	$(CXX) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

##  Compile .c files  ##
$(OBJDIR)/%$(OBJSUFFIX): %.c
	@$(ECHO) Building $@
	$(CC) $(CCFLAGS) $(USERCCFLAGS) -I$(INCLUDEDIR) -o $@  $<

%$(CLHEADER): %.cl
	xxd -i $(CLPROGRAM).cl | sed 's/\([0-9a-f]\)$$/\0, 0x00/' > $(CLPROGRAM).clh

##  Link  ##
$(BINDIR)/$(BINARY): makedirs $(CLSOURCE) $(OBJS)
	@$(ECHO) Building $@
	$(AR)  $(ARFLAGS) $(AROUTPUT) $@ $(OBJS) $(USERARFLAGS) $(EXTRALIB)

##  Documentation  ##
doc: makedocdir
	$(DOXY) $(DOXYCFG)

## Clean
clean:
	rm -rf $(OBJS) $(OBJDIR) $(BINDIR) $(DOCDIR)
ifdef SUBDIRS
	for i in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$i clean; done
endif

.PHONY: regression clean arch makedirs dep $(SUBDIRS)

.DEFAULT_GOAL := sse
