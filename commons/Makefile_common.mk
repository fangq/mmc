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
CCFLAGS    += -c -Wall -g -fno-strict-aliasing#-pedantic -std=c99 -mfpmath=sse -ffast-math -mtune=core2
INCLUDEDIR := $(MMCDIR)/src
EXTRALIB   += -lm
AROUTPUT   += -o
MAKE       := make

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
MKMEXOPT    =CC='$(CC)' CXX='$(CXX)' CXXFLAGS='$(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='-L$$TMW_ROOT$$MATLABROOT/sys/os/$$ARCH $$LDFLAGS $(MEXLINKOPT)' $(FASTMATH) -cxx -outdir $(BINDIR) -liomp5
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
endif

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

ARFLAGS    := 

OBJSUFFIX  := .o
BINSUFFIX  := 

OBJS       := $(addprefix $(OBJDIR)/, $(FILES))
OBJS       := $(subst $(OBJDIR)/$(MMCSRC)/,$(MMCSRC)/,$(OBJS))
OBJS       := $(addsuffix $(OBJSUFFIX), $(OBJS))

release:   CCFLAGS+= -O3
sse ssemath mex oct mexsse octsse: CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse4
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
oct octsse:     AR=CC=$(CC) CXX=$(CXX) LFLAGS='$(OPENMP)' LDFLAGS='$(LFLAGS)' CPPFLAGS='$(CCFLAGS) $(USERCCFLAGS) -std=c++11' $(USEROCTOPT) $(MKOCT)
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

all release sse ssemath prof omp mex oct mexsse octsse pnacl: $(SUBDIRS) $(BINDIR)/$(BINARY)

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

##  Link  ##
$(BINDIR)/$(BINARY): makedirs $(OBJS)
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
