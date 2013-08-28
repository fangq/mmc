########################################################
#  MMC: Mesh-based Monte Carlo
#  Copyright (C) 2009 Qianqian Fang 
#                    <fangq at nmr.mgh.harvard.edu>
#
#  $Id$
########################################################

########################################################
# Base Makefile for all example/tests and main program
#
# This file specifies the compiler options for compiling
# and linking
########################################################

ifndef ROOTDIR
	ROOTDIR  := .
endif

ifndef MMCDIR
	MMCDIR   := $(ROOTDIR)
endif

BXDSRC :=$(MMCDIR)/src

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
FASTMATH   := #-ffast-math

ECHO	   := echo
MKDIR      := mkdir

MKMEX      :=mex
MKMEXOPT    =CC='$(CC)' CXX='$(CXX)' CXXFLAGS='$(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='$$LDFLAGS $(OPENMP)' $(FASTMATH) -cxx
MKOCT      :=mkoctfile

DLLFLAG=-fPIC

ARCH = $(shell uname -m)
PLATFORM = $(shell uname -s)
ifeq ($(findstring MINGW32,$(PLATFORM)), MINGW32)
    MKMEX      :=cmd //c mex.bat
    MKMEXOPT    =COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LINKFLAGS='$$LINKFLAGS $(OPENMP)' $(FASTMATH)
    DLLFLAG     =
endif

DOXY       := doxygen
DOCDIR     := $(MMCDIR)/doc

ifeq ($(CC),icc)
	OPENMP   := -openmp
	FASTMATH :=
	EXTRALIB :=
endif

ARFLAGS    := 

OBJSUFFIX  := .o
BINSUFFIX  := 

OBJS       := $(addprefix $(OBJDIR)/, $(FILES))
OBJS       := $(subst $(OBJDIR)/$(BXDSRC)/,$(BXDSRC)/,$(OBJS))
OBJS       := $(addsuffix $(OBJSUFFIX), $(OBJS))

TARGETSUFFIX:=$(suffix $(BINARY))

release:   CCFLAGS+= -O3
sse ssemath mexsse octsse: CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse4
sse ssemath omp mex oct mexsse octsse:   CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
sse ssemath omp:   ARFLAGS+= $(OPENMP) $(FASTMATH)
ssemath mexsse octsse:   CCFLAGS+= -DUSE_SSE2 -DMMC_USE_SSE_MATH
mex mexsse:        ARFLAGS+=$(MKMEXOPT)
prof:      CCFLAGS+= -O3 -pg
prof:      ARFLAGS+= -O3 -g -pg

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
oct octsse:     AR=CC=$(CC) CXX=$(CXX) LDFLAGS='$(OPENMP)' CPPFLAGS='$(CCFLAGS) $(USERCCFLAGS)' $(MKOCT)


ifeq ($(TARGETSUFFIX),.so)
	CCFLAGS+= $(DLLFLAG) 
	ARFLAGS+= -shared -Wl,-soname,$(BINARY).1 
endif

ifeq ($(TARGETSUFFIX),.a)
        CCFLAGS+=
	AR         := ar
        ARFLAGS    := r
	AROUTPUT   :=
endif

all release sse ssemath prof omp mex oct mexsse octsse: $(SUBDIRS) $(BINDIR)/$(BINARY)

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

