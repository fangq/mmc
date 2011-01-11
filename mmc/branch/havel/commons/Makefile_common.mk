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
CCFLAGS    := -c -Wall -g -fno-strict-aliasing #-mfpmath=sse -ffast-math -mtune=core2
INCLUDEDIR := $(MMCDIR)/src
EXTRALIB   := -lm
AROUTPUT   := -o
MAKE       := make

OPENMP     := -fopenmp
FASTMATH   := #-ffast-math

ECHO	   := echo
MKDIR      := mkdir

DOXY       := doxygen
DOCDIR     := $(MMCDIR)/doc

ifeq ($(CC),icc)
	OPENMP   := -openmp
	FASTMATH :=
	EXTRALIB :=
endif

ARFLAGS    := $(EXTRALIB)

OBJSUFFIX  := .o
BINSUFFIX  := 

OBJS       := $(addprefix $(OBJDIR)/, $(FILES))
OBJS       := $(subst $(OBJDIR)/$(BXDSRC)/,$(BXDSRC)/,$(OBJS))
OBJS       := $(addsuffix $(OBJSUFFIX), $(OBJS))

TARGETSUFFIX:=$(suffix $(BINARY))

release:   CCFLAGS+= -O3
sse:       CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse4
sse omp:   CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
sse omp:   ARFLAGS+= $(OPENMP) $(FASTMATH)
prof:      CCFLAGS+= -O3 -pg
prof:      ARFLAGS+= -O3 -g -pg
dp dpomp:  CCFLAGS+= -mdouble-float -O3
dpomp:     CCFLAGS+= $(OPENMP)
dpomp:     ARFLAGS+= $(OPENMP) 

ifeq ($(TARGETSUFFIX),.so)
	CCFLAGS+= -fPIC 
	ARFLAGS+= -shared -Wl,-soname,$(BINARY).1 
endif

ifeq ($(TARGETSUFFIX),.a)
        CCFLAGS+=
	AR         := ar
        ARFLAGS    := r
	AROUTPUT   :=
endif

all release sse prof omp dp icc: $(SUBDIRS) $(BINDIR)/$(BINARY)

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
	$(AR)  $(ARFLAGS) $(AROUTPUT) $@ $(OBJS) $(USERARFLAGS)

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

