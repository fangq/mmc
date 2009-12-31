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
ROOTDIR := .
endif

ifndef BXDDIR
BXDDIR := $(ROOTDIR)
endif

BXDSRC :=$(BXDDIR)/src

CXX        := g++
CC         := gcc
AR         := g++
BIN        := bin
BUILT      := built
BINDIR     := $(BIN)
OBJDIR 	   := $(BUILT)
CCFLAGS    := -c -Wall -g -O3
INCLUDEDIR := $(BXDDIR)/src
ARFLAGS    :=
AROUTPUT   := -o
MAKE       :=make

ECHO	 := echo
MKDIR    := mkdir

OBJSUFFIX        := .o
BINSUFFIX        := 

OBJS      := $(addprefix $(OBJDIR)/, $(FILES))
OBJS      := $(subst $(OBJDIR)/$(BXDSRC)/,$(BXDSRC)/,$(OBJS))
OBJS      := $(addsuffix $(OBJSUFFIX), $(OBJS))

TARGETSUFFIX:=$(suffix $(BINARY))

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

all: $(SUBDIRS) makedirs $(BINDIR)/$(BINARY)

$(SUBDIRS):
	$(MAKE) -C $@ --no-print-directory

makedirs:
	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); fi
	@if test ! -d $(BINDIR); then $(MKDIR) $(BINDIR); fi

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
$(BINDIR)/$(BINARY): $(OBJS)
	@$(ECHO) Building $@
	$(AR)  $(ARFLAGS) $(AROUTPUT) $@ $(OBJS) $(USERARFLAGS)

## Clean
clean:
	rm -rf $(OBJS) $(OBJDIR) $(BINDIR)
ifdef SUBDIRS
	for i in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$i clean; done
endif

.PHONY: regression clean arch makedirs dep $(SUBDIRS)

