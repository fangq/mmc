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
CUDACC     :=nvcc
BIN        := bin
BUILT      := built
BINDIR     := $(BIN)
OBJDIR 	   := $(BUILT)
CCFLAGS    += -c -Wall -g -DMCX_EMBED_CL -fno-strict-aliasing#-pedantic -std=c99 -mfpmath=sse -ffast-math -mtune=core2
INCLUDEDIR := $(MMCDIR)/src -I$(MMCDIR)/src/zmat/easylzma -I$(MMCDIR)/src/ubj
AROUTPUT   += -o
MAKE       ?= make

ZMATLIB    :=libzmat.a
USERARFLAGS?=$(ZMATLIB) -lz

LIBOPENCLDIR ?= /usr/local/cuda/lib64
LIBOPENCL  ?=-lOpenCL
EXTRALIB   += -lm -lstdc++ -L$(LIBOPENCLDIR)

OPENMP     := -fopenmp
OPENMPLIB  := -fopenmp
FASTMATH   := #-ffast-math
CUCCOPT    +=-Xcompiler $(OPENMP) -use_fast_math -Xptxas -O3,-v
CUDA_STATIC=--cudart static -Xcompiler "-static-libgcc -static-libstdc++"

ECHO	   := echo
MKDIR      := mkdir

ARCH = $(shell uname -m)
ifeq ($(findstring x86_64,$(ARCH)), x86_64)
     CCFLAGS+=-m64
endif

MEXLINKOPT +=$(OPENMPLIB)
MKMEX      :=mex
MKMEXOPT    =CC='$(CC)' CXX='$(CXX)' CXXLIBS='$$CXXLIBS $(LIBOPENCL)' CXXFLAGS='$(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='-L$$TMW_ROOT$$MATLABROOT/sys/os/$$ARCH $$LDFLAGS $(MEXLINKOPT)' $(FASTMATH) -cxx -outdir $(BINDIR)
MKOCT      :=mkoctfile -v

DLLFLAG=-fPIC

PLATFORM = $(shell uname -s)
ifeq ($(findstring MINGW64,$(PLATFORM)), MINGW64)
    MW_MINGW64_LOC=/c/msys64/usr/
    MKMEX      :=cmd //c mex.bat
    INCLUDEDIRS+=-I"./mingw64/include"
    LIBOPENCL   ="c:\Windows\System32\OpenCL.dll"
    MKMEXOPT    =-f mexopts_msys2_gcc.xml COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='$$LDFLAGS -static $(OPENMPLIB) $(LIBOPENCL) $(MEXLINKOPT)' $(FASTMATH) -outdir ../mmclab
    EXTRALIB   +=-static
    CCFLAGS    +=-D__USE_MINGW_ANSI_STDIO=1
    DLLFLAG    =
else ifeq ($(findstring MSYS,$(PLATFORM)), MSYS)
    MKMEX      :=cmd //c mex.bat
    INCLUDEDIRS+=-I"./mingw64/include"
    LIBOPENCL   ="c:\Windows\System32\OpenCL.dll"
    MKMEXOPT    =-f mexopts_msys2_gcc.xml COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='$$LDFLAGS -static $(OPENMPLIB) $(LIBOPENCL) $(MEXLINKOPT)' $(FASTMATH) -outdir ../mmclab
    EXTRALIB   +=-static
    CCFLAGS    +=-D__USE_MINGW_ANSI_STDIO=1
    DLLFLAG    =
else ifeq ($(findstring CYGWIN,$(PLATFORM)), CYGWIN)
    MKMEX      :=cmd /c mex.bat
    MKMEXOPT    =-f mexopts_msys2_gcc.xml COMPFLAGS='$$COMPFLAGS $(CCFLAGS) $(USERCCFLAGS)' LDFLAGS='$$LDFLAGS -static $(OPENMPLIB) $(LIBOPENCL) $(MEXLINKOPT)' $(FASTMATH) -outdir ../mmclab
    LIBOPENCL   ="c:\Windows\System32\OpenCL.dll"
    INCLUDEDIRS+=-I"./mingw64/include"
    EXTRALIB   +=-static
    CCFLAGS    +=-D__USE_MINGW_ANSI_STDIO=1
    DLLFLAG     =
else ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
    INCLUDEDIRS=-I/System/Library/Frameworks/OpenCL.framework/Headers
    LIBOPENCL=-framework OpenCL
    LIBOPENCLDIR=/System/Library/Frameworks/OpenCL.framework/Versions/A
    OPENMPLIB=-static-libgcc /usr/local/lib/libgomp.a
    CUDA_STATIC=--cudart static
endif

ifeq ($(BACKEND),ocelot)
  LINKOPT=-L/usr/local/lib `OcelotConfig -l` -ltinfo
  CUCCOPT=-D__STRICT_ANSI__ -g #--maxrregcount 32
else ifeq ($(BACKEND),cudastatic)
  ifeq ($(findstring Darwin,$(PLATFORM)), Darwin)
      CUDART=-lcudadevrt -lcudart_static -ldl -static-libgcc -static-libstdc++
  else
      CUDART=-lcudadevrt -lcudart_static -ldl -lrt -static-libgcc -static-libstdc++
  endif
  LINKOPT=-L/usr/local/cuda/lib -lm $(CUDART)
  CUCCOPT+=-g -lineinfo -Xcompiler -Wall#-arch compute_20 #--maxrregcount 32
else
  LINKOPT=-L/usr/local/cuda/lib -lm $(CUDART)
  CUCCOPT+=-g -lineinfo -Xcompiler -Wall#-arch compute_20 #--maxrregcount 32
endif

CUGENCODE?=-arch=sm_30

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
sse ssemath mex oct: CCFLAGS+= -DMMC_USE_SSE -DHAVE_SSE2 -msse -msse2 -msse3 -mssse3 -msse4.1
sse ssemath omp mex oct mexomp octomp:   CCFLAGS+= -O3 $(OPENMP) $(FASTMATH)
sse ssemath omp:   ARFLAGS+= $(OPENMPLIB) $(FASTMATH)
ssemath:   CCFLAGS+=-DUSE_SSE2 -DMMC_USE_SSE_MATH
mex mexomp:        ARFLAGS+=$(MKMEXOPT)
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

mex oct mexomp octomp:   EXTRALIB=
mex oct mexomp octomp:   CCFLAGS+=$(DLLFLAG) -DMCX_CONTAINER
mex oct mexomp octomp:   CPPFLAGS+=-g $(DLLFLAG) -DMCX_CONTAINER
mex oct mexomp octomp:   BINDIR=../mmclab
mex mexomp:     AR=$(MKMEX)
mex mexomp:     AROUTPUT=-output
mex mexomp:     ARFLAGS+=mmclab.cpp -I$(INCLUDEDIR)

oct:            BINARY=mmc.mex
oct octomp:     ARFLAGS+=--mex -DMATLAB_MEX_FILE mmclab.cpp -I$(INCLUDEDIR)
oct octomp:     AR=CC=$(CC) CXX=$(CXX) LFLAGS='$(LFLAGS) $(OPENMPLIB) $(LIBOPENCL) $(MEXLINKOPT)' CPPFLAGS='$(CCFLAGS) $(USERCCFLAGS) -std=c++11' $(USEROCTOPT) $(MKOCT)
oct octomp:     USERARFLAGS+=-o $(BINDIR)/mmc

debug:     sse
debug:     CUCCOPT+=-DMCX_DEBUG

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

cuda: sse

all release sse ssemath prof omp mex oct mexomp octomp pnacl web debug cuda: $(SUBDIRS) $(BINDIR)/$(BINARY)

$(SUBDIRS):
	$(MAKE) -C $@ --no-print-directory

makedirs:
	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); fi
	@if test ! -d $(OBJDIR)/ubj; then $(MKDIR) $(OBJDIR)/ubj; fi
	@if test ! -d $(BINDIR); then $(MKDIR) $(BINDIR); fi

makedocdir:
	@if test ! -d $(DOCDIR); then $(MKDIR) $(DOCDIR); fi

.SUFFIXES : $(OBJSUFFIX) .cpp


##  Compile .cu files ##
$(OBJDIR)/%$(OBJSUFFIX): %.cu
	$(CUDACC) -c $(CUCCOPT) -o $@  $<

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
$(BINDIR)/$(BINARY): makedirs $(CLSOURCE) $(ZMATLIB) $(OBJS)
	@$(ECHO) Building $@
	$(AR)  $(ARFLAGS) $(AROUTPUT) $@ $(OBJS) $(USERARFLAGS) $(EXTRALIB)


$(ZMATLIB):
	-$(MAKE) -C zmat lib AR=ar CPPOPT="$(DLLFLAG)" CCOPT="$(DLLFLAG)" USERLINKOPT=

##  Documentation  ##
doc: makedocdir
	$(DOXY) $(DOXYCFG)

## Clean
clean:
	-$(MAKE) -C zmat clean
	rm -rf $(OBJS) $(OBJDIR) $(BINDIR) $(DOCDIR)
ifdef SUBDIRS
	for i in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$i clean; done
endif

.PHONY: regression clean arch makedirs dep $(SUBDIRS)

.DEFAULT_GOAL := sse
