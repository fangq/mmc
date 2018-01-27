# Copyright (c) 2013 The Native Client Authors. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

#
# GNU Make based build file.  For details on GNU Make see:
# http://www.gnu.org/software/make/manual/make.html
#

#
# Get pepper directory for toolchain and includes.
#
# If NACL_SDK_ROOT is not set, then assume it can be found three directories up.
#
THIS_MAKEFILE := $(abspath $(lastword $(MAKEFILE_LIST)))
NACL_SDK_ROOT ?= $(abspath $(dir $(THIS_MAKEFILE))../..)
MMC_ROOT ?= ..

# Project Build flags
WARNINGS := -Wno-long-long -Wall -Wswitch-enum -pedantic -Werror
CXXFLAGS := -pthread -std=gnu++98 $(WARNINGS)

#
# Compute tool paths
#
GETOS := python $(NACL_SDK_ROOT)/tools/getos.py
OSHELPERS = python $(NACL_SDK_ROOT)/tools/oshelpers.py
OSNAME := $(shell $(GETOS))
RM := $(OSHELPERS) rm
MAKE := make

PNACL_TC_PATH := $(abspath $(NACL_SDK_ROOT)/toolchain/$(OSNAME)_pnacl)
PNACL_CXX := $(PNACL_TC_PATH)/bin/pnacl-clang++
PNACL_FINALIZE := $(PNACL_TC_PATH)/bin/pnacl-finalize
CXXFLAGS := -I$(NACL_SDK_ROOT)/include -I$(MMC_ROOT)/src -DMMC_SFMT
LDFLAGS := -L$(NACL_SDK_ROOT)/lib/pnacl/Release -lppapi_cpp -lppapi -L$(MMC_ROOT)/src/bin -lmmc-pnacl

#
# Disable DOS PATH warning when using Cygwin based tools Windows
#
CYGWIN ?= nodosfilewarning
export CYGWIN


# Declare the ALL target first, to make the 'all' target the default build
all: webmmc.pexe

clean:
	$(RM) webmmc.pexe webmmc.bc

webmmc.bc: webmmc.cc
	$(PNACL_CXX) -o $@ $< -O2 $(CXXFLAGS) $(LDFLAGS)

webmmc.pexe: webmmc.bc $(MMC_ROOT)/src/bin/libmmc-pnacl.a
	$(PNACL_FINALIZE) -o $@ $<

$(MMC_ROOT)/src/bin/libmmc-pnacl.a:
	$(MAKE) -C $(MMC_ROOT)/src clean
	$(MAKE) -C $(MMC_ROOT)/src -f makefile_sfmt pnacl BINARY=libmmc-pnacl.a NACL_SDK_ROOT=$(PNACL_TC_PATH)

#
# Makefile target to run the SDK's simple HTTP server and serve this example.
#
HTTPD_PY := python $(NACL_SDK_ROOT)/tools/httpd.py --no_dir_check

.PHONY: serve
serve: all
	$(HTTPD_PY) -C $(CURDIR)
