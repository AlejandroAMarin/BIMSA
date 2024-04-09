#
# Copyright (c) 2002-2006 Hewlett-Packard Development Company, L.P.
# Contributed by Stephane Eranian <eranian@hpl.hp.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.  
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# This file is part of libpfm, a performance monitoring support library for
# applications on Linux.
#

#
# This file defines the global compilation settings.
# It is included by every Makefile
#
#
SYS  := $(shell uname -s)
ARCH := $(shell uname -m)
ifeq (i686,$(findstring i686,$(ARCH)))
override ARCH=ia32
endif
ifeq (i586,$(findstring i586,$(ARCH)))
override ARCH=ia32
endif
ifeq (i486,$(findstring i486,$(ARCH)))
override ARCH=ia32
endif
ifeq (i386,$(findstring i386,$(ARCH)))
override ARCH=ia32
endif
ifeq (ppc,$(findstring ppc,$(ARCH)))
override ARCH=powerpc
endif
ifeq (sparc64,$(findstring sparc64,$(ARCH)))
override ARCH=sparc
endif

#
# CONFIG_PFMLIB_SHARED: y=compile static and shared versions, n=static only
# CONFIG_PFMLIB_OLD_PFMV2: enable old ( 2.x, x <=4) perfmon2 (mutually exclusive with v3 support)
CONFIG_PFMLIB_SHARED?=y
CONFIG_PFMLIB_OLD_PFMV2?=n

#
# Cray-X2 is cross-compiled. Check the programming environment
#
PE := $(shell echo $${CRAY_PE_TARGET})
ifeq (cray-x2,$(PE))
override ARCH=crayx2
endif

#
# Cell Broadband Engine is reported as PPC but needs special handling.
#
ifeq ($(SYS),Linux)
MACHINE := $(shell grep -q 'Cell Broadband Engine' /proc/cpuinfo && echo cell)
ifeq (cell,$(MACHINE))
override ARCH=cell
endif
endif

#
# Library version
#
VERSION=3
REVISION=10
AGE=0

#
# Where should things (lib, headers, man) go in the end.
#
install_prefix?=/usr/local
PREFIX?=$(install_prefix)
LIBDIR=$(PREFIX)/lib
INCDIR=$(PREFIX)/include
MANDIR=$(PREFIX)/share/man
EXAMPLESDIR=$(PREFIX)/share/doc/libpfm-$(VERSION).$(REVISION).$(AGE)/examples

#
# Configuration Paramaters for libpfm library
#
ifeq ($(ARCH),ia64)
CONFIG_PFMLIB_ARCH_IA64=y
endif

ifeq ($(ARCH),x86_64)
CONFIG_PFMLIB_ARCH_X86_64=y
endif

ifeq ($(ARCH),ia32)
CONFIG_PFMLIB_ARCH_I386=y
endif

ifeq ($(ARCH),mips64)
CONFIG_PFMLIB_ARCH_MIPS64=y
#
# SiCortex/Linux
#
MACHINE := $(shell test -f /etc/sicortex-release && echo sicortex)
ifeq (sicortex,$(MACHINE))
CONFIG_PFMLIB_ARCH_SICORTEX=y
endif
endif

ifeq ($(ARCH),powerpc)
CONFIG_PFMLIB_ARCH_POWERPC=y
endif

ifeq ($(ARCH),sparc)
CONFIG_PFMLIB_ARCH_SPARC=y
endif

ifeq ($(XTPE_COMPILE_TARGET),linux)
CONFIG_PFMLIB_ARCH_CRAYXT=y
CONFIG_PFMLIB_SHARED=n
CONFIG_PFMLIB_OLD_PFMV2=y
endif

ifeq ($(XTPE_COMPILE_TARGET),catamount)
CONFIG_PFMLIB_ARCH_CRAYXT=y
CONFIG_PFMLIB_SHARED=n
CONFIG_PFMLIB_OLD_PFMV2=y
endif

ifeq ($(ARCH),crayx2)
CONFIG_PFMLIB_ARCH_CRAYX2=y
CONFIG_PFMLIB_SHARED=n
CONFIG_PFMLIB_OLD_PFMV2=y
endif

ifeq ($(ARCH),cell)
CONFIG_PFMLIB_CELL=y
endif

# handle special cases for 64-bit builds
ifeq ($(BITMODE),64)
ifeq ($(ARCH),powerpc)
CONFIG_PFMLIB_ARCH_POWERPC64=y
endif
endif

#
# you shouldn't have to touch anything beyond this point
#

#
# The entire package can be compiled using 
# icc the Intel Itanium Compiler (7.x,8.x, 9.x)
# or GNU C
#CC=icc
CC?=gcc
LIBS=
INSTALL=install
LN?=ln -sf
PFMINCDIR=$(TOPDIR)/include
PFMLIBDIR=$(TOPDIR)/lib
DBG?=-g -Wall -Werror
# gcc/mips64 bug
ifeq ($(CONFIG_PFMLIB_ARCH_SICORTEX),y)
OPTIM?=-O
else
OPTIM?=-O2
endif
CFLAGS+=$(OPTIM) $(DBG) -I$(PFMINCDIR)
MKDEP=makedepend
PFMLIB=$(PFMLIBDIR)/libpfm.a

# Reset options for Cray XT
ifeq ($(CONFIG_PFMLIB_ARCH_CRAYXT),y)
LDFLAGS+=-static
CONFIG_PFMLIB_OLD_PFMV2=y
endif

# Reset the compiler for Cray-X2 (load x2-gcc module)
ifeq ($(CONFIG_PFMLIB_ARCH_CRAYX2),y)
CC=craynv-cray-linux-gnu-gcc
LDFLAGS+=-static
CONFIG_PFMLIB_OLD_PFMV2=y
endif

ifeq ($(CONFIG_PFMLIB_ARCH_SICORTEX),y)
CONFIG_PFMLIB_OLD_PFMV2=y
endif

ifeq ($(CONFIG_PFMLIB_ARCH_POWERPC64),y)
CFLAGS+= -m64
LDFLAGS+= -m64
LIBDIR=$(PREFIX)/lib64
endif

ifeq ($(CONFIG_PFMLIB_OLD_PFMV2),y)
CFLAGS +=-DPFMLIB_OLD_PFMV2
endif
