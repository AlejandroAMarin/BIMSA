DPU_DIR := dpu
HOST_DIR := host
BUILDDIR ?= bin
NR_TASKLETS ?= 1
NR_DPUS ?= 1
BL ?= 8
BLI ?= 5
SZ ?= 1024
VERSION ?= iterative
PERF_CYCLES ?= 0
PERF_INSTRUCTIONS ?= 0
PERF ?= 0
ADAPTIVE ?=0
BANDED ?=0
#SZ ?=8192

define conf_filename
	${BUILDDIR}/.NR_DPUS_$(1)_NR_TASKLETS_$(2).conf
endef
CONF := $(call conf_filename,${NR_DPUS},${NR_TASKLETS})

COMMON_INCLUDES := support
HOST_TARGET := ${BUILDDIR}/ulsapim_host
DPU_TARGET := ${BUILDDIR}/ulsapim_${VERSION}

HOST_SOURCES := $(wildcard ${HOST_DIR}/ulsapim.c)
DPU_SOURCES := $(wildcard ${DPU_DIR}/ulsapim_${VERSION}.c)

.PHONY: all clean test

__dirs := $(shell mkdir -p ${BUILDDIR})

COMMON_FLAGS := -Wall -Wextra  -g -I${COMMON_INCLUDES}
HOST_FLAGS := ${COMMON_FLAGS} -std=c11 -O3 `dpu-pkg-config --cflags --libs dpu` -DNR_TASKLETS=${NR_TASKLETS} -DNR_DPUS=${NR_DPUS} -DBLI=${BLI} -DBL=${BL} -DSTACK_SIZE_DEFAULT=${SZ} -DPERF=${PERF} -DPERF_INSTRUCTIONS=${PERF_INSTRUCTIONS} -DPERF_CYCLES=${PERF_CYCLES} -DDPU_BINARY=\"./bin/ulsapim_${VERSION}\" -lm -fopenmp 

DPU_FLAGS := ${COMMON_FLAGS} -O2 -DNR_TASKLETS=${NR_TASKLETS} -DBL=${BL} -DBLI=${BLI} -DSTACK_SIZE_DEFAULT=${SZ} -DPERF=${PERF} -DADAPTIVE=${ADAPTIVE} -DBANDED=${BANDED}

all: ${HOST_TARGET} ${DPU_TARGET}

${CONF}:
	$(RM) $(call conf_filename,*,*)
	touch ${CONF}

${HOST_TARGET}: ${HOST_SOURCES} ${COMMON_INCLUDES} ${CONF}
	$(CC) -o $@ ${HOST_SOURCES} ${HOST_FLAGS}

${DPU_TARGET}: ${DPU_SOURCES} ${COMMON_INCLUDES} ${CONF}
	dpu-upmem-dpurte-clang ${DPU_FLAGS} -o $@ ${DPU_SOURCES}

clean:
	$(RM) -r $(BUILDDIR)

test: all
	./${HOST_TARGET} 