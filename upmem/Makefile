DPU_DIR := dpu
HOST_DIR := host
BUILDDIR ?= bin
NR_TASKLETS ?= 1
NR_DPUS ?= 1
WFT ?= 9
SEQT ?= 4
CIGART ?= 9
LENT ?= 3
TASKT ?= 9
SZ ?= 1024
ADAPTIVE ?=0
BANDED ?=0
PRINT ?= 0
DEBUG ?= 0
MAX_DISTANCE_THRESHOLD ?= 50
MAX_ERROR ?= 0.20
DYNAMIC ?= 0
BATCH_SIZE ?= 0
#SZ ?=8192
PAPI_DIR=../papi/src/install

SRC_WFA_CPU := $(wildcard ${HOST_DIR}/wfa_cpu.c)
ARGS=-I . -Ilib/
ARGS_WFA_CPU=-L../WFA2-lib/lib/ $(ARGS) -I../WFA2-lib/  -I/${PAPI_DIR}/include -L/${PAPI_DIR}/lib -lwfa -lpapi


define conf_filename
	${BUILDDIR}/.NR_DPUS_$(1)_NR_TASKLETS_$(2).conf
endef
CONF := $(call conf_filename,${NR_DPUS},${NR_TASKLETS})

COMMON_INCLUDES := support ../WFA2-lib
HOST_TARGET := ${BUILDDIR}/bimsa_host
DPU_TARGET := ${BUILDDIR}/bimsa

HOST_SOURCES := $(wildcard ${HOST_DIR}/bimsa_host.c) $(SRC_WFA_CPU) ${HOST_DIR}/timer.c
DPU_SOURCES := $(wildcard ${DPU_DIR}/bimsa.c)

.PHONY: all clean test

__dirs := $(shell mkdir -p ${BUILDDIR})

COMMON_FLAGS := -Wall -Wextra  -g $(addprefix -I, $(COMMON_INCLUDES))
HOST_FLAGS := ${COMMON_FLAGS} -L../WFA2-lib/lib/ -lwfa -std=c11 -O3 `dpu-pkg-config --cflags --libs dpu` -DNR_TASKLETS=${NR_TASKLETS} -DNR_DPUS=${NR_DPUS} -DWFT=${WFT} -DSEQT=${SEQT} -DLEN=${LENT} -DCIGART=${CIGART} -DTASKT=${TASKT} -DSTACK_SIZE_DEFAULT=${SZ} -DPRINT=${PRINT} -DDEBUG=${DEBUG} -DMAX_DISTANCE_THRESHOLD=${MAX_DISTANCE_THRESHOLD} -DMAX_ERROR=${MAX_ERROR} -DDYNAMIC=${DYNAMIC} -DBATCH_SIZE=${BATCH_SIZE} -lm -fopenmp 

DPU_FLAGS := ${COMMON_FLAGS} -O2 -DNR_TASKLETS=${NR_TASKLETS} -DWFT=${WFT} -DSEQT=${SEQT} -DLEN=${LENT} -DCIGART=${CIGART} -DNR_DPUS=${NR_DPUS} -DTASKT=${TASKT} -DSTACK_SIZE_DEFAULT=${SZ} -DADAPTIVE=${ADAPTIVE} -DBANDED=${BANDED} -DPRINT=${PRINT} -DDEBUG=${DEBUG} -DMAX_DISTANCE_THRESHOLD=${MAX_DISTANCE_THRESHOLD} -DMAX_ERROR=${MAX_ERROR} -DDYNAMIC=${DYNAMIC} -DBATCH_SIZE=${BATCH_SIZE}

all: wfa-cpu ${HOST_TARGET} ${DPU_TARGET} 

${CONF}:
	$(RM) $(call conf_filename,*,*)
	touch ${CONF}

${HOST_TARGET}: ${HOST_SOURCES} ${COMMON_INCLUDES} ${CONF}
	$(CC) -o $@ ${HOST_SOURCES} ${HOST_FLAGS}

${DPU_TARGET}: ${DPU_SOURCES} ${COMMON_INCLUDES} ${CONF}
	dpu-upmem-dpurte-clang ${DPU_FLAGS} -o $@ ${DPU_SOURCES}

clean:
	$(RM) -r $(BUILDDIR)

wfa-cpu: $(SRC_WFA_CPU) 
	$(MAKE) -C ../WFA2-lib
	mkdir -p build
	gcc $(ARGS) $(ARGS_WFA_CPU) -Isupport -O3 -fopenmp -Wall -fPIC -c $(SRC_WFA_CPU)
	mv *.o build/

test: all
	./${HOST_TARGET} 
