/*
 * ita2_dear.c - example of how use the D-EAR with the Itanium 2 PMU
 *
 * Copyright (c) 2003-2006 Hewlett-Packard Development Company, L.P.
 * Contributed by Stephane Eranian <eranian@hpl.hp.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * This file is part of libpfm, a performance monitoring support library for
 * applications on Linux/ia64.
 */
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>
#include <fcntl.h>

#include <perfmon/perfmon.h>
#include <perfmon/perfmon_default_smpl.h>
#include <perfmon/pfmlib_itanium2.h>

#define NUM_PMCS PFMLIB_MAX_PMCS
#define NUM_PMDS PFMLIB_MAX_PMDS

#define MAX_EVT_NAME_LEN	128
#define MAX_PMU_NAME_LEN	32

#define SMPL_PERIOD	(40)

#define EVENT_NAME	"data_ear_cache_lat4"

#define M_PMD(x)		(1UL<<(x))
#define DEAR_REGS_MASK		(M_PMD(2)|M_PMD(3)|M_PMD(17))

typedef pfm_default_smpl_hdr_t		dear_hdr_t;
typedef pfm_default_smpl_entry_t	dear_entry_t;
typedef pfm_default_smpl_ctx_arg_t	dear_ctx_t;
#define DEAR_FMT_UUID	        	PFM_DEFAULT_SMPL_UUID

static pfm_uuid_t buf_fmt_id = DEAR_FMT_UUID;


static void *smpl_vaddr;
static unsigned long entry_size;
static int id;

#if defined(__ECC) && defined(__INTEL_COMPILER)

/* if you do not have this file, your compiler is too old */
#include <ia64intrin.h>

#define hweight64(x)	_m64_popcnt(x)

#elif defined(__GNUC__)

static __inline__ int
hweight64 (unsigned long x)
{
	unsigned long result;
	__asm__ ("popcnt %0=%1" : "=r" (result) : "r" (x));
	return (int)result;
}

#else
#error "you need to provide inline assembly from your compiler"
#endif



long
do_test(unsigned long size)
{
    unsigned long i, sum  = 0;
    int *array;

    printf("buffer size %.1fMB\n", (size*sizeof(int))/1024.0);
    array = (int *)malloc(size * sizeof(int));
    if (array == NULL ) {
        printf("line = %d No memory available!\n", __LINE__);
        exit(1);
    }
    for(i=0; i<size; i++) {
        array[i]=1;
    }
    return sum;
}

static void fatal_error(char *fmt,...) __attribute__((noreturn));

static void
fatal_error(char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);

	exit(1);
}

/*
 * print content of sampling buffer
 *
 * XXX: using stdio to print from a signal handler is not safe with multi-threaded
 * applications
 */
#define safe_printf	printf

static void
process_smpl_buffer(void)
{
	dear_hdr_t *hdr;
	dear_entry_t *ent;
	unsigned long pos;
	unsigned long smpl_entry = 0;
	pfm_ita2_pmd_reg_t *reg;
	int ret;
	unsigned int i;
	static unsigned long last_ovfl = ~0UL;

	hdr = (dear_hdr_t *)smpl_vaddr;

	/*
	 * check that we are not diplaying the previous set of samples again.
	 * Required to take care of the last batch of samples.
	 */
	if (hdr->hdr_overflows <= last_ovfl && last_ovfl != ~0UL) {
		printf("skipping identical set of samples %lu <= %lu\n", hdr->hdr_overflows, last_ovfl);
		return;
	}

	pos = (unsigned long)(hdr+1);

	/*
	 * walk through all the entries recored in the buffer
	 */
	for(i=0; i < hdr->hdr_count; i++) {

		ret = 0;

		ent = (dear_entry_t *)pos;
		/*
		 * print entry header
		 */
		safe_printf("Entry %ld PID:%d CPU:%d STAMP:0x%lx IIP:0x%016lx\n",
			smpl_entry++,
			ent->pid,
			ent->cpu,
			ent->tstamp,
			ent->ip);

		/*
		 * point to first recorded register (always contiguous with entry header)
		 */
		reg = (pfm_ita2_pmd_reg_t*)(ent+1);

		safe_printf("PMD2 : 0x%016lx\n", reg->pmd_val);

		reg++;

		safe_printf("PMD3 : 0x%016lx, latency %u\n",
			    reg->pmd_val,
			    reg->pmd3_ita2_reg.dear_latency);

		reg++;

		safe_printf("PMD17: 0x%016lx, valid %c, address 0x%016lx\n",
				reg->pmd_val,
				reg->pmd17_ita2_reg.dear_vl ? 'Y': 'N',
				(reg->pmd17_ita2_reg.dear_iaddr << 4) |
				(unsigned long)reg->pmd17_ita2_reg.dear_slot);

		/*
		 * move to next entry
		 */
		pos += entry_size;
	}
}

static void
overflow_handler(int n, struct siginfo *info, struct sigcontext *sc)
{
	/* dangerous */
	printf("Notification received\n");

	process_smpl_buffer();
	/*
	 * And resume monitoring
	 */
	if (perfmonctl(id, PFM_RESTART,NULL, 0) == -1) {
		perror("PFM_RESTART");
		exit(1);
	}
}

int
main(void)
{
	pfarg_reg_t pd[NUM_PMDS];
	pfarg_reg_t pc[NUM_PMCS];
	pfmlib_input_param_t inp;
	pfmlib_output_param_t outp;
	pfmlib_event_t ev;
	dear_ctx_t ctx[1];
	pfarg_load_t load_args;
	pfmlib_options_t pfmlib_options;
	struct sigaction act;
	unsigned int i;
	int ret, type = 0;

	/*
	 * Initialize pfm library (required before we can use it)
	 */
	if (pfm_initialize() != PFMLIB_SUCCESS) {
		fatal_error("Can't initialize library\n");
	}

	/*
	 * Let's make sure we run this on the right CPU
	 */
	pfm_get_pmu_type(&type);
	if (type != PFMLIB_ITANIUM2_PMU) {
		char model[MAX_PMU_NAME_LEN];
		pfm_get_pmu_name(model, MAX_PMU_NAME_LEN);
		fatal_error("this program does not work with %s PMU\n", model);
	}

	/*
	 * Install the overflow handler (SIGIO)
	 */
	memset(&act, 0, sizeof(act));
	act.sa_handler = (sig_t)overflow_handler;
	sigaction (SIGIO, &act, 0);


	/*
	 * pass options to library (optional)
	 */
	memset(&pfmlib_options, 0, sizeof(pfmlib_options));
	pfmlib_options.pfm_debug = 0; /* set to 1 for debug */
	pfmlib_options.pfm_verbose = 1; /* set to 1 for debug */
	pfm_set_options(&pfmlib_options);



	memset(pd, 0, sizeof(pd));
	memset(pc, 0, sizeof(pc));
	memset(pc, 0, sizeof(pc));
	memset(ctx, 0, sizeof(ctx));
	memset(&load_args, 0, sizeof(load_args));

	/*
	 * prepare parameters to library. we don't use any Itanium
	 * specific features here. so the pfp_model is NULL.
	 */
	memset(&inp,0, sizeof(inp));
	memset(&outp,0, sizeof(outp));

	/*
	 * To count the number of occurence of this instruction, we must
	 * program a counting monitor with the IA64_TAGGED_INST_RETIRED_PMC8
	 * event.
	 */
	if (pfm_find_full_event(EVENT_NAME, &ev) != PFMLIB_SUCCESS) {
		fatal_error("cannot find event %s\n", EVENT_NAME);
	}

	/*
	 * set the (global) privilege mode:
	 * 	PFM_PLM0 : kernel level only
	 */
	inp.pfp_dfl_plm   = PFM_PLM3|PFM_PLM0;

	/*
	 * how many counters we use
	 */
	inp.pfp_event_count = 1;

	/*
	 * propagate the event descriptor
	 */
	inp.pfp_events[0] = ev;
	
	/*
	 * let the library figure out the values for the PMCS
	 *
	 * We use all global settings for this EAR.
	 */
	if ((ret=pfm_dispatch_events(&inp, NULL, &outp, NULL)) != PFMLIB_SUCCESS) {
		fatal_error("cannot configure events: %s\n", pfm_strerror(ret));
	}

	/*
	 * prepare context structure.
	 *
	 * format specific parameters MUST be concatenated to the regular
	 * pfarg_context_t structure. For convenience, the default sampling
	 * format provides a data structure that already combines the pfarg_context_t
	 * with what is needed fot this format.
	 */

	 /*
	  * We initialize the format specific information.
	  * The format is identified by its UUID which must be copied
	  * into the ctx_buf_fmt_id field.
	  */
	memcpy(ctx[0].ctx_arg.ctx_smpl_buf_id, buf_fmt_id, sizeof(pfm_uuid_t));

	/*
	 * the size of the buffer is indicated in bytes (not entries).
	 *
	 * The kernel will record into the buffer up to a certain point.
	 * No partial samples are ever recorded.
	 */
	ctx[0].buf_arg.buf_size = 4096;


	/*
	 * now create the context for self monitoring/per-task
	 */
	if (perfmonctl(0, PFM_CREATE_CONTEXT, ctx, 1) == -1 ) {
		if (errno == ENOSYS) {
			fatal_error("Your kernel does not have performance monitoring support!\n");
		}
		fatal_error("Can't create PFM context %s\n", strerror(errno));
	}
	/*
	 * extract the file descriptor we will use to
	 * identify this newly created context
	 */
	id = ctx[0].ctx_arg.ctx_fd;

	printf("Sampling buffer mapped at %p\n", ctx[0].ctx_arg.ctx_smpl_vaddr);

	smpl_vaddr = ctx[0].ctx_arg.ctx_smpl_vaddr;

	/*
	 * Now prepare the argument to initialize the PMDs and PMCS.
	 * We must pfp_pmc_count to determine the number of PMC to intialize.
	 * We must use pfp_event_count to determine the number of PMD to initialize.
	 * Some events causes extra PMCs to be used, so  pfp_pmc_count may be >= pfp_event_count.
	 *
	 * This step is new compared to libpfm-2.x. It is necessary because the library no
	 * longer knows about the kernel data structures.
	 */

	for (i=0; i < outp.pfp_pmc_count; i++) {
		pc[i].reg_num   = outp.pfp_pmcs[i].reg_num;
		pc[i].reg_value = outp.pfp_pmcs[i].reg_value;
	}

	/*
	 * the PMC controlling the event ALWAYS come first, that's why this loop
	 * is safe even when extra PMC are needed to support a particular event.
	 */
	for (i=0; i < inp.pfp_event_count; i++) {
		pd[i].reg_num   = pc[i].reg_num;
	}

	/*
	 * indicate we want notification when buffer is full
	 */
	pc[0].reg_flags |= PFM_REGFL_OVFL_NOTIFY;

	/*
	 * indicate which PMD to include in the sample
	 */
	pc[0].reg_smpl_pmds[0] = DEAR_REGS_MASK;

	/*
	 * compute size of each sample: fixed-size header + all our DEAR regs
	 */
	entry_size = sizeof(dear_entry_t)+(hweight64(DEAR_REGS_MASK)<<3);

	/*
	 * initialize the PMD and the sampling period
	 */
	pd[0].reg_value       = (~0UL) - SMPL_PERIOD +1;
	pd[0].reg_long_reset  = (~0UL) - SMPL_PERIOD +1;
	pd[0].reg_short_reset = (~0UL) - SMPL_PERIOD +1;

	/*
	 * Now program the registers
	 *
	 * We don't use the save variable to indicate the number of elements passed to
	 * the kernel because, as we said earlier, pc may contain more elements than
	 * the number of events we specified, i.e., contains more thann coutning monitors.
	 */
	if (perfmonctl(id, PFM_WRITE_PMCS, pc, outp.pfp_pmc_count) == -1) {
		fatal_error("perfmonctl error PFM_WRITE_PMCS errno %d\n",errno);
	}
	if (perfmonctl(id, PFM_WRITE_PMDS, pd, inp.pfp_event_count) == -1) {
		fatal_error("perfmonctl error PFM_WRITE_PMDS errno %d\n",errno);
	}
	/*
	 * attach context to stopped task
	 */
	load_args.load_pid = getpid();
	if (perfmonctl(id, PFM_LOAD_CONTEXT, &load_args, 1) == -1) {
		fatal_error("perfmonctl error PFM_LOAD_CONTEXT errno %d\n",errno);
	}
	/*
	 * setup asynchronous notification on the file descriptor
	 */
	ret = fcntl(id, F_SETFL, fcntl(id, F_GETFL, 0) | O_ASYNC);
	if (ret == -1) {
		fatal_error("cannot set ASYNC: %s\n", strerror(errno));
	}

	/*
	 * get ownership of the descriptor
	 */
	ret = fcntl(id, F_SETOWN, getpid());
	if (ret == -1) {
		fatal_error("cannot setown: %s\n", strerror(errno));
	}

	/*
	 * Let's roll now.
	 */
	pfm_self_start(id);

	do_test(10000);

	pfm_self_stop(id);

	/*
	 * We must call the processing routine to cover the last entries recorded
	 * in the sampling buffer, i.e. which may not be full
	 */
	process_smpl_buffer();

	/*
	 * let's stop this now
	 */
	close(id);
	return 0;
}
