/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.in by autoheader.  */

/* cpu type */
#define CPU x86

/* POSIX 1b clock */
#define HAVE_CLOCK_GETTIME 1

/* POSIX 1b realtime clock */
#define HAVE_CLOCK_GETTIME_REALTIME 1

/* POSIX 1b realtime HR clock */
/* #undef HAVE_CLOCK_GETTIME_REALTIME_HR */

/* POSIX 1b per-thread clock */
#define HAVE_CLOCK_GETTIME_THREAD CLOCK_THREAD_CPUTIME_ID

/* Native access to a hardware cycle counter */
/* #undef HAVE_CYCLE */

/* Define to 1 if you have the <c_asm.h> header file. */
/* #undef HAVE_C_ASM_H */

/* This platform has the ffsll() function */
#define HAVE_FFSLL 1

/* Define to 1 if you have the `gethrtime' function. */
/* #undef HAVE_GETHRTIME */

/* Full gettid function */
#define HAVE_GETTID 1

/* Normal gettimeofday timer */
/* #undef HAVE_GETTIMEOFDAY */

/* Define if hrtime_t is defined in <sys/time.h> */
/* #undef HAVE_HRTIME_T */

/* Define to 1 if you have the <intrinsics.h> header file. */
/* #undef HAVE_INTRINSICS_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `cpc' library (-lcpc). */
/* #undef HAVE_LIBCPC */

/* perfctr header file */
/* #undef HAVE_LIBPERFCTR_H */

/* Define to 1 if you have the `mach_absolute_time' function. */
/* #undef HAVE_MACH_ABSOLUTE_TIME */

/* Define to 1 if you have the <mach/mach_time.h> header file. */
/* #undef HAVE_MACH_MACH_TIME_H */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Altix memory mapped global cycle counter */
/* #undef HAVE_MMTIMER */

/* Define to 1 if you have the <perfmon/pfmlib.h> header file. */
/* #undef HAVE_PERFMON_PFMLIB_H */

/* Montecito headers */
#define HAVE_PERFMON_PFMLIB_MONTECITO_H 1

/* Working per thread getrusage */
/* #undef HAVE_PER_THREAD_GETRUSAGE */

/* Working per thread timer */
/* #undef HAVE_PER_THREAD_TIMES */

/* new pfmlib_output_param_t */
#define HAVE_PFMLIB_OUTPUT_PFP_PMD_COUNT 1

/* event description function */
#define HAVE_PFM_GET_EVENT_DESCRIPTION 1

/* new pfm_msg_t */
/* #undef HAVE_PFM_MSG_TYPE */

/* old reg_evt_idx */
/* #undef HAVE_PFM_REG_EVT_IDX */

/* Define to 1 if you have the `read_real_time' function. */
/* #undef HAVE_READ_REAL_TIME */

/* Define to 1 if you have the `sched_getcpu' function. */
#define HAVE_SCHED_GETCPU 1

/* Define to 1 if you have the <sched.h> header file. */
#define HAVE_SCHED_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* gettid syscall function */
/* #undef HAVE_SYSCALL_GETTID */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Keyword for per-thread variables */
#define HAVE_THREAD_LOCAL_STORAGE __thread

/* Define to 1 if you have the `time_base_to_time' function. */
/* #undef HAVE_TIME_BASE_TO_TIME */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define for _rtc() intrinsic. */
/* #undef HAVE__RTC */

/* Define if _rtc() is not found. */
/* #undef NO_RTC_INTRINSIC */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ptools-perfapi@icl.utk.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "PAPI"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "PAPI 7.1.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "papi"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "7.1.0.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Use the perfctr virtual TSC for per-thread times */
/* #undef USE_PERFCTR_PTTIMER */

/* Use /proc for per-thread times */
/* #undef USE_PROC_PTTIMER */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
