#ifndef _COMMON_H_
#define _COMMON_H_

// Transfer size between MRAM and WRAM
#ifdef WFT
#define WF_TRANSFER_LOG2 WFT
#define WF_TRANSFER (1 << WF_TRANSFER_LOG2)
#else
#define WF_TRANSFER_LOG2 9
#define WF_TRANSFER (1 << WF_TRANSFER_LOG2)
#endif

#ifdef TASKT
#define TASK_TRANSFER_LOG2 TASKT
#define TASK_TRANSFER (1 << TASK_TRANSFER_LOG2)
#else
#define TASK_TRANSFER_LOG2 3
#define TASK_TRANSFER (1 << TASK_TRANSFER_LOG2)
#endif

#ifdef LENT
#define LEN_TRANSFER_LOG2 LENT
#define LEN_TRANSFER (1 << LEN_TRANSFER_LOG2)
#else
#define LEN_TRANSFER_LOG2 3
#define LEN_TRANSFER (1 << LEN_TRANSFER_LOG2)
#endif

#ifdef CIGART
#define CIGAR_TRANSFER_LOG2 CIGART
#define CIGAR_TRANSFER (1 << CIGAR_TRANSFER_LOG2)
#else
#define CIGAR_TRANSFER_LOG2 3
#define CIGAR_TRANSFER (1 << CIGAR_TRANSFER_LOG2)
#endif

#ifdef SEQT
#define SEQ_TRANSFER_LOG2 SEQT
#define SEQ_TRANSFER (1 << SEQ_TRANSFER_LOG2)
#else
#define SEQ_TRANSFER_LOG2 3
#define SEQ_TRANSFER (1 << SEQ_TRANSFER_LOG2)
#endif

#define MAX_ERRORS 1048576
/*
 * Translate k and offset to coordinates h,v
 */
#define EWAVEFRONT_V(k,offset) ((offset)-(k))
#define EWAVEFRONT_H(k,offset) (offset)

#define EWAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define EWAVEFRONT_OFFSET(h,v)   (h)

#define MAX(a,b) (((a)>=(b))?(a):(b))
#define MIN(a,b) (((a)<=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))

#define SWAP(a,b) do {__typeof__(a) aux = a; a = b; b = aux;} while (0)

#define WAVEFRONT_K_INVERSE(k,plen,tlen)  ((tlen)-(plen)-(k))

#define TYPE_BYTES_SIZE 4

#define REPETITIONS 1
#define PADDING     64
#define MIN_VAL INT16_MIN

#define MRAM_LIMIT 64000000
#define WRAM_LIMIT 64000

#define QUEUE_SZ 32768

char op_to_str[5] = {'M', 'X', 'I', 'D', '?'};

typedef int32_t ewf_offset_t;  // Edit Wavefront Offset

typedef struct  {
    uint32_t longest_seq;
    uint32_t num_pairs_per_dpu;
    uint32_t batch_pairs_per_dpu;
    uint32_t batch_idx;
    ewf_offset_t threshold;
    enum kernels {
		kernel1 = 0,
		nr_kernels = 1,
	} kernel;
}dpu_arguments_t;

typedef enum {
    M = 0,
    X,
    I,
    D,
    NOOP
} ewf_op_t;

typedef struct {
    ewf_offset_t offset;
    unsigned int k;
    ewf_op_t op;
} ewf_breakpoint_t;

/*
 * Wavefront
 */
typedef struct {
  int lo;                      // Effective lowest diagonal (inclusive)
  int hi;                      // Effective highest diagonal (inclusive)
  ewf_offset_t* offsets;       // Offsets
  ewf_offset_t* offsets_mem;   // Offsets memory
  ewf_op_t* ops;           // Operations
  ewf_op_t* ops_mem;       // Operations memory
} edit_wavefront_t;

/*
 * Edit Wavefronts
 */
typedef struct {
  // Dimensions
  int pattern_length;
  int text_length;
  int max_distance;
  // Wavefronts forward and reverse
  edit_wavefront_t* wavefronts_forward[2];
  edit_wavefront_t* wavefronts_reverse[2];
  // Breakpoints
  ewf_breakpoint_t* breakpoints;
  unsigned int num_breakpoints;
  // CIGAR
  char* edit_cigar;
  int edit_cigar_length;
} edit_wavefronts_t;

/*
 * pair metadata for a tasklet job
 */
typedef struct {
  int32_t p_len;
  int32_t t_len;
  int32_t ma_p_start;
  int32_t ma_t_start;
} pair_meta_t; // 64 bytes

typedef struct {
    uint64_t overlap;
    uint64_t base_case;
    uint64_t compute;
    uint64_t extend;
    uint64_t main;
} dpu_results_t;

#ifndef ENERGY
#define ENERGY 0
#endif
#define PRINT 0
#define DEBUG 0
#define WAVEFRONT_BAND 30
#define BASE_BAND 10
#define MIN_WFA_LEN 10
#define MAX_DISTANCE_ADAPTIVE 50
#define MAX_DISTANCE_THRESHOLD 50
#define MAX_ERROR 0.20
#define BATCH_SIZE 0

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#endif
