#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "common.h"

// Params ---------------------------------------------------------------------
typedef struct Params {
  char *   input_file;
  uint32_t  n_warmup;
  uint32_t  n_reps;
  uint32_t  max_pairs;
  int* sets;
  int num_sets;
}Params;

void usage() {
  fprintf(stderr,
    "\nUsage:  ./program [options]"
    "\n"
    "\nGeneral options:"
    "\n    -h        help"
    "\n    -w <W>    # of untimed warmup iterations (default=0)"
    "\n    -e <E>    # of timed repetition iterations (default=1)"
    "\n"
    "\nBenchmark-specific options:"
    "\n    -i <I>    Input file name"
    "\n    -s <S>    Sequence sets and size e.g. '100 1000 1000' "
    "\n    -n <N>    Limits number of pairs to compute"
    "\n");
  }

  struct Params input_params(int argc, char **argv) {
    struct Params p;
    p.input_file    =  malloc(1024 * sizeof(char));
    p.input_file    =  "/home/upmem0046/aalonso/inputs/PacBio.CSS.1M.seq";
    //p.input_file    =  "UNSET";
    p.n_warmup      = 0;
    p.n_reps        = 1;
    p.max_pairs     = UINT32_MAX;
    //p.max_pairs = 70;
    //p.max_pairs     = 2500;
    p.num_sets = 1;
    p.sets = (int*)realloc(p.sets, p.num_sets * sizeof(int));
    p.sets[0] = 18319;
    // p.sets[1] = 10000;
    // p.sets[2] = 12000;
    // p.sets[3] = 18318;
    char* token;
    int opt;
    while((opt = getopt(argc, argv, "hw:e:i:n:s:")) >= 0) {
      switch(opt) {
        case 'h':
        usage();
        exit(0);
        break;
        case 'w': p.n_warmup        = atoi(optarg); break;
        case 'e': p.n_reps          = atoi(optarg); break;
        case 'i': p.input_file      = optarg;       break;
        case 'n': p.max_pairs       = atoi(optarg); break;
        case 's': 
          // Parse the list of integers and store them in an array
          token = strtok(optarg, " ");
          while (token != NULL) {
              p.num_sets++;
              p.sets = (int*)realloc(p.sets, p.num_sets * sizeof(int));
              p.sets[p.num_sets - 1] = atoi(token);
              token = strtok(NULL, " ");
          }        
        break;
        default:
        fprintf(stderr, "\nUnrecognized option!\n");
        usage();
        exit(0);
      }
    }
    assert(NR_DPUS > 0 && "Invalid # of dpus!");
    return p;
  }
#endif
