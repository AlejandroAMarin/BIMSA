#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "timer.h"

void cpu_recovery(char* patterns, char* texts, uint32_t* pattern_lengths, 
                    uint32_t* text_lengths, uint32_t* cpu_pairs_idx, uint64_t n_sequences, uint64_t longest_seq, Timer* time);