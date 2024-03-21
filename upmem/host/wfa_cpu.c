#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "wavefront/wavefront_align.h"
#include "wfa_cpu.h"
#include <omp.h>

void cpu_recovery(char* patterns, char* texts, uint32_t* pattern_lengths, 
                    uint32_t* text_lengths, uint32_t* cpu_pairs_idx , uint64_t n_sequences, uint64_t longest_seq, Timer* time){
    //start(time, 4, 1);
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = edit;
    attributes.memory_mode = wavefront_memory_ultralow;
    attributes.heuristic.strategy = wf_heuristic_none;
    //attributes.system.verbose = 2;
    wavefront_aligner_t* wf_aligner = wavefront_aligner_new(&attributes);
    #pragma omp for schedule(static,1) nowait
    for(uint32_t seq = 0; seq < n_sequences; seq++){
        const char* text = &texts[cpu_pairs_idx[seq]*longest_seq];
        const char* pattern = &patterns[cpu_pairs_idx[seq]*longest_seq];
        size_t tlen = pattern_lengths[cpu_pairs_idx[seq]];
        size_t plen = text_lengths[cpu_pairs_idx[seq]];        
        wavefront_align(wf_aligner, pattern, plen, text, tlen);
        //printf("[WFA2lib] output: %d\n", wf_aligner->cigar->score);
        //printf("%s\n", wf_aligner->cigar->operations);
    }    
    //stop(time, 4);
}