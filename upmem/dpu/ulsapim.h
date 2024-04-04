#include <stdint.h>
#include <stdio.h>
#include <defs.h>
#include <mram.h>
#include <alloc.h>
#include <mram.h>
#include <seqread.h>
#include <barrier.h>
#include "common.h"
#include <profiling.h>
#include <perfcounter.h>
#include "counterperf.h"
//#include <built_ins.h>

// PROFILING_INIT(mem_read);
// PROFILING_INIT(mem_write);

// PROFILING_INIT(compute);
// PROFILING_INIT(extend);
//PROFILING_INIT(overlap);
//PROFILING_INIT(base_case);

inline void memory_align(uint32_t* ma, uint16_t* read_idx, int bytes_size){
	/* 
	 * Align the memory address to 8 bytes aligned and shift the index to acces the desired data
	 */
	// Initialize pattern lenghts cache
	int reminder = *ma % 8;
	if(reminder == 0)
	{
			*read_idx = 0;
	}
	else
	{
			*read_idx = reminder / bytes_size;
			*ma -= reminder; // -4 bytes because we use int32
	}
}

inline void mram_read_aligned(uint32_t* ma, void* cache_pointer, uint16_t* read_idx, int block_size, int bytes_size){
	//profiling_start&mem_read);
	/*
	 * Align the data to read, then read the data
	 */
	memory_align(ma, read_idx, bytes_size);
	mram_read((__mram_ptr void *) *ma, cache_pointer, block_size);
	//profiling_stop&mem_read);
}

inline void mram_read_aligned_reverse(uint32_t* ma, void* cache_pointer, int16_t* read_idx, int block_size, int bytes_size){
	//profiling_start&mem_read);
	/*
	 * When reading from right to left an array, the readings must be aligned backwards
	 */
	// Initialize pattern lenghts cache
	int reminder = *ma % 8;
	if(reminder == 0)
	{
			*read_idx = 0;
	}
	else
	{
			*read_idx = - ((8 - reminder) / bytes_size);
			*ma += 8 - reminder; // +4 bytes because we use int32
	}
	mram_read((__mram_ptr void *) *ma, cache_pointer, block_size);
	//profiling_stop&mem_read);
	
}

// static void compare_4bytes(uint32_t *cmp, const uint8_t *av, const uint8_t *bv)
// {
//     __asm__("cmpb4 %[cmp], %[av], %[bv]"
//             : [cmp] "=r"(*cmp)
//             : [av] "r"(*(uint32_t *)av), [bv] "r"(*(uint32_t *)bv));
// }

// void opt3_umem_wfa_extend(int k_min, int k_max,
// 	uint32_t loop_limit, ewf_offset_t max_distance, 
// 	int32_t pattern_length, int32_t text_length, 
// 	uint32_t ma_pattern_start, uint32_t ma_text_start, 
// 	uint32_t ma_wf_offsets, 
// 	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets, dpu_results_t* result){
// 	//profiling_start(&extend);
// 	#if PERF_SECTIONS
// 	perfcounter_cycles counter;
// 	timer_start(&counter);
// 	#endif
// 	int32_t k, v, h;
// 	uint16_t read_idx_patt_char;
// 	uint16_t read_idx_text_char;
// 	uint16_t read_idx_wf;
// 	uint32_t cma_wf_offsets;
// 	uint32_t cma_pattern;
// 	uint32_t cma_text;
// 	// Extend diagonally each wavefront point
// 	// Initialize wavefront cache
// 	cma_wf_offsets = (ma_wf_offsets + ((k_min + max_distance) * sizeof(ewf_offset_t))); 
// 	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

// 	// Iterate between k_min and k_max
// 	for (k = k_min; k <= k_max; ++k)
// 	{
// 		const int bases_to_cmp = 4;
// 		int acc = 0;
		
// 		// Init v and h
// 		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
// 		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

// 		while (v < pattern_length && h < text_length) {

// 			// read and get the displacement for both words (displacement = read_idx_xx_char)
// 			cma_pattern = ma_pattern_start + v * sizeof(char);
// 			mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
// 			cma_text  = ma_text_start + h * sizeof(char);
// 			mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);

// 			const uint8_t* word_p_ptr = (uint8_t*)&cache_pattern[read_idx_patt_char];
// 			const uint8_t* word_t_ptr = (uint8_t*)&cache_text[read_idx_text_char];
// 			uint32_t cmp4out;


// 			compare_4bytes(&cmp4out, word_p_ptr, word_t_ptr);

// 			printf("[DEBUG] cmp4out %d pattern %s text %s \n", cmp4out, &cache_pattern[read_idx_patt_char], &cache_text[read_idx_text_char]);

// 			// Iterate over each bit position
// 			for (int i = 7; i >= 0; i--) {
// 				// Extract the bit at position i
// 				uint8_t bit = (cmp4out >> i) & 0x01;
// 				// Print the bit (0 or 1)
// 				printf("%u", bit);
// 			}
// 			printf("\n");


// 			//int eq_elements = __builtin_clo(cmp4out);

// 			//printf("[DEBUG] COUNT LEADING ONES %d\n", eq_elements);

// 			int eq_elements = 3;

// 			// Check pattern and text boundaries and count the zeros to know the mismatching offset
// 			const int next_v = v + bases_to_cmp;
// 			const int next_h = h + bases_to_cmp;
// 			int max_pattern_values = pattern_length - v;
// 			int max_text_values = text_length - h;
// 			if(next_v > pattern_length){
// 				eq_elements = MIN(max_pattern_values, eq_elements);
// 			}
// 			if(next_h > text_length){
// 				eq_elements = MIN(max_text_values, eq_elements);
// 			}

// 			acc += eq_elements;
// 			// printf("[DPU] equivalent elements %d\n", eq_elements);
// 			if (eq_elements < bases_to_cmp) {
// 				break;
// 			}

// 			v += bases_to_cmp;
// 			h += bases_to_cmp;
// 		}

// 		cache_wavefront_offsets[read_idx_wf] += acc;

// 		// Increase index of wavefront cache and write/read if neccesary
// 		read_idx_wf++;
// 		if(read_idx_wf == loop_limit)
// 		{
// 			//profiling_start(&mem_write);
// 			mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_write);
// 			cma_wf_offsets += WF_TRANSFER;
// 			//profiling_start&mem_read);
// 			mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_read);
// 			read_idx_wf = 0;
// 		}
		
// 	}

// 	mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);

// 	#if PERF_SECTIONS
// 	result->extend += timer_stop(&counter);
// 	#endif
// }

// void opt2_umem_wfa_extend(int k_min, int k_max,
// 	uint32_t loop_limit, ewf_offset_t max_distance, 
// 	int32_t pattern_length, int32_t text_length, 
// 	uint32_t ma_pattern_start, uint32_t ma_text_start, 
// 	uint32_t ma_wf_offsets, 
// 	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets, dpu_results_t* result){
// 	//profiling_start(&extend);
// 	#if PERF_SECTIONS
// 	perfcounter_cycles counter;
// 	timer_start(&counter);
// 	#endif
// 	int32_t k, v, h;
// 	uint16_t read_idx_patt_char;
// 	uint16_t read_idx_text_char;
// 	uint16_t read_idx_wf;
// 	uint32_t cma_wf_offsets;
// 	uint32_t cma_pattern;
// 	uint32_t cma_text;
// 	// Extend diagonally each wavefront point
// 	// Initialize wavefront cache
// 	cma_wf_offsets = (ma_wf_offsets + ((k_min + max_distance) * sizeof(ewf_offset_t))); 
// 	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

// 	// Iterate between k_min and k_max
// 	for (k = k_min; k <= k_max; ++k)
// 	{
// 		const int bases_to_cmp = 4;
// 		int acc = 0;
		
// 		// Init v and h
// 		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
// 		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

// 		while (v < pattern_length && h < text_length) {

// 			// read and get the displacement for both words (displacement = read_idx_xx_char)
// 			cma_pattern = ma_pattern_start + v * sizeof(char);
// 			mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
// 			cma_text  = ma_text_start + h * sizeof(char);
// 			mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
// 			uint32_t cache_pattern_offset = 0;
// 			uint32_t cache_text_offset = 0;
// 			if(read_idx_patt_char >= 4) cache_pattern_offset = 4;
// 			if(read_idx_text_char >= 4) cache_text_offset = 4;

// 			//printf("[DEBUG] read idx pattern %d read idx text %d\n", read_idx_patt_char, read_idx_text_char);
// 			// Separate the read pattern and text into 2 words, word and next_word
// 			uint32_t* word_p_ptr = (uint32_t*)&cache_pattern[cache_pattern_offset];
// 			uint32_t* word_p_ptr_2 = word_p_ptr + 1;
// 			uint32_t* word_t_ptr = (uint32_t*)&cache_text[cache_text_offset];
// 			uint32_t* word_t_ptr_2 = word_t_ptr + 1;
	
// 			// for (int i = (sizeof(uint32_t) * 8) - 1; i >= 0; i--) {
// 			// 	printf("%d", (*word_p_ptr >> i) & 1);
// 			// 	if (i % 8 == 0)
// 			// 		printf(" ");
// 			// }
// 			// printf("\n");
// 			// for (int i = (sizeof(uint32_t) * 8) - 1; i >= 0; i--) {
// 			// 	printf("%d", (*word_p_ptr_2 >> i) & 1);
// 			// 	if (i % 8 == 0)
// 			// 		printf(" ");
// 			// }
// 			// printf("\n");

// 			uint32_t sub_word_p_1 = *word_p_ptr >> ((read_idx_patt_char - cache_pattern_offset) * 8);
// 			uint32_t sub_word_p_2 = ((uint64_t) *word_p_ptr_2) << ((bases_to_cmp - (read_idx_patt_char - cache_pattern_offset)) * 8);
// 			uint32_t sub_word_t_1 = *word_t_ptr >> ((read_idx_text_char - cache_text_offset) * 8);
// 			uint32_t sub_word_t_2 = ((uint64_t) *word_t_ptr_2) << ((bases_to_cmp - (read_idx_text_char - cache_text_offset)) * 8);
// 			//printf("[DEBUG] pasa la rola sub_word_p_1 %d sub_word_p_2 %d\n", sub_word_p_1, sub_word_p_2);

// 			// for (int i = (sizeof(uint32_t) * 8) - 1; i >= 0; i--) {
// 			// 	printf("%d", (sub_word_p_1 >> i) & 1);
// 			// 	if (i % 8 == 0)
// 			// 		printf(" ");
// 			// }
// 			// printf("\n");
// 			// for (int i = (sizeof(uint32_t) * 8) - 1; i >= 0; i--) {
// 			// 	printf("%d", (sub_word_p_2 >> i) & 1);
// 			// 	if (i % 8 == 0)
// 			// 		printf(" ");
// 			// }
// 			// printf("\n");
// 			// or of the sub words to eliminate the trash values obtained from unaligned reads
// 			const uint32_t word_p = sub_word_p_1 | sub_word_p_2;
// 			const uint32_t word_t = sub_word_t_1 | sub_word_t_2;
// 			// printf("[DEBUG] target chars %s\n", &cache_pattern[cache_pattern_offset]);
// 			// printf("[DEBUG] output chars %s\n", &word_p);
// 			// xor comparison of the bases
// 			uint32_t diff = word_p ^ word_t;
// 			int eq_elements = __builtin_ctzl(diff) / 8;

// 			// Check pattern and text boundaries and count the zeros to know the mismatching offset
// 			const int next_v = v + bases_to_cmp;
// 			const int next_h = h + bases_to_cmp;
// 			int max_pattern_values = pattern_length - v;
// 			int max_text_values = text_length - h;
// 			if(next_v > pattern_length){
// 				eq_elements = MIN(max_pattern_values, eq_elements);
// 			}
// 			if(next_h > text_length){
// 				eq_elements = MIN(max_text_values, eq_elements);
// 			}

// 			acc += eq_elements;
// 			// printf("[DPU] equivalent elements %d\n", eq_elements);
// 			if (eq_elements < bases_to_cmp) {
// 				break;
// 			}

// 			v += bases_to_cmp;
// 			h += bases_to_cmp;
// 		}

// 		cache_wavefront_offsets[read_idx_wf] += acc;

// 		// Increase index of wavefront cache and write/read if neccesary
// 		read_idx_wf++;
// 		if(read_idx_wf == loop_limit)
// 		{
// 			//profiling_start(&mem_write);
// 			mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_write);
// 			cma_wf_offsets += WF_TRANSFER;
// 			//profiling_start&mem_read);
// 			mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_read);
// 			read_idx_wf = 0;
// 		}
		
// 	}

// 	mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);

// 	#if PERF_SECTIONS
// 	result->extend += timer_stop(&counter);
// 	#endif
// }


// void opt_umem_wfa_extend(int k_min, int k_max,
// 	uint32_t loop_limit, ewf_offset_t max_distance, 
// 	int32_t pattern_length, int32_t text_length, 
// 	uint32_t ma_pattern_start, uint32_t ma_text_start, 
// 	uint32_t ma_wf_offsets, 
// 	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets, dpu_results_t* result){
// 	//profiling_start(&extend);
// 	#if PERF_SECTIONS
// 	perfcounter_cycles counter;
// 	timer_start(&counter);
// 	#endif
// 	int32_t k, v, h;
// 	uint16_t read_idx_patt_char;
// 	uint16_t read_idx_text_char;
// 	uint16_t read_idx_wf;
// 	uint32_t cma_wf_offsets;
// 	uint32_t cma_pattern;
// 	uint32_t cma_text;
// 	// Extend diagonally each wavefront point
// 	// Initialize wavefront cache
// 	cma_wf_offsets = (ma_wf_offsets + ((k_min + max_distance) * sizeof(ewf_offset_t))); 
// 	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

// 	// Iterate between k_min and k_max
// 	for (k = k_min; k <= k_max; ++k)
// 	{
// 		const int bases_to_cmp = 8;
// 		int acc = 0;
		
// 		// Init v and h
// 		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
// 		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

// 		while (v < pattern_length && h < text_length) {

// 			// read and get the displacement for both words (displacement = read_idx_xx_char)
// 			cma_pattern = ma_pattern_start + v * sizeof(char);
// 			mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
// 			cma_text  = ma_text_start + h * sizeof(char);
// 			mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);

// 			// Separate the read pattern and text into 2 words, word and next_word
// 			uint32_t* word_p_ptr = (uint32_t*)cache_pattern;
// 			uint32_t* word_p_ptr_2 = word_p_ptr + 1;
// 			uint32_t* word_p_ptr_3 = word_p_ptr_2 + 1;
// 			uint32_t* word_p_ptr_4 = word_p_ptr_3 + 1;
// 			uint32_t* word_t_ptr = (uint32_t*)cache_text;
// 			uint32_t* word_t_ptr_2 = word_t_ptr + 1;
// 			uint32_t* word_t_ptr_3 = word_t_ptr_2 + 1;
// 			uint32_t* word_t_ptr_4 = word_t_ptr_3 + 1;

// 			// Move word and next word to sub words so they align with 
// 			// the beggining of sub_word_p1 and the end of sub_word_p2 respectively 
// 			// * 8 because each element is 8 bits
// 			uint32_t sub_word_p_1;
// 			uint32_t sub_word_p_2;
// 			uint32_t sub_word_p_3;
// 			uint32_t sub_word_p_4;
// 			uint32_t sub_word_t_1;
// 			uint32_t sub_word_t_2;
// 			uint32_t sub_word_t_3;
// 			uint32_t sub_word_t_4;
// 			if (read_idx_patt_char > 0) {
// 			if (read_idx_patt_char < 4)
// 			{
// 				sub_word_p_1 = *word_p_ptr >> (read_idx_patt_char * 8);
// 				sub_word_p_2 = *word_p_ptr_2 << ((4 - read_idx_patt_char) * 8);
// 				sub_word_p_3 = *word_p_ptr_2 >> (read_idx_patt_char * 8);
// 				sub_word_p_4 = *word_p_ptr_3 >> ((4 - read_idx_patt_char) * 8);
// 			}else{
// 				sub_word_p_1 = *word_p_ptr_2 >> ((read_idx_patt_char - 4) * 8);
// 				sub_word_p_2 = *word_p_ptr_3 << ((8 - read_idx_patt_char) * 8);
// 				sub_word_p_3 = *word_p_ptr_3 >> ((read_idx_patt_char - 4) * 8);
// 				sub_word_p_4 = *word_p_ptr_4 << ((8 - read_idx_patt_char) * 8);
// 			}
// 			}else{
// 				sub_word_p_1 = *word_p_ptr;
// 				sub_word_p_2 = 0;
// 				sub_word_p_3 = *word_p_ptr_2;
// 				sub_word_p_4 = 0;
// 			}

// 			if (read_idx_text_char > 0) {
// 			if (read_idx_text_char < 4)
// 			{
// 				sub_word_t_1 = *word_t_ptr >> (read_idx_text_char * 8);
// 				sub_word_t_2 = *word_t_ptr_2 << ((4 - read_idx_text_char) * 8);
// 				sub_word_t_3 = *word_t_ptr_2 >> (read_idx_text_char * 8);
// 				sub_word_t_4 = *word_t_ptr_3 >> ((4 - read_idx_text_char) * 8);
// 			}else{
// 				sub_word_t_1 = *word_t_ptr_2 >> ((read_idx_text_char - 4) * 8);
// 				sub_word_t_2 = *word_t_ptr_3 << ((8 - read_idx_text_char) * 8);
// 				sub_word_t_3 = *word_t_ptr_3 >> ((read_idx_text_char - 4) * 8);
// 				sub_word_t_4 = *word_t_ptr_4 << ((8 - read_idx_text_char) * 8);
// 			}
// 			}else{
// 				sub_word_t_1 = *word_t_ptr;
// 				sub_word_t_2 = 0;
// 				sub_word_t_3 = *word_t_ptr_2;
// 				sub_word_t_4 = 0;
// 			}

// 			// or of the sub words to eliminate the trash values obtained from unaligned reads
// 			const uint32_t word_p_1 = sub_word_p_1 | sub_word_p_2;
// 			const uint32_t word_p_2 = sub_word_p_3 | sub_word_p_4;
// 			const uint32_t word_t_1 = sub_word_t_1 | sub_word_t_2;
// 			const uint32_t word_t_2 = sub_word_t_3 | sub_word_t_4;
			
// 			// xor comparison of the bases
// 			uint32_t diff_1 = word_p_1 ^ word_t_1;
// 			uint32_t diff_2;
// 			int eq_elements = __builtin_ctzl(diff_1) / 8;
// 			if (eq_elements == 4){
// 				diff_2 = word_p_2 ^ word_t_2;
// 				eq_elements += __builtin_ctzl(diff_2) / 8;
// 			}

// 			// Check pattern and text boundaries and count the zeros to know the mismatching offset
// 			const int next_v = v + bases_to_cmp;
// 			const int next_h = h + bases_to_cmp;
// 			int max_pattern_values = pattern_length - v;
// 			int max_text_values = text_length - h;
// 			if(next_v > pattern_length){
// 				eq_elements = MIN(max_pattern_values, eq_elements);
// 			}
// 			if(next_h > text_length){
// 				eq_elements = MIN(max_text_values, eq_elements);
// 			}

// 			acc += eq_elements;
// 			// printf("[DPU] equivalent elements %d\n", eq_elements);
// 			if (eq_elements < bases_to_cmp) {
// 				break;
// 			}

// 			v += bases_to_cmp;
// 			h += bases_to_cmp;
// 		}

// 		cache_wavefront_offsets[read_idx_wf] += acc;

// 		// Increase index of wavefront cache and write/read if neccesary
// 		read_idx_wf++;
// 		if(read_idx_wf == loop_limit)
// 		{
// 			//profiling_start(&mem_write);
// 			mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_write);
// 			cma_wf_offsets += WF_TRANSFER;
// 			//profiling_start&mem_read);
// 			mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, WF_TRANSFER);
// 			//profiling_stop&mem_read);
// 			read_idx_wf = 0;
// 		}
		
// 	}

// 	mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);

// 	#if PERF_SECTIONS
// 	result->extend += timer_stop(&counter);
// 	#endif
// }

void umem_wfa_extend(int k_min, int k_max,
	uint32_t loop_limit, ewf_offset_t max_distance, 
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	uint32_t ma_wf_offsets, 
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets, dpu_results_t* result){
	//profiling_start(&extend);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	int32_t k, v, h;
	uint16_t read_idx_patt_char;
	uint16_t read_idx_text_char;
	uint16_t read_idx_wf;
	uint32_t cma_wf_offsets;
	uint32_t cma_pattern;
	uint32_t cma_text;
	// Extend diagonally each wavefront point
	// Initialize wavefront cache
	cma_wf_offsets = (ma_wf_offsets + ((k_min + max_distance) * sizeof(ewf_offset_t))); 
	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

	// Iterate between k_min and k_max
	for (k = k_min; k <= k_max; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + v * sizeof(char);
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
		cma_pattern += SEQ_TRANSFER;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
		cma_text += SEQ_TRANSFER;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_pattern += SEQ_TRANSFER;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_text, cache_text, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_text += SEQ_TRANSFER;
				read_idx_text_char = 0;
			}

		}

		// Increase index of wavefront cache and write/read if neccesary
		read_idx_wf++;
		if(read_idx_wf == loop_limit)
		{
			//profiling_start(&mem_write);
			mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
			//profiling_stop&mem_write);
			cma_wf_offsets += WF_TRANSFER;
			//profiling_start&mem_read);
			mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, WF_TRANSFER);
			//profiling_stop&mem_read);
			read_idx_wf = 0;
		}
	}
	// Write back the latest updates to wavefront cache
	//profiling_start(&mem_write);
	mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
	//profiling_stop&mem_write);
	//profiling_stop(&extend);
	#if PERF_SECTIONS
	result->extend += timer_stop(&counter);
	#endif
}

void umem_wfa_extend_rv(int k_min, int k_max,
	uint32_t loop_limit, ewf_offset_t max_distance, 
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	uint32_t ma_wf_offsets, 
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets, dpu_results_t* result){
	//profiling_start(&extend);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	int32_t k, v, h;
	int16_t read_idx_patt_char;
	int16_t read_idx_text_char;
	uint16_t read_idx_wf;
	uint32_t cma_wf_offsets;
	uint32_t cma_pattern;
	uint32_t cma_text;
	// Extend diagonally each wavefront point
	// Initialize wavefront cache
	cma_wf_offsets = (ma_wf_offsets + ((k_min + max_distance) * sizeof(ewf_offset_t)));
	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

	// Iterate between k_min and k_max
	for (k = k_min; k <= k_max; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + ((pattern_length)-v) * sizeof(char) - SEQ_TRANSFER;
		mram_read_aligned_reverse(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
		read_idx_patt_char = SEQ_TRANSFER -1 + read_idx_patt_char;
		cma_pattern -= SEQ_TRANSFER;

		// Init text chars cache
		cma_text  = ma_text_start + ((text_length)-h) * sizeof(char) - SEQ_TRANSFER;
		mram_read_aligned_reverse(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
		read_idx_text_char = SEQ_TRANSFER -1 + read_idx_text_char;
		cma_text -= SEQ_TRANSFER;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char--;
			read_idx_text_char--;

			if(read_idx_patt_char < 0)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_pattern -= SEQ_TRANSFER;
				read_idx_patt_char = SEQ_TRANSFER - 1;
			}


			if(read_idx_text_char < 0)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_text, cache_text, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_text -= SEQ_TRANSFER;
				read_idx_text_char = SEQ_TRANSFER - 1;
			}

		}
		
	}

	// Increase index of wavefront cache and write/read if neccesary
	read_idx_wf++;
	if(read_idx_wf == loop_limit)
	{
		//profiling_start(&mem_write);
		mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
		//profiling_stop&mem_write);
		cma_wf_offsets += WF_TRANSFER;
		//profiling_start&mem_read);
		mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, WF_TRANSFER);
		//profiling_stop&mem_read);
		read_idx_wf = 0;
	}
		//profiling_start(&mem_write);
		mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, WF_TRANSFER);
		//profiling_stop&mem_write);
	//profiling_stop(&extend);
	#if PERF_SECTIONS
	result->extend += timer_stop(&counter);
	#endif
}

ewf_offset_t umem_wfa_extend_after_compute(int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start,  
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets,
	uint16_t read_idx_nextwf, int32_t starting_k, int32_t ending_k, dpu_results_t* result){
	//profiling_stop(&compute);
	//profiling_start(&extend);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	int32_t k, v, h;
	uint16_t read_idx_patt_char;
	uint16_t read_idx_text_char;
	uint32_t cma_pattern;
	uint32_t cma_text;
	ewf_offset_t max = 0;

	// Iterate between k_min and k_max
	for (k = starting_k; k <= ending_k; ++k)
	{

		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_nextwf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_nextwf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + v * sizeof(char);
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
		cma_pattern += SEQ_TRANSFER;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
		cma_text += SEQ_TRANSFER;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_nextwf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_pattern += SEQ_TRANSFER;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_text, cache_text, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_text += SEQ_TRANSFER;
				read_idx_text_char = 0;
			}

		}
		max = (cache_wavefront_offsets[read_idx_nextwf] > max) ? cache_wavefront_offsets[read_idx_nextwf] : max;
		// Increase index of wavefront cache and write/read if neccesary
		read_idx_nextwf++;
	}
	//profiling_stop(&extend);
	//profiling_start(&compute);
	#if PERF_SECTIONS
	result->extend += timer_stop(&counter);
	#endif
	return max;
}

ewf_offset_t umem_wfa_extend_after_compute_rv(int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start,  
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets,
	uint16_t read_idx_nextwf, int32_t starting_k, int32_t ending_k, dpu_results_t* result){
	//profiling_stop(&compute);
	//profiling_start(&extend);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	int32_t k, v, h;
	int16_t read_idx_patt_char;
	int16_t read_idx_text_char;
	uint32_t cma_pattern;
	uint32_t cma_text;
	ewf_offset_t max = 0;

	// Iterate between k_min and k_max
	for (k = starting_k; k <= ending_k; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_nextwf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_nextwf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + ((pattern_length)-v) * sizeof(char) - SEQ_TRANSFER;
		mram_read_aligned_reverse(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
		read_idx_patt_char = SEQ_TRANSFER -1 + read_idx_patt_char;
		cma_pattern -= SEQ_TRANSFER;

		// Init text chars cache
		cma_text  = ma_text_start + ((text_length)-h) * sizeof(char) - SEQ_TRANSFER;
		mram_read_aligned_reverse(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
		read_idx_text_char = SEQ_TRANSFER -1 + read_idx_text_char;
		cma_text -= SEQ_TRANSFER;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_nextwf]);
			v++;
			h++;
			read_idx_patt_char--;
			read_idx_text_char--;

			if(read_idx_patt_char < 0)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_pattern -= SEQ_TRANSFER;
				read_idx_patt_char = SEQ_TRANSFER - 1;
			}


			if(read_idx_text_char < 0)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_text, cache_text, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_text -= SEQ_TRANSFER;
				read_idx_text_char = SEQ_TRANSFER - 1;
			}

		}
		max = (cache_wavefront_offsets[read_idx_nextwf] > max) ? cache_wavefront_offsets[read_idx_nextwf] : max;
		// Increase index of wavefront cache and write/read if neccesary
		read_idx_nextwf++;
	}
	//profiling_stop(&extend);
	//profiling_start(&compute);
	#if PERF_SECTIONS
	result->extend += timer_stop(&counter);
	#endif
	return max;
}

ewf_offset_t umem_wfa_compute_and_extend(int32_t lo, int32_t hi,
	int32_t* nextwavefront_lo, int32_t* nextwavefront_hi, 
	uint32_t loop_limit, ewf_offset_t max_distance, 
	uint32_t ma_wf_offsets, uint32_t ma_wf_next_offsets, 
	ewf_offset_t* cache_wavefront_offsets, ewf_offset_t* cache_nextwavefront_offsets,
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	char* cache_pattern, char* cache_text, uint16_t wf_transfer_size, dpu_results_t* result){
	//profiling_start(&compute);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	uint16_t read_idx_wf;
	uint16_t read_idx_nextwf, starting_idx;
	int32_t k, starting_k;
	uint32_t cma_wf_offsets, cma_nextwf_offsets;
    // Fetch wavefronts
    *nextwavefront_hi = hi + 1;
    *nextwavefront_lo = lo - 1;

    int ins;
    int sub;
    int del;
	ewf_offset_t max = 0;
	ewf_offset_t extend_max = 0;

    // lo
    cma_wf_offsets = (ma_wf_offsets + ((lo - 2 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, wf_transfer_size, TYPE_BYTES_SIZE);

	cma_nextwf_offsets = (ma_wf_next_offsets + ((lo - 1 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_nextwf_offsets, cache_nextwavefront_offsets, &read_idx_nextwf, wf_transfer_size, TYPE_BYTES_SIZE);

	starting_idx = read_idx_nextwf;
	// ins contains offsets[lo]
	ins = cache_wavefront_offsets[read_idx_wf];

	// sub contains offsets[lo+1]
	read_idx_wf++;
	sub = cache_wavefront_offsets[read_idx_wf];

	// del contains offsets[lo+2]
	read_idx_wf++;
	del = cache_wavefront_offsets[read_idx_wf];

    // Compute next wavefront starting point
    for (k = lo - 1; k <= hi + 1; ++k)
    {	
        cache_nextwavefront_offsets[read_idx_nextwf] = MAX((sub + 1), MAX((ins + 1), del)); // MAX
		ins = sub; // Lower
        sub = del; // Mid
        read_idx_wf++;
        if(read_idx_wf == loop_limit)
        {
            cma_wf_offsets += wf_transfer_size;
			//profiling_start&mem_read);
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, wf_transfer_size);
			//profiling_stop&mem_read);
            read_idx_wf = 0;
        }
		del = cache_wavefront_offsets[read_idx_wf]; // Upper
        read_idx_nextwf++;
        if(read_idx_nextwf == loop_limit)
        {
			starting_k = k + (int16_t)(starting_idx - loop_limit + 1);
			#if PERF_SECTIONS
			result->compute += timer_stop(&counter);
			#endif
			extend_max = umem_wfa_extend_after_compute(pattern_length, text_length, 
										ma_pattern_start, ma_text_start, 
										cache_pattern, cache_text, cache_nextwavefront_offsets,
										starting_idx, starting_k, k, result);
			#if PERF_SECTIONS
			timer_start(&counter);
			#endif
			max = (extend_max > max) ? extend_max : max;
			//profiling_start(&mem_write);
            mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, wf_transfer_size);
			//profiling_stop&mem_write);
            cma_nextwf_offsets += wf_transfer_size;
			//profiling_start&mem_read);
            mram_read((__mram_ptr void *) cma_nextwf_offsets, cache_nextwavefront_offsets, wf_transfer_size);
			//profiling_stop&mem_read);
            read_idx_nextwf = 0;
			starting_idx = 0;
        }
    }
	starting_k = hi + 2 - (int16_t)(read_idx_nextwf - starting_idx);
	#if PERF_SECTIONS
	result->compute += timer_stop(&counter);
	#endif
	extend_max = umem_wfa_extend_after_compute(pattern_length, text_length, 
								ma_pattern_start, ma_text_start,  
								cache_pattern, cache_text, cache_nextwavefront_offsets,
								starting_idx, starting_k, hi+1, result);
	#if PERF_SECTIONS
	timer_start(&counter);
	#endif
	max = (extend_max > max) ? extend_max : max;
	//profiling_start(&mem_write);
    mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, wf_transfer_size);
	//profiling_stop(&compute);
	//profiling_stop&mem_write);
	#if PERF_SECTIONS
	result->compute += timer_stop(&counter);
	#endif
	return max;
}

ewf_offset_t umem_wfa_compute_and_extend_rv(int32_t lo, int32_t hi,
	int32_t* nextwavefront_lo, int32_t* nextwavefront_hi, 
	uint32_t loop_limit, ewf_offset_t max_distance, 
	uint32_t ma_wf_offsets, uint32_t ma_wf_next_offsets, 
	ewf_offset_t* cache_wavefront_offsets, ewf_offset_t* cache_nextwavefront_offsets,
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	char* cache_pattern, char* cache_text, uint16_t wf_transfer_size, dpu_results_t* result){
	//profiling_start(&compute);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	uint16_t read_idx_wf;
	uint16_t read_idx_nextwf, starting_idx;
	int32_t k, starting_k;
	uint32_t cma_wf_offsets, cma_nextwf_offsets;
    // Fetch wavefronts
    *nextwavefront_hi = hi + 1;
    *nextwavefront_lo = lo - 1;

    int ins;
    int sub;
    int del;
	ewf_offset_t max = 0;

    // lo
    cma_wf_offsets = (ma_wf_offsets + ((lo - 2 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, wf_transfer_size, TYPE_BYTES_SIZE);

	cma_nextwf_offsets = (ma_wf_next_offsets + ((lo - 1 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_nextwf_offsets, cache_nextwavefront_offsets, &read_idx_nextwf, wf_transfer_size, TYPE_BYTES_SIZE);
	starting_idx = read_idx_nextwf;
	// ins contains offsets[lo]
	ins = cache_wavefront_offsets[read_idx_wf];

	// sub contains offsets[lo+1]
	read_idx_wf++;
	sub = cache_wavefront_offsets[read_idx_wf];

	// del contains offsets[lo+2]
	read_idx_wf++;
	del = cache_wavefront_offsets[read_idx_wf];

    // Compute next wavefront starting point
    for (k = lo - 1; k <= hi + 1; ++k)
    {	
        cache_nextwavefront_offsets[read_idx_nextwf] = MAX((sub + 1), MAX((ins + 1), del)); // MAX
		max = (cache_nextwavefront_offsets[read_idx_nextwf] > max) ? cache_nextwavefront_offsets[read_idx_nextwf] : max;
		ins = sub; // Lower
        sub = del; // Mid
        read_idx_wf++;
        if(read_idx_wf == loop_limit)
        {
            cma_wf_offsets += wf_transfer_size;
			//profiling_start&mem_read);
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, wf_transfer_size);
			//profiling_stop&mem_read);
            read_idx_wf = 0;
        }
		del = cache_wavefront_offsets[read_idx_wf]; // Upper
        read_idx_nextwf++;
        if(read_idx_nextwf == loop_limit)
        {
			starting_k = k + (int16_t)(starting_idx - loop_limit + 1);
			#if PERF_SECTIONS
			result->compute += timer_stop(&counter);
			#endif
			umem_wfa_extend_after_compute_rv(pattern_length, text_length, 
										ma_pattern_start, ma_text_start, 
										cache_pattern, cache_text, cache_nextwavefront_offsets,
										starting_idx, starting_k, k, result);
			#if PERF_SECTIONS
			timer_start(&counter);
			#endif
			//profiling_start(&mem_write);
            mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, wf_transfer_size);
			//profiling_stop&mem_write);
            cma_nextwf_offsets += wf_transfer_size;
			//profiling_start&mem_read);
            mram_read((__mram_ptr void *) cma_nextwf_offsets, cache_nextwavefront_offsets, wf_transfer_size);
			//profiling_stop&mem_read);
            read_idx_nextwf = 0;
			starting_idx = 0;
        }
    }
	starting_k = hi + 2 - (int16_t)(read_idx_nextwf - starting_idx);
	#if PERF_SECTIONS
	result->compute += timer_stop(&counter);
	#endif
	umem_wfa_extend_after_compute_rv(pattern_length, text_length, 
								ma_pattern_start, ma_text_start,  
								cache_pattern, cache_text, cache_nextwavefront_offsets,
								starting_idx, starting_k, hi+1, result);
	#if PERF_SECTIONS
	timer_start(&counter);
	#endif
	//profiling_start(&mem_write);
    mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, wf_transfer_size);
	//profiling_stop(&compute);
	//profiling_stop&mem_write);
	#if PERF_SECTIONS
	result->compute += timer_stop(&counter);
	#endif
	return max;
}

int check_overlap(int32_t wf_fw_lo, int32_t wf_fw_hi,
				   int32_t wf_rv_lo, int32_t wf_rv_hi,
				   ewf_offset_t* cache_wf_fw, ewf_offset_t* cache_wf_rv,
				   uint32_t ma_wf_fw, uint32_t ma_wf_rv,
				   ewf_offset_t max_distance, uint32_t loop_limit,
				   int32_t pattern_length, int32_t text_length,
                   ewf_offset_t* const out_offset,
                   int* const out_k){
	//profiling_start(&overlap);
	uint16_t read_idx_wf_fw;
	int16_t read_idx_wf_rv;
	uint32_t cma_wf_fw, cma_wf_rv;
	const int lo_rev = WAVEFRONT_K_INVERSE(wf_rv_hi,pattern_length,text_length);
  	const int hi_rev = WAVEFRONT_K_INVERSE(wf_rv_lo,pattern_length,text_length);
	const int kmin = MAX(wf_fw_lo, lo_rev);
	const int kmax = MIN(wf_fw_hi, hi_rev);
	// DEBUG diagonals
	// printf("DPU fw_hi %d fw_lo %d rv_hi %d rv_lo %d wf_rv_hi %d wf_rv_lo %d plen %d tlen %d \n", wf_fw_hi, wf_fw_lo, hi_rev, lo_rev, wf_rv_hi, wf_rv_lo, pattern_length, text_length);
	// printf("DPU kmin %d kmax %d \n", kmin, kmax);
	// printf("DPU kmin reverse %d kmax reverse %d \n", WAVEFRONT_K_INVERSE(kmin,pattern_length,text_length), WAVEFRONT_K_INVERSE(kmax,pattern_length,text_length));

	cma_wf_fw = (ma_wf_fw + ((kmin + max_distance) * sizeof(ewf_offset_t)));
	mram_read_aligned(&cma_wf_fw, cache_wf_fw, &read_idx_wf_fw, WF_TRANSFER, TYPE_BYTES_SIZE);
	// Reverse is read backwards to emulate the k index inversion using k inversion macro and subtracting block_size
	// It is needed to subtract 1 size of offset since subtracting block_size will fall out of the first position we want
	cma_wf_rv = (ma_wf_rv + ((WAVEFRONT_K_INVERSE(kmin,pattern_length,text_length) + max_distance +1) * sizeof(ewf_offset_t))- WF_TRANSFER);
	mram_read_aligned_reverse(&cma_wf_rv, cache_wf_rv, &read_idx_wf_rv, WF_TRANSFER, TYPE_BYTES_SIZE);
	// The index will go backwards, so we start at loop ending -1 (since index ends at 0) and align it with read_idx_wf_rv
	read_idx_wf_rv = loop_limit -1 + read_idx_wf_rv;

	for (int k = kmin; k <= kmax; k++) {
		const ewf_offset_t fw_offset  = cache_wf_fw[read_idx_wf_fw];
		const ewf_offset_t rv_offset  = cache_wf_rv[read_idx_wf_rv];
		read_idx_wf_fw++;
		read_idx_wf_rv--;

		if(read_idx_wf_fw == loop_limit)
		{
			cma_wf_fw += WF_TRANSFER;
			//profiling_start&mem_read);
			mram_read((__mram_ptr void *) cma_wf_fw, cache_wf_fw, WF_TRANSFER);
			//profiling_stop&mem_read);
			read_idx_wf_fw = 0;
		}
		if(read_idx_wf_rv < 0) // less than 0 because the 0 index contains information
		{
			cma_wf_rv -= WF_TRANSFER;
			//profiling_start&mem_read);
			mram_read((__mram_ptr void *) cma_wf_rv, cache_wf_rv, WF_TRANSFER);
			//profiling_stop&mem_read);
			read_idx_wf_rv = loop_limit-1;
		}
		if ((fw_offset + rv_offset) >= text_length) 
		{
			// DEBUG breakpoint
			//printf("DPU breakpoint: fw_offset=%d, rev_offset=%d, text_length=%d\n", fw_offset, rv_offset, text_length);

			*out_offset = fw_offset;
			*out_k = k;
			// DEBUG breakpoint
			// printf("DPU Breakpoint at (h, v, k, offset) = (%d, %d)\n", 
			// EWAVEFRONT_H(k, fw_offset), EWAVEFRONT_V(k, fw_offset));
			//profiling_stop(&overlap);
			return 1;
		}
	}
	// DEBUG breakpoint
	// printf("\n");
	//profiling_stop(&overlap);
	return 0;
}

inline void reset_wavefronts_upmem(uint32_t ma_wf_fw, uint32_t ma_wf_fw_next, 
							uint32_t ma_wf_rv, uint32_t ma_wf_rv_next,
							ewf_offset_t* cache_wf_fw, ewf_offset_t* cache_wf_rv,
							uint32_t loop_limit, uint32_t offsets_size_per_tl,
							ewf_offset_t max_distance){

	uint16_t read_idx_wf;
	// Initialize wavefronts to the minimum int16 value,
	// since the reverse wavefront will access uncomputed regions
	for(uint32_t i = 0; i < loop_limit; i++){
		cache_wf_fw[i] = MIN_VAL;
		cache_wf_rv[i] = MIN_VAL;
	} 

	for(uint32_t i = 0; i < offsets_size_per_tl; i+=WF_TRANSFER)
	{
		//profiling_start(&mem_write);
		mram_write(cache_wf_fw, (__mram_ptr void *) (ma_wf_fw      + i), WF_TRANSFER);
		mram_write(cache_wf_fw, (__mram_ptr void *) (ma_wf_fw_next + i), WF_TRANSFER);
		mram_write(cache_wf_rv, (__mram_ptr void *) (ma_wf_rv      + i), WF_TRANSFER);
		mram_write(cache_wf_rv, (__mram_ptr void *) (ma_wf_rv_next + i), WF_TRANSFER);
		//profiling_stop&mem_write);
	}

	// Align the center of the wavefront and write a 0, 
	// since the first position we access must be extended from 0
	uint32_t cma_wf_fw = (ma_wf_fw + (max_distance * sizeof(ewf_offset_t)));
	uint32_t cma_wf_rv = (ma_wf_rv + (max_distance * sizeof(ewf_offset_t)));

	
	memory_align(&cma_wf_fw, &read_idx_wf, TYPE_BYTES_SIZE);
	cache_wf_fw[read_idx_wf] = 0;
	//profiling_start(&mem_write);
	mram_write(cache_wf_fw, (__mram_ptr void *) cma_wf_fw, WF_TRANSFER);
	//profiling_stop&mem_write);
	memory_align(&cma_wf_rv, &read_idx_wf, TYPE_BYTES_SIZE);
	cache_wf_rv[read_idx_wf] = 0;
	//profiling_start(&mem_write);
	mram_write(cache_wf_rv, (__mram_ptr void *) cma_wf_rv, WF_TRANSFER);
	//profiling_stop&mem_write);
}

ewf_offset_t* umem_wfa_base_extend(int k_min, int k_max, 
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets){

	int32_t k, v, h;
	uint16_t read_idx_patt_char;
	uint16_t read_idx_text_char;
	uint16_t read_idx_wf = 0;
	uint32_t cma_pattern;
	uint32_t cma_text;
	// Extend diagonally each wavefront point
	/* ** DEBUG sequences **
	//printf("DPU Pattern: ");
	// for (int i = 0; i < pattern_length; ++i)
	// {
	// 	// v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
	// 	// h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);
	// 	*ma_pattern = ma_pattern_start + i * sizeof(char);
	// 	mram_read_aligned(ma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
	// 	*ma_pattern += BLOCK_SIZE_INPUTS;
	// 	printf("%c",cache_pattern[read_idx_patt_char]);
	// 	read_idx_patt_char++;
	// 	if(read_idx_patt_char == BLOCK_SIZE_INPUTS)
	// 	{
	// 		mram_read((__mram_ptr void *) *ma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
	// 		*ma_pattern += BLOCK_SIZE_INPUTS;
	// 		read_idx_patt_char = 0;
	// 	}
	// 	read_idx_wf++;
	// }
	// printf("\n");
	// read_idx_wf = 0;

	// //printf("DPU text: ");
	// for (int i = 0; i < text_length; ++i)
	// {
	// 	// v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
	// 	// h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);
	// 	*ma_text  = ma_text_start + i * sizeof(char);
	// 	mram_read_aligned(ma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
	// 	*ma_text += BLOCK_SIZE_INPUTS;
	// 	printf("%c",cache_text[read_idx_text_char]);
	// 	read_idx_text_char++;
	// 	if(read_idx_text_char == BLOCK_SIZE_INPUTS)
	// 	{
	// 		mram_read((__mram_ptr void *) *ma_text, cache_text, BLOCK_SIZE_INPUTS);
	// 		*ma_text += BLOCK_SIZE_INPUTS;
	// 		read_idx_text_char = 0;
	// 	}
	// 	read_idx_wf++;
	// }
	// printf("\n");
	// read_idx_wf = 0;
	**/

	// Iterate between k_min and k_max
	for (k = k_min; k <= k_max; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + v * sizeof(char);
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, SEQ_TRANSFER, 1);
		cma_pattern += SEQ_TRANSFER;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, SEQ_TRANSFER, 1);
		cma_text += SEQ_TRANSFER;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_pattern += SEQ_TRANSFER;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == SEQ_TRANSFER)
			{
				//profiling_start&mem_read);
				mram_read((__mram_ptr void *) cma_text, cache_text, SEQ_TRANSFER);
				//profiling_stop&mem_read);
				cma_text += SEQ_TRANSFER;
				read_idx_text_char = 0;
			}

		}

		// Print all the wavefronts DEBUG
		//printf("%d ",cache_wavefront_offsets[read_idx_wf]); 

		// Increase index of wavefront cache and write/read if neccesary
		read_idx_wf++;
	}
	return &cache_wavefront_offsets[read_idx_wf];
}

void umem_wfa_base_compute_next(int32_t lo, int32_t hi,
	int32_t* nextwavefront_lo, int32_t* nextwavefront_hi, 
	ewf_offset_t* cache_wavefront_offsets, ewf_offset_t* cache_nextwavefront_offsets){

	uint16_t read_idx_wf = 0;
	uint16_t read_idx_nextwf = 0;
	int32_t k;
    // Fetch wavefronts
    *nextwavefront_hi = hi + 1;
    *nextwavefront_lo = lo - 1;

    int ins;
    int sub;
    int del;

    // ins contains offsets[lo]
    ins = cache_wavefront_offsets[read_idx_wf];

    // sub contains offsets[lo+1]
    read_idx_wf++;
    sub = cache_wavefront_offsets[read_idx_wf];

    // del contains offsets[lo+2]
    read_idx_wf++;
    del = cache_wavefront_offsets[read_idx_wf];
    
    cache_nextwavefront_offsets[read_idx_nextwf] = ins; // next[lo - 1] = current[lo]

    read_idx_nextwf++;
    // Loop peeling (k=lo)
    cache_nextwavefront_offsets[read_idx_nextwf] = MAX(ins + 1, sub); //lo

    read_idx_nextwf++;
    // Compute next wavefront starting point
    for (k = lo + 1; k <= hi - 1; ++k)
    {
        cache_nextwavefront_offsets[read_idx_nextwf] = MAX((sub + 1), MAX((ins + 1), del)); // MAX
        ins = sub; // Lower
        sub = del; // Mid
        read_idx_wf++;

        del = cache_wavefront_offsets[read_idx_wf]; // Upper

        read_idx_nextwf++;
    }

    read_idx_wf++;

    cache_nextwavefront_offsets[read_idx_nextwf] = MAX(sub, ins) + 1;
    // Loop peeling (k=hi+1)
    read_idx_nextwf++;

    cache_nextwavefront_offsets[read_idx_nextwf] = sub + 1;
}

void edit_wavefronts_backtrace(    
    ewf_offset_t* cache_wf,
	int32_t middle,
    int target_k,
    int distance,
	uint32_t* cma_cigar,
	char* cigar,
	char* cache_cigar_aux,
	uint32_t* cma_cigar_aux,
	int* j){

	int i = CIGAR_TRANSFER-1;
	int cigar_len = 0;
	int del, sub, ins, next_target_k;
	char operation;
	ewf_offset_t max_op;
	ewf_offset_t* cache_wf_curr = &cache_wf[(distance)*(distance)];
	ewf_offset_t* cache_wf_prev_curr;
	middle --;

  	for(;distance>0;--distance){
		cache_wf_prev_curr = &cache_wf[(distance-1)*(distance-1)];
		if(middle + target_k+1 <= middle + middle){
			del = cache_wf_prev_curr[middle + target_k+1]; // Upper
		}else{
			del = MIN_VAL;
		}
		
		if (middle + target_k <= middle + middle && middle + target_k >= 0)
		{
			sub = cache_wf_prev_curr[middle + target_k]+1; // Mid
		}else{
			sub = MIN_VAL;
		}

		if(middle + target_k-1 >= 0){
			ins = cache_wf_prev_curr[middle + target_k-1]+1; // Lower
		}else{
			ins = MIN_VAL;
		}

		if (sub >= ins && sub >= del)
		{
		operation = 'X';
		max_op = sub;
		next_target_k = target_k;
		}else if(ins >= sub && ins >= del){
		operation = 'I';
		max_op = ins;
		next_target_k = target_k -1;
		}else if(del >= ins && del >= sub){
		operation = 'D';
		max_op = del;
		next_target_k = target_k +1;
		}else{
		operation = 'D';
		max_op = del;
		next_target_k = target_k +1;
		}

		while(cache_wf_curr[middle+1 + target_k] > max_op){
			cache_wf_curr[middle+1 + target_k]--;

			cigar[i] = 'M';
			i--;
			if (i < 0){
				*cma_cigar_aux -= CIGAR_TRANSFER;
				//profiling_start(&mem_write);
				mram_write(cigar, (__mram_ptr void *) *cma_cigar_aux, CIGAR_TRANSFER);
				//profiling_stop&mem_write);
				cigar_len += CIGAR_TRANSFER;
				i = CIGAR_TRANSFER-1;
			}
		}

		cigar[i] = operation;
		i--;
		if (i < 0){
			*cma_cigar_aux -= CIGAR_TRANSFER;
			//profiling_start(&mem_write);
			mram_write(cigar, (__mram_ptr void *) *cma_cigar_aux, CIGAR_TRANSFER);
			//profiling_stop&mem_write);
			cigar_len += CIGAR_TRANSFER;
			i = CIGAR_TRANSFER-1;
		}
		target_k = next_target_k;
		cache_wf_curr = cache_wf_prev_curr;
		middle--;

	}
	while(cache_wf_curr[middle+1 + target_k] > 0){
		cache_wf_curr[middle+1 + target_k]--;
		cigar[i] = 'M';
		i--;

		if (i < 0){
			*cma_cigar_aux -= CIGAR_TRANSFER;
			//profiling_start(&mem_write);
			mram_write(cigar, (__mram_ptr void *) *cma_cigar_aux, CIGAR_TRANSFER);
			//profiling_stop&mem_write);
			cigar_len += CIGAR_TRANSFER;
			i = CIGAR_TRANSFER-1;
		}
	}
	i++;
	cigar_len += CIGAR_TRANSFER -i;
	if(i== CIGAR_TRANSFER){
		i = 0;
		*cma_cigar_aux += CIGAR_TRANSFER;
	}

	while (cigar_len > 0)
	{
		cache_cigar_aux[(*j)] = cigar[i];
		(*j)++;
		i++;
		if ((*j)>=CIGAR_TRANSFER)
		{
			//profiling_start(&mem_write);
			mram_write(cache_cigar_aux, (__mram_ptr void *) *cma_cigar, CIGAR_TRANSFER);
			//profiling_stop&mem_write);
			(*j)=0;
			*cma_cigar += CIGAR_TRANSFER;
		}
		if (i>=CIGAR_TRANSFER)
		{
			//profiling_start(&mem_write);
			mram_read((__mram_ptr void *) *cma_cigar_aux, cigar, CIGAR_TRANSFER);
			//profiling_stop&mem_write);
			*cma_cigar_aux += CIGAR_TRANSFER;
			i=0;
		}
		cigar_len--;
		
	}
	cache_cigar_aux[(*j)] = '\0';
}

/*
 * Edit distance alignment using wavefronts
 */
int base_wavefront(uint32_t ma_pattern_start, uint32_t ma_text_start, 
		uint32_t ma_wf, uint32_t ma_wf_next,
		char* cache_pattern, char* cache_text, 
		ewf_offset_t* cache_wf,
		const int pattern_length,
		const int text_length,
		const ewf_offset_t threshold, uint32_t* cma_cigar, char* cache_cigar, char* cache_cigar_aux,
		uint32_t* cma_cigar_aux,
		int* cigar_aux_pos, const uint32_t base_limit, dpu_results_t* result
		){
	//profiling_start(&base_case);
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	timer_start(&counter);
	#endif
	// Parameters
	int32_t wf_lo = 0; 
	int32_t wf_hi = 0;
	int32_t wf_next_lo = 0; 
	int32_t wf_next_hi = 0;
	int longest_seq;
	const int target_k = EWAVEFRONT_DIAGONAL(text_length,pattern_length);
	const int target_k_abs = ABS(target_k);
	const ewf_offset_t target_offset = EWAVEFRONT_OFFSET(text_length,pattern_length);
	ewf_offset_t* cache_wf_curr = cache_wf;
	ewf_offset_t* cache_wf_next;
	int distance = 0;
	#if BANDED
	int32_t diff;
	#endif

	if (pattern_length == 0 && text_length == 0) return distance; 

	for(uint16_t i = 1; i < base_limit; i++){
		cache_wf[i] = MIN_VAL;
	}
	cache_wf[0] = 0;

	// Init wavefronts
	// Compute wavefronts for increasing distance
	for (;distance<threshold;++distance)
	{
		// Extend diagonally each wavefront point
		cache_wf_next = umem_wfa_base_extend(wf_lo, wf_hi,
			pattern_length, text_length, 
			ma_pattern_start, ma_text_start,  
			cache_pattern, cache_text, cache_wf_curr);

		//Exit condition
		if (target_k_abs <= distance && cache_wf_curr[ABS(wf_lo)+target_k] == target_offset) break;

		//Compute next wavefront starting point
		umem_wfa_base_compute_next(wf_lo, wf_hi,
					&wf_next_lo, &wf_next_hi, 
					cache_wf_curr, cache_wf_next);

		#if BANDED
		diff = ABS(wf_next_hi - wf_next_lo);
		if(diff > BASE_BAND){
			wf_next_lo = wf_lo;
			wf_next_hi = wf_hi;
		}
		#endif

		//Swap
		SWAP(wf_lo, wf_next_lo);
		SWAP(wf_hi, wf_next_hi);
		SWAP(ma_wf, ma_wf_next);
		cache_wf_curr = cache_wf_next;
	}

	longest_seq = pattern_length > text_length ? pattern_length : text_length;
	edit_wavefronts_backtrace(    
							cache_wf,
							ABS(wf_lo),
							target_k,
							distance,
							cma_cigar,
							cache_cigar,
							cache_cigar_aux,
							cma_cigar_aux,
							cigar_aux_pos);
	//profiling_stop(&base_case);
	#if PERF_SECTIONS
	result->base_case += timer_stop(&counter);
	#endif
	// Return distance
	return distance;
}

void wfa_adaptive(int32_t* wf_lo, int32_t* wf_hi, int32_t ma_wf, ewf_offset_t* cache_wf, 
				int32_t pattern_length, int32_t text_length,
				ewf_offset_t max_distance, uint32_t loop_limit)
{
    int alignment_k = EWAVEFRONT_DIAGONAL(text_length, pattern_length);
	uint32_t cma_wf_offsets;
	uint16_t read_idx_wf;
	int16_t read_idx_wf_rv;
	ewf_offset_t offset;
	int32_t wf_lo_init = *wf_lo;
	int32_t wf_hi_init = *wf_hi;

    if ((*wf_hi - *wf_lo + 1) < MIN_WFA_LEN)
        return;

    int min_distance = MAX(pattern_length, text_length);

    int klo = *wf_lo;
    int khi = *wf_hi;

	cma_wf_offsets = (ma_wf + ((klo + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wf, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

    for (int k = klo; k <= khi; ++k)
    {
        offset = cache_wf[read_idx_wf];
        int v = EWAVEFRONT_V(k, offset);
        int h = EWAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);
        min_distance = MIN(distance, min_distance);

		read_idx_wf++;
        if(read_idx_wf == loop_limit)
        {			
            cma_wf_offsets += WF_TRANSFER;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, WF_TRANSFER);
            read_idx_wf = 0;
        }
    }

    // reduce m
    // reduce from bottom
    int top_limit = MIN(alignment_k - 1, khi);

	cma_wf_offsets = (ma_wf + ((klo + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wf, &read_idx_wf, WF_TRANSFER, TYPE_BYTES_SIZE);

    for (int k = *wf_lo; k < top_limit; ++k)
    {
        offset = cache_wf[read_idx_wf];
        int v = EWAVEFRONT_V(k, offset);
        int h = EWAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);
        if ((distance - min_distance) <= MAX_DISTANCE_ADAPTIVE)
            break;
        *wf_lo = *wf_lo + 1;

		read_idx_wf++;
        if(read_idx_wf == loop_limit)
        {			
            cma_wf_offsets += WF_TRANSFER;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, WF_TRANSFER);
            read_idx_wf = 0;
        }
    }

    // //reduce from top
    int bottom_limit = MAX(alignment_k + 1, *wf_lo);

	cma_wf_offsets = (ma_wf - WF_TRANSFER + ((khi + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned_reverse(&cma_wf_offsets, cache_wf, &read_idx_wf_rv, WF_TRANSFER, TYPE_BYTES_SIZE);
	read_idx_wf_rv = loop_limit -1 + read_idx_wf_rv;

    for (int k = khi; k > bottom_limit; --k)
    {
        offset = cache_wf[read_idx_wf_rv];
        int v = EWAVEFRONT_V(k, offset);
        int h = EWAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);

        if (distance - min_distance <= MAX_DISTANCE_ADAPTIVE)
            break;
        *wf_hi = *wf_hi - 1;

		read_idx_wf_rv--;
        if(read_idx_wf_rv < 0)
        {			
            cma_wf_offsets -= WF_TRANSFER;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, WF_TRANSFER);
            read_idx_wf_rv = loop_limit-1;
        }
    }

    if (*wf_lo > *wf_hi)
    {
		*wf_hi = wf_hi_init;
		*wf_lo = wf_lo_init;
        return;
    }
}

int find_breakpoint_iterative(int32_t pattern_length, int32_t text_length,
						ewf_offset_t* cache_wf,
						ewf_offset_t* cache_wf_fw, ewf_offset_t* cache_wf_rv,
						ewf_offset_t* cache_wf_fw_next, ewf_offset_t* cache_wf_rv_next,
						char* cache_pattern, char* cache_text,
						uint32_t ma_wf_fw, uint32_t ma_wf_fw_next, 
						uint32_t ma_wf_rv, uint32_t ma_wf_rv_next,
						uint32_t ma_pattern_start, uint32_t ma_text_start,
						const uint32_t loop_limit, const uint32_t base_limit, uint32_t offsets_size_per_tl, const ewf_offset_t threshold,
						uint32_t* cma_cigar, char* cache_cigar, 
						char* cache_cigar_aux,
						uint32_t* cma_cigar_aux, int* cigar_aux_pos,
						pair_meta_t* cache_tasks, int16_t* task_idx, uint32_t* cma_tasks,
						int32_t task_limit, ewf_offset_t max_distance, dpu_results_t *result, perfcounter_cycles *main_counter
						){

	if(max_distance == 0) return 0;
	uint32_t cma_wf_fw;
	uint16_t read_idx_wf;
	uint32_t cma_wf_rv;
	uint16_t read_idx_wf_rv;
	//ewf_offset_t max_distance = (pattern_length + text_length) * MAX_ERROR;
	//ewf_offset_t limit_distance = (pattern_length + text_length) * 0.2;
	// If max_distance equal to 0, then this pair does not have to be computed
	ewf_offset_t distance_fw   = 0;
	ewf_offset_t distance_rv   = 0;
	ewf_offset_t max_fw   	   = 0;
	ewf_offset_t max_rv   	   = 0;
	ewf_offset_t overlap_threshold = (text_length < pattern_length) ? text_length : pattern_length;
	int32_t wf_fw_lo = 0;
	int32_t wf_fw_hi = 0;
	int32_t wf_fw_next_lo = 0;
	int32_t wf_fw_next_hi = 0;
	int32_t wf_rv_lo = 0;
	int32_t wf_rv_hi = 0;
	int32_t wf_rv_next_lo = 0;
	int32_t wf_rv_next_hi = 0;
	uint16_t wffw_transfer_size = 8;
	uint16_t wfrv_transfer_size = 8;
	#if PERF_SECTIONS
	perfcounter_cycles counter;
	#endif
	pair_meta_t pair_metadata;
	#if BANDED
	int32_t diff;
	#endif

	reset_wavefronts_upmem(ma_wf_fw, ma_wf_fw_next, 
							ma_wf_rv, ma_wf_rv_next,
							cache_wf_fw, cache_wf_rv,
							loop_limit, offsets_size_per_tl,
							max_distance);
	#if PERF_MAIN
	result->main += timer_stop(main_counter);
	timer_start(main_counter);
	#endif
	// main wavefront computation loop
	ewf_offset_t offset;
	int k;
	int overlap;

	// extend forward
	umem_wfa_extend(wf_fw_lo, wf_fw_hi,
			loop_limit, max_distance, 
			pattern_length, text_length, 
			ma_pattern_start, ma_text_start, 
			ma_wf_fw, 
			cache_pattern, cache_text, cache_wf_fw, result);	
	#if PERF_MAIN
	result->main += timer_stop(main_counter);
	timer_start(main_counter);
	#endif
	// extend reverse
	umem_wfa_extend_rv(wf_rv_lo, wf_rv_hi,
			loop_limit, max_distance, 
			pattern_length, text_length, 
			ma_pattern_start, ma_text_start, 
			ma_wf_rv,
			cache_pattern, cache_text, cache_wf_rv, result);
	#if PERF_MAIN
	result->main += timer_stop(main_counter);
	timer_start(main_counter);
	#endif
	// 1. Check if both are equal (same offset) (only on first iteration)
	// 2. Check if some of them has extended until the end
	cma_wf_fw = (ma_wf_fw + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_fw, &read_idx_wf, TYPE_BYTES_SIZE);
	cma_wf_rv = (ma_wf_rv + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_rv, &read_idx_wf_rv, TYPE_BYTES_SIZE);
	if ((text_length==pattern_length) && (cache_wf_fw[read_idx_wf] + cache_wf_rv[read_idx_wf_rv]) >= text_length){
		int cigar_len = text_length;
		while (cigar_len > 0)
		{
			cache_cigar_aux[(*cigar_aux_pos)] = 'M';
			(*cigar_aux_pos)++;
			if ((*cigar_aux_pos)>=CIGAR_TRANSFER)
			{
				//profiling_start(&mem_write);
				mram_write(cache_cigar_aux, (__mram_ptr void *) *cma_cigar, CIGAR_TRANSFER);
				//profiling_stop&mem_write);
				(*cigar_aux_pos)=0;
				*cma_cigar += CIGAR_TRANSFER;
			}
			cigar_len--;
			
		}
		cache_cigar_aux[(*cigar_aux_pos)] = '\0';
		
		(*task_idx)--;
		return 0;
	}

	while(1)
	{
		#if ADAPTIVE
		wfa_adaptive(&wf_fw_lo, &wf_fw_hi, ma_wf_fw, cache_wf_fw, 
			pattern_length, text_length, max_distance, loop_limit);
		#endif
		// compute and extend forward
		distance_fw++;
		if(distance_fw >= MAX_DISTANCE_THRESHOLD){
			#if PRINT
			printf("[DPU] distance exceeded %d limit %d plen %d tlen %d\n", distance_fw, max_distance, pattern_length, text_length);
			#endif
			return -1;
		}

		if((((wf_fw_hi+1) - (wf_fw_lo-1)) + 1)*4 > wffw_transfer_size && wffw_transfer_size < WF_TRANSFER){
			wffw_transfer_size = wffw_transfer_size << 1;
			//wffw_transfer_size += 16;
		}
		//wffw_transfer_size = WF_TRANSFER;
		//printf("WE HERE distance_fw %d distance_rv %d\n", distance_fw, distance_rv);

		max_fw = umem_wfa_compute_and_extend(wf_fw_lo, wf_fw_hi,
				&wf_fw_next_lo, &wf_fw_next_hi, 
				loop_limit, max_distance, 
				ma_wf_fw, ma_wf_fw_next,
				cache_wf_fw, cache_wf_fw_next,
				pattern_length, text_length, 
				ma_pattern_start, ma_text_start, 
				cache_pattern, cache_text, wffw_transfer_size, result);
		#if PERF_MAIN
		result->main += timer_stop(main_counter);
		timer_start(main_counter);
		#endif
		#if BANDED
		diff = ABS(wf_fw_next_hi - wf_fw_next_lo);
		if(diff > WAVEFRONT_BAND){
			wf_fw_next_lo = wf_fw_lo;
			wf_fw_next_hi = wf_fw_hi;
		}
		#endif

		// check forward overlap
		if(max_fw + max_rv >= overlap_threshold){
			#if PERF_SECTIONS
			timer_start(&counter);
			#endif
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_lo, wf_rv_hi,
							cache_wf_fw_next, cache_wf_rv,
							ma_wf_fw_next, ma_wf_rv,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);
			#if PERF_MAIN
			result->main += timer_stop(main_counter);
			timer_start(main_counter);
			#endif
			#if PERF_SECTIONS
			result->overlap += timer_stop(&counter);
			#endif
		}
		if(overlap) {
			if(distance_fw + distance_rv + 4 >= threshold){
				//Write 2 tasks meta data
				// task right
				pair_metadata.p_len = pattern_length - EWAVEFRONT_V(k,offset); //pattern_length;
				pair_metadata.t_len = text_length - EWAVEFRONT_H(k,offset); //text_length;
				pair_metadata.ma_p_start = ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char);
				pair_metadata.ma_t_start = ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char);

				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", pair_metadata.ma_p_start, 
				// pair_metadata.ma_t_start, pair_metadata.ma_rev_p_start, pair_metadata.ma_rev_t_start, pair_metadata.p_len, pair_metadata.t_len);

				if (*task_idx >= task_limit){
					//profiling_start(&mem_write);
					mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, TASK_TRANSFER);
					//profiling_stop&mem_write);
					*cma_tasks += TASK_TRANSFER;
					*task_idx = 0;
				}
				cache_tasks[*task_idx] = pair_metadata;

				(*task_idx)++;

				// task left
				pair_metadata.p_len = EWAVEFRONT_V(k,offset); //pattern_length;
				pair_metadata.t_len = EWAVEFRONT_H(k,offset); //text_length;
				pair_metadata.ma_p_start = ma_pattern_start;
				pair_metadata.ma_t_start = ma_text_start;

				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", pair_metadata.ma_p_start, 
				// pair_metadata.ma_t_start, pair_metadata.ma_rev_p_start, pair_metadata.ma_rev_t_start, pair_metadata.p_len, pair_metadata.t_len);

				if (*task_idx >= task_limit){
					//profiling_start(&mem_write);
					mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, TASK_TRANSFER);
					//profiling_stop&mem_write);
					*cma_tasks += TASK_TRANSFER;
					*task_idx = 0;
				}
				cache_tasks[*task_idx] = pair_metadata;

				return distance_fw + distance_rv;
			}else{
				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d plen %d tlen %d\n", ma_pattern_start, 
				// ma_text_start, EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset));
				
				base_wavefront(ma_pattern_start, ma_text_start, 
										ma_wf_fw, ma_wf_fw_next,
										cache_pattern, cache_text, 
										cache_wf,
										EWAVEFRONT_V(k,offset),
										EWAVEFRONT_H(k,offset),
										threshold, cma_cigar, cache_cigar, 
										cache_cigar_aux,
										cma_cigar_aux, cigar_aux_pos, base_limit, result);
				#if PERF_MAIN
				result->main += timer_stop(main_counter);
				timer_start(main_counter);
				#endif
				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d plen %d tlen %d\n", ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
				// ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset));

				base_wavefront(ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
										ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), 
										ma_wf_fw, ma_wf_fw_next,
										cache_pattern, cache_text, 
										cache_wf,
										pattern_length - EWAVEFRONT_V(k,offset),
										text_length - EWAVEFRONT_H(k,offset),
										threshold, cma_cigar, cache_cigar, 
										cache_cigar_aux,
										cma_cigar_aux, cigar_aux_pos, base_limit, result);
				#if PERF_MAIN
				result->main += timer_stop(main_counter);
				timer_start(main_counter);
				#endif
				// decrement task_idx
				(*task_idx)--;
				return distance_fw + distance_rv;
			}
			break;
		}
		#if ADAPTIVE
		wfa_adaptive(&wf_fw_lo, &wf_fw_hi, ma_wf_fw, cache_wf_fw, 
			pattern_length, text_length, max_distance, loop_limit);
		#endif
		distance_rv++;
		if(distance_rv >= MAX_DISTANCE_THRESHOLD){			
			#if PRINT	
			printf("[DPU] distance exceeded %d limit %d plen %d tlen %d\n", distance_rv, max_distance, pattern_length, text_length);
			#endif
			return -1;
		}
		
		if((((wf_fw_hi+1) - (wf_fw_lo-1)) + 1)*4 > wfrv_transfer_size && wfrv_transfer_size < WF_TRANSFER){
			wfrv_transfer_size = wfrv_transfer_size << 1;
			//wfrv_transfer_size += 16;
		}
		//wfrv_transfer_size = WF_TRANSFER;

		// compute and extend reverse
		max_rv = umem_wfa_compute_and_extend_rv(wf_rv_lo, wf_rv_hi,
				&wf_rv_next_lo, &wf_rv_next_hi, 
				loop_limit, max_distance, 
				ma_wf_rv, ma_wf_rv_next, 
				cache_wf_rv, cache_wf_rv_next,
				pattern_length, text_length, 
				ma_pattern_start, ma_text_start, 
				cache_pattern, cache_text, wfrv_transfer_size, result);
		#if PERF_MAIN
		result->main += timer_stop(main_counter);
		timer_start(main_counter);
		#endif

		#if BANDED
		diff = ABS(wf_rv_next_hi - wf_rv_next_lo);
		if(diff > WAVEFRONT_BAND){
			wf_fw_next_lo = wf_fw_lo;
			wf_fw_next_hi = wf_fw_hi;
		}
		#endif

		// check reverse overlap
		if(max_fw + max_rv >= overlap_threshold){
			#if PERF_SECTIONS
			timer_start(&counter);
			#endif
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_next_lo, wf_rv_next_hi,
							cache_wf_fw_next, cache_wf_rv_next,
							ma_wf_fw_next, ma_wf_rv_next,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);
			#if PERF_MAIN
			result->main += timer_stop(main_counter);
			timer_start(main_counter);
			#endif
			#if PERF_SECTIONS
			result->overlap += timer_stop(&counter);
			#endif
		}

		if(overlap) {
			
			if(distance_fw + distance_rv + 4 >= threshold){
				//Write 2 tasks meta data
				// task right
				pair_metadata.p_len = pattern_length - EWAVEFRONT_V(k,offset); // pattern_length;
				pair_metadata.t_len = text_length - EWAVEFRONT_H(k,offset); //text_length;
				pair_metadata.ma_p_start = ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char);
				pair_metadata.ma_t_start = ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char);

				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", pair_metadata.ma_p_start, 
				// pair_metadata.ma_t_start, pair_metadata.ma_rev_p_start, pair_metadata.ma_rev_t_start, pair_metadata.p_len, pair_metadata.t_len);

				if (*task_idx >= task_limit){
					//profiling_start(&mem_write);
					mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, TASK_TRANSFER);
					//profiling_stop&mem_write);
					*cma_tasks += TASK_TRANSFER;
					*task_idx = 0;
				}
				cache_tasks[*task_idx] = pair_metadata;

				(*task_idx)++;

				// task left
				pair_metadata.p_len = EWAVEFRONT_V(k,offset); // pattern_length;
				pair_metadata.t_len = EWAVEFRONT_H(k,offset); //text_length;
				pair_metadata.ma_p_start = ma_pattern_start;
				pair_metadata.ma_t_start = ma_text_start;

				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", pair_metadata.ma_p_start, 
				// pair_metadata.ma_t_start, pair_metadata.ma_rev_p_start, pair_metadata.ma_rev_t_start, pair_metadata.p_len, pair_metadata.t_len);

				if (*task_idx >= task_limit){
					//profiling_start(&mem_write);
					mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, TASK_TRANSFER);
					//profiling_stop&mem_write);
					*cma_tasks += TASK_TRANSFER;
					*task_idx = 0;
				}
				cache_tasks[*task_idx] = pair_metadata;

				return distance_fw + distance_rv;
			}else{
				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d plen %d tlen %d\n", ma_pattern_start, 
				// ma_text_start, EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset));

				base_wavefront(ma_pattern_start, ma_text_start, 
										ma_wf_fw, ma_wf_fw_next,
										cache_pattern, cache_text, 
										cache_wf,
										EWAVEFRONT_V(k,offset),
										EWAVEFRONT_H(k,offset),
										threshold, cma_cigar, cache_cigar, 
										cache_cigar_aux,
										cma_cigar_aux, cigar_aux_pos, base_limit, result);
				#if PERF_MAIN
				result->main += timer_stop(main_counter);
				timer_start(main_counter);
				#endif
				// DEBUG task
				// printf("DPU ma_p_s %d ma_t_s %d plen %d tlen %d\n", ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
				// ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset));

				base_wavefront(ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
										ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), 
										ma_wf_fw, ma_wf_fw_next,
										cache_pattern, cache_text, 
										cache_wf,
										pattern_length - EWAVEFRONT_V(k,offset),
										text_length - EWAVEFRONT_H(k,offset),
										threshold, cma_cigar, cache_cigar, 
										cache_cigar_aux,
										cma_cigar_aux, cigar_aux_pos, base_limit, result);
				#if PERF_MAIN
				result->main += timer_stop(main_counter);
				timer_start(main_counter);
				#endif
				// decrement task_idx
				(*task_idx)--;
				return distance_fw + distance_rv;
			}
			break;
		}
		// ################################## SWAP ##################################
		SWAP(wf_fw_lo, wf_fw_next_lo);
		SWAP(wf_fw_hi, wf_fw_next_hi);
		SWAP(wf_rv_lo, wf_rv_next_lo);
		SWAP(wf_rv_hi, wf_rv_next_hi);
		SWAP(ma_wf_fw, ma_wf_fw_next);
		SWAP(ma_wf_rv, ma_wf_rv_next);
		// ############################# SWAP ENDS HERE #############################
	}
}

