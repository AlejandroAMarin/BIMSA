#include <stdint.h>
#include <stdio.h>
#include <defs.h>
#include <mram.h>
#include <alloc.h>
#include <mram.h>
#include <seqread.h>
#include <barrier.h>
#include "common.h"

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
	/*
	 * Align the data to read, then read the data
	 */
	memory_align(ma, read_idx, bytes_size);
	mram_read((__mram_ptr void *) *ma, cache_pointer, block_size);
}

inline void mram_read_aligned_reverse(uint32_t* ma, void* cache_pointer, int16_t* read_idx, int block_size, int bytes_size){
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
}

void umem_wfa_extend(int k_min, int k_max,
	uint32_t loop_limit, ewf_offset_t max_distance, 
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	uint32_t ma_wf_offsets, 
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets){

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
	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

	// Iterate between k_min and k_max
	for (k = k_min; k <= k_max; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + v * sizeof(char);
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
		cma_pattern += BLOCK_SIZE_INPUTS;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
		cma_text += BLOCK_SIZE_INPUTS;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
				cma_pattern += BLOCK_SIZE_INPUTS;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_text, cache_text, BLOCK_SIZE_INPUTS);
				cma_text += BLOCK_SIZE_INPUTS;
				read_idx_text_char = 0;
			}

		}

		// Increase index of wavefront cache and write/read if neccesary
		read_idx_wf++;
		if(read_idx_wf == loop_limit)
		{
			mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, BLOCK_SIZE);
			cma_wf_offsets += BLOCK_SIZE;
			mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, BLOCK_SIZE);
			read_idx_wf = 0;
		}
	}
	// Write back the latest updates to wavefront cache
	mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, BLOCK_SIZE);
}

void umem_wfa_extend_rv(int k_min, int k_max,
	uint32_t loop_limit, ewf_offset_t max_distance, 
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	uint32_t ma_wf_offsets, 
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets){

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
	mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

	// Iterate between k_min and k_max
	for (k = k_min; k <= k_max; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_wf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_wf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + ((pattern_length)-v) * sizeof(char) - BLOCK_SIZE_INPUTS;
		mram_read_aligned_reverse(&cma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
		read_idx_patt_char = BLOCK_SIZE_INPUTS -1 + read_idx_patt_char;
		cma_pattern -= BLOCK_SIZE_INPUTS;

		// Init text chars cache
		cma_text  = ma_text_start + ((text_length)-h) * sizeof(char) - BLOCK_SIZE_INPUTS;
		mram_read_aligned_reverse(&cma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
		read_idx_text_char = BLOCK_SIZE_INPUTS -1 + read_idx_text_char;
		cma_text -= BLOCK_SIZE_INPUTS;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char--;
			read_idx_text_char--;

			if(read_idx_patt_char < 0)
			{
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
				cma_pattern -= BLOCK_SIZE_INPUTS;
				read_idx_patt_char = BLOCK_SIZE_INPUTS - 1;
			}


			if(read_idx_text_char < 0)
			{
				mram_read((__mram_ptr void *) cma_text, cache_text, BLOCK_SIZE_INPUTS);
				cma_text -= BLOCK_SIZE_INPUTS;
				read_idx_text_char = BLOCK_SIZE_INPUTS - 1;
			}

		}
	}

	// Increase index of wavefront cache and write/read if neccesary
	read_idx_wf++;
	if(read_idx_wf == loop_limit)
	{
		mram_write(cache_wavefront_offsets, (__mram_ptr void *) cma_wf_offsets, BLOCK_SIZE);
		cma_wf_offsets += BLOCK_SIZE;
		mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, BLOCK_SIZE);
		read_idx_wf = 0;
	}
		mram_write(cache_wavefront_offsets,(__mram_ptr void *) cma_wf_offsets, BLOCK_SIZE);
}

void umem_wfa_extend_after_compute(int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start,  
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets,
	uint16_t read_idx_nextwf, int32_t starting_k, int32_t ending_k){

	int32_t k, v, h;
	uint16_t read_idx_patt_char;
	uint16_t read_idx_text_char;
	uint32_t cma_pattern;
	uint32_t cma_text;

	// Iterate between k_min and k_max
	for (k = starting_k; k <= ending_k; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_nextwf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_nextwf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + v * sizeof(char);
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
		cma_pattern += BLOCK_SIZE_INPUTS;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
		cma_text += BLOCK_SIZE_INPUTS;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_nextwf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
				cma_pattern += BLOCK_SIZE_INPUTS;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_text, cache_text, BLOCK_SIZE_INPUTS);
				cma_text += BLOCK_SIZE_INPUTS;
				read_idx_text_char = 0;
			}

		}

		// Increase index of wavefront cache and write/read if neccesary
		read_idx_nextwf++;
	}
}

void umem_wfa_extend_after_compute_rv(int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start,  
	char* cache_pattern, char* cache_text, ewf_offset_t* cache_wavefront_offsets,
	uint16_t read_idx_nextwf, int32_t starting_k, int32_t ending_k){

	int32_t k, v, h;
	int16_t read_idx_patt_char;
	int16_t read_idx_text_char;
	uint32_t cma_pattern;
	uint32_t cma_text;

	// Iterate between k_min and k_max
	for (k = starting_k; k <= ending_k; ++k)
	{
		// Init v and h
		v = EWAVEFRONT_V(k, cache_wavefront_offsets[read_idx_nextwf]);
		h = EWAVEFRONT_H(k, cache_wavefront_offsets[read_idx_nextwf]);

		// Init pattern chars cache
		cma_pattern = ma_pattern_start + ((pattern_length)-v) * sizeof(char) - BLOCK_SIZE_INPUTS;
		mram_read_aligned_reverse(&cma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
		read_idx_patt_char = BLOCK_SIZE_INPUTS -1 + read_idx_patt_char;
		cma_pattern -= BLOCK_SIZE_INPUTS;

		// Init text chars cache
		cma_text  = ma_text_start + ((text_length)-h) * sizeof(char) - BLOCK_SIZE_INPUTS;
		mram_read_aligned_reverse(&cma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
		read_idx_text_char = BLOCK_SIZE_INPUTS -1 + read_idx_text_char;
		cma_text -= BLOCK_SIZE_INPUTS;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_nextwf]);
			v++;
			h++;
			read_idx_patt_char--;
			read_idx_text_char--;

			if(read_idx_patt_char < 0)
			{
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
				cma_pattern -= BLOCK_SIZE_INPUTS;
				read_idx_patt_char = BLOCK_SIZE_INPUTS - 1;
			}


			if(read_idx_text_char < 0)
			{
				mram_read((__mram_ptr void *) cma_text, cache_text, BLOCK_SIZE_INPUTS);
				cma_text -= BLOCK_SIZE_INPUTS;
				read_idx_text_char = BLOCK_SIZE_INPUTS - 1;
			}

		}
		// Increase index of wavefront cache and write/read if neccesary
		read_idx_nextwf++;
	}
}

void umem_wfa_compute_and_extend(int32_t lo, int32_t hi,
	int32_t* nextwavefront_lo, int32_t* nextwavefront_hi, 
	uint32_t loop_limit, ewf_offset_t max_distance, 
	uint32_t ma_wf_offsets, uint32_t ma_wf_next_offsets, 
	ewf_offset_t* cache_wavefront_offsets, ewf_offset_t* cache_nextwavefront_offsets,
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	char* cache_pattern, char* cache_text){

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

    // lo
    cma_wf_offsets = (ma_wf_offsets + ((lo - 2 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

	cma_nextwf_offsets = (ma_wf_next_offsets + ((lo - 1 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_nextwf_offsets, cache_nextwavefront_offsets, &read_idx_nextwf, BLOCK_SIZE, TYPE_BYTES_SIZE);
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
            cma_wf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, BLOCK_SIZE);
            read_idx_wf = 0;
        }
		del = cache_wavefront_offsets[read_idx_wf]; // Upper
        read_idx_nextwf++;
        if(read_idx_nextwf == loop_limit)
        {
			starting_k = k + (int16_t)(starting_idx - loop_limit + 1);
			umem_wfa_extend_after_compute(pattern_length, text_length, 
										ma_pattern_start, ma_text_start, 
										cache_pattern, cache_text, cache_nextwavefront_offsets,
										starting_idx, starting_k, k);
			
            mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, BLOCK_SIZE);
            cma_nextwf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_nextwf_offsets, cache_nextwavefront_offsets, BLOCK_SIZE);
            read_idx_nextwf = 0;
			starting_idx = 0;
        }
    }
	starting_k = hi + 2 - (int16_t)(read_idx_nextwf - starting_idx);
	umem_wfa_extend_after_compute(pattern_length, text_length, 
								ma_pattern_start, ma_text_start,  
								cache_pattern, cache_text, cache_nextwavefront_offsets,
								starting_idx, starting_k, hi+1);
    mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, BLOCK_SIZE);
}

void umem_wfa_compute_and_extend_rv(int32_t lo, int32_t hi,
	int32_t* nextwavefront_lo, int32_t* nextwavefront_hi, 
	uint32_t loop_limit, ewf_offset_t max_distance, 
	uint32_t ma_wf_offsets, uint32_t ma_wf_next_offsets, 
	ewf_offset_t* cache_wavefront_offsets, ewf_offset_t* cache_nextwavefront_offsets,
	int32_t pattern_length, int32_t text_length, 
	uint32_t ma_pattern_start, uint32_t ma_text_start, 
	char* cache_pattern, char* cache_text){

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

    // lo
    cma_wf_offsets = (ma_wf_offsets + ((lo - 2 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wavefront_offsets, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

	cma_nextwf_offsets = (ma_wf_next_offsets + ((lo - 1 + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_nextwf_offsets, cache_nextwavefront_offsets, &read_idx_nextwf, BLOCK_SIZE, TYPE_BYTES_SIZE);
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
            cma_wf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wavefront_offsets, BLOCK_SIZE);
            read_idx_wf = 0;
        }
		del = cache_wavefront_offsets[read_idx_wf]; // Upper
        read_idx_nextwf++;
        if(read_idx_nextwf == loop_limit)
        {
			starting_k = k + (int16_t)(starting_idx - loop_limit + 1);
			umem_wfa_extend_after_compute_rv(pattern_length, text_length, 
										ma_pattern_start, ma_text_start, 
										cache_pattern, cache_text, cache_nextwavefront_offsets,
										starting_idx, starting_k, k);
			
            mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, BLOCK_SIZE);
            cma_nextwf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_nextwf_offsets, cache_nextwavefront_offsets, BLOCK_SIZE);
            read_idx_nextwf = 0;
			starting_idx = 0;
        }
    }
	starting_k = hi + 2 - (int16_t)(read_idx_nextwf - starting_idx);
	umem_wfa_extend_after_compute_rv(pattern_length, text_length, 
								ma_pattern_start, ma_text_start,  
								cache_pattern, cache_text, cache_nextwavefront_offsets,
								starting_idx, starting_k, hi+1);
    mram_write(cache_nextwavefront_offsets, (__mram_ptr void *) cma_nextwf_offsets, BLOCK_SIZE);
}

int check_overlap(int32_t wf_fw_lo, int32_t wf_fw_hi,
				   int32_t wf_rv_lo, int32_t wf_rv_hi,
				   ewf_offset_t* cache_wf_fw, ewf_offset_t* cache_wf_rv,
				   uint32_t ma_wf_fw, uint32_t ma_wf_rv,
				   ewf_offset_t max_distance, uint32_t loop_limit,
				   int32_t pattern_length, int32_t text_length,
                   ewf_offset_t* const out_offset,
                   int* const out_k){

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
	mram_read_aligned(&cma_wf_fw, cache_wf_fw, &read_idx_wf_fw, BLOCK_SIZE, TYPE_BYTES_SIZE);
	// Reverse is read backwards to emulate the k index inversion using k inversion macro and subtracting block_size
	// It is needed to subtract 1 size of offset since subtracting block_size will fall out of the first position we want
	cma_wf_rv = (ma_wf_rv + ((WAVEFRONT_K_INVERSE(kmin,pattern_length,text_length) + max_distance +1) * sizeof(ewf_offset_t))- BLOCK_SIZE);
	mram_read_aligned_reverse(&cma_wf_rv, cache_wf_rv, &read_idx_wf_rv, BLOCK_SIZE, TYPE_BYTES_SIZE);
	// The index will go backwards, so we start at loop ending -1 (since index ends at 0) and align it with read_idx_wf_rv
	read_idx_wf_rv = loop_limit -1 + read_idx_wf_rv;

	for (int k = kmin; k <= kmax; k++) {
		const ewf_offset_t fw_offset  = cache_wf_fw[read_idx_wf_fw];
		const ewf_offset_t rv_offset  = cache_wf_rv[read_idx_wf_rv];
		read_idx_wf_fw++;
		read_idx_wf_rv--;

		if(read_idx_wf_fw == loop_limit)
		{
			cma_wf_fw += BLOCK_SIZE;
			mram_read((__mram_ptr void *) cma_wf_fw, cache_wf_fw, BLOCK_SIZE);
			read_idx_wf_fw = 0;
		}
		if(read_idx_wf_rv < 0) // less than 0 because the 0 index contains information
		{
			cma_wf_rv -= BLOCK_SIZE;
			mram_read((__mram_ptr void *) cma_wf_rv, cache_wf_rv, BLOCK_SIZE);
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

			return 1;
		}
	}
	// DEBUG breakpoint
	// printf("\n");

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

	for(uint32_t i = 0; i < offsets_size_per_tl; i+=BLOCK_SIZE)
	{
		mram_write(cache_wf_fw, (__mram_ptr void *) (ma_wf_fw      + i), BLOCK_SIZE);
		mram_write(cache_wf_fw, (__mram_ptr void *) (ma_wf_fw_next + i), BLOCK_SIZE);
		mram_write(cache_wf_rv, (__mram_ptr void *) (ma_wf_rv      + i), BLOCK_SIZE);
		mram_write(cache_wf_rv, (__mram_ptr void *) (ma_wf_rv_next + i), BLOCK_SIZE);
	}

	// Align the center of the wavefront and write a 0, 
	// since the first position we access must be extended from 0
	uint32_t cma_wf_fw = (ma_wf_fw + (max_distance * sizeof(ewf_offset_t)));
	uint32_t cma_wf_rv = (ma_wf_rv + (max_distance * sizeof(ewf_offset_t)));

	
	memory_align(&cma_wf_fw, &read_idx_wf, TYPE_BYTES_SIZE);
	cache_wf_fw[read_idx_wf] = 0;
	mram_write(cache_wf_fw, (__mram_ptr void *) cma_wf_fw, BLOCK_SIZE);
	memory_align(&cma_wf_rv, &read_idx_wf, TYPE_BYTES_SIZE);
	cache_wf_rv[read_idx_wf] = 0;
	mram_write(cache_wf_rv, (__mram_ptr void *) cma_wf_rv, BLOCK_SIZE);
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
		mram_read_aligned(&cma_pattern, cache_pattern, &read_idx_patt_char, BLOCK_SIZE_INPUTS, 1);
		cma_pattern += BLOCK_SIZE_INPUTS;

		// Init text chars cache
		cma_text  = ma_text_start + h * sizeof(char);
		mram_read_aligned(&cma_text, cache_text, &read_idx_text_char, BLOCK_SIZE_INPUTS, 1);
		cma_text += BLOCK_SIZE_INPUTS;

		while (v < pattern_length && h < text_length && cache_pattern[read_idx_patt_char]==cache_text[read_idx_text_char])
		{
			++(cache_wavefront_offsets[read_idx_wf]);
			v++;
			h++;
			read_idx_patt_char++;
			read_idx_text_char++;

			if(read_idx_patt_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_pattern, cache_pattern, BLOCK_SIZE_INPUTS);
				cma_pattern += BLOCK_SIZE_INPUTS;
				read_idx_patt_char = 0;
			}


			if(read_idx_text_char == BLOCK_SIZE_INPUTS)
			{
				mram_read((__mram_ptr void *) cma_text, cache_text, BLOCK_SIZE_INPUTS);
				cma_text += BLOCK_SIZE_INPUTS;
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

	int i = BLOCK_SIZE_INPUTS-1;
	int cigar_len = 0;
	int del, sub, ins, next_target_k;
	uint32_t cma_init_cigar = *cma_cigar_aux;
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
				cma_init_cigar -= BLOCK_SIZE_INPUTS;
				mram_write(cigar, (__mram_ptr void *) cma_init_cigar, BLOCK_SIZE_INPUTS);
				cigar_len += BLOCK_SIZE_INPUTS;
				i = BLOCK_SIZE_INPUTS-1;
			}
		}

		cigar[i] = operation;
		i--;
		if (i < 0){
			cma_init_cigar -= BLOCK_SIZE_INPUTS;
			mram_write(cigar, (__mram_ptr void *) cma_init_cigar, BLOCK_SIZE_INPUTS);
			cigar_len += BLOCK_SIZE_INPUTS;
			i = BLOCK_SIZE_INPUTS-1;
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
			cma_init_cigar -= BLOCK_SIZE_INPUTS;
			mram_write(cigar, (__mram_ptr void *) cma_init_cigar, BLOCK_SIZE_INPUTS);
			cigar_len += BLOCK_SIZE_INPUTS;
			i = BLOCK_SIZE_INPUTS-1;
		}
	}
	i++;
	cigar_len += BLOCK_SIZE_INPUTS -i;
	if(i== BLOCK_SIZE_INPUTS){
		i = 0;
		cma_init_cigar += BLOCK_SIZE_INPUTS;
	}

	while (cigar_len > 0)
	{
		cache_cigar_aux[(*j)] = cigar[i];
		(*j)++;
		i++;
		if ((*j)>=BLOCK_SIZE_INPUTS)
		{
			mram_write(cache_cigar_aux, (__mram_ptr void *) *cma_cigar, BLOCK_SIZE_INPUTS);
			(*j)=0;
			*cma_cigar += BLOCK_SIZE_INPUTS;
		}
		if (i>=BLOCK_SIZE_INPUTS)
		{
			mram_read((__mram_ptr void *) cma_init_cigar, cigar, BLOCK_SIZE_INPUTS);
			cma_init_cigar += BLOCK_SIZE_INPUTS;
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
		int threshold, uint32_t* cma_cigar, char* cache_cigar, char* cache_cigar_aux,
		uint32_t* cma_cigar_aux,
		int* cigar_aux_pos
		){
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

	for(uint16_t i = 1; i < BLOCK_SIZE; i++){
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
    mram_read_aligned(&cma_wf_offsets, cache_wf, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

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
            cma_wf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, BLOCK_SIZE);
            read_idx_wf = 0;
        }
    }

    // reduce m
    // reduce from bottom
    int top_limit = MIN(alignment_k - 1, khi);

	cma_wf_offsets = (ma_wf + ((klo + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned(&cma_wf_offsets, cache_wf, &read_idx_wf, BLOCK_SIZE, TYPE_BYTES_SIZE);

    for (int k = *wf_lo; k < top_limit; ++k)
    {
        offset = cache_wf[read_idx_wf];
        int v = EWAVEFRONT_V(k, offset);
        int h = EWAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);
        if ((distance - min_distance) <= MAX_DISTANCE_THRESHOLD)
            break;
        *wf_lo = *wf_lo + 1;

		read_idx_wf++;
        if(read_idx_wf == loop_limit)
        {			
            cma_wf_offsets += BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, BLOCK_SIZE);
            read_idx_wf = 0;
        }
    }

    // //reduce from top
    int bottom_limit = MAX(alignment_k + 1, *wf_lo);

	cma_wf_offsets = (ma_wf - BLOCK_SIZE + ((khi + max_distance) * sizeof(ewf_offset_t)));
    mram_read_aligned_reverse(&cma_wf_offsets, cache_wf, &read_idx_wf_rv, BLOCK_SIZE, TYPE_BYTES_SIZE);
	read_idx_wf_rv = loop_limit -1 + read_idx_wf_rv;

    for (int k = khi; k > bottom_limit; --k)
    {
        offset = cache_wf[read_idx_wf_rv];
        int v = EWAVEFRONT_V(k, offset);
        int h = EWAVEFRONT_H(k, offset);
        int left_v = pattern_length - v;
        int left_h = text_length - h;
        int distance = MAX(left_v, left_h);

        if (distance - min_distance <= MAX_DISTANCE_THRESHOLD)
            break;
        *wf_hi = *wf_hi - 1;

		read_idx_wf_rv--;
        if(read_idx_wf_rv < 0)
        {			
            cma_wf_offsets -= BLOCK_SIZE;
            mram_read((__mram_ptr void *) cma_wf_offsets, cache_wf, BLOCK_SIZE);
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
						uint32_t loop_limit, uint32_t offsets_size_per_tl, ewf_offset_t threshold,
						uint32_t* cma_cigar, char* cache_cigar, 
						char* cache_cigar_aux,
						uint32_t* cma_cigar_aux, int* cigar_aux_pos,
						pair_meta_t* cache_tasks, int16_t* task_idx, uint32_t* cma_tasks,
						int32_t task_limit, ewf_offset_t max_distance
						){

	uint32_t cma_wf_fw;
	uint16_t read_idx_wf;
	uint32_t cma_wf_rv;
	uint16_t read_idx_wf_rv;
	//ewf_offset_t max_distance = (pattern_length + text_length) * MAX_ERROR;
	//ewf_offset_t limit_distance = (pattern_length + text_length) * 0.2;
	// If max_distance equal to 0, then this pair does not have to be computed
	if(max_distance == 0)
		return 0;
	ewf_offset_t distance_fw   = 0;
	ewf_offset_t distance_rv   = 0;
	int32_t wf_fw_lo = 0;
	int32_t wf_fw_hi = 0;
	int32_t wf_fw_next_lo = 0;
	int32_t wf_fw_next_hi = 0;
	int32_t wf_rv_lo = 0;
	int32_t wf_rv_hi = 0;
	int32_t wf_rv_next_lo = 0;
	int32_t wf_rv_next_hi = 0;
	pair_meta_t pair_metadata;
	#if BANDED
	int32_t diff;
	#endif

	reset_wavefronts_upmem(ma_wf_fw, ma_wf_fw_next, 
							ma_wf_rv, ma_wf_rv_next,
							cache_wf_fw, cache_wf_rv,
							loop_limit, offsets_size_per_tl,
							max_distance);

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
			cache_pattern, cache_text, cache_wf_fw);

	// extend reverse
	umem_wfa_extend_rv(wf_rv_lo, wf_rv_hi,
			loop_limit, max_distance, 
			pattern_length, text_length, 
			ma_pattern_start, ma_text_start, 
			ma_wf_rv,
			cache_pattern, cache_text, cache_wf_rv);

	// 1. Check if both are equal (same offset) (only on first iteration)
	// 2. Check if some of them has extended until the end
	cma_wf_fw = (ma_wf_fw + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_fw, &read_idx_wf, TYPE_BYTES_SIZE);
	cma_wf_rv = (ma_wf_rv + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_rv, &read_idx_wf_rv, TYPE_BYTES_SIZE);
	if ((text_length==pattern_length) && (cache_wf_fw[read_idx_wf] + cache_wf_rv[read_idx_wf_rv]) >= text_length){
		distance_fw = 0;
		distance_rv = 0;
		int cigar_len = text_length;

		while (cigar_len > 0)
		{
			cache_cigar_aux[(*cigar_aux_pos)] = 'M';
			(*cigar_aux_pos)++;
			if ((*cigar_aux_pos)>=BLOCK_SIZE_INPUTS && cache_cigar_aux[(*cigar_aux_pos)-1] != '\0')
			{
				mram_write(cache_cigar_aux, (__mram_ptr void *) *cma_cigar, BLOCK_SIZE_INPUTS);
				(*cigar_aux_pos)=0;
				*cma_cigar += BLOCK_SIZE_INPUTS;
			}
			cigar_len--;
			
		}
		if((*cigar_aux_pos) > 0)(*cigar_aux_pos)--;
		
		(*task_idx)--;

	}else{
		while(1)
		{
			#if ADAPTIVE
			wfa_adaptive(&wf_fw_lo, &wf_fw_hi, ma_wf_fw, cache_wf_fw, 
				pattern_length, text_length, max_distance, loop_limit);
			#endif
			// compute and extend forward
			distance_fw++;
			if(distance_fw >= max_distance){
				printf("[DPU] distance exceeded %d limit %d plen %d tlen %d\n", distance_fw, max_distance, pattern_length, text_length);
				return -1;
			}
			umem_wfa_compute_and_extend(wf_fw_lo, wf_fw_hi,
					&wf_fw_next_lo, &wf_fw_next_hi, 
					loop_limit, max_distance, 
					ma_wf_fw, ma_wf_fw_next,
					cache_wf_fw, cache_wf_fw_next,
					pattern_length, text_length, 
					ma_pattern_start, ma_text_start, 
					cache_pattern, cache_text);

			#if BANDED
			diff = ABS(wf_fw_next_hi - wf_fw_next_lo);
			if(diff > WAVEFRONT_BAND){
				wf_fw_next_lo = wf_fw_lo;
				wf_fw_next_hi = wf_fw_hi;
			}
			#endif

			// check forward overlap
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_lo, wf_rv_hi,
							cache_wf_fw_next, cache_wf_rv,
							ma_wf_fw_next, ma_wf_rv,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);
			if(overlap) {
				if(distance_fw + distance_rv > threshold){
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
						mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, BLOCK_SIZE);
						*cma_tasks += BLOCK_SIZE;
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
						mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, BLOCK_SIZE);
						*cma_tasks += BLOCK_SIZE;
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
											cma_cigar_aux, cigar_aux_pos);
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
											cma_cigar_aux, cigar_aux_pos);
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
			if(distance_rv >= max_distance){				
				printf("[DPU] distance exceeded %d limit %d plen %d tlen %d\n", distance_rv, max_distance, pattern_length, text_length);
				return -1;
			}
			// compute and extend reverse
			umem_wfa_compute_and_extend_rv(wf_rv_lo, wf_rv_hi,
					&wf_rv_next_lo, &wf_rv_next_hi, 
					loop_limit, max_distance, 
					ma_wf_rv, ma_wf_rv_next, 
					cache_wf_rv, cache_wf_rv_next,
					pattern_length, text_length, 
					ma_pattern_start, ma_text_start, 
					cache_pattern, cache_text);

			#if BANDED
			diff = ABS(wf_rv_next_hi - wf_rv_next_lo);
			if(diff > WAVEFRONT_BAND){
				wf_fw_next_lo = wf_fw_lo;
				wf_fw_next_hi = wf_fw_hi;
			}
			#endif

			// check reverse overlap
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_next_lo, wf_rv_next_hi,
							cache_wf_fw_next, cache_wf_rv_next,
							ma_wf_fw_next, ma_wf_rv_next,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);

			if(overlap) {
				if(distance_fw + distance_rv > threshold){
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
						mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, BLOCK_SIZE);
						*cma_tasks += BLOCK_SIZE;
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
						mram_write(cache_tasks, (__mram_ptr void *) *cma_tasks, BLOCK_SIZE);
						*cma_tasks += BLOCK_SIZE;
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
											cma_cigar_aux, cigar_aux_pos);
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
											cma_cigar_aux, cigar_aux_pos);
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
	return distance_fw + distance_rv;
}

int find_breakpoint_recursive(int32_t pattern_length, int32_t text_length,
						ewf_offset_t* cache_wf,
						ewf_offset_t* cache_wf_fw, ewf_offset_t* cache_wf_rv,
						ewf_offset_t* cache_wf_fw_next, ewf_offset_t* cache_wf_rv_next,
						char* cache_pattern, char* cache_text,
						uint32_t ma_wf_fw, uint32_t ma_wf_fw_next, 
						uint32_t ma_wf_rv, uint32_t ma_wf_rv_next,
						uint32_t ma_pattern_start, uint32_t ma_text_start,
						uint32_t loop_limit, uint32_t offsets_size_per_tl, ewf_offset_t threshold,
						uint32_t* cma_cigar, char* cache_cigar, 
						char* cache_cigar_aux,
						uint32_t* cma_cigar_aux, int* cigar_aux_pos
						){

	ewf_offset_t max_distance  = (pattern_length + text_length) * MAX_ERROR;
	uint32_t cma_wf_fw;
	uint16_t read_idx_wf;
	uint32_t cma_wf_rv;
	uint16_t read_idx_wf_rv;
	// If max_distance equal to 0, then this pair does not have to be computed
	if(max_distance == 0)
		return 0;
	ewf_offset_t distance_fw   = 0;
	ewf_offset_t distance_rv   = 0;
	int32_t wf_fw_lo = 0;
	int32_t wf_fw_hi = 0;
	int32_t wf_fw_next_lo = 0;
	int32_t wf_fw_next_hi = 0;
	int32_t wf_rv_lo = 0;
	int32_t wf_rv_hi = 0;
	int32_t wf_rv_next_lo = 0;
	int32_t wf_rv_next_hi = 0;
	#if BANDED
	int32_t diff;
	#endif

	reset_wavefronts_upmem(ma_wf_fw, ma_wf_fw_next, 
							ma_wf_rv, ma_wf_rv_next,
							cache_wf_fw, cache_wf_rv,
							loop_limit, offsets_size_per_tl,
							max_distance);

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
			cache_pattern, cache_text, cache_wf_fw);

	// extend reverse
	umem_wfa_extend_rv(wf_rv_lo, wf_rv_hi,
			loop_limit, max_distance, 
			pattern_length, text_length, 
			ma_pattern_start, ma_text_start, 
			ma_wf_rv,
			cache_pattern, cache_text, cache_wf_rv);

	// // extend reverse
	// umem_wfa_extend(wf_rv_lo, wf_rv_hi,
	// 		loop_limit, max_distance, 
	// 		pattern_length, text_length, 
	// 		ma_rv_pattern_start, ma_rv_text_start, 
	// 		ma_wf_rv,
	// 		cache_pattern, cache_text, cache_wf_rv);

	// 1. Check if both are equal (same offset) (only on first iteration)
	// 2. Check if some of them has extended until the end
	cma_wf_fw = (ma_wf_fw + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_fw, &read_idx_wf, TYPE_BYTES_SIZE);
	cma_wf_rv = (ma_wf_rv + (max_distance * sizeof(ewf_offset_t)));
	memory_align(&cma_wf_rv, &read_idx_wf_rv, TYPE_BYTES_SIZE);
	if ((text_length==pattern_length) && (cache_wf_fw[read_idx_wf] + cache_wf_rv[read_idx_wf_rv]) >= text_length){
		distance_fw = 0;
		distance_rv = 0;
		int cigar_len = text_length;

		while (cigar_len > 0)
		{
			cache_cigar_aux[(*cigar_aux_pos)] = 'M';
			(*cigar_aux_pos)++;
			if ((*cigar_aux_pos)>=BLOCK_SIZE_INPUTS && cache_cigar_aux[(*cigar_aux_pos)-1] != '\0')
			{
				mram_write(cache_cigar_aux, (__mram_ptr void *) *cma_cigar, BLOCK_SIZE_INPUTS);
				(*cigar_aux_pos)=0;
				*cma_cigar += BLOCK_SIZE_INPUTS;
			}
			cigar_len--;
			
		}
		if((*cigar_aux_pos) > 0)(*cigar_aux_pos)--;

	}else{
		while(1)
		{
			#if ADAPTIVE
			wfa_adaptive(&wf_fw_lo, &wf_fw_hi, ma_wf_fw, cache_wf_fw, 
				pattern_length, text_length, max_distance, loop_limit);
			#endif
			// compute and extend forward
			distance_fw++;
			umem_wfa_compute_and_extend(wf_fw_lo, wf_fw_hi,
					&wf_fw_next_lo, &wf_fw_next_hi, 
					loop_limit, max_distance, 
					ma_wf_fw, ma_wf_fw_next,
					cache_wf_fw, cache_wf_fw_next,
					pattern_length, text_length, 
					ma_pattern_start, ma_text_start, 
					cache_pattern, cache_text);

			#if BANDED
			diff = ABS(wf_fw_next_hi - wf_fw_next_lo);
			if(diff > WAVEFRONT_BAND){
				wf_fw_next_lo = wf_fw_lo;
				wf_fw_next_hi = wf_fw_hi;
			}
			#endif

			// check forward overlap
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_lo, wf_rv_hi,
							cache_wf_fw_next, cache_wf_rv,
							ma_wf_fw_next, ma_wf_rv,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);
			if(overlap) {
				if(distance_fw + distance_rv > threshold){
					// DEBUG recursion
					//printf("DPU doing recursion curr distance %d\n", distance_fw + distance_rv);
					//printf("DPU left flank length %d\n", EWAVEFRONT_H(k,offset));
					// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", ma_pattern_start, 
					// ma_text_start, ma_rv_pattern_start + (pattern_length - EWAVEFRONT_V(k,offset)) * sizeof(char), 
					// ma_rv_text_start + (text_length - EWAVEFRONT_H(k,offset)) * sizeof(char), EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset));

					find_breakpoint_recursive(EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset),
									cache_wf,
									cache_wf_fw, cache_wf_rv,
									cache_wf_fw_next, cache_wf_rv_next,
									cache_pattern, cache_text,
									ma_wf_fw, ma_wf_fw_next, 
									ma_wf_rv, ma_wf_rv_next,
									ma_pattern_start, ma_text_start,
									loop_limit, offsets_size_per_tl, threshold,
									cma_cigar, cache_cigar, cache_cigar_aux, 
									cma_cigar_aux, cigar_aux_pos
									);
					// DEBUG recursion
					//printf("DPU right flank length %d\n", text_length - EWAVEFRONT_H(k,offset));
					// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
					// ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), ma_rv_pattern_start, ma_rv_text_start, pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset));

					find_breakpoint_recursive(pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset),
									cache_wf,
									cache_wf_fw, cache_wf_rv,
									cache_wf_fw_next, cache_wf_rv_next,
									cache_pattern, cache_text,
									ma_wf_fw, ma_wf_fw_next, 
									ma_wf_rv, ma_wf_rv_next,
									ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
									ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char),
									loop_limit, offsets_size_per_tl, threshold,
									cma_cigar, cache_cigar, cache_cigar_aux, 
									cma_cigar_aux, cigar_aux_pos
									);
				}else{
					// DEBUG recursion
					//printf("DPU base case found curr distance %d\n", distance_fw + distance_rv);
					//max_distance = EWAVEFRONT_V(k,offset) + EWAVEFRONT_H(k,offset);
					// printf("DPU printing max distance for first case: %d\n", max_distance);
					//printf("DPU left flank length %d\n", EWAVEFRONT_V(k,offset));
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
											cma_cigar_aux, cigar_aux_pos);
					// DEBUG recursion
					//max_distance = pattern_length - EWAVEFRONT_V(k,offset) + text_length - EWAVEFRONT_H(k,offset);
					// printf("DPU printing max distance: %d\n", max_distance);
					//printf("DPU Solution distance %d\n", distance);
					//printf("DPU right flank length %d\n", pattern_length - EWAVEFRONT_V(k,offset));
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
											cma_cigar_aux, cigar_aux_pos);

				}
				break;
			}
			#if ADAPTIVE
			wfa_adaptive(&wf_fw_lo, &wf_fw_hi, ma_wf_fw, cache_wf_fw, 
				pattern_length, text_length, max_distance, loop_limit);
			#endif
			distance_rv++;
			
			// compute and extend reverse
			umem_wfa_compute_and_extend_rv(wf_rv_lo, wf_rv_hi,
					&wf_rv_next_lo, &wf_rv_next_hi, 
					loop_limit, max_distance, 
					ma_wf_rv, ma_wf_rv_next, 
					cache_wf_rv, cache_wf_rv_next,
					pattern_length, text_length, 
					ma_pattern_start, ma_text_start, 
					cache_pattern, cache_text);
			
			#if BANDED
			diff = ABS(wf_rv_next_hi - wf_rv_next_lo);
			if(diff > WAVEFRONT_BAND){
				wf_fw_next_lo = wf_fw_lo;
				wf_fw_next_hi = wf_fw_hi;
			}
			#endif

			// check reverse overlap
			overlap = check_overlap(wf_fw_next_lo, wf_fw_next_hi, 
							wf_rv_next_lo, wf_rv_next_hi,
							cache_wf_fw_next, cache_wf_rv_next,
							ma_wf_fw_next, ma_wf_rv_next,
							max_distance, loop_limit,
							pattern_length, text_length,
							&offset,
							&k);

			if(overlap) {
				if(distance_fw + distance_rv > threshold){
					// DEBUG recursion
					//printf("DPU doing recursion curr distance %d\n", distance_fw + distance_rv);
					//printf("DPU left flank length %d\n", EWAVEFRONT_H(k,offset));
					// printf("DPU ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", ma_pattern_start, 
					// ma_text_start, ma_rv_pattern_start + (pattern_length - EWAVEFRONT_V(k,offset)) * sizeof(char), 
					// ma_rv_text_start + (text_length - EWAVEFRONT_H(k,offset)) * sizeof(char), EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset));

					find_breakpoint_recursive(EWAVEFRONT_V(k,offset), EWAVEFRONT_H(k,offset),
									cache_wf,
									cache_wf_fw, cache_wf_rv,
									cache_wf_fw_next, cache_wf_rv_next,
									cache_pattern, cache_text,
									ma_wf_fw, ma_wf_fw_next, 
									ma_wf_rv, ma_wf_rv_next,
									ma_pattern_start, ma_text_start,
									loop_limit, offsets_size_per_tl, threshold,
									cma_cigar, cache_cigar, 
									cache_cigar_aux, cma_cigar_aux, cigar_aux_pos
									);
					// DEBUG recursion
					//printf("DPU right flank length %d\n", text_length - EWAVEFRONT_H(k,offset));
					// printf("DPU ma_p_s %d ma_t_s %d plen %d tlen %d\n", ma_pattern_start + EWAVEFRONT_V(k,offset) * sizeof(char), 
					// ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char), pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset));

					find_breakpoint_recursive(pattern_length - EWAVEFRONT_V(k,offset), text_length - EWAVEFRONT_H(k,offset),
									cache_wf,
									cache_wf_fw, cache_wf_rv,
									cache_wf_fw_next, cache_wf_rv_next,
									cache_pattern, cache_text,
									ma_wf_fw, ma_wf_fw_next, 
									ma_wf_rv, ma_wf_rv_next,
									ma_pattern_start + EWAVEFRONT_V(k,offset) * 	sizeof(char), 
									ma_text_start + EWAVEFRONT_H(k,offset) * sizeof(char),
									loop_limit, offsets_size_per_tl, threshold,
									cma_cigar, cache_cigar, 
									cache_cigar_aux, cma_cigar_aux, cigar_aux_pos
									);
				}else{
					// DEBUG recursion
					//printf("DPU base case found curr distance %d\n", distance_fw + distance_rv);
					//max_distance = EWAVEFRONT_V(k,offset) + EWAVEFRONT_H(k,offset);
					// printf("DPU printing max distance for first case: %d\n", max_distance);
					//printf("DPU left flank length %d\n", EWAVEFRONT_V(k,offset));
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
											cma_cigar_aux, cigar_aux_pos);
					// DEBUG recursion
					//max_distance = pattern_length - EWAVEFRONT_V(k,offset) + text_length - EWAVEFRONT_H(k,offset);
					// printf("DPU printing max distance: %d\n", max_distance);
					//printf("DPU Solution distance %d\n", distance);
					//printf("DPU right flank length %d\n", pattern_length - EWAVEFRONT_V(k,offset));
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
											cma_cigar_aux, cigar_aux_pos);
				}
				break;
			}
			// DEBUG recursivity
			//printf("DPU loop number %d\n", distance_fw + distance_rv);

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
	return distance_fw + distance_rv;
}