/*
 * WFA Distance implementation with multiple tasklets
 *
 */
#include <stdint.h>
#include <stdio.h>
#include <defs.h>
#include <mram.h>
#include <alloc.h>
#include <mram.h>
#include <seqread.h>
#include <barrier.h>
#include "common.h"
#include "ulsapim.h"
#include <perfcounter.h>
#include "counterperf.h"


__host dpu_arguments_t DPU_INPUT_ARGUMENTS;

__host dpu_results_t DPU_RESULTS[NR_TASKLETS];

BARRIER_INIT(my_barrier, NR_TASKLETS);

extern int main_kernel1(void);

int(*kernels[nr_kernels])(void) = {main_kernel1};


int main(void){
	// Kernel
	return kernels[DPU_INPUT_ARGUMENTS.kernel]();
}


// main_kernel1
int main_kernel1() {
	uint32_t tasklet_id = me();
#if PRINT
	//printf("tasklet_id = %u\n", tasklet_id);
#endif
	if(tasklet_id == 0){
		mem_reset(); // Reset the heap
		#if PERF_CYCLES
		perfcounter_config(COUNT_CYCLES, true);
		#endif
		#if PERF_INSTRUCTIONS
		perfcounter_config(COUNT_INSTRUCTIONS, true);
		#endif
	}
	// Barrier
	barrier_wait(&my_barrier);

	//printf("DPU STARTING STACK SIZE %d\n", check_stack());

	#if PERF
		perfcounter_cycles counter;
		dpu_results_t *result = &DPU_RESULTS[tasklet_id];
		result->counter = 0;
		timer_start(&counter);
	#endif

	// Arguments
	const uint32_t longest_seq           = DPU_INPUT_ARGUMENTS.longest_seq;
	const uint32_t num_pairs_per_dpu     = DPU_INPUT_ARGUMENTS.num_pairs_per_dpu;
	const ewf_offset_t threshold     	 = DPU_INPUT_ARGUMENTS.threshold;

	// Parameters
	const uint32_t num_pairs_per_tasklet = num_pairs_per_dpu / NR_TASKLETS;
	ewf_offset_t   max_distance = 2 * longest_seq;
	ewf_offset_t cigar_length = max_distance;
	ewf_offset_t cigar_length_adjustment = cigar_length % BLOCK_SIZE_INPUTS;
	if(cigar_length_adjustment != 0)
		cigar_length += (BLOCK_SIZE_INPUTS - cigar_length_adjustment);

	// Offsets per tl size calculation
	uint32_t offsets_size_per_tl      = 2 * max_distance * sizeof(ewf_offset_t); //TODO double check if cigar distance or max_length
	uint32_t offsets_size_adjustement = offsets_size_per_tl % BLOCK_SIZE;
	if(offsets_size_adjustement != 0)
		offsets_size_per_tl += (BLOCK_SIZE - offsets_size_adjustement);

	// Wavefront variables
	int32_t pattern_length, text_length;
	int32_t target_k, target_k_abs;
	ewf_offset_t target_offset, distance;
	uint32_t my_first_pair = tasklet_id * num_pairs_per_tasklet;
		
	int length_size = num_pairs_per_tasklet * sizeof(uint32_t);
	int reminder = length_size % 8;
	length_size = length_size + reminder;

	// Memory addresses
	const uint32_t ma_start = (uint32_t) DPU_MRAM_HEAP_POINTER + BLOCK_SIZE_INPUTS;
	uint32_t       ma_pattern_start, ma_text_start;
	uint32_t       ma_results         = ma_start +
											 (num_pairs_per_dpu * (2 * longest_seq) * sizeof(char)) +  			 // Sequences (patterns and texts)
											 (2 * length_size) + 				         // Sequences' lengths
											 (tasklet_id * num_pairs_per_tasklet * sizeof(int64_t));		 // My offset inside results
	uint32_t       ma_wf_fw      = (ma_results - (tasklet_id * num_pairs_per_tasklet * sizeof(int64_t))) + 	 // Sequences' lengths starting position
	                                         (num_pairs_per_dpu * sizeof(int64_t)) + 				       	 // Results
					                         (tasklet_id * offsets_size_per_tl);                                 // My offset inside wavefronts
	uint32_t       ma_wf_fw_next = (ma_wf_fw - (tasklet_id * offsets_size_per_tl)) + 		 			 // wavefronts forward starting position
	                                       	 (NR_TASKLETS * offsets_size_per_tl) + 			       				 // forward_next wavefronts
					       					 (tasklet_id * offsets_size_per_tl);                             	 // My offset inside forward_next wavefronts
	uint32_t       ma_wf_rv 	  = (ma_wf_fw_next - (tasklet_id * offsets_size_per_tl)) + 		 	 // wavefronts forward_next starting position
	                                       	 (NR_TASKLETS * offsets_size_per_tl) + 			       				 // reverse wavefronts
					       					 (tasklet_id * offsets_size_per_tl);                             	 // My offset inside reverse wavefronts
	uint32_t       ma_wf_rv_next = (ma_wf_rv - (tasklet_id * offsets_size_per_tl)) + 		 	 			 // wavefronts reverse starting position
	                                       	 (NR_TASKLETS * offsets_size_per_tl) + 			       				 // reverse_next wavefronts
					       					 (tasklet_id * offsets_size_per_tl);                             	 // My offset inside reverse_next wavefronts
	uint32_t       ma_cigar = (ma_wf_rv_next - (tasklet_id * offsets_size_per_tl)) + 		 	 			 // wavefronts reverse starting position
											 (NR_TASKLETS * offsets_size_per_tl) + 			       				 // cigars 
											 (tasklet_id * cigar_length * num_pairs_per_tasklet); // My offset inside ma_cigar wavefronts
	uint32_t       ma_cigar_aux = (ma_cigar - (tasklet_id * cigar_length * num_pairs_per_tasklet)) +
											 (NR_TASKLETS * cigar_length * num_pairs_per_tasklet) + 	// starting position
								  			 (tasklet_id * (cigar_length + 8)); // My offset	

	// DEBUG mem pos
	//if(tasklet_id == 0)	printf("DPU ma cigar mem pos: %d\n", ma_cigar - ma_start);
	//	printf("DPU results mem pos: %d\n", ma_results - ma_start);
	//	printf("DPU threshold: %d\n", threshold);

	uint32_t       cma_cigar = ma_cigar;
	uint32_t       cma_cigar_aux = ma_cigar_aux + cigar_length*sizeof(char)+8;

  	uint32_t ma_pattern_lengths = ma_start + num_pairs_per_dpu * (2 * longest_seq * sizeof(char)) + 
								  (tasklet_id * length_size);
	uint32_t ma_text_lengths    = ma_pattern_lengths + (NR_TASKLETS * length_size);

	//printf("DPU ma pattern len %d ma text len %d\n", ma_pattern_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER, ma_text_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER);

	const uint32_t loop_limit            = BLOCK_SIZE 	     / sizeof(ewf_offset_t);

#if PRINT
 	if(tasklet_id == 0)
     {
 		//printf("Memory usage per DPU: %u bytes.\n", ma_cigar_aux + (NR_TASKLETS * cigar_length * num_pairs_per_tasklet)); 	// starting position
		printf("[INFO] Memory usage per DPU: %u/%d bytes %d remaining.\n", (ma_cigar_aux + (NR_TASKLETS * cigar_length * num_pairs_per_tasklet)), MRAM_LIMIT, MRAM_LIMIT - (ma_cigar_aux + (NR_TASKLETS * cigar_length * num_pairs_per_tasklet)));      // starting position
 	}
#endif

	if(ma_cigar_aux + (NR_TASKLETS * cigar_length * num_pairs_per_tasklet) > MRAM_LIMIT){ 
		printf("[WARNING] MRAM memory can be exceeded\n");
		return 1;
	}


	// Allocate local caches to store the MRAM data // WRAM = 64k for all tasklets, in reality it is 40k beacuse of tasklets stack
	ewf_offset_t *cache_wf      		= (ewf_offset_t *) mem_alloc(BLOCK_SIZE*4); // multiplied by 4 because elements are 
	ewf_offset_t *cache_wf_fw 			= (ewf_offset_t *) &cache_wf[0];
	ewf_offset_t *cache_wf_fw_next 		= (ewf_offset_t *) &cache_wf[loop_limit];
	ewf_offset_t *cache_wf_rv      		= (ewf_offset_t *) &cache_wf[loop_limit*2]; // TODO unify cache reverse and forward (only when compute --> extend)
	ewf_offset_t *cache_wf_rv_next 		= (ewf_offset_t *) &cache_wf[loop_limit*3];
	char         *cache_pattern         = (char *)         mem_alloc(BLOCK_SIZE_INPUTS);
	char         *cache_text            = (char *)         mem_alloc(BLOCK_SIZE_INPUTS);
	char         *cache_cigar           = (char *)         mem_alloc(BLOCK_SIZE_INPUTS);
	int32_t      *cache_pattern_lengths = (int32_t *)      mem_alloc(BLOCK_SIZE_INPUTS);
	int32_t      *cache_text_lengths    = (int32_t *)      mem_alloc(BLOCK_SIZE_INPUTS);
	char 	 	 *cache_cigar_aux 		= (char *) 	       mem_alloc(BLOCK_SIZE_INPUTS);
	int32_t      *cache_distances       = (int32_t *)      mem_alloc(8);

	// Memory access helper variables
	uint16_t       read_idx_patt_length;
	uint16_t       read_idx_text_length;
	uint16_t       read_idx_distances;
	int 		   cigar_aux_pos;
	
	const uint32_t num_ints_inputs_block = BLOCK_SIZE_INPUTS / sizeof(int32_t);
	read_idx_distances = 0;

	// Initialize pattern lenghts cache
	mram_read_aligned(&ma_pattern_lengths, cache_pattern_lengths, &read_idx_patt_length, BLOCK_SIZE_INPUTS, TYPE_BYTES_SIZE);
	ma_pattern_lengths += BLOCK_SIZE_INPUTS;

	// Initialize text lenghts cache
	mram_read_aligned(&ma_text_lengths, cache_text_lengths, &read_idx_text_length, BLOCK_SIZE_INPUTS, TYPE_BYTES_SIZE);
	ma_text_lengths += BLOCK_SIZE_INPUTS;
	// Iterate over pairs
	for(uint32_t pair = 0; pair < num_pairs_per_tasklet; pair++)
	{
		// Initialize memory addresses
		ma_pattern_start 	= ma_start         	  + ((my_first_pair + pair) * longest_seq * sizeof(char));
		ma_text_start   	= ma_pattern_start + (num_pairs_per_dpu  * longest_seq * sizeof(char));
		cma_cigar 			= ma_cigar 			  + (pair * cigar_length * sizeof(char));
		cigar_aux_pos = 0;

		// Obtain current pattern and text lengths
		pattern_length = cache_pattern_lengths [read_idx_patt_length];
		text_length    = cache_text_lengths    [read_idx_text_length];

		// Increment pattern lengths index inside caches and fetch from MRAM if necessary
		read_idx_patt_length++;
		if(read_idx_patt_length == num_ints_inputs_block)
		{
			mram_read((__mram_ptr void *) ma_pattern_lengths, cache_pattern_lengths, BLOCK_SIZE_INPUTS);
			ma_pattern_lengths += BLOCK_SIZE_INPUTS;
			read_idx_patt_length = 0;
		}

		// Increment text lengths index inside caches and fetch from MRAM if necessary
		read_idx_text_length++;
		if(read_idx_text_length == num_ints_inputs_block)
		{
			mram_read((__mram_ptr void *) ma_text_lengths, cache_text_lengths, BLOCK_SIZE_INPUTS);
			ma_text_lengths += BLOCK_SIZE_INPUTS;
			read_idx_text_length = 0;
		}

		// Initialize some parameters and variables
		target_k      = EWAVEFRONT_DIAGONAL(text_length, pattern_length);
		target_k_abs  = ABS(target_k);
		target_offset = EWAVEFRONT_OFFSET(text_length, pattern_length);

		distance = find_breakpoint_recursive(pattern_length, text_length,
						cache_wf,
						cache_wf_fw, cache_wf_rv,
						cache_wf_fw_next, cache_wf_rv_next,
						cache_pattern, cache_text,
						ma_wf_fw, ma_wf_fw_next, 
						ma_wf_rv, ma_wf_rv_next,
						ma_pattern_start, ma_text_start,
						loop_limit, offsets_size_per_tl, threshold,
						&cma_cigar, cache_cigar, cache_cigar_aux,
						&cma_cigar_aux, &cigar_aux_pos
						);
		
		if ((cigar_aux_pos)==BLOCK_SIZE_INPUTS-1)
		{
			mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, BLOCK_SIZE_INPUTS);
			(cigar_aux_pos)=-1;
			cma_cigar += BLOCK_SIZE_INPUTS;
		}
		cache_cigar_aux[cigar_aux_pos+1] = '\0';

		mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, BLOCK_SIZE_INPUTS);	

		// DEBUG cigar
		// for (int i = 0; i < BLOCK_SIZE_INPUTS; i++)
		// {
		// 	printf("%c\n", cache_cigar_aux[i]);
		// }
			
		// Write back te results to HOST
		// cache_distances[read_idx_distances] = distance;
		// read_idx_distances++;
		// if(read_idx_distances == 2) // Min numpairs per tasklet to compute
		// {
		// 	read_idx_distances = 0;
		// 	mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
		// 	ma_results += 8;
		// }
		cache_distances[0] = (int64_t) distance;
		//read_idx_distances++;
		// if(read_idx_distances == 2) // Min numpairs per tasklet to compute
		// {
			//read_idx_distances = 0;
		mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
		ma_results += 8;
		//}
	}
	#if PERF
	result->counter = timer_stop(&counter);
	#endif
	return 0;
}
