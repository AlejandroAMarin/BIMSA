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
	const uint64_t longest_seq           = DPU_INPUT_ARGUMENTS.longest_seq;
	const uint64_t num_pairs_per_dpu     = DPU_INPUT_ARGUMENTS.num_pairs_per_dpu;
	const ewf_offset_t threshold     	 = DPU_INPUT_ARGUMENTS.threshold;

	// Parameters
	const uint64_t num_pairs_per_tasklet = num_pairs_per_dpu / NR_TASKLETS;	
	size_t max_distance = 2 * longest_seq * MAX_ERROR;
	size_t cigar_length = 2 * longest_seq;
	ewf_offset_t cigar_length_adjustment = cigar_length % BLOCK_SIZE_INPUTS;
	if(cigar_length_adjustment != 0)
		cigar_length += (BLOCK_SIZE_INPUTS - cigar_length_adjustment);

	// Offsets per tl size calculation
	size_t offsets_size_per_tl      = 2 * max_distance * sizeof(ewf_offset_t);
	size_t offsets_size_adjustement = offsets_size_per_tl % BLOCK_SIZE;
	if(offsets_size_adjustement != 0)
		offsets_size_per_tl += (BLOCK_SIZE - offsets_size_adjustement);

	// Wavefront variables
	int32_t pattern_length, text_length;
	int32_t target_k, target_k_abs;
	ewf_offset_t target_offset, distance;
	uint32_t my_first_pair = tasklet_id * num_pairs_per_tasklet;

	// int length_size = num_pairs_per_dpu * sizeof(uint32_t);
	// int reminder = length_size % 8;
	// length_size = length_size + reminder;

	// Memory addresses
	const uint32_t ma_start = (uint32_t) DPU_MRAM_HEAP_POINTER + BLOCK_SIZE_INPUTS;
	uint32_t       ma_pattern_start, ma_text_start;
	uint32_t 	   ma_pattern_lengths = ma_start + num_pairs_per_dpu * 2 * longest_seq * sizeof(char) + // patterns and texts
								  		(tasklet_id * num_pairs_per_tasklet * sizeof(uint32_t)); // my offset inside pattern lengths
	uint32_t 	   ma_text_lengths    = ma_pattern_lengths + (NR_TASKLETS * num_pairs_per_tasklet * sizeof(uint32_t));
	//uint32_t       ma_results         = ma_text_lengths + (NR_TASKLETS * num_pairs_per_tasklet * sizeof(uint32_t));
	uint32_t       ma_results         = ma_start +
											 (num_pairs_per_dpu * (2 * longest_seq) * sizeof(char)) +  			 // Sequences (patterns and texts)
											 (2 * num_pairs_per_dpu * sizeof(uint32_t)) + 				         // Sequences' lengths
											 (tasklet_id * num_pairs_per_tasklet * sizeof(int32_t));		 // My offset inside results
	uint32_t       ma_wf_fw      = (ma_results - (tasklet_id * num_pairs_per_tasklet * sizeof(int32_t))) + 	 // Sequences' lengths starting position
	                                         (num_pairs_per_dpu * sizeof(int32_t)) + 				       	 // Results
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
	uint32_t       ma_tasks = (ma_cigar_aux - (tasklet_id * (cigar_length + 8))) +
											 (NR_TASKLETS * (cigar_length + 8)) + 	// starting position
								  			 (tasklet_id * QUEUE_SZ); // My offset
	

	// DEBUG mem pos
	//if(tasklet_id == 0)	printf("DPU ma cigar mem pos: %d\n", ma_cigar -(uint32_t)DPU_MRAM_HEAP_POINTER);
	//	printf("DPU results mem pos: %d\n", ma_results - ma_start);
	//	printf("DPU threshold: %d\n", threshold);

	uint32_t       cma_cigar = ma_cigar;
	uint32_t       cma_cigar_aux = ma_cigar_aux + cigar_length*sizeof(char)+8;
	uint32_t       cma_tasks = ma_tasks;

	//printf("DPU ma pattern len %d ma text len %d\n", ma_pattern_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER, ma_text_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER);
	//printf("DPU ma  %d\n", ma_results -(uint32_t)DPU_MRAM_HEAP_POINTER);


	const uint32_t loop_limit            = BLOCK_SIZE 	     / sizeof(ewf_offset_t);
	const int32_t task_limit            = BLOCK_SIZE 	     / sizeof(pair_meta_t);

#if PRINT
 	if(tasklet_id == 0)
    {
 		//printf("[INFO] Memory usage per DPU: %u/%d bytes %d remaining.\n", (ma_tasks + (NR_TASKLETS * QUEUE_SZ)), MRAM_LIMIT, MRAM_LIMIT - (ma_tasks + (NR_TASKLETS * QUEUE_SZ))); 	// starting position
		// printf("LOOP LIMIT %u BLOCK_SIZE %d offset size %d\n", loop_limit, BLOCK_SIZE, sizeof(ewf_offset_t));
		// printf("TASK LIMIT %u BLOCK_SIZE %d offset size %d\n", task_limit, BLOCK_SIZE, sizeof(pair_meta_t));
 	}
#endif

	if(ma_tasks + (NR_TASKLETS * QUEUE_SZ) > MRAM_LIMIT){
		printf("[ERROR] MRAM memory exceeded\n");
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
	pair_meta_t  *cache_tasks    		= (pair_meta_t *)  mem_alloc(BLOCK_SIZE);

	// Memory access helper variables
	uint16_t       read_idx_patt_length;
	uint16_t       read_idx_text_length;
	uint16_t       read_idx_distances;
	int 		   cigar_aux_pos;
	int16_t        task_idx;
	pair_meta_t	   pair_metadata;
	int 		   loop;
	//int			   max_loop;
	
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
		task_idx = 0;
		cigar_aux_pos = 0;
		loop = 1;
		//max_loop = 0;

		// if(tasklet_id == 0){
		// 	printf("DPU ma patterns %d\n", ma_pattern_start -(uint32_t)DPU_MRAM_HEAP_POINTER);
		// 	printf("DPU ma texts %d\n", ma_text_start -(uint32_t)DPU_MRAM_HEAP_POINTER);
		// }

		// Obtain current pattern and text lengths
		pattern_length = cache_pattern_lengths [read_idx_patt_length];
		text_length    = cache_text_lengths    [read_idx_text_length];

		if (pattern_length + text_length == 0){
			//printf("Plen and tlen %d and %d \n", pattern_length, text_length);
			if(read_idx_distances > 0) // Min numpairs per tasklet to compute
			{
				mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
			}
			break;
		}

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

		distance = find_breakpoint_iterative(pattern_length, text_length,
					cache_wf,
					cache_wf_fw, cache_wf_rv,
					cache_wf_fw_next, cache_wf_rv_next,
					cache_pattern, cache_text,
					ma_wf_fw, ma_wf_fw_next, 
					ma_wf_rv, ma_wf_rv_next,
					ma_pattern_start, ma_text_start,
					loop_limit, offsets_size_per_tl, threshold,
					&cma_cigar, cache_cigar, cache_cigar_aux,
					&cma_cigar_aux, &cigar_aux_pos, cache_tasks, &task_idx, &cma_tasks,
					task_limit, max_distance
					);
		if(distance == -1){
			loop = 0;
			printf("[DPU] first breakpoint exceeds mas distance");
		} 
		if(task_idx<0 && cma_tasks == ma_tasks) loop = 0;
		// if(distance < 7){
		 	//printf("DPU distance %d\n", distance);
		// 	printf("Plen and tlen %d and %d \n", pattern_length, text_length);
		// }

		while(loop){
			// read a task
			if(task_idx<0 && cma_tasks == ma_tasks) {
				break;
			}
			
			pair_metadata = cache_tasks[task_idx];
			
			text_length = pair_metadata.t_len;
			pattern_length = pair_metadata.p_len;
			ma_pattern_start = pair_metadata.ma_p_start;
			ma_text_start = pair_metadata.ma_t_start;

			// DEBUG task
			//printf("\nDPU looping task %d in block %d tasklet id %d\n", task_idx, (ma_tasks - cma_tasks)/BLOCK_SIZE, tasklet_id);
			// printf("DPU CONSUMING TASK ma_p_s %d ma_t_s %d ma_p_rev_s %d ma_t_rev_s %d plen %d tlen %d\n", pair_metadata.ma_p_start, 
			// 		pair_metadata.ma_t_start, pair_metadata.ma_rev_p_start, pair_metadata.ma_rev_t_start, pair_metadata.p_len, pair_metadata.t_len);

			find_breakpoint_iterative(pattern_length, text_length,
								cache_wf,
								cache_wf_fw, cache_wf_rv,
								cache_wf_fw_next, cache_wf_rv_next,
								cache_pattern, cache_text,
								ma_wf_fw, ma_wf_fw_next, 
								ma_wf_rv, ma_wf_rv_next,
								ma_pattern_start, ma_text_start,
								loop_limit, offsets_size_per_tl, threshold,
								&cma_cigar, cache_cigar, cache_cigar_aux,
								&cma_cigar_aux, &cigar_aux_pos, cache_tasks, &task_idx, &cma_tasks,
								task_limit, max_distance
								);
			if(task_idx<0 && cma_tasks == ma_tasks) {
				break;
			}

			// read if cache empty
			if(task_idx<0){
				cma_tasks -= BLOCK_SIZE;
				mram_read((__mram_ptr void *) cma_tasks, cache_tasks, BLOCK_SIZE);
				task_idx = (task_limit)-1;
			}

		}
		
		if ((cigar_aux_pos)==BLOCK_SIZE_INPUTS-1)
		{
			mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, BLOCK_SIZE_INPUTS);
			(cigar_aux_pos)=-1;
			cma_cigar += BLOCK_SIZE_INPUTS;
		}
		cache_cigar_aux[cigar_aux_pos+1] = '\0';

		mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, BLOCK_SIZE_INPUTS);				
		
		// Write back te results to HOST
		// cache_distances[0] = (int64_t) distance;
		// mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
		// ma_results += 8;
		// Write back te results to HOST
		cache_distances[read_idx_distances] = distance;
		//printf("DPU distance = %d\n", distance);
		read_idx_distances++;
		if(read_idx_distances == 2) // Min numpairs per tasklet to compute
		{
			read_idx_distances = 0;
			mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
			ma_results += 8;
		}

	}

	#if PERF
	result->counter = timer_stop(&counter);
	#endif
	return 0;
}
