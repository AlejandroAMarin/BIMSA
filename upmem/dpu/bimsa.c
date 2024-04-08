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
#include "bimsa.h"
#include "mutex.h"

__host dpu_arguments_t DPU_INPUT_ARGUMENTS;

__host dpu_results_t DPU_RESULTS[NR_TASKLETS];

BARRIER_INIT(my_barrier, NR_TASKLETS);

extern int main_kernel1(void);

MUTEX_INIT(my_mutex);

uint32_t global_pair = 0;

int(*kernels[nr_kernels])(void) = {main_kernel1};

uint32_t get_next_pair(const uint32_t offset){
	uint32_t pair;
	mutex_lock(my_mutex);
		pair = global_pair;
		global_pair += offset;
	mutex_unlock(my_mutex);

	return pair;
}


int main(void){
	// Kernel
	return kernels[DPU_INPUT_ARGUMENTS.kernel]();
}


// main_kernel1
int main_kernel1() {
	uint32_t tasklet_id = me();
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
		result->main = 0;	
	#endif
	#if PERF_MAIN
		timer_start(&counter);
	#endif

	// Arguments
	const uint32_t longest_seq            	   = DPU_INPUT_ARGUMENTS.longest_seq;
	const uint32_t num_pairs_per_dpu       	   = DPU_INPUT_ARGUMENTS.num_pairs_per_dpu;
	const uint32_t batch_pairs_per_dpu 		   = DPU_INPUT_ARGUMENTS.batch_pairs_per_dpu;
	const uint32_t batch_idx     		   	   = DPU_INPUT_ARGUMENTS.batch_idx;
	const ewf_offset_t threshold     	   	   = DPU_INPUT_ARGUMENTS.threshold;

	uint32_t spare_tasklet_pairs = 0;
	uint32_t spare_dpu_pairs = 0;
	if(batch_pairs_per_dpu * (batch_idx+1) > num_pairs_per_dpu){
		spare_dpu_pairs = ((batch_pairs_per_dpu * (batch_idx+1)) - (num_pairs_per_dpu));
		spare_tasklet_pairs =  spare_dpu_pairs / NR_TASKLETS;
	}
	// Parameters
	//const uint32_t num_pairs_per_tasklet = num_pairs_per_dpu / NR_TASKLETS;
	//const uint32_t batch_pairs_per_tasklet = batch_pairs_per_dpu / NR_TASKLETS;
	#if DYNAMIC
		const uint32_t offset = 2;
	#else
		const uint32_t offset = batch_pairs_per_dpu / NR_TASKLETS;
	#endif
	//uint32_t my_first_pair = tasklet_id * (batch_pairs_per_tasklet - spare_tasklet_pairs); // first pair for all tasklets in first batch
	// uint32_t pair = batch_pairs_per_dpu * batch_idx + // What batch I'm in
	// 				(batch_pairs_per_tasklet - spare_tasklet_pairs) * tasklet_id; // Offset pair inside batch
	//uint32_t ending_pair = batch_pairs_per_dpu * (batch_idx+1) - // Where does the next batch start and this one ends
	//					   (batch_pairs_per_tasklet - spare_tasklet_pairs) * tasklet_id; // Offset pair inside batch
	//uint32_t ending_pair = pair + (batch_pairs_per_tasklet - spare_tasklet_pairs);

	int32_t max_distance = 2 * longest_seq * MAX_ERROR;
	//int64_t max_distance = 2 * 100;
	int32_t cigar_length = 2 * longest_seq;
	ewf_offset_t cigar_length_adjustment = cigar_length % CIGAR_TRANSFER;
	if(cigar_length_adjustment != 0)
		cigar_length += (CIGAR_TRANSFER - cigar_length_adjustment);

	// Offsets per tl size calculation
	int32_t offsets_size_per_tl      = 2 * max_distance * sizeof(ewf_offset_t);
	int32_t offsets_size_adjustement = offsets_size_per_tl % WF_TRANSFER;
	if(offsets_size_adjustement != 0)
		offsets_size_per_tl += (WF_TRANSFER - offsets_size_adjustement);

	// Wavefront variables
	int32_t pattern_length, text_length;
	int32_t target_k, target_k_abs;
	ewf_offset_t target_offset, distance;

	// Memory addresses
	const uint32_t ma_start = (uint32_t) DPU_MRAM_HEAP_POINTER + SEQ_TRANSFER;
	uint32_t       ma_pattern_start, ma_text_start;
	uint32_t	   ma_pattern_lengths_start = ma_start + num_pairs_per_dpu * 2 * longest_seq * sizeof(char);
	//uint32_t 	   ma_pattern_lengths = ma_pattern_lengths_start + // patterns and texts
	//							  		(tasklet_id * num_pairs_per_tasklet * sizeof(uint32_t)); // my offset inside pattern lengths
	//uint32_t 	   ma_text_lengths    = ma_pattern_lengths + (NR_TASKLETS * num_pairs_per_tasklet * sizeof(uint32_t));
	uint32_t       ma_results_start   = ma_pattern_lengths_start + (num_pairs_per_dpu * sizeof(uint32_t))*2;
	// uint32_t       ma_uint32_t       ma_results        = ma_start +
	// 										 (num_pairs_per_dpu * (2 * longest_seq) * sizeof(char)) +  			 // Sequences (patterns and texts)
	// 										 (2 * num_pairs_per_dpu * sizeof(uint32_t)) + 				         // Sequences' lengths
	// 										 (tasklet_id * num_pairs_per_tasklet * sizeof(int32_t));		 // My offset inside results
	uint32_t       ma_wf_fw      = ma_results_start + (num_pairs_per_dpu * sizeof(int32_t)) + 				       	 // Results
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
	uint32_t	   ma_cigar_start = (ma_wf_rv_next - (tasklet_id * offsets_size_per_tl)) + 		 	 			 // wavefronts reverse starting position
											 (NR_TASKLETS * offsets_size_per_tl);
	uint32_t       ma_cigar_aux = ma_cigar_start + (num_pairs_per_dpu * cigar_length) + 	// starting position
								  			 (tasklet_id * (cigar_length + 8)); // My offset
	uint32_t       ma_tasks = (ma_cigar_aux - (tasklet_id * (cigar_length + 8))) +
											 (NR_TASKLETS * (cigar_length + 8)) + 	// starting position
								  			 (tasklet_id * QUEUE_SZ); // My offset
	

	// DEBUG mem pos
	//if(tasklet_id == 0)	printf("DPU ma cigar mem pos: %d\n", ma_cigar -(uint32_t)DPU_MRAM_HEAP_POINTER);
	//	printf("DPU results mem pos: %d\n", ma_results - ma_start);
	//	printf("DPU threshold: %d\n", threshold);

	uint32_t cma_cigar;
	uint32_t cma_cigar_aux;
	uint32_t cma_tasks;
	uint32_t ma_results;
	uint32_t ma_pattern_lengths;
	uint32_t ma_text_lengths;
	
	//printf("DPU ma pattern len %d ma text len %d\n", ma_pattern_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER, ma_text_lengths -(uint32_t)DPU_MRAM_HEAP_POINTER);
	//printf("DPU ma  %d\n", ma_results -(uint32_t)DPU_MRAM_HEAP_POINTER);

	const uint32_t wf_limit              = WF_TRANSFER / sizeof(ewf_offset_t);
	const int32_t task_limit             = TASK_TRANSFER / sizeof(pair_meta_t);
	const uint32_t num_ints_inputs_block = LEN_TRANSFER / sizeof(int32_t);
	const uint32_t base_limit 			 = (4*WF_TRANSFER)/sizeof(uint32_t);

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
	ewf_offset_t *cache_wf      		= (ewf_offset_t *) mem_alloc(WF_TRANSFER*4); // multiplied by 4 because elements are 
	ewf_offset_t *cache_wf_fw 			= (ewf_offset_t *) &cache_wf[0];
	ewf_offset_t *cache_wf_fw_next 		= (ewf_offset_t *) &cache_wf[wf_limit];
	ewf_offset_t *cache_wf_rv      		= (ewf_offset_t *) &cache_wf[wf_limit*2]; // TODO unify cache reverse and forward (only when compute --> extend)
	ewf_offset_t *cache_wf_rv_next 		= (ewf_offset_t *) &cache_wf[wf_limit*3];
	char         *cache_pattern         = (char *)         mem_alloc(SEQ_TRANSFER);
	char         *cache_text            = (char *)         mem_alloc(SEQ_TRANSFER);
	int32_t      *cache_pattern_lengths = (int32_t *)      mem_alloc(LEN_TRANSFER);
	int32_t      *cache_text_lengths    = (int32_t *)      mem_alloc(LEN_TRANSFER);
	char         *cache_cigar           = (char *)         mem_alloc(CIGAR_TRANSFER);
	char 	 	 *cache_cigar_aux 		= (char *) 	       mem_alloc(CIGAR_TRANSFER);
	int32_t      *cache_distances       = (int32_t *)      mem_alloc(8);
	pair_meta_t  *cache_tasks    		= (pair_meta_t *)  mem_alloc(TASK_TRANSFER);

	// Memory access helper variables
	uint16_t       read_idx_length = 0;
	uint16_t       read_idx_distances = 0;
	int 		   cigar_aux_pos;
	int16_t        task_idx;
	pair_meta_t	   pair_metadata;
	int 		   loop;
	uint32_t pair;
	//int			   max_loop;

	//printf("[DPU] batch pairs per tl %d batch pairs per dpu %lu total pairs in dpus %lu reminder %d\n", batch_pairs_per_tasklet, batch_pairs_per_dpu * (batch_idx), num_pairs_per_dpu * NR_DPUS, spare_tasklet_pairs);
	// printf("[DPU] pair %lu last pair %lu batch idx %d reminder_pairs %lu\n", batch_pairs_per_tasklet * batch_idx, batch_pairs_per_tasklet * (batch_idx+1), batch_idx, reminder_pairs);
	// Iterate over pairs
	while(1)
	{
		pair = get_next_pair(offset);
		if(pair >= batch_pairs_per_dpu) break;

	for(uint32_t i = 0; i < offset; i++)
	{
		printf("[DPU] TL %d pair %d\n", me(), pair);
		cma_cigar_aux = ma_cigar_aux + cigar_length*sizeof(char)+8;
		cma_tasks = ma_tasks;
		//printf("[DPU] OG ma_plen %d ma_tlen %d\n", ma_pattern_lengths, ma_text_lengths);
		// Initialize memory addresses
		ma_pattern_start 	= ma_start         	 	   + (pair * longest_seq * sizeof(char));
		ma_text_start   	= ma_pattern_start 	 	   + (num_pairs_per_dpu  * longest_seq * sizeof(char));
		ma_pattern_lengths 	= ma_pattern_lengths_start + (pair * sizeof(int32_t));
		ma_text_lengths		= ma_pattern_lengths 	   + (num_pairs_per_dpu * sizeof(int32_t));
		cma_cigar 			= ma_cigar_start 		   + (pair * cigar_length * sizeof(char));
		ma_results 			= ma_results_start 		   + (pair * sizeof(int32_t));
		task_idx = 0;
		cigar_aux_pos = 0;
		loop = 1;
		//max_loop = 0;
		//printf("[DPU] NEW ma_plen %d ma_tlen %d\n", ma_pattern_lengths, ma_text_lengths);

		// if(tasklet_id == 0){
		// 	printf("DPU ma patterns %d\n", ma_pattern_start -(uint32_t)DPU_MRAM_HEAP_POINTER);
		// 	printf("DPU ma texts %d\n", ma_text_start -(uint32_t)DPU_MRAM_HEAP_POINTER);
		// }

		// Increment lengths index inside caches and fetch from MRAM if necessary
		if(read_idx_length == 0)
		{
			mram_read_aligned(&ma_pattern_lengths, cache_pattern_lengths, &read_idx_length, LEN_TRANSFER, TYPE_BYTES_SIZE);
			mram_read_aligned(&ma_text_lengths, cache_text_lengths, &read_idx_length, LEN_TRANSFER, TYPE_BYTES_SIZE);
			// ma_pattern_lengths += LEN_TRANSFER;
			// ma_text_lengths += LEN_TRANSFER;
		}
		// Obtain current pattern and text lengths
		pattern_length = cache_pattern_lengths [read_idx_length];
		text_length    = cache_text_lengths    [read_idx_length];
		//printf("[DPU] tasklet %d processing pair %d/%d total pairs per tl %d total pairs per dpu %lu | batch idx %d \n", me(), pair, ending_pair, (batch_pairs_per_tasklet - spare_tasklet_pairs), batch_pairs_per_dpu, batch_idx);
		//printf("[DPU] processing pair %d/%lu plen = %d tlen = %d idx %d\n", pair, batch_pairs_per_tasklet * tasklet_id + batch_pairs * (batch_idx+1) - spare_tasklet_pairs, pattern_length, text_length, read_idx_length);
		read_idx_length++;
		if(read_idx_length >= num_ints_inputs_block) read_idx_length = 0;

		if (pattern_length + text_length <= 0){
			//printf("Plen and tlen %d and %d \n", pattern_length, text_length);
			if(read_idx_distances > 0) // Min numpairs per tasklet to compute
			{
				mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
			}
			#if PRINT
			printf("[DPU] plen and tlen are 0\n");
			#endif
			break;
		}

		// Initialize some parameters and variables
		target_k      = EWAVEFRONT_DIAGONAL(text_length, pattern_length);
		target_k_abs  = ABS(target_k);
		target_offset = EWAVEFRONT_OFFSET(text_length, pattern_length);

		// profiling_start(&first_breakpoint);
		#if PERF_MAIN
		result->main += timer_stop(&counter);
		timer_start(&counter);
		#endif
		distance = find_breakpoint_iterative(pattern_length, text_length,
					cache_wf,
					cache_wf_fw, cache_wf_rv,
					cache_wf_fw_next, cache_wf_rv_next,
					cache_pattern, cache_text,
					ma_wf_fw, ma_wf_fw_next, 
					ma_wf_rv, ma_wf_rv_next,
					ma_pattern_start, ma_text_start,
					wf_limit, base_limit, offsets_size_per_tl, threshold,
					&cma_cigar, cache_cigar, cache_cigar_aux,
					&cma_cigar_aux, &cigar_aux_pos, cache_tasks, &task_idx, &cma_tasks,
					task_limit, max_distance
					);
		// profiling_stop(&first_breakpoint);

		if(distance == -1){
			loop = 0;
			#if PRINT
			printf("[DPU] first breakpoint exceeds max distance\n");
			#endif
		} 
		//if(task_idx<0 && cma_tasks == ma_tasks || distance == 0) loop = 0;
		if(distance == 0) loop = 0;
		// profiling_start(&other_breakpoints);
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
			#if PERF_MAIN
			result->main += timer_stop(&counter);
			timer_start(&counter);
			#endif
			find_breakpoint_iterative(pattern_length, text_length,
								cache_wf,
								cache_wf_fw, cache_wf_rv,
								cache_wf_fw_next, cache_wf_rv_next,
								cache_pattern, cache_text,
								ma_wf_fw, ma_wf_fw_next, 
								ma_wf_rv, ma_wf_rv_next,
								ma_pattern_start, ma_text_start,
								wf_limit, base_limit, offsets_size_per_tl, threshold,
								&cma_cigar, cache_cigar, cache_cigar_aux,
								&cma_cigar_aux, &cigar_aux_pos, cache_tasks, &task_idx, &cma_tasks,
								task_limit, max_distance
								);
			if(task_idx<0 && cma_tasks == ma_tasks) {
				break;
			}

			// read if cache empty
			if(task_idx<0){
				cma_tasks -= TASK_TRANSFER;
				mram_read((__mram_ptr void *) cma_tasks, cache_tasks, TASK_TRANSFER);
				task_idx = (task_limit)-1;
			}

		}
		// profiling_stop(&other_breakpoints);
		if ((cigar_aux_pos)==CIGAR_TRANSFER-1)
		{
			mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, CIGAR_TRANSFER);
			cigar_aux_pos = 0;
			cma_cigar += CIGAR_TRANSFER;
		}
		if(cigar_aux_pos > 0) cigar_aux_pos++;
		cache_cigar_aux[cigar_aux_pos] = '\0';
			//if(distance == 0) printf("pair %d cache cigar aux 2 %s\n", pair, cache_cigar_aux);
		mram_write(cache_cigar_aux, (__mram_ptr void *) cma_cigar, CIGAR_TRANSFER);				
		
		// Write back te results to HOST
		// cache_distances[0] = (int64_t) distance;
		// mram_write(cache_distances, (__mram_ptr void *) ma_results, 8);
		// ma_results += 8;
		// Write back te results to HOST
		cache_distances[read_idx_distances] = distance;
		//printf("[DPU] TL %d pair %d distance = %d\n", me(), pair, distance);
		read_idx_distances++;
		if(read_idx_distances == 2) // Min numpairs per tasklet to compute
		{
			read_idx_distances = 0;
			mram_write(cache_distances, (__mram_ptr void *) (ma_results - sizeof(int32_t)), 8);
		}

		pair++;
	}
	}
	
	#if PERF_MAIN
	result->main += timer_stop(&counter);
	#endif
	return 0;
}
