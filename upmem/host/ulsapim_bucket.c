/**
 * app.c
 * WFA Host Application Source File
 *
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <time.h> 
#include <omp.h>

#include <dpu.h>
#include <dpu_log.h>

#include "wfa_cpu.h"

#if ENERGY
#include <dpu_probe.h>
#endif

#include "params.h"
//#include "timer.h"
#include "common.h"

// Define the DPU Binary path as DPU_BINARY here
//#define DPU_BINARY "./bin/ulsapim_iterative"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1500000



static bool read_sequence_data(char* input_filename, char** patterns, char** texts,
								uint32_t** pattern_lengths, uint32_t** text_lengths,
								uint32_t num_sets, uint64_t* num_pairs, int64_t* longest_seq, uint32_t max_pairs)
{
	FILE *file;
	int pattern_length;
	int text_length;
	uint32_t total_pairs = 0;
	size_t max_line_size = longest_seq[num_sets-1] + 1500000;
	int32_t set;
	uint64_t discarded = 0;


    // Open the file
    file = fopen(input_filename, "r");
    if (file == NULL) {
        printf("[ERROR] Input file not found!\n");
        return false;
    }

	char* line1 = malloc(max_line_size * sizeof(char));
	char* line2 = malloc(max_line_size * sizeof(char));

	while(1)
	{
		// Read queries
		pattern_length = getline(&line1, &max_line_size, file);
		text_length    = getline(&line2, &max_line_size, file);
		if (pattern_length == -1 || text_length == -1) break;

		total_pairs++;
		if(total_pairs > max_pairs) break;

		pattern_length -= 2;
		text_length -= 2;

		if(line1[0] != '>'){
			printf("[ERROR] Input file not valid\n");
        	return false;
		}
		if(line2[0] != '<'){
			printf("[ERROR] Input file not valid\n");
        	return false;
		}

		set = -1;

		for(uint32_t i=0; i<num_sets; i++){
			if(pattern_length < longest_seq[i] && text_length < longest_seq[i]){
				set = i;
				break;
			}
		}
		if (set == -1){
			printf("A pair has been discarded, with length %d pattern %d text\n", pattern_length, text_length);
			discarded++;
			continue;
		}

		strncpy(patterns[set] + (longest_seq[set] * num_pairs[set]), line1 + 1, pattern_length);
		strncpy(texts[set]    + (longest_seq[set] * num_pairs[set]), line2 + 1, text_length);

		patterns[set][(longest_seq[set] * num_pairs[set]) + pattern_length] = '\0';
		texts[set][(longest_seq[set] * num_pairs[set]) + text_length] = '\0';
		// Save sizes after removing start and end of line characters
		pattern_lengths[set][num_pairs[set]] = pattern_length;
		text_lengths[set][num_pairs[set]]    = text_length;
		num_pairs[set]++;

		if (num_pairs[set] % 1000 == 0)
		{
			pattern_lengths[set] = (uint32_t *) realloc(pattern_lengths[set],  (size_t)(num_pairs[set] + 1000) * (size_t)sizeof(uint32_t));
			text_lengths[set] = (uint32_t *) realloc(text_lengths[set],  (size_t)(num_pairs[set] + 1000) * (size_t)sizeof(uint32_t));
			patterns[set] = (char *) realloc(patterns[set],  (size_t)(num_pairs[set] + 1000) * (size_t)sizeof(char) * (size_t)longest_seq[set]);
			texts[set] = (char *) realloc(texts[set],  (size_t)(num_pairs[set] + 1000) * (size_t)sizeof(char) * (size_t)longest_seq[set]);
			memset(&pattern_lengths[set][num_pairs[set]], 0, 1000 * sizeof(uint32_t));
			memset(&text_lengths[set][num_pairs[set]],    0, 1000 * sizeof(uint32_t));
		}

	}
	if(discarded > 0)printf("[WARNING] %ld pairs have been discarded for beeing to large for the specified sets\n", discarded);
	free(line1);
	free(line2);
    return true;
}

char* cigar_pprint(char* input, FILE *output) {
    int length = strlen(input);
    char* result = (char*) malloc((2 * length + 1) * sizeof(char));

	int count = 1;
    int index = 0;
    int i;

	for (i = 1; i <= length; i++) {
        if (input[i] == input[i-1]) {
            count++;
        }
        else {
			if (count > 0){
				char countStr[12];
				sprintf(countStr, "%d", count);
				int j;
				for (j = 0; countStr[j] != '\0'; j++)
				{
					result[index++] = countStr[j];
				}

			}
            result[index++] = input[i-1];
            count = 1;
        }
    }
	result[index] = '\0'; // add null terminator
	fprintf(output, "%s\n", result);
	return result;
	//printf("%s\n", result);
}

int cigar_distance(char* input) {
    int length = strlen(input);

	int distance = 0;

	for (int i = 0; i < length; i++) {
		if (input[i] != 'M')
		{
			distance++;
		}
    }
	return distance;
}

int check_cigar(int plen, int tlen,
				char* pattern,
				char* text,
				const char* curr_cigar) {
        int text_pos = 0, pattern_pos = 0;

       int cigar_len = strlen(curr_cigar);
	   char curr_cigar_element;

        for (int i=0; i<cigar_len; i++) {
            curr_cigar_element = curr_cigar[i];
			// printf("Current cigar %c ", curr_cigar[i]);
			// printf("pattern %c and text %c\n", pattern[pattern_pos], text[text_pos]);
            if (curr_cigar_element == 'M') {
				if (pattern[pattern_pos] != text[text_pos]) {
					// printf("CIGAR: %s\n", curr_cigar);
					// printf("Alignment not matching at CCIGAR index %d"
					// 		" (pattern[%d] = %c != text[%d] = %c)\n)",
					// 		i, pattern_pos,	pattern[pattern_pos], text_pos, text[text_pos]);
					return 1;
				}
				++pattern_pos;
				++text_pos;
            } else {
                switch (curr_cigar_element) {
                    case 'I':
                        ++text_pos;
                        break;
                    case 'D':
                        ++pattern_pos;
                        break;
                    case 'X':
                        if (pattern[pattern_pos] == text[text_pos]) {
                            // printf("CIGAR: %s\n", curr_cigar);
                            // printf("Alignment not mismatching at CCIGAR index %d"
                            //       " (pattern[%d] = %c == text[%d] = %c)\n)",
                            //       i, pattern_pos, pattern[pattern_pos], text_pos, text[text_pos]);
                            return 1;
                        }
                        ++pattern_pos;
                        ++text_pos;
                        break;
                    default:
                        // printf("Invalid CIGAR generated.\n");
                        break;
                }
            }
        }

        if (pattern_pos != plen) {
            // printf("Alignment incorrect length, pattern-aligned: %d, "
            //       "pattern-length: %d.", pattern_pos, plen);
            // printf("TEXT: %s", text);
            // printf("PATTERN: %s", pattern);
            return 1;
        }

        if (text_pos != tlen) {
            // printf("Alignment incorrect length, text-aligned: %d, "
            //       "text-length: %d", text_pos, tlen);
            return 1;
        }

    return 0;
}

void run_set(uint64_t num_pairs_per_dpu, uint64_t batch_pairs_per_dpu,
			int64_t longest_seq, struct dpu_set_t* dpu_set,
			ewf_offset_t threshold, uint32_t rep,
			uint32_t n_warmup, Timer* timer, uint32_t batch_idx){

	unsigned int kernel = 0;
	uint32_t i = 0;
	struct dpu_set_t dpu;
	dpu_arguments_t input_arguments = {longest_seq, num_pairs_per_dpu, batch_pairs_per_dpu, batch_idx, threshold, kernel};
	// Copy input arguments to dpu
	DPU_FOREACH(*dpu_set, dpu) {
		DPU_ASSERT(dpu_copy_to(dpu, "DPU_INPUT_ARGUMENTS", 0, (const void *) &input_arguments, sizeof(input_arguments)));
		i++;
	}

	// Run kernel on DPUs
	if (rep >= n_warmup)
	{
		start(timer, 2);
	#if ENERGY
		DPU_ASSERT(dpu_probe_start(&probe));
	#endif
	}

	DPU_ASSERT(dpu_launch(*dpu_set, DPU_ASYNCHRONOUS)); // DPU can be used asyncronously (have into account for future improvements)
}

void send_data(char* patterns, char* texts,
			uint32_t* pattern_lengths, uint32_t* text_lengths,
			uint64_t num_pairs, uint64_t num_pairs_per_dpu,
			int64_t longest_seq, struct dpu_set_t* dpu_set,
			uint32_t* nr_of_dpus, uint32_t rep,
			uint32_t n_warmup, Timer* timer, uint64_t* mem_offset){

		printf("[INFO] " ANSI_COLOR_BLUE "num_dpus required: \t %d" ANSI_COLOR_RESET "\n", NR_DPUS);
		printf("[INFO] " ANSI_COLOR_BLUE "num_pairs: \t %ld" ANSI_COLOR_RESET "\n", num_pairs);
		printf("[INFO] " ANSI_COLOR_BLUE "longest_sequence: %ld" ANSI_COLOR_RESET "\n", longest_seq);
		printf("[INFO] Pairs per DPU %ld (adjusted): (%ld)\n", num_pairs/NR_DPUS, num_pairs_per_dpu);
		printf("[INFO] Total WRAM usage: %.2f%% \n", ((float)NR_TASKLETS * (WF_TRANSFER*4 + SEQ_TRANSFER*2 + CIGAR_TRANSFER*2 + LEN_TRANSFER*2 + TASK_TRANSFER + 8 + STACK_SIZE_DEFAULT)/(float)WRAM_LIMIT)*100);
		printf("[INFO] Reading input... "ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
		printf("[INFO] Allocating DPUs... ");
		fflush(stdout);

		DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, dpu_set));
		DPU_ASSERT(dpu_load(*dpu_set, DPU_BINARY, NULL));
		DPU_ASSERT(dpu_get_nr_dpus(*dpu_set, nr_of_dpus));

		printf("\n[INFO] " ANSI_COLOR_BLUE "num_dpus allocated: %d" ANSI_COLOR_RESET "\n", *nr_of_dpus);

		uint32_t i = 0;
		struct dpu_set_t dpu;
		printf("[INFO] Transferring inputs to DPUs... ");
		fflush(stdout);

		if (rep >= n_warmup)
			start(timer, 1);

		*mem_offset = SEQ_TRANSFER; // offset points to first MRAM address
		// printf("[VERBOSE] HOST ma pattern %d\n", mem_offset[set]);
		// Patterns transfer
		DPU_FOREACH(*dpu_set, dpu, i)
		{
			DPU_ASSERT(dpu_prepare_xfer(dpu, patterns + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		//printf("[VERBOSE] HOST pattern sent %ld bytes\n", longest_seq * num_pairs_per_dpu * sizeof(char));

    	i = 0;
		*mem_offset += (longest_seq * num_pairs_per_dpu * sizeof(char)); // of
		// printf("[VERBOSE] HOST ma texts %d\n", mem_offset[set]);
		// Texts transfer
		DPU_FOREACH(*dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, texts + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		//printf("[VERBOSE] HOST text sent %ld bytes\n", longest_seq * num_pairs_per_dpu * sizeof(char));

		i = 0;
		*mem_offset += (longest_seq * num_pairs_per_dpu * sizeof(char)); // offset points to end of texts array
		// printf("[VERBOSE] HOST ma pattern lenghts %d\n", mem_offset[set]);
		// int length_size = num_pairs_per_dpu[set] * sizeof(uint32_t);
		// reminder = fmod(length_size, 8);
		// printf("reminder %f\n", reminder);
		// length_size = length_size + reminder;

		// Pattern lengths array transfer
		DPU_FOREACH(*dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, pattern_lengths + i * num_pairs_per_dpu));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, num_pairs_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));
		// printf("[VERBOSE] HOST pattern lengths sent %ld bytes\n", num_pairs_per_dpu[set] * sizeof(uint32_t));

		i = 0;
		*mem_offset += num_pairs_per_dpu * sizeof(uint32_t); // offset points to end of pattern lengths array
		// printf("[VERBOSE] HOST ma text lengths %d\n", mem_offset[set]);
		// Texts lengths array transfer
		DPU_FOREACH(*dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, text_lengths + i * num_pairs_per_dpu));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, num_pairs_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));
		// printf("[VERBOSE] HOST text lengths sent %ld bytes\n", num_pairs_per_dpu[set] * sizeof(uint32_t));

		*mem_offset += num_pairs_per_dpu * sizeof(uint32_t); // offset points to end of texts lengths array (start of results array)
		if (rep >= n_warmup)
			stop(timer, 1);
		printf(""ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
}

void retrieve_cigars(struct dpu_set_t dpu_set, uint64_t num_pairs_per_dpu, 
					uint64_t offsets_size_per_tl, uint64_t mem_offset, 
					char* cigar, uint64_t cigar_length,	Timer* timer, 
					uint32_t rep, uint32_t n_warmup){

	printf("\n[INFO] Retrieving CIGARS to host... ");
	fflush(stdout);
	struct dpu_set_t dpu;
	if (rep >= n_warmup)
		start(timer, 5);
	
	mem_offset += NR_TASKLETS * (4 * offsets_size_per_tl);
	int i = 0;
	// printf("\n[VERBOSE] HOST ma CIGAR %d\n", mem_offset[set]);
	DPU_FOREACH(dpu_set, dpu, i) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, cigar + i * num_pairs_per_dpu * cigar_length));
	}
	DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
				mem_offset, num_pairs_per_dpu * cigar_length * sizeof(char), DPU_XFER_DEFAULT));
	//printf("[VERBOSE] HOST CIGAR recieved %ld bytes\n", num_pairs_per_dpu * cigar_length * sizeof(char));

	if(rep >= n_warmup)
		stop(timer, 5);

	printf(""ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
	// #if PRINT
	// printf("LOGS\n");
	// DPU_FOREACH(dpu_set, dpu) {
	// 	DPU_ASSERT(dpu_log_read(dpu, stdout));
	// }
	// #endif
	DPU_ASSERT(dpu_free(dpu_set));
}

void retrieve_and_recover(struct dpu_set_t dpu_set,	uint32_t batch_pairs_per_dpu, uint32_t num_pairs_per_dpu,
					uint64_t* mem_offset, char* patterns, char* texts,
					uint32_t* pattern_lengths, uint32_t* text_lengths,
					uint32_t spare_dpu_pairs, int64_t longest_seq,
					int32_t* dpu_distances, Timer* timer, uint32_t rep, uint32_t n_warmup, 
					uint32_t batch_idx, uint64_t num_pairs, uint32_t* cpu_pairs_idx, 
					dpu_results_t* results, dpu_results_t* results_min, uint64_t* perf_main,
					uint64_t* perf_sum, int* perf_valid, uint64_t* over_distance){
	
	#pragma omp master
	{
	printf("\n[INFO] Retrieving distances for batch %d... ", (batch_idx));
	//printf("[DEBUG] spare dpu pairs: %d batch offset: %d\n", spare_dpu_pairs, (batch_pairs_per_dpu * batch_idx) - spare_dpu_pairs);
	fflush(stdout);
	struct dpu_set_t dpu;
	if (rep >= n_warmup)
		start(timer, 3);

	int i = 0;
	// Copy back the DPU results (mem_offset[set] already points to the results starting memory position)
	// printf("\n[VERBOSE] HOST ma distances %d\n", mem_offset[set]);
	DPU_FOREACH(dpu_set, dpu, i) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_distances + i * num_pairs_per_dpu + (batch_pairs_per_dpu * batch_idx)));
		//printf("\n[DEBUG] dpu %d writing on position %ld n elements %d\n", i,  + i * batch_pairs_per_dpu, batch_pairs_per_dpu);
		#if PERF
		dpu_results_t result;
		// Retrieve tasklet timings
		for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
			DPU_ASSERT(dpu_copy_from(dpu, "DPU_RESULTS", each_tasklet * sizeof(dpu_results_t), &result, sizeof(dpu_results_t)));
			//printf("\n[PERF] DPU kernel perf main dpu %d tasklet %d: %g\n", i, each_tasklet, (double)result.main);
			#if PERF_MAIN
			if (result.main > 100000){
				*perf_sum += result.main;
				perf_main[i * each_tasklet] += result.main;
				*perf_valid++;
			}
			if (result.main > results->main && result.main > 100000)
			 	results[i].main += result.main;
			if (result.main < results_min->main && result.main > 100000)
			 	results_min[i].main += result.main;
			#endif
			#if PERF_SECTIONS
			if (result.overlap > results[i].overlap)
				results[i].overlap += result.overlap;
			if (result.base_case > results[i].base_case)
				results[i].base_case += result.base_case;
			if (result.compute > results->compute)
				results[i].compute += result.compute;
			if (result.extend > results[i].extend)
				results[i].extend += result.extend;
			#endif
		}
		#endif
	}

	DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
				*mem_offset, (batch_pairs_per_dpu - spare_dpu_pairs) * sizeof(int32_t), DPU_XFER_DEFAULT));
	//printf("[VERBOSE] HOST distances recieved %ld bytes\n", (batch_pairs_per_dpu - spare_dpu_pairs) * sizeof(int32_t) * NR_DPUS);

	if(rep >= n_warmup)
		stop(timer, 3);

	printf(""ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");

	*mem_offset += (batch_pairs_per_dpu - spare_dpu_pairs) * sizeof(int32_t);
	//printf("[DEBUG] batch pairs adjusted %d \n", (batch_pairs_per_dpu - spare_dpu_pairs));

	}

	#pragma omp barrier

	/// CPU RECOVERY ///
	#pragma omp single
	{
	(*over_distance) = 0;
	for(int i = 0; i<NR_DPUS; i++){
		for (uint32_t j = (batch_pairs_per_dpu*(batch_idx)); j < (batch_pairs_per_dpu*(batch_idx+1)) - spare_dpu_pairs; j++)
		{
			//printf("\n[DEBUG] bad distance %d for batch %d\n", dpu_distances[j], batch_idx);
			if(dpu_distances[i * num_pairs_per_dpu + j] == -1){
				cpu_pairs_idx[(*over_distance)] = i * num_pairs_per_dpu + j;
				(*over_distance)++;
			}
		}
	}

	printf("\n[INFO] Batch %d executing CPU recovery for %ld pairs (%.2f%%)... ", batch_idx, *over_distance, ((double)*over_distance/(double)num_pairs)*100);
	fflush(stdout);
	}

	#pragma omp barrier

	if( omp_get_thread_num() > 0){
	cpu_recovery(patterns, texts, pattern_lengths, 
                    text_lengths, cpu_pairs_idx, (*over_distance), longest_seq);
	// #pragma omp single
	// {
	// printf(""ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
	// }
	/// CPU RECOVERY ///
	}
}

void validate_results(char* patterns, char* texts,
			uint32_t* pattern_lengths, uint32_t* text_lengths,
			uint64_t num_pairs, int64_t longest_seq,
			char* cigar, uint64_t cigar_length,
			int32_t* dpu_distances, Timer timer, uint32_t n_reps){

	// Print timing results
	printf("[TIME] CPU-DPU pattern & texts (ms): ");
	print(&timer, 1, n_reps);
	printf("\n[TIME] DPU Kernel (ms): ");
	print(&timer, 2, n_reps);
	printf("\n[TIME] DISTANCES retrieve Time (ms): ");
	print(&timer, 3, n_reps);
	printf("\n[TIME] CIGAR retrieve Time (ms): ");
	print(&timer, 5, n_reps);
	printf("\n[TIME] FULL EXECUTION Time (distances + kernel + cpu) (ms): ");
	print(&timer, 4, n_reps);
	printf("\n[TIME] Added transfer Time (ms): ");
	print(&timer, 6, n_reps);
	printf("\n");
	printf("[INFO] Validating Results...\n");

	FILE * output;
	// Result validation
	int valid = 0;
	int failed = 0;
	int cigar_failed = 0;
	int distance_failed = 0;
	int over_distance = 0;
	output = fopen("dpu_output", "w");

	#pragma omp parallel
	{
	char* cigar_debug = (char*)malloc(cigar_length * sizeof(char));
	int cigar_d;
	int invalid;
	#pragma omp for reduction(+:failed) reduction(+:over_distance) reduction(+:cigar_failed) reduction(+:distance_failed)
	for (uint32_t j = 0; j < num_pairs; j++)
	{
		if(dpu_distances[j] == -1){
			over_distance++;

			#if DEBUG
			//printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] pair number %d has higher error than the maximum permited\n", j);
			#endif
			continue;
		}

		fprintf(output,"%d\t",dpu_distances[j]);

		// printf("Printing cigar %s", cigar[set][j*cigar_length[set]]);
		// printf("Printing distances %ld\n", dpu_distances[set][j]);
		cigar_debug = cigar_pprint(&cigar[j*cigar_length], output);
		cigar_d = cigar_distance(&cigar[j*cigar_length]);
		if(cigar_d != dpu_distances[j])
		{
			#if DEBUG
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED"Mismatch result at index %d: dpu_distance=%d cigar_distance=%d\n" ANSI_COLOR_RESET, j,dpu_distances[j], cigar_d);
			printf("Printing CIGAR: %s\n",cigar_debug);
			#endif
			distance_failed++;
		}
		invalid = check_cigar(pattern_lengths[j], text_lengths[j],
							&patterns[j*longest_seq],
							&texts[j*longest_seq],
							&cigar[j*cigar_length]);
		if(invalid){
			#if DEBUG
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED"CIGAR IS NOT VALID at index %d: dpu_distance=%d cigar_distance=%d\n"ANSI_COLOR_RESET, j,dpu_distances[j], cigar_d);
			printf("Printing INVALID CIGAR: %s\n",cigar_debug);
			#endif
			cigar_failed++;
		}else{
			valid++;
		}
		//if(j % 100000 == 0) printf("[INFO] %d pairs validated...\n", j);
	}
	}

	failed = cigar_failed;
	if(failed < distance_failed) failed = distance_failed;

	if (!failed) {
		printf("[" ANSI_COLOR_GREEN "-OK-" ANSI_COLOR_RESET "] Results are equal.\n");
	} else {
		printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] Results differ! %d/%ld failed %.2f%% failure...\n", failed, num_pairs, (float)((float)failed/(float)num_pairs)*100);
		//printf("[" ANSI_COLOR_YELLOW "Valid" ANSI_COLOR_RESET "] Valid pairs %d.\n", valid);
	}
	//if(over_distance>0) printf("[" ANSI_COLOR_YELLOW "ERR-" ANSI_COLOR_RESET "] %ld pairs have a higher error than the maximum stipulated %.2f%% pairs discarded.\n", over_distance,  (float)((float)over_distance/(float)num_pairs)*100);
}

int main(int argc, char **argv){
	struct Params p = input_params(argc, argv);

	printf(ANSI_COLOR_YELLOW "###############################################\n" ANSI_COLOR_RESET);
	printf(ANSI_COLOR_YELLOW "|------------------ ULSAPIM ------------------|\n" ANSI_COLOR_RESET);
	printf(ANSI_COLOR_YELLOW "###############################################\n" ANSI_COLOR_RESET);

  if(strcmp(p.input_file, "UNSET") == 0)
	{
    printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED" You need to specify an input filename!\n"ANSI_COLOR_RESET);
    printf("[INFO] Please run %s -h for help.\n", argv[0]);

    return 0;
	}
  if(p.num_sets == 0)
	{
    printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED" You need to specify at least one set \n"ANSI_COLOR_RESET);
    printf("[INFO] Please run %s -h for help.\n", argv[0]);

    return 0;
	}
	ewf_offset_t threshold = (ewf_offset_t)sqrt(4*WF_TRANSFER);
	printf("[INFO] Threshold %d\n",threshold);
	bool file_status;
	printf("[INFO] sets %d\n",p.sets[0]);
	if(NR_TASKLETS * (WF_TRANSFER*4 + SEQ_TRANSFER*2 + CIGAR_TRANSFER*2 + LEN_TRANSFER*2 + TASK_TRANSFER + 8 + STACK_SIZE_DEFAULT) > WRAM_LIMIT){
		printf("[ERROR] Exceeded WRAM capacity, reconfigure BL/BLI/NR_TASKLETS\n");
		return 0;
	}

	printf("[INFO] Number of sets %d\n", p.num_sets);
	printf("[INFO] " ANSI_COLOR_BLUE "Input file: %s" ANSI_COLOR_RESET "\n", p.input_file);
	printf("[INFO] Reading input file...\n");
	fflush(stdout);
	// Timer declaration
	Timer timer;
	timer_init(&timer, 1);
	timer_init(&timer, 2);
	timer_init(&timer, 3);
	timer_init(&timer, 4);
	timer_init(&timer, 5);
	timer_init(&timer, 6);
	Timer initialization;
	timer_init(&initialization, 1);
	//start(&initialization, 1, 1);
	struct dpu_set_t dpu_set;
	uint32_t nr_of_dpus;
	int set = 0;

	// Allocate and create input data
	// All patterns and text are allocated of the maximum size, then read only the required size
	char** patterns_set = (char**)malloc(p.num_sets * sizeof(char*));
	char** texts_set = (char**)malloc(p.num_sets * sizeof(char*));
	uint32_t** pattern_lengths_set = (uint32_t**)malloc(p.num_sets * sizeof(uint32_t*));
	uint32_t** text_lengths_set = (uint32_t**)malloc(p.num_sets * sizeof(uint32_t*));
	uint64_t num_pairs_set[p.num_sets];
	int64_t longest_seq_set[p.num_sets];
	int64_t longest_seq = 0;

	for(set=0; set<p.num_sets; set++){
		num_pairs_set[set]=0;
		if(p.sets != NULL){
			longest_seq_set[set] = p.sets[set];
			longest_seq_set[set] = (round(longest_seq_set[set] / 8) * 8)+8;
			longest_seq = (longest_seq_set[set] > longest_seq) ? longest_seq_set[set] : longest_seq;
		}
	}

	for (set = 0; set < p.num_sets; set++) {
        patterns_set[set] = (char*)malloc((uint64_t)1000 * (uint64_t) longest_seq_set[set] * (uint64_t) sizeof(char)); // Allocate memory for each dynamic array
        if (patterns_set[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		texts_set[set] = (char*)malloc((uint64_t)1000 * (uint64_t) longest_seq_set[set] * (uint64_t) sizeof(char)); // Allocate memory for each dynamic array
        if (texts_set[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }

		pattern_lengths_set[set] = (uint32_t*)malloc(1000 * sizeof(uint32_t)); // Allocate memory for each dynamic array
        if (pattern_lengths_set[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		text_lengths_set[set] = (uint32_t*)malloc(1000 * sizeof(uint32_t)); // Allocate memory for each dynamic array
        if (text_lengths_set[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		memset(pattern_lengths_set[set], 0, 1000 * sizeof(uint32_t));
		memset(text_lengths_set[set],    0, 1000 * sizeof(uint32_t));

	}
	//stop(&initialization, 1);
	start(&initialization, 1);
  	file_status = read_sequence_data(p.input_file, patterns_set, texts_set,
						pattern_lengths_set, text_lengths_set,
						p.num_sets, num_pairs_set, longest_seq_set, p.max_pairs);
	assert(file_status);
	stop(&initialization, 1);
	printf("[TIME] Reading input file (ms): ");
	print(&initialization, 1, 1);
	printf("\n");
	// printf("[VERBOSE] num pairs %ld and longest seq %ld\n", num_pairs[set], longest_seq[set]);

	start(&initialization, 1);
	float flt_num_pairs_per_dpu;
	uint64_t num_pairs_per_dpu;
	uint64_t num_pairs_alloc;
	float reminder;
	uint64_t cigar_length;
	uint64_t cigar_length_adjustment;
	uint64_t offsets_size_per_tl;
	uint64_t offsets_size_adjustement;
	uint64_t mem_offset;
	size_t mram_usage;
	uint64_t num_pairs = 0;
	uint64_t max_set_pairs = 0;
	int32_t* dpu_distances;
	char* cigar;
	char* patterns;
	char* texts;
	uint32_t* pattern_lengths;
	uint32_t* text_lengths;
	// Allocate DPUs
	for(set=0; set<p.num_sets; set++){
		num_pairs += num_pairs_set[set];
		max_set_pairs = (num_pairs_set[set] > max_set_pairs) ? num_pairs_set[set] : max_set_pairs;
		printf("[INFO] set %d has %ld pairs of size %ld\n", set, num_pairs_set[set], longest_seq_set[set]);
	}

	flt_num_pairs_per_dpu = (float) num_pairs / (float) NR_DPUS;
	num_pairs_per_dpu = (uint64_t) flt_num_pairs_per_dpu;
	reminder = fmod(flt_num_pairs_per_dpu, (float) NR_TASKLETS*2);
	if(reminder != 0) num_pairs_per_dpu = (uint64_t) ((flt_num_pairs_per_dpu + NR_TASKLETS*2) - reminder);
	num_pairs_alloc = num_pairs_per_dpu * NR_DPUS;
	//if(num_pairs_alloc[set] != num_pairs[set])printf("[WARNING] num pairs is not adjusted\n");
	// printf("[VERBOSE] num pairs alloc %ld vs num pairs %ld\n", num_pairs_alloc[set], num_pairs[set]);

	cigar_length = 2 * longest_seq;
	cigar_length_adjustment = cigar_length % CIGAR_TRANSFER;
	if(cigar_length_adjustment != 0) cigar_length += (CIGAR_TRANSFER - cigar_length_adjustment);
	// printf("[VERBOSE] CIGAR length %ld\n", cigar_length[set]);

	offsets_size_per_tl = 4*longest_seq* MAX_ERROR * sizeof(ewf_offset_t);
	//offsets_size_per_tl = 4*100 * sizeof(ewf_offset_t);
	offsets_size_adjustement = offsets_size_per_tl % WF_TRANSFER;
	if(offsets_size_adjustement != 0) offsets_size_per_tl += (WF_TRANSFER - offsets_size_adjustement);
	// printf("[VERBOSE] offset_size length %d\n", offsets_size_per_tl[set]);

	mram_usage = SEQ_TRANSFER // starting offset
	+ 2 * (longest_seq * num_pairs_per_dpu * sizeof(char)) // patterns and texts
	+ 2 * (num_pairs_per_dpu * sizeof(uint32_t)) // pattern and text lengths
	+ NR_TASKLETS * (4 * offsets_size_per_tl) // wavefronts
	+ (num_pairs_per_dpu * sizeof(int32_t)) // distances
	+ num_pairs_per_dpu * cigar_length * sizeof(char) // cigars
	+ NR_TASKLETS * (cigar_length + 8) // cigar tmp
	+ NR_TASKLETS * QUEUE_SZ; // task queue
	
	printf("[INFO] mram usage %.2f%%\n", ((float)mram_usage/(float)(MRAM_LIMIT))*100);
	if (mram_usage > MRAM_LIMIT){
		printf("[ERROR] set is way too large for the DPUS to process, increase DPUs or reduce tasklets \n");
		return 0;
	}

	// Allocating result arrays
	dpu_distances     = (int32_t*)malloc((size_t)num_pairs_alloc * (size_t)sizeof(int32_t));
	cigar        		= (char*)malloc((size_t)num_pairs_alloc *((size_t)cigar_length * (size_t)sizeof(char)));
	pattern_lengths = (uint32_t *) malloc((size_t)(num_pairs_alloc) * (size_t)sizeof(uint32_t));
	text_lengths = (uint32_t *) malloc((size_t)(num_pairs_alloc) * (size_t)sizeof(uint32_t));
	patterns = (char *) malloc((size_t)(num_pairs_alloc) * (size_t)sizeof(char) * (size_t)longest_seq);
	texts = (char *) malloc((size_t)(num_pairs_alloc) * (size_t)sizeof(char) * (size_t)longest_seq);
	uint32_t pair = 0;

	for(uint32_t set_pair = 0; set_pair<max_set_pairs; set_pair++){
		for(set=0; set<p.num_sets; set++){
			if(set_pair<num_pairs_set[set]){
				strcpy(&patterns[pair*longest_seq], &patterns_set[set][set_pair*longest_seq_set[set]]);
				pattern_lengths[pair] = pattern_lengths_set[set][set_pair];
				strcpy(&texts[pair*longest_seq], &texts_set[set][set_pair*longest_seq_set[set]]);
				text_lengths[pair] = text_lengths_set[set][set_pair];
				pair++;
			}
		}
	}

	memset(&pattern_lengths[pair], 0, (num_pairs_alloc - num_pairs) * sizeof(uint32_t));
	memset(&text_lengths[pair],    0, (num_pairs_alloc - num_pairs)  * sizeof(uint32_t));

	free(patterns_set);
	free(pattern_lengths_set);
	free(texts_set);
	free(text_lengths_set);
	
	dpu_results_t results[NR_DPUS];
	dpu_results_t results_min[NR_DPUS];
	uint64_t perf_main[NR_DPUS * NR_TASKLETS];
	uint64_t perf_sum = 0;
	int perf_valid = 0;

	#if PERF
	for(uint32_t i = 0; i<NR_DPUS; i++){
		for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
			perf_main[i * each_tasklet] = 0;
		}
		results[i].main = 0;
		results_min[i].main = 0xFFFFFFFFFFFFFFFF;
		results[i].overlap = 0;
		results[i].base_case = 0;
		results[i].compute = 0;
		results[i].extend = 0;

	}
	#endif

	uint32_t n_batches;
	uint32_t batch_pairs_per_dpu;
	uint64_t batch_pairs;
	uint32_t spare_dpu_pairs = 0;
	uint64_t over_distance = 0;
	bool terminate = 0;
	if(num_pairs_alloc <= BATCH_SIZE || BATCH_SIZE == 0){
		n_batches = 1;
		batch_pairs_per_dpu = num_pairs_alloc / NR_DPUS;
		batch_pairs = num_pairs_alloc;
	}else{
		flt_num_pairs_per_dpu = (float) BATCH_SIZE / (float) NR_DPUS;
		batch_pairs_per_dpu = (uint64_t) flt_num_pairs_per_dpu;
		reminder = fmod(flt_num_pairs_per_dpu, (float) NR_TASKLETS*2);
		if(reminder != 0) batch_pairs_per_dpu = (uint64_t) ((flt_num_pairs_per_dpu + NR_TASKLETS*2) - reminder);
		batch_pairs =  batch_pairs_per_dpu * NR_DPUS; 
		n_batches = (uint32_t) ceil((float)num_pairs_per_dpu / (float) batch_pairs_per_dpu);
		
	}
	printf("[INFO] batch pairs per dpu %d total batch pairs %ld\n", batch_pairs_per_dpu, batch_pairs);
	uint32_t* cpu_pairs_idx = (uint32_t *) malloc((size_t)(batch_pairs) * (size_t)sizeof(uint32_t));
	omp_set_num_threads(17);

	for (uint32_t rep = 0; rep < p.n_warmup + p.n_reps; rep++) { // Loop adding compute weight and warmup
		start(&timer, 6);
		send_data(patterns, texts,
					pattern_lengths, text_lengths,
					num_pairs, num_pairs_per_dpu,
					longest_seq, &dpu_set, &nr_of_dpus,
					rep, p.n_warmup, &timer, &mem_offset);
		
		start(&timer, 4);
		#pragma omp parallel shared(over_distance, cpu_pairs_idx)
		{
		for(uint32_t batch_idx = 0; batch_idx < n_batches; batch_idx++){
			#pragma omp master
			{
			run_set(num_pairs_per_dpu, batch_pairs_per_dpu,
					longest_seq, &dpu_set, threshold, 
					rep, p.n_warmup, &timer,batch_idx);
			
			if(batch_pairs_per_dpu * (batch_idx+1) > num_pairs_per_dpu){
				spare_dpu_pairs = ((batch_pairs_per_dpu * (batch_idx+1)) - (num_pairs_per_dpu));
				//printf("[DEBUG] BATCH PAIRS VS DPU PAIRS %d/%ld | spare dpu pairs %d\n", batch_pairs_per_dpu * (batch_idx+1), num_pairs_per_dpu, spare_dpu_pairs);
			}

			printf("\n[INFO] DPU execution launched batch %d/%d\n", batch_idx, n_batches-1);
			fflush(stdout);
			bool set_finished = false;
			bool set_failed;
			#if PRINT
			struct dpu_set_t dpu_fail;
			#endif

			while(!set_finished){
				dpu_status(dpu_set, &set_finished, &set_failed);
				if(set_failed){
					printf("[ERROR] DPU EXECUTION FAILED\n");
					#if PRINT
					DPU_FOREACH(dpu_set, dpu_fail) {
						DPU_ASSERT(dpu_log_read(dpu_fail, stdout));
					}
					#endif
					DPU_ASSERT(dpu_free(dpu_set));
					#pragma omp atomic write
					terminate = 1;
					break;
				}
				if (set_finished)
				{
					if (rep >= p.n_warmup)
					{
						stop(&timer, 2);
					#if ENERGY
						DPU_ASSERT(dpu_probe_stop(&probe));
					#endif
					}
					printf("[INFO] DPU execution finished succesfully\n");
					#if PRINT
					DPU_FOREACH(dpu_set, dpu_fail) {
						DPU_ASSERT(dpu_log_read(dpu_fail, stdout));
					}
					#endif
				}
			}
			}
			//printf("\n[DEBUG] thread %d ready and waiting for batch %d\n", omp_get_thread_num(), batch_idx);
			#pragma omp barrier

			// #pragma omp single
			// {
			// while(over_distance<100){
			// 	over_distance++;
			// }
			// }

			// printf("\n[DEBUG] OVER DISTANCE VALUE %ld\n", over_distance);

			if(terminate) break;
			retrieve_and_recover(dpu_set, batch_pairs_per_dpu, num_pairs_per_dpu, 
					&mem_offset, patterns, texts,
					pattern_lengths, text_lengths,
					spare_dpu_pairs, longest_seq,
					dpu_distances, &timer, rep, p.n_warmup, batch_idx, 
					num_pairs, cpu_pairs_idx, results, results_min, perf_main, 
					&perf_sum, &perf_valid, &over_distance);
		}
		}
		stop(&timer, 4);
		//printf("[DEBUG] reached end of parellel reagion...\n");
		if(terminate) return 0;
		retrieve_cigars(dpu_set, num_pairs_per_dpu, offsets_size_per_tl,
						mem_offset, cigar, cigar_length, 
						&timer, rep, p.n_warmup);
		stop(&timer, 6);

		#if PERF_MAIN
		uint64_t max_cycles = 0;
		uint64_t min_cycles = 0xFFFFFFFFFFFFFFFF;
		double perf_mean = (double)perf_sum/perf_valid;
		double sum_squared_diff = 0;
		double perf_diff;
		// Print performance results
		if(rep >= p.n_warmup){
			for(uint32_t i = 0; i<NR_DPUS; i++){
				if(results[i].main > max_cycles)
					max_cycles = results[i].main;
				if(results_min[i].main < min_cycles && results_min[i].main > 0)
					min_cycles = results_min[i].main;
				for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
					if(perf_main[i * each_tasklet] > 0){
						perf_diff = ((double)perf_main[i * each_tasklet]) - perf_mean;
						sum_squared_diff += perf_diff * perf_diff;
					}
				}
				i++;
			}
		}
		printf("\n[PERF] DPU kernel perf main MAX: %g\n", (double)max_cycles);
		printf("\n[PERF] DPU kernel standard deviation: %g\n", sqrt(sum_squared_diff / perf_valid));
		printf("\n[PERF] DPU kernel mean: %g\n", perf_mean);
		printf("\n[PERF] DPU kernel CV: %.2f\n", sqrt(sum_squared_diff / perf_valid)/perf_mean);
		printf("\n[PERF] DPU kernel perf main MIN: %g\n", (double)min_cycles);
		#endif
		#if PERF_SECTIONS
		uint64_t max_cycles = 0;
		uint64_t min_cycles = 0xFFFFFFFFFFFFFFFF;
		max_cycles = 0;
		min_cycles = 0xFFFFFFFFFFFFFFFF;
		// Print performance results
		if(rep >= p.n_warmup){
			for(uint32_t i = 0; i<NR_DPUS; i++){
				if(results[i].overlap > max_cycles)
					max_cycles = results[i].overlap;
				if(results[i].overlap < min_cycles)
					min_cycles = results[i].overlap;
				i++;
			}
		}
		printf("[PERF] DPU kernel perf overlap MAX: %g\n", (double)max_cycles);
		printf("[PERF] DPU kernel perf overlap MIN: %g\n", (double)min_cycles);
		max_cycles = 0;
		min_cycles = 0xFFFFFFFFFFFFFFFF;
		// Print performance results
		if(rep >= p.n_warmup){
			DPU_FOREACH(dpu_set, dpu) {
				if(results[i].base_case > max_cycles)
					max_cycles = results[i].base_case;
				if(results[i].base_case < min_cycles)
					min_cycles = results[i].base_case;
				i++;
			}
		}
		printf("[PERF] DPU kernel perf base_case MAX: %g\n", (double)max_cycles);
		printf("[PERF] DPU kernel perf base_case MIN: %g\n", (double)min_cycles);
		max_cycles = 0;
		min_cycles = 0xFFFFFFFFFFFFFFFF;
		// Print performance results
		if(rep >= p.n_warmup){
			for(uint32_t i = 0; i<NR_DPUS; i++){
				if(results[i].compute > max_cycles)
					max_cycles = results[i].compute;
				if(results[i].compute < min_cycles)
					min_cycles = results[i].compute;
				i++;
			}
		}
		printf("[PERF] DPU kernel perf compute MAX: %g\n", (double)max_cycles);
		printf("[PERF] DPU kernel perf compute MIN: %g\n", (double)min_cycles);
		max_cycles = 0;
		min_cycles = 0xFFFFFFFFFFFFFFFF;
		// Print performance results
		if(rep >= p.n_warmup){
			for(uint32_t i = 0; i<NR_DPUS; i++){
				if(results[i].extend > max_cycles)
					max_cycles = results[i].extend;
				if(results[i].extend < min_cycles)
					min_cycles = results[i].extend;
				i++;
			}
		}
		printf("[PERF] DPU kernel perf extend MAX: %g\n", (double)max_cycles);
		printf("[PERF] DPU kernel perf extend MIN: %g\n", (double)min_cycles);
		#endif
	}


	validate_results(patterns, texts,
			pattern_lengths, text_lengths,
			num_pairs, longest_seq,
			cigar, cigar_length,
			dpu_distances, timer, p.n_reps);

	free(patterns);
	free(texts);
	free(pattern_lengths);
	free(text_lengths);
	free(cigar);
	free(dpu_distances);
	free(cpu_pairs_idx);

	return 0;
}
