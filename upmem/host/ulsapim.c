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
#include <math.h>

#include <dpu.h>
#include <dpu_log.h>

#if ENERGY
#include <dpu_probe.h>
#endif

#include "params.h"
#include "timer.h"
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

		for(int i=0; i<num_sets; i++){
			if(pattern_length < longest_seq[i] && text_length < longest_seq[i]){
				set = i;
				break;
			}
		}
		if (set == -1){
			// printf("A pair has been discarded, with length %d pattern %d text\n", pattern_length, text_length);
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

void run_set(char* patterns, char* texts,
			uint32_t* pattern_lengths, uint32_t* text_lengths,
			uint64_t num_pairs, uint64_t num_pairs_per_dpu,
			int64_t longest_seq, struct dpu_set_t* dpu_set,
			int32_t nr_dpus_set, uint32_t* nr_of_dpus, ewf_offset_t threshold, uint32_t rep,
			uint32_t n_warmup, Timer* timer, int set, int num_sets, uint64_t* mem_offset){


		printf("[INFO] Transferring inputs from to DPUs... ");
		fflush(stdout);

		if (rep >= n_warmup)
			start(timer, 1, rep - n_warmup);

		printf("\n[INFO] " ANSI_COLOR_BLUE "set number: \t %d" ANSI_COLOR_RESET "\n", set);
		printf("[INFO] " ANSI_COLOR_BLUE "num_dpus required: \t %d" ANSI_COLOR_RESET "\n", nr_dpus_set);
		printf("[INFO] " ANSI_COLOR_BLUE "num_pairs: \t %ld" ANSI_COLOR_RESET "\n", num_pairs);
		printf("[INFO] " ANSI_COLOR_BLUE "longest_sequence: %ld" ANSI_COLOR_RESET "\n", longest_seq);
		printf("[INFO] Pairs per DPU (adjusted): (%ld)\n", num_pairs_per_dpu);
		printf("[INFO] Total WRAM usage: (%f) bytes \n", ((float)NR_TASKLETS * (WF_TRANSFER*4 + SEQ_TRANSFER*2 + CIGAR_TRANSFER*2 + LEN_TRANSFER*2 + TASK_TRANSFER + 8 + STACK_SIZE_DEFAULT)/(float)WRAM_LIMIT)*100);
		printf("[INFO] Reading input... "ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
		printf("[INFO] Allocating DPUs... ");
		fflush(stdout);

		DPU_ASSERT(dpu_alloc(nr_dpus_set, NULL, dpu_set));
		DPU_ASSERT(dpu_load(*dpu_set, DPU_BINARY, NULL));
		DPU_ASSERT(dpu_get_nr_dpus(*dpu_set, nr_of_dpus));

		printf("\n[INFO] " ANSI_COLOR_BLUE "num_dpus allocated: %d" ANSI_COLOR_RESET "\n", *nr_of_dpus);

		unsigned int kernel = 0;
		dpu_arguments_t input_arguments = {longest_seq, num_pairs_per_dpu, threshold, kernel};
		uint32_t i = 0;
		struct dpu_set_t dpu;
		printf("[INFO] Transferring inputs to DPUs... ");

		// Copy input arguments to dpu
		DPU_FOREACH(*dpu_set, dpu) {
			DPU_ASSERT(dpu_copy_to(dpu, "DPU_INPUT_ARGUMENTS", 0, (const void *) &input_arguments, sizeof(input_arguments)));
			i++;
		}

		i = 0;
		*mem_offset = SEQ_TRANSFER; // offset points to first MRAM address
		// printf("[VERBOSE] HOST ma pattern %d\n", mem_offset[set]);
		// Patterns transfer
		DPU_FOREACH(*dpu_set, dpu, i)
		{
			DPU_ASSERT(dpu_prepare_xfer(dpu, patterns + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		// printf("[VERBOSE] HOST pattern sent %ld bytes\n", longest_seq[set] * num_pairs_per_dpu[set] * sizeof(char));

    	i = 0;
		*mem_offset += (longest_seq * num_pairs_per_dpu * sizeof(char)); // of
		// printf("[VERBOSE] HOST ma texts %d\n", mem_offset[set]);
		// Texts transfer
		DPU_FOREACH(*dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, texts + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(*dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, *mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		// printf("[VERBOSE] HOST text sent %ld bytes\n", longest_seq[set] * num_pairs_per_dpu[set] * sizeof(char));

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
		// Run kernel on DPUs

		if (rep >= n_warmup)
		{
			start(timer, 2, rep - n_warmup);
		#if ENERGY
			DPU_ASSERT(dpu_probe_start(&probe));
		#endif
		}
		if (NR_DPUS >= num_sets){
			DPU_ASSERT(dpu_launch(*dpu_set, DPU_ASYNCHRONOUS)); // DPU can be used asyncronously (have into account for future improvements)
		}else{
			DPU_ASSERT(dpu_launch(*dpu_set, DPU_SYNCHRONOUS));
		}
}

void retrieve_results(struct dpu_set_t dpu_set, int32_t* dpu_distances,
					uint64_t num_pairs_per_dpu, uint64_t offsets_size_per_tl,
					uint64_t mem_offset, char* cigar, uint64_t cigar_length,
					Timer* timer, uint32_t rep, uint32_t n_warmup, uint32_t nr_of_dpus){

	printf("[INFO] Retrieving results to host... ");
	fflush(stdout);
	struct dpu_set_t dpu;
	if (rep >= n_warmup)
		start(timer, 3, rep - n_warmup);
	
	#if PERF
	double cc;
    double cc_min;
	cc = 0;
    cc_min = 0;
	dpu_results_t results[nr_of_dpus];
	#endif
	int i = 0;
	// Copy back the DPU results (mem_offset[set] already points to the results starting memory position)
	// printf("\n[VERBOSE] HOST ma distances %d\n", mem_offset[set]);
	DPU_FOREACH(dpu_set, dpu, i) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_distances + i * num_pairs_per_dpu));

		#if PERF
		dpu_results_t result;
		results[i].main = 0;
		results[i].overlap = 0;
		results[i].base_case = 0;
		results[i].compute = 0;
		results[i].extend = 0;
		// Retrieve tasklet timings
		for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
			// results[i].main = 0;
			// results[i].overlap = 0;
			// results[i].base_case = 0;
			// results[i].compute = 0;
			// results[i].extend = 0;
			DPU_ASSERT(dpu_copy_from(dpu, "DPU_RESULTS", each_tasklet * sizeof(dpu_results_t), &result, sizeof(dpu_results_t)));
			if (result.main > results->main)
			 	results[i].main = result.main;
			if (result.overlap > results[i].overlap)
				results[i].overlap = result.overlap;
			if (result.base_case > results[i].base_case)
				results[i].base_case = result.base_case;
			if (result.compute > results->compute)
				results[i].compute = result.compute;
			if (result.extend > results[i].extend)
				results[i].extend = result.extend;
		}
		#endif
	}

	DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
				mem_offset, num_pairs_per_dpu * sizeof(int32_t), DPU_XFER_DEFAULT));
	// printf("[VERBOSE] HOST distances recieved %ld bytes\n", num_pairs_per_dpu[set] * sizeof(int32_t));

	#if PERF
	uint64_t max_cycles = 0;
	uint64_t min_cycles = 0xFFFFFFFFFFFFFFFF;
	// Print performance results
	if(rep >= n_warmup){
		i = 0;
		DPU_FOREACH(dpu_set, dpu) {
			if(results[i].main > max_cycles)
				max_cycles = results[i].main;
			if(results[i].main < min_cycles)
				min_cycles = results[i].main;
			i++;
		}
		cc = (double)max_cycles;
		cc_min = (double)min_cycles;
	}
	printf("\n[PERF] DPU kernel perf main: %g\n", cc);
	max_cycles = 0;
	min_cycles = 0xFFFFFFFFFFFFFFFF;
	cc = 0;
	cc_min = 0;
	// Print performance results
	if(rep >= n_warmup){
		i = 0;
		DPU_FOREACH(dpu_set, dpu) {
			if(results[i].overlap > max_cycles)
				max_cycles = results[i].overlap;
			if(results[i].overlap < min_cycles)
				min_cycles = results[i].overlap;
			i++;
		}
		cc = (double)max_cycles;
		cc_min = (double)min_cycles;
	}
	printf("[PERF] DPU kernel perf overlap: %g\n", cc);
	max_cycles = 0;
	min_cycles = 0xFFFFFFFFFFFFFFFF;
	cc = 0;
	cc_min = 0;
	// Print performance results
	if(rep >= n_warmup){
		i = 0;
		DPU_FOREACH(dpu_set, dpu) {
			if(results[i].base_case > max_cycles)
				max_cycles = results[i].base_case;
			if(results[i].base_case < min_cycles)
				min_cycles = results[i].base_case;
			i++;
		}
		cc = (double)max_cycles;
		cc_min = (double)min_cycles;
	}
	printf("[PERF] DPU kernel perf base_case: %g\n", cc);
	max_cycles = 0;
	min_cycles = 0xFFFFFFFFFFFFFFFF;
	cc = 0;
	cc_min = 0;
	// Print performance results
	if(rep >= n_warmup){
		i = 0;
		DPU_FOREACH(dpu_set, dpu) {
			if(results[i].compute > max_cycles)
				max_cycles = results[i].compute;
			if(results[i].compute < min_cycles)
				min_cycles = results[i].compute;
			i++;
		}
		cc = (double)max_cycles;
		cc_min = (double)min_cycles;
	}
	printf("[PERF] DPU kernel perf compute: %g\n", cc);
	max_cycles = 0;
	min_cycles = 0xFFFFFFFFFFFFFFFF;
	cc = 0;
	cc_min = 0;
	// Print performance results
	if(rep >= n_warmup){
		i = 0;
		DPU_FOREACH(dpu_set, dpu) {
			if(results[i].extend > max_cycles)
				max_cycles = results[i].extend;
			if(results[i].extend < min_cycles)
				min_cycles = results[i].extend;
			i++;
		}
		cc = (double)max_cycles;
		cc_min = (double)min_cycles;
	}
	printf("[PERF] DPU kernel perf extend: %g\n", cc);
	#endif

	// move mem offset to the end of the 4 wavefronts
	// Offsets per tl size calculatio
	mem_offset += NR_TASKLETS * (4 * offsets_size_per_tl) + (num_pairs_per_dpu * sizeof(int32_t));
	i = 0;
	// printf("\n[VERBOSE] HOST ma CIGAR %d\n", mem_offset[set]);
	DPU_FOREACH(dpu_set, dpu, i) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, cigar + i * num_pairs_per_dpu * cigar_length));
	}
	DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
				mem_offset, num_pairs_per_dpu * cigar_length * sizeof(char), DPU_XFER_DEFAULT));
	// printf("[VERBOSE] HOST CIGAR recieved %ld bytes\n", num_pairs_per_dpu[set] * cigar_length[set] * sizeof(char));

	if(rep >= n_warmup)
		stop(timer, 3);

	printf(""ANSI_COLOR_GREEN "OK" ANSI_COLOR_RESET".\n");
	#if PRINT
	printf("LOGS\n");
	DPU_FOREACH(dpu_set, dpu) {
		DPU_ASSERT(dpu_log_read(dpu, stdout));
	}
	#endif
	DPU_ASSERT(dpu_free(dpu_set));
}

void validate_results(char* patterns, char* texts,
			uint32_t* pattern_lengths, uint32_t* text_lengths,
			uint64_t num_pairs, int64_t longest_seq, struct dpu_set_t dpu_set,
			char* cigar, uint64_t cigar_length,
			int32_t* dpu_distances, Timer timer, int set, uint32_t n_reps){

	// Print timing results
	printf("\n[TIME] CPU-DPU (ms): ");
	print(&timer, 1, n_reps);
	printf("\n[TIME] DPU Kernel (ms): ");
	print(&timer, 2, n_reps);
	printf("\n[TIME] DPU-CPU Time (ms): ");
	print(&timer, 3, n_reps);
	printf("\n");
	printf("[INFO] Validating Results for set %d ...\n", set);

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
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] pair number %d has higher error than the maximum permited\n", j);
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
			printf("Printing wrong CIGAR: %s\n",cigar_debug);
			#endif
			distance_failed++;
		}
		invalid = check_cigar(pattern_lengths[j], text_lengths[j],
							&patterns[j*longest_seq],
							&texts[j*longest_seq],
							&cigar[j*cigar_length]);
		if(invalid && cigar_d == dpu_distances[j]){
			#if DEBUG
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED"CIGAR IS NOT VALID at index %d: dpu_distance=%d cigar_distance=%d\n"ANSI_COLOR_RESET, j,dpu_distances[j], cigar_d);
			printf("Printing wrowng CIGAR: %s\n",cigar_debug);
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
		printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] Results differ! %d/%ld failed...\n", failed, num_pairs);
		//printf("[" ANSI_COLOR_YELLOW "Valid" ANSI_COLOR_RESET "] Valid pairs %d.\n", valid);
	}
	if(over_distance>0) printf("[" ANSI_COLOR_YELLOW "ERR-" ANSI_COLOR_RESET "] %d pairs have a higher error than the maximum stipulated.\n", over_distance);
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
	ewf_offset_t threshold = (ewf_offset_t)sqrt(3*WF_TRANSFER);
	printf("[INFO] Threshold %d\n",threshold);
	bool file_status;
	printf("[INFO] sets %d\n",p.sets[0]);
	if(NR_TASKLETS * (WF_TRANSFER*4 + SEQ_TRANSFER*2 + CIGAR_TRANSFER*2 + LEN_TRANSFER*2 + TASK_TRANSFER + 8 + STACK_SIZE_DEFAULT) > WRAM_LIMIT){
		printf("[ERROR] Exceeded WRAM capacity, reconfigure BL/BLI/NR_TASKLETS");
		return 0;
	}

	printf("[INFO] Number of sets %d\n", p.num_sets);
	printf("[INFO] " ANSI_COLOR_BLUE "Input file: %s" ANSI_COLOR_RESET "\n", p.input_file);
	printf("[INFO] Reading input file...\n");
	fflush(stdout);
	// Timer declaration
	Timer timer[p.num_sets];
	Timer initialization;
	//start(&initialization, 1, 1);
	struct dpu_set_t* dpu_set = (struct dpu_set_t*)malloc(p.num_sets * sizeof(struct dpu_set_t));
	uint32_t nr_of_dpus[p.num_sets];
	int set=0;

	// Allocate and create input data
	// All patterns and text are allocated of the maximum size, then read only the required size
	char** patterns = (char**)malloc(p.num_sets * sizeof(char*));
	char** texts = (char**)malloc(p.num_sets * sizeof(char*));
	uint32_t** pattern_lengths = (uint32_t**)malloc(p.num_sets * sizeof(uint32_t*));
	uint32_t** text_lengths = (uint32_t**)malloc(p.num_sets * sizeof(uint32_t*));
	uint64_t num_pairs[p.num_sets];
	int64_t longest_seq[p.num_sets];

	for(set=0; set<p.num_sets; set++){
		num_pairs[set]=0;
		if(p.sets != NULL){
			longest_seq[set] = p.sets[set];
			longest_seq[set] = (round(longest_seq[set] / 8) * 8)+8;
		}
	}

	for (set = 0; set < p.num_sets; set++) {
        patterns[set] = (char*)malloc((uint64_t)1000 * (uint64_t) longest_seq[set] * (uint64_t) sizeof(char)); // Allocate memory for each dynamic array
        if (patterns[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		texts[set] = (char*)malloc((uint64_t)1000 * (uint64_t) longest_seq[set] * (uint64_t) sizeof(char)); // Allocate memory for each dynamic array
        if (texts[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }

		pattern_lengths[set] = (uint32_t*)malloc(1000 * sizeof(uint32_t)); // Allocate memory for each dynamic array
        if (pattern_lengths[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		text_lengths[set] = (uint32_t*)malloc(1000 * sizeof(uint32_t)); // Allocate memory for each dynamic array
        if (text_lengths[set] == NULL) {
            printf("[ERROR] Memory allocation failed. Exiting...\n");
            // You may need to free memory for previously allocated arrays before exiting.
            return 1;
        }
		memset(pattern_lengths[set], 0, 1000 * sizeof(uint32_t));
		memset(text_lengths[set],    0, 1000 * sizeof(uint32_t));

	}
	//stop(&initialization, 1);
	start(&initialization, 1, 1);
  	file_status = read_sequence_data(p.input_file, patterns, texts,
						pattern_lengths, text_lengths,
						p.num_sets, num_pairs, longest_seq, p.max_pairs);
	assert(file_status);
	stop(&initialization, 1);
	printf("[TIME] Reading input file (ms): ");
	print(&initialization, 1, 1);
	printf("\n");
	// printf("[VERBOSE] num pairs %ld and longest seq %ld\n", num_pairs[set], longest_seq[set]);

	start(&initialization, 1, 1);
	int32_t unused_dpus = 0;
	int32_t current_dpus = NR_DPUS;
	int64_t total_bytes_size = 0;
	int32_t nr_dpus_set[p.num_sets];
	float flt_num_pairs_per_dpu;
	uint64_t num_pairs_per_dpu[p.num_sets];
	uint64_t num_pairs_alloc[p.num_sets];
	float reminder;
	uint64_t cigar_length[p.num_sets];
	uint64_t cigar_length_adjustment;
	uint64_t offsets_size_per_tl[p.num_sets];
	uint64_t offsets_size_adjustement;
	uint64_t mem_offset[p.num_sets];
	int32_t** 	dpu_distances        = (int32_t**)malloc(p.num_sets * sizeof(int32_t*));
	char** 		cigar        		 = (char**)malloc(p.num_sets * sizeof(char*));
	//char** 		cigar_debug        	 = (char**)malloc(p.num_sets * sizeof(char*));
	size_t mram_usage[p.num_sets];
	uint32_t queued_sets = 0;
	uint32_t set_partitions[p.num_sets];
	bool finished_set[p.num_sets];

	// Allocate DPUs
	for(set=0; set<p.num_sets; set++){
		if(num_pairs[set]==0) continue;
		total_bytes_size += num_pairs[set] * longest_seq[set];
	}
	for(set=0; set<p.num_sets; set++){
		if(num_pairs[set]==0) continue;
		if (NR_DPUS >= p.num_sets)
		{
			nr_dpus_set[set] = current_dpus / ((float)total_bytes_size/(float)(num_pairs[set] * longest_seq[set]));
			printf("[INFO] set %d has %d dpus asigned total bytes %ld set bytes %ld division %f \n", set, nr_dpus_set[set],total_bytes_size, (num_pairs[set] * longest_seq[set]), ((float)total_bytes_size/(float)(num_pairs[set] * longest_seq[set])));
			if(nr_dpus_set[set]==0){
				nr_dpus_set[set] = 1;
				current_dpus--;
			}
		}else{
			nr_dpus_set[set] = NR_DPUS;
		}
	}
	// for(set=p.num_sets; set>0; set--){
	// 	if(nr_dpus_set[set] > num_pairs[set]){
	// 		unused_dpus += nr_dpus_set[set] - num_pairs[set];
	// 	}else{
	// 		if(unused_dpus - (num_pairs[set] - nr_dpus_set[set]) > 0){
	// 			nr_dpus_set[set] += num_pairs[set] - unused_dpus;
	// 			unused_dpus -= (num_pairs[set] - nr_dpus_set[set]);
	// 			if(unused_dpus<0) unused_dpus = 0;
	// 		}else{
	// 			nr_dpus_set[set] += unused_dpus;
	// 			unused_dpus = 0;
	// 		}
	// 	}
	// }


	for(set=0; set<p.num_sets; set++){
		if(num_pairs[set]==0) continue;
		set_partitions[set] = 1;
		finished_set[set] = false;
		flt_num_pairs_per_dpu = (float) num_pairs[set] / (float) nr_dpus_set[set];
		num_pairs_per_dpu[set] = (uint64_t) flt_num_pairs_per_dpu;
		reminder = fmod(flt_num_pairs_per_dpu, (float) NR_TASKLETS*2);
		if(reminder != 0) num_pairs_per_dpu[set] = (uint64_t) ((flt_num_pairs_per_dpu + NR_TASKLETS*2) - reminder);
		num_pairs_alloc[set] = num_pairs_per_dpu[set] * nr_dpus_set[set];
		//if(num_pairs_alloc[set] != num_pairs[set])printf("[WARNING] num pairs is not adjusted\n");
		// printf("[VERBOSE] num pairs alloc %ld vs num pairs %ld\n", num_pairs_alloc[set], num_pairs[set]);

		cigar_length[set] = 2 * longest_seq[set];
		cigar_length_adjustment = cigar_length[set] % CIGAR_TRANSFER;
		if(cigar_length_adjustment != 0) cigar_length[set] += (CIGAR_TRANSFER - cigar_length_adjustment);
		// printf("[VERBOSE] CIGAR length %ld\n", cigar_length[set]);

		offsets_size_per_tl[set] = 4*longest_seq[set]* MAX_ERROR * sizeof(ewf_offset_t);
		offsets_size_adjustement = offsets_size_per_tl[set] % WF_TRANSFER;
		if(offsets_size_adjustement != 0) offsets_size_per_tl[set] += (WF_TRANSFER - offsets_size_adjustement);
		// printf("[VERBOSE] offset_size length %d\n", offsets_size_per_tl[set]);

		mram_usage[set] = SEQ_TRANSFER // starting offset
		+ 2 * (longest_seq[set] * num_pairs_per_dpu[set] * sizeof(char)) // patterns and texts
		+ 2 * (num_pairs_per_dpu[set] * sizeof(uint32_t)) // pattern and text lengths
		+ NR_TASKLETS * (4 * offsets_size_per_tl[set]) // wavefronts
		+ (num_pairs_per_dpu[set] * sizeof(int32_t)) // distances
		+ num_pairs_per_dpu[set] * cigar_length[set] * sizeof(char) // cigars
		+ NR_TASKLETS * (cigar_length[set] + 8) // cigar tmp
		+ NR_TASKLETS*QUEUE_SZ; // task queue
		
		printf("[INFO] mram usage vs MRAM_LIMIT %ld/%d\n", mram_usage[set], MRAM_LIMIT);

		while (mram_usage[set] > MRAM_LIMIT)
		{
			set_partitions[set]++;
			mram_usage[set] = SEQ_TRANSFER // starting offset
			+ 2 * (longest_seq[set] * (num_pairs_per_dpu[set] / set_partitions[set])  * sizeof(char)) // patterns and texts
			+ 2 * ((num_pairs_per_dpu[set] / set_partitions[set]) * sizeof(uint32_t)) // pattern and text lengths
			+ NR_TASKLETS * (4 * offsets_size_per_tl[set]) // wavefronts
			+ ((num_pairs_per_dpu[set] / set_partitions[set]) * sizeof(int32_t)) // distances
			+ (num_pairs_per_dpu[set] / set_partitions[set]) * cigar_length[set] * sizeof(char) // cigars
			+ NR_TASKLETS * (cigar_length[set] + 8) // cigar tmp
			+ NR_TASKLETS*QUEUE_SZ; // task queue
			queued_sets++;
			if (set_partitions[set] > 10){
				printf("[ERROR] set is way too large for the DPUS to process, consider creating different sets or reducing file size \n");
				return 0;
			}
		}

		if(set_partitions[set] > 1){
			printf("[INFO] divided set %d into %d partitions \n", set, set_partitions[set]);
		}

		// Allocating result arrays
		dpu_distances[set]      = (int32_t*)malloc((size_t)num_pairs_alloc[set] * (size_t)sizeof(int32_t));
		cigar[set]        		= (char*)malloc((size_t)num_pairs_alloc[set] *((size_t)cigar_length[set] * (size_t)sizeof(char)));
		//cigar_debug[set]        = (char*)malloc(cigar_length[set] * sizeof(char));
		pattern_lengths[set] = (uint32_t *) realloc(pattern_lengths[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(uint32_t));
		text_lengths[set] = (uint32_t *) realloc(text_lengths[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(uint32_t));
		patterns[set] = (char *) realloc(patterns[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(char) * (size_t)longest_seq[set]);
		texts[set] = (char *) realloc(texts[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(char) * (size_t)longest_seq[set]);
		memset(&pattern_lengths[set][num_pairs[set]], 0, (num_pairs_alloc[set] - num_pairs[set]) * sizeof(uint32_t));
		memset(&text_lengths[set][num_pairs[set]],    0, (num_pairs_alloc[set] - num_pairs[set])  * sizeof(uint32_t));
		// printf("[VERBOSE] Distance allocated %ld bytes\n", (uint64_t)num_pairs_alloc[set] * (uint64_t)sizeof(int32_t));
		// printf("[VERBOSE] CIGAR allocated %ld bytes\n", (uint64_t)num_pairs_alloc[set] *((uint64_t)cigar_length[set] * (uint64_t)sizeof(char)));
		// printf("[VERBOSE] CIGAR debug allocated %ld bytes\n", cigar_length[set] * sizeof(char));
		// printf("[VERBOSE] Pattern lengths allocated %ld bytes\n", (uint64_t)(num_pairs_alloc[set]) * (uint64_t)sizeof(uint32_t));
		// printf("[VERBOSE] Text lengths allocated %ld bytes\n", (uint64_t)(num_pairs_alloc[set]) * (uint64_t)sizeof(uint32_t));
		// printf("[VERBOSE] Patterns allocated %ld bytes\n", (uint64_t)(num_pairs_alloc[set]) * (uint64_t)sizeof(char) * (uint64_t)longest_seq[set]);
		// printf("[VERBOSE] Texts allocated %ld bytes\n", (uint64_t)(num_pairs_alloc[set]) * (uint64_t)sizeof(char) * (uint64_t)longest_seq[set]);

	}


	char** Qpatterns;
	char** Qtexts;
	uint32_t** Qpattern_lengths;
	uint32_t** Qtext_lengths;
	int32_t** 	Qdpu_distances;
	char** 		Qcigar;
	char** 		Qcigar_debug;
	struct dpu_set_t* Qdpu_set;
	uint64_t Qnum_pairs[queued_sets];
	uint32_t queue = 0;
	uint32_t Qnr_dpus_set[queued_sets];
	Timer Qtimer[queued_sets];
	uint64_t Qnum_pairs_per_dpu[queued_sets];
	uint64_t Qcigar_length[queued_sets];
	uint64_t Qoffsets_size_per_tl[queued_sets];
	uint64_t Qnum_pairs_alloc[queued_sets];
	uint64_t Qlongest_seq[queued_sets];
	uint32_t Qnr_of_dpus[queued_sets];
	uint64_t Qmem_offset[queued_sets];
	uint32_t partition;
	uint32_t pair_reminder[p.num_sets];
	bool finished_queue[queued_sets];
	bool execution_queue[queued_sets];

	Qpatterns = (char**)malloc(queued_sets * sizeof(char*));
	Qtexts = (char**)malloc(queued_sets * sizeof(char*));
	Qpattern_lengths = (uint32_t**)malloc(queued_sets * sizeof(uint32_t*));
	Qtext_lengths = (uint32_t**)malloc(queued_sets * sizeof(uint32_t*));
	Qdpu_distances      = (int32_t**)malloc(queued_sets * sizeof(int32_t*));
	Qcigar        		= (char**)malloc(queued_sets * sizeof(char*));
	Qcigar_debug        = (char**)malloc(queued_sets * sizeof(char*));
	Qdpu_set = (struct dpu_set_t*)malloc(queued_sets * sizeof(struct dpu_set_t));

		// create queued sets
		for(set=0; set<p.num_sets; set++){
			if(num_pairs[set]==0) continue;
			// TODO check si las divisiones son 8 byte aligned
			if(queued_sets <= 0) continue;
			pair_reminder[set] = num_pairs[set] % (int64_t) set_partitions[set];
			num_pairs[set] = num_pairs[set] / (int64_t) set_partitions[set];
			flt_num_pairs_per_dpu = (float) num_pairs[set] / (float) nr_dpus_set[set];
			num_pairs_per_dpu[set] = (uint64_t) flt_num_pairs_per_dpu;
			reminder = fmod(flt_num_pairs_per_dpu, (float) NR_TASKLETS*2);
			if(reminder != 0) num_pairs_per_dpu[set] = (uint64_t) ((flt_num_pairs_per_dpu + NR_TASKLETS*2) - reminder);
			num_pairs_alloc[set] = num_pairs_per_dpu[set] * nr_dpus_set[set];

			dpu_distances[set]      = (int32_t*)realloc(dpu_distances[set], (size_t)num_pairs_alloc[set] * (size_t)sizeof(int32_t) * set_partitions[set]);
			cigar[set]        		= (char*)realloc(cigar[set], (size_t)num_pairs_alloc[set] *((size_t)cigar_length[set] * (size_t)sizeof(char) * set_partitions[set]));
			//cigar_debug[set]        = (char*)malloc(cigar_length[set] * sizeof(char));
			pattern_lengths[set] = (uint32_t *) realloc(pattern_lengths[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(uint32_t) * set_partitions[set]);
			text_lengths[set] = (uint32_t *) realloc(text_lengths[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(uint32_t) * set_partitions[set]);
			patterns[set] = (char *) realloc(patterns[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(char) * (size_t)longest_seq[set] * set_partitions[set]);
			texts[set] = (char *) realloc(texts[set],  (size_t)(num_pairs_alloc[set]) * (size_t)sizeof(char) * (size_t)longest_seq[set] * set_partitions[set]);
			memset(&pattern_lengths[set][num_pairs[set] * set_partitions[set]], 0, (num_pairs_alloc[set] * set_partitions[set] - num_pairs[set]* set_partitions[set]) * sizeof(uint32_t) );
			memset(&text_lengths[set][num_pairs[set] * set_partitions[set]],    0, (num_pairs_alloc[set] * set_partitions[set] - num_pairs[set]* set_partitions[set])  * sizeof(uint32_t) );

			// num_pairs_alloc[set] = num_pairs_alloc[set] / set_partitions[set];
			// num_pairs_per_dpu[set] = num_pairs_per_dpu[set] / set_partitions[set];
			for (partition = 1; partition < set_partitions[set]; partition++){
				// Slice sizes in half
				if(partition == set_partitions[set]-1 && pair_reminder[set]) num_pairs[set]++;
				Qnum_pairs[queue] = num_pairs[set];
				Qnum_pairs_alloc[queue] = num_pairs_alloc[set];
				Qnum_pairs_per_dpu[queue] = num_pairs_per_dpu[set];
				Qcigar_length[queue] = cigar_length[set];
				Qoffsets_size_per_tl[queue] = offsets_size_per_tl[set];
				Qlongest_seq[queue] = longest_seq[set];
				finished_queue[queue] = false;
				execution_queue[queue] = false;
				Qnr_dpus_set[queue] = nr_dpus_set[set];
				printf("[INFO] Number of dpus required %d by set %d in queue %d\n",Qnr_dpus_set[queue], set, queue);

				// Alloc queued structures
				Qdpu_distances[queue]      = (int32_t*)malloc((uint64_t)Qnum_pairs_alloc[queue] * (uint64_t)sizeof(int32_t));
				Qcigar[queue]        		= (char*)malloc((uint64_t)Qnum_pairs_alloc[queue] *((uint64_t)Qcigar_length[queue] * (uint64_t)sizeof(char)));
				Qcigar_debug[queue]        = (char*)malloc(Qcigar_length[queue] * sizeof(char));
				// Qpattern_lengths[queue] = (uint32_t *) malloc((uint64_t)(Qnum_pairs_alloc[queue]) * (uint64_t)sizeof(uint32_t));
				// Qtext_lengths[queue] = (uint32_t *) malloc((uint64_t)(Qnum_pairs_alloc[queue]) * (uint64_t)sizeof(uint32_t));
				// Qpatterns[queue] = (char *) malloc((uint64_t)(Qnum_pairs_alloc[queue]) * (uint64_t)sizeof(char) * (uint64_t)Qlongest_seq[queue]);
				// Qtexts[queue] = (char *) malloc((uint64_t)(Qnum_pairs_alloc[queue]) * (uint64_t)sizeof(char) * (uint64_t)Qlongest_seq[queue]);
				// memset(Qpattern_lengths[queue], 0, Qnum_pairs_alloc[queue] * sizeof(uint32_t));
				// memset(Qtext_lengths[queue],    0, Qnum_pairs_alloc[queue] * sizeof(uint32_t));

				// // move data from old structures to queue structures
				// Qpattern_lengths[queue] = &pattern_lengths[set][num_pairs_alloc[set]*partition];
				// Qtext_lengths[queue] = &text_lengths[set][num_pairs_alloc[set]*partition];
				// Qpatterns[queue] = &patterns[set][num_pairs_alloc[set]*partition];
				// Qtexts[queue] = &texts[set][num_pairs_alloc[set]*partition];
				Qpattern_lengths[queue] = &pattern_lengths[set][num_pairs[set]*partition];
				Qtext_lengths[queue] = &text_lengths[set][num_pairs[set]*partition];
				Qpatterns[queue] = &patterns[set][num_pairs[set]*partition];
				Qtexts[queue] = &texts[set][num_pairs[set]*partition];
				queue++;
			}
	}
	// stop(&initialization, 1);
	// printf("[TIME] Initializing variables (ms): ");
	// print(&initialization, 1, 1);
	// printf("\n");

for (uint32_t rep = 0; rep < p.n_warmup + p.n_reps; rep++) { // Loop adding compute weight and warmup
	for(set=0; set<p.num_sets; set++){
		if(num_pairs[set]==0)continue;
		run_set(patterns[set], texts[set],
				pattern_lengths[set], text_lengths[set],
				num_pairs[set], num_pairs_per_dpu[set],
				longest_seq[set], &dpu_set[set], nr_dpus_set[set], &nr_of_dpus[set],
				threshold, rep, p.n_warmup, &timer[set], set, p.num_sets, &mem_offset[set]);
		printf("[INFO] DPU set number %d launched\n", set);
		fflush(stdout);
	}

	bool set_finished = false;
	bool set_failed;
	uint32_t finished_sets = 0;
	uint32_t dpus_free = 0;
	uint64_t Qmram_usage;
	uint32_t queue;
	uint64_t debug_running=1;
	struct dpu_set_t dpu_fail;

	if (NR_DPUS >= p.num_sets){
		while(finished_sets<p.num_sets + queued_sets){
			for(set=0; set<p.num_sets; set++){
				if(!finished_set[set]) dpu_status(dpu_set[set], &set_finished, &set_failed);
				if(set_failed){
					printf("[ERROR] SET %d FAILED\n", set);
					DPU_FOREACH(dpu_set[set], dpu_fail) {
						DPU_ASSERT(dpu_log_read(dpu_fail, stdout));
					}
					DPU_ASSERT(dpu_free(dpu_set[set]));
					return 0;
				}
				if (set_finished)
				{
					finished_sets++;
					set_finished = false;
					finished_set[set] = true;
					dpus_free += nr_dpus_set[set];
					if (rep >= p.n_warmup)
					{
						stop(&timer[set], 2);
					#if ENERGY
						DPU_ASSERT(dpu_probe_stop(&probe));
					#endif
					}
					printf("[INFO] DPU set number %d finished\n", set);
					retrieve_results(dpu_set[set], dpu_distances[set],
						num_pairs_per_dpu[set], offsets_size_per_tl[set],
						mem_offset[set], cigar[set], cigar_length[set],
						&timer[set], rep, p.n_warmup, nr_of_dpus[set]);
				}
			}
			if(queued_sets<=0) continue;
			for(queue=0; queue<queued_sets; queue++){
				//if(queue == 1) continue;
				if(finished_queue[queue]) continue;
				if(execution_queue[queue]) dpu_status(Qdpu_set[queue], &set_finished, &set_failed);
				if(set_failed){
					printf("[ERROR] QUEUED SET %d FAILED\n", queue);
					DPU_ASSERT(dpu_free(Qdpu_set[queue]));
					return 1;
				}
				if (set_finished)
				{
					finished_sets++;
					finished_queue[queue] = true;
					dpus_free += Qnr_dpus_set[queue];
					execution_queue[queue] = false;
					set_finished = false;
					if (rep >= p.n_warmup)
					{
						stop(&Qtimer[queue], 2);
					#if ENERGY
						DPU_ASSERT(dpu_probe_stop(&probe));
					#endif
					}
					printf("[INFO] DPU QUEUED set %d finished\n", queue);
					retrieve_results(Qdpu_set[queue], Qdpu_distances[queue],
						Qnum_pairs_per_dpu[queue], Qoffsets_size_per_tl[queue],
						Qmem_offset[queue], Qcigar[queue], Qcigar_length[queue],
						&Qtimer[queue], rep, p.n_warmup, nr_of_dpus[set]);
					continue;
				}
				if(dpus_free < Qnr_dpus_set[queue]) continue;
				flt_num_pairs_per_dpu = (float) Qnum_pairs[queue] / (float) Qnr_dpus_set[queue];
				Qnum_pairs_per_dpu[queue] = (uint64_t) flt_num_pairs_per_dpu;
				reminder = fmod(flt_num_pairs_per_dpu, (float) NR_TASKLETS*2);
				if(reminder != 0) Qnum_pairs_per_dpu[queue] = (uint64_t) ((flt_num_pairs_per_dpu + NR_TASKLETS*2) - reminder);

				dpus_free = dpus_free - Qnr_dpus_set[queue];
				//printf("About to run a queue %d\n", queue);
				// printf("hello \n");
				run_set(Qpatterns[queue], Qtexts[queue],
					Qpattern_lengths[queue], Qtext_lengths[queue],
					Qnum_pairs[queue], Qnum_pairs_per_dpu[queue],
					Qlongest_seq[queue], &Qdpu_set[queue], Qnr_dpus_set[queue], &Qnr_of_dpus[queue],
					threshold, rep, p.n_warmup, &Qtimer[queue],queue, queued_sets, &Qmem_offset[queue]);
				execution_queue[queue] = true;
				printf("[INFO] DPU QUEUED set number %d launched\n", queue);
				fflush(stdout);
			}
		}
	}
}
for(set=0; set<p.num_sets; set++){
	if(num_pairs[set]==0) continue;
	validate_results(patterns[set], texts[set],
			pattern_lengths[set], text_lengths[set],
			num_pairs[set], longest_seq[set], dpu_set[set],
			cigar[set], cigar_length[set],
			dpu_distances[set], timer[set], set, p.n_reps);
}

	if(queued_sets>0){
		for(queue=0; queue<queued_sets; queue++){
			if(num_pairs[queue]==0) continue;
			validate_results(Qpatterns[queue], Qtexts[queue],
					Qpattern_lengths[queue], Qtext_lengths[queue],
					Qnum_pairs[queue], Qlongest_seq[queue], Qdpu_set[queue],
					Qcigar[queue], Qcigar_length[queue],
					Qdpu_distances[queue], Qtimer[queue], queue, p.n_reps);
		}
	}
	free(patterns);
	free(texts);
	free(Qpatterns);
	free(Qtexts);
	return 0;
}