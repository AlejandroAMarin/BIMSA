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
#define DPU_BINARY "./bin/ulsapim_iterative"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1500000

static bool read_sequence_metadata(char* input_filename, uint64_t* num_pairs, uint64_t* longest_seq)
{
    FILE *file;
    char line[MAX_LINE_LENGTH];
	uint64_t length;
	bool eof;

    // Open the file
    file = fopen(input_filename, "r");
    if (file == NULL) {
        printf("[ERROR] Input file not found!\n");
        return false;
    }
	if(*num_pairs==0){
    // Read each line and find the line with the maximum number of characters
    while (fgets(line, sizeof(line), file)) {
        length = strlen(line);
        if (length > *longest_seq) {
            *longest_seq = length;
        }
        (*num_pairs)++;
    }
	}else{
		for(uint64_t i=0; i<=*num_pairs; i++){
			eof = fgets(line, sizeof(line), file);
			if(!eof) break;
			length = strlen(line);
			if (length > *longest_seq) {
				*longest_seq = length;
			}
		}
	}

	(*longest_seq)--;
    // Close the file
    fclose(file);

    printf("[INFO] Number of lines: %ld\n", *num_pairs);
    printf("[INFO] Line with maximum characters (%ld characters):\n", *longest_seq);

    return true;
}

static bool read_sequence_data(char* input_filename, char* patterns, char* texts,
								uint32_t* pattern_lengths, uint32_t* text_lengths,
								uint64_t* num_pairs, uint64_t longest_seq)
{
	FILE *file;
	int pattern_length;
	int text_length;
	uint32_t seq_number = 0;
	size_t max_line_size = longest_seq + 16;

    // Open the file
    file = fopen(input_filename, "r");
    if (file == NULL) {
        printf("[ERROR] Input file not found!\n");
        return false;
    }

	char* line1 = malloc(max_line_size * sizeof(char));
	char* line2 = malloc(max_line_size * sizeof(char));

	while(seq_number < *num_pairs)
	{
		// Read queries
		pattern_length = getline(&line1, &max_line_size, file);
		text_length    = getline(&line2, &max_line_size, file);
		if (pattern_length == -1 || text_length == -1) break;

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

		strncpy(patterns + (longest_seq * seq_number), line1 + 1, pattern_length);
		strncpy(texts    + (longest_seq * seq_number), line2 + 1, text_length);

		patterns[(longest_seq * seq_number) + pattern_length] = '\0';
		texts[(longest_seq * seq_number) + text_length] = '\0';

		// Save sizes after removing start and end of line characters
		pattern_lengths[seq_number] = pattern_length;
		text_lengths[seq_number]    = text_length;
		seq_number++;
	}

	*num_pairs = seq_number;
	free(line1);
	free(line2);
    return true;
}


// Fill input arrays
// static bool read_inputs(char * input_filename, char* patterns, char* texts, uint64_t * num_pairs, uint32_t * pattern_lengths, uint32_t * text_lengths, uint64_t longest_seq)
// {
// 	FILE * input_file = fopen(input_filename, "r");

// 	if(input_file == NULL)
// 	{
// 		printf("[ERROR] Input file not found!\n");
// 		return false;
// 	} else {
// 		// Parameters
// 		int pattern_length;
// 		int text_length;
// 		uint32_t seq_number       = 0;
// 		size_t max_line_size = longest_seq + 16;

// 		char* line1 = malloc(max_line_size * sizeof(char));
// 		char* line2 = malloc(max_line_size * sizeof(char));

// 		// Discard first line
// 		pattern_length = getline(&line1, &max_line_size, input_file);

// 		while(seq_number < *num_pairs)
// 		{
// 			// Read queries
// 			pattern_length = getline(&line1, &max_line_size, input_file);
// 			text_length    = getline(&line2, &max_line_size, input_file);
// 			if (pattern_length == -1 || text_length == -1) break;

// 			pattern_length -= 2;
// 			text_length -= 2;

// 			strncpy(patterns + (longest_seq * seq_number), line1 + 1, pattern_length);
// 			strncpy(texts    + (longest_seq * seq_number), line2 + 1, text_length);

// 			patterns[(longest_seq * seq_number) + pattern_length] = '\0';
// 			texts[(longest_seq * seq_number) + text_length] = '\0';

// 			// Save sizes after removing start and end of line characters
// 			pattern_lengths[seq_number] = pattern_length;
// 			text_lengths[seq_number]    = text_length;
// 			seq_number++;
// 		}

// 		*num_pairs = seq_number;
// 		free(line1);
// 		free(line2);
// 		return true;
// 	}
// }

char* cigar_pprint2(char* input, FILE *output) {
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

int main(int argc, char **argv){
  // Timer declaration
	Timer timer;
	struct Params p = input_params(argc, argv);
	struct dpu_set_t dpu_set, dpu;
	uint32_t nr_of_dpus;

	printf(ANSI_COLOR_YELLOW "###############################################\n" ANSI_COLOR_RESET);
	printf(ANSI_COLOR_YELLOW "|---------------- BI-WFA-UPMEM ---------------|\n" ANSI_COLOR_RESET);
	printf(ANSI_COLOR_YELLOW "###############################################\n" ANSI_COLOR_RESET);

  if(strcmp(p.input_file, "UNSET") == 0)
	{
    printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED" You need to specify an input filename!\n"ANSI_COLOR_RESET);
    printf("[INFO] Please run %s -h for help.\n", argv[0]);

    return 0;
	}
	printf("[INFO] Reading input file...\n");
	printf("[INFO] " ANSI_COLOR_BLUE "Input file: %s" ANSI_COLOR_RESET "\n", p.input_file);
	fflush(stdout);
	uint64_t longest_seq = 0;
	uint64_t num_pairs = 0;
	ewf_offset_t threshold = (ewf_offset_t)sqrt(4*BLOCK_SIZE);
	bool file_status;
	#if PERF
	double cc = 0;
    double cc_min = 0;
	#endif

	if(NR_TASKLETS * (BLOCK_SIZE*5 + BLOCK_SIZE_INPUTS*6 + 8 + STACK_SIZE_DEFAULT) > WRAM_LIMIT){
		printf("[ERROR] Exceeded WRAM capacity");
		return 0;
	}
	//printf("Printing blocksize: %d and threshold %d\n", BLOCK_SIZE, threshold);

	//read_inputs(p.input_file, patterns, texts, &num_pairs, pattern_lengths, text_lengths, longest_seq);

  	//FILE * fp = fopen(p.input_file, "r");
	FILE * output;

	// Allocate and create input data
	// All patterns and text are allocated of the maximum size, then read only the required size
	char* patterns;
	char* texts;
  	char* rev_patterns;
  	char* rev_texts;

	// assert(fscanf(fp, "%ld%ld", &lines, &longest_seq) == 2);
	if(p.max_pairs != 0)
		num_pairs = p.max_pairs;
	// fclose(fp);
	file_status = read_sequence_metadata(p.input_file, &num_pairs, &longest_seq);
	assert(file_status);

	if(longest_seq % 8 != 0) longest_seq = ((longest_seq + 8 - 1) / 8) * 8;

	num_pairs = num_pairs / 2;

	if(p.max_pairs != 0)
		num_pairs = p.max_pairs;

	printf("[VERBOSE] num pairs %ld and longest seq %ld\n", num_pairs, longest_seq);

	// Adjust the number of pairs to the number of dpus
	float flt_num_pairs_per_dpu = (float) num_pairs / (float) NR_DPUS;

	uint64_t num_pairs_per_dpu = (uint64_t) flt_num_pairs_per_dpu;

	float reminder = fmod(flt_num_pairs_per_dpu, (float) (NR_TASKLETS*2));

	if(reminder != 0) num_pairs_per_dpu = (uint64_t) ((flt_num_pairs_per_dpu + (NR_TASKLETS*2)) - reminder);

	uint64_t num_pairs_alloc = num_pairs_per_dpu * NR_DPUS;
	printf("[VERBOSE] num pairs alloc %ld vs num pairs %ld\n", num_pairs_alloc, num_pairs);

	// Allocating arrays of sequences
	patterns = malloc(num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));
	texts    = malloc(num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));
  	rev_patterns = malloc(num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));
	rev_texts    = malloc(num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));
	// Allocating arrays of pattern and text lengths
	uint32_t * pattern_lengths = malloc(num_pairs_alloc * sizeof(uint32_t));
	uint32_t * text_lengths    = malloc(num_pairs_alloc * sizeof(uint32_t));
	memset(pattern_lengths, 0, num_pairs_alloc * sizeof(uint32_t));
	memset(text_lengths,    0, num_pairs_alloc * sizeof(uint32_t));
	uint64_t cigar_length = 2 * longest_seq;
	uint64_t cigar_length_adjustment = cigar_length % BLOCK_SIZE_INPUTS;
	if(cigar_length_adjustment != 0)
		cigar_length += (BLOCK_SIZE_INPUTS - cigar_length_adjustment);
	printf("[VERBOSE] CIGAR length %ld\n", cigar_length);
	// Allocating result arrays
	ewf_offset_t * host_distances       = malloc(num_pairs * sizeof(ewf_offset_t));
	int32_t * dpu_distances        		= malloc(num_pairs_alloc * sizeof(int32_t));
	char 		 * cigar        		= malloc(num_pairs_alloc * ((uint64_t)cigar_length * (uint64_t)sizeof(char)));
	char 		 * cigar_debug        	= malloc(cigar_length * sizeof(char));
	printf("[VERBOSE] Distance allocated %ld bytes\n", num_pairs_alloc * sizeof(int32_t));
	printf("[VERBOSE] CIGAR allocated %ld bytes\n", num_pairs_alloc * ((uint64_t)cigar_length * (uint64_t)sizeof(char)));
	printf("[VERBOSE] CIGAR debug allocated %ld bytes\n", cigar_length * sizeof(char));
	printf("[VERBOSE] Pattern lengths allocated %ld bytes\n", num_pairs_alloc * sizeof(uint32_t));
	printf("[VERBOSE] Text lengths allocated %ld bytes\n", num_pairs_alloc * sizeof(uint32_t));
	printf("[VERBOSE] Patterns allocated %ld bytes\n", num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));
	printf("[VERBOSE] Texts allocated %ld bytes\n", num_pairs_alloc * longest_seq * (uint64_t) sizeof(char));

	uint32_t offsets_size_per_tl      = 4*longest_seq * sizeof(ewf_offset_t) * MAX_ERROR;
	uint32_t offsets_size_adjustement = offsets_size_per_tl % BLOCK_SIZE;
	if(offsets_size_adjustement != 0)
		offsets_size_per_tl += (BLOCK_SIZE - offsets_size_adjustement);
	printf("[VERBOSE] offset_size length %d\n", offsets_size_per_tl);
	

	printf("[INFO] " ANSI_COLOR_BLUE "num_pairs: \t %ld" ANSI_COLOR_RESET "\n", num_pairs);
	printf("[INFO] " ANSI_COLOR_BLUE "longest_sequence: %ld" ANSI_COLOR_RESET "\n", longest_seq);
	printf("[INFO] Pairs per DPU (adjusted): (%ld)\n", num_pairs_per_dpu);
	printf("[INFO] Total WRAM usage: (%d) bytes \n", NR_TASKLETS * (BLOCK_SIZE*5 + BLOCK_SIZE_INPUTS*6 + 8 + STACK_SIZE_DEFAULT));



	if(patterns == NULL)
	{
		printf("Error when allocating patterns memory\n");
		return 0;
	}

	if(texts == NULL)
	{
		printf("Error when allocating text memory\n");
		return 0;
	}

	if(rev_patterns == NULL)
	{
		printf("Error when allocating patterns memory\n");
		return 0;
	}

	if(rev_texts == NULL)
	{
		printf("Error when allocating text memory\n");
		return 0;
	}

  //read_inputs(p.input_file, patterns, texts, &num_pairs, pattern_lengths, text_lengths, longest_seq);
  	file_status = read_sequence_data(p.input_file, patterns, texts,
						pattern_lengths, text_lengths,
						&num_pairs, longest_seq);
	assert(file_status);

	printf("[INFO] Reading input OK.\n");
	printf("[INFO] Allocating DPUs...");
	fflush(stdout);

  // Reading inputs and alocating memory finished

  // Allocate DPUs and load binary
	DPU_ASSERT(dpu_alloc(NR_DPUS, NULL, &dpu_set));
	DPU_ASSERT(dpu_load(dpu_set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(dpu_set, &nr_of_dpus));

  printf("OK.\n");
	printf("[INFO] " ANSI_COLOR_BLUE "num_dpus:%d" ANSI_COLOR_RESET "\n", nr_of_dpus);

  unsigned int kernel = 0;
	dpu_arguments_t input_arguments = {longest_seq, num_pairs_per_dpu, threshold, kernel};
	uint32_t mem_offset;
	uint32_t final_offset = BLOCK_SIZE_INPUTS + 3* (longest_seq * num_pairs_per_dpu * sizeof(char))
	+ 2*(num_pairs_per_dpu * sizeof(uint32_t)) + NR_TASKLETS * (4 * offsets_size_per_tl) + (num_pairs_per_dpu * sizeof(int64_t))
	+ num_pairs_per_dpu * cigar_length * sizeof(char) + NR_TASKLETS * (cigar_length + 8) + NR_TASKLETS*QUEUE_SZ;

	if(final_offset> MRAM_LIMIT){
	       	printf("[ERROR] MRAM limit surpassed\n");
		return 1;
	}

for (uint32_t rep = 0; rep < p.n_warmup + p.n_reps; rep++) { // Loop adding compute weight and warmup

		printf("[INFO] Transferring inputs from to DPUs...");
		fflush(stdout);
		if (rep >= p.n_warmup)
			start(&timer, 1, rep - p.n_warmup);
		uint32_t i = 0;

		// Copy input arguments to dpu
		DPU_FOREACH(dpu_set, dpu) {
			DPU_ASSERT(dpu_copy_to(dpu, "DPU_INPUT_ARGUMENTS", 0, (const void *) &input_arguments, sizeof(input_arguments)));
			i++;
		}

		i = 0;
		mem_offset = BLOCK_SIZE_INPUTS; // offset points to first MRAM address
		printf("[VERBOSE] HOST ma pattern %d\n", mem_offset);
		// Patterns transfer
		DPU_FOREACH(dpu_set, dpu, i)
		{
			DPU_ASSERT(dpu_prepare_xfer(dpu, patterns + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST pattern sent %ld bytes\n", longest_seq * num_pairs_per_dpu * sizeof(char));

    	i = 0;
		mem_offset += (longest_seq * num_pairs_per_dpu * sizeof(char)); // offset points to end of patterns array
		// Texts transfer
		printf("[VERBOSE] HOST ma texts %d\n", mem_offset);
		DPU_FOREACH(dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, texts + i * num_pairs_per_dpu * longest_seq));
		}

		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, mem_offset, longest_seq * num_pairs_per_dpu * sizeof(char), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST text sent %ld bytes\n", longest_seq * num_pairs_per_dpu * sizeof(char));

		i = 0;
		mem_offset += (longest_seq * num_pairs_per_dpu * sizeof(char)); // offset points to end of texts array
		// Pattern lengths array transfer
		printf("[VERBOSE] HOST ma pattern lenghts %d\n", mem_offset);
		DPU_FOREACH(dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, pattern_lengths + i * num_pairs_per_dpu));
		}

		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, mem_offset, num_pairs_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST pattern lengths sent %ld bytes\n", num_pairs_per_dpu * sizeof(uint32_t));

		i = 0;
		mem_offset += (num_pairs_per_dpu * sizeof(uint32_t)); // offset points to end of pattern lengths array
		// Texts lengths array transfer
		printf("[VERBOSE] HOST ma text lengths %d\n", mem_offset);
		DPU_FOREACH(dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, text_lengths + i * num_pairs_per_dpu));
		}

		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, DPU_MRAM_HEAP_POINTER_NAME, mem_offset, num_pairs_per_dpu * sizeof(uint32_t), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST text lengths sent %ld bytes\n", num_pairs_per_dpu * sizeof(uint32_t));

		mem_offset += (num_pairs_per_dpu * sizeof(uint32_t)); // offset points to end of texts lengths array (start of results array)

		if (rep >= p.n_warmup)
			stop(&timer, 1);
		printf("OK.\n");
		// Run kernel on DPUs

		printf("[INFO] Running kernel in DPUs...");
		fflush(stdout);
		if (rep >= p.n_warmup)
		{
			start(&timer, 2, rep - p.n_warmup);
#if ENERGY
			DPU_ASSERT(dpu_probe_start(&probe));
#endif
		}

		DPU_ASSERT(dpu_launch(dpu_set, DPU_SYNCHRONOUS)); // DPU can be used asyncronously (have into account for future improvements)

		if (rep >= p.n_warmup)
		{
			stop(&timer, 2);
#if ENERGY
			DPU_ASSERT(dpu_probe_stop(&probe));
#endif
		}

		printf("OK.\n");

		printf("[INFO] Retrieving results to host...\n");
		fflush(stdout);
		if (rep >= p.n_warmup)
			start(&timer, 3, rep - p.n_warmup);
		#if PERF
		dpu_results_t results[nr_of_dpus];
		#endif
		i = 0;
		// Copy back the DPU results (mem_offset already points to the results starting memory position)
		printf("\n[VERBOSE] HOST ma distances %d\n", mem_offset);
		DPU_FOREACH(dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_distances + i * num_pairs_per_dpu));

			#if PERF
			results[i].counter = 0;
			// Retrieve tasklet timings
			for (unsigned int each_tasklet = 0; each_tasklet < NR_TASKLETS; each_tasklet++) {
				dpu_results_t result;
				result.counter = 0;
				DPU_ASSERT(dpu_copy_from(dpu, "DPU_RESULTS", each_tasklet * sizeof(dpu_results_t), &result, sizeof(dpu_results_t)));
				if (result.counter > results[i].counter)
					results[i].counter = result.counter;
			}
			#endif
		}
		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
					mem_offset, num_pairs_per_dpu * sizeof(int32_t), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST distances recieved %ld bytes\n", num_pairs_per_dpu * sizeof(int32_t));


		#if PERF
        uint64_t max_cycles = 0;
        uint64_t min_cycles = 0xFFFFFFFFFFFFFFFF;
        // Print performance results
        if(rep >= p.n_warmup){
            i = 0;
            DPU_FOREACH(dpu_set, dpu) {
                if(results[i].counter > max_cycles)
                    max_cycles = results[i].counter;
                if(results[i].counter < min_cycles)
                    min_cycles = results[i].counter;
                i++;
            }
            cc += (double)max_cycles;
            cc_min += (double)min_cycles;
        }
		#endif

		// move mem offset to the end of the 4 wavefronts
		// Offsets per tl size calculatio
		mem_offset += NR_TASKLETS * (4 * offsets_size_per_tl) + (num_pairs_per_dpu * sizeof(int32_t));
		i = 0;
		printf("\n[VERBOSE] HOST ma CIGAR %d\n", mem_offset);
		DPU_FOREACH(dpu_set, dpu, i) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, cigar + i * num_pairs_per_dpu * cigar_length));
		}
		DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, DPU_MRAM_HEAP_POINTER_NAME,
					mem_offset, num_pairs_per_dpu * cigar_length * sizeof(char), DPU_XFER_DEFAULT));
		printf("[VERBOSE] HOST CIGAR recieved %ld bytes\n", num_pairs_per_dpu * cigar_length * sizeof(char));
		if(rep >= p.n_warmup)
			stop(&timer, 3);

		printf("OK.\n");
#if PRINT
		printf("LOGS\n");
		DPU_FOREACH(dpu_set, dpu) {
			DPU_ASSERT(dpu_log_read(dpu, stdout));
		}
#endif

}

	printf("OK.\n");

	// Print timing results
	printf("\n[TIME] CPU-DPU (ms): ");
	print(&timer, 1, p.n_reps);
	printf("\n[TIME] DPU Kernel (ms): ");
	print(&timer, 2, p.n_reps);
	printf("\n[TIME] DPU-CPU Time (ms): ");
	print(&timer, 3, p.n_reps);
	printf("\n");
	#if PERF
	printf("[PERF] DPU kernel perf data: %g\n", cc);
	#endif

	int cigar_d;
	int invalid;
	int failed = 0;
	int cigar_failed = 0;
	int distance_failed = 0;
	output = fopen("dpu_output", "w");
	printf("[INFO] Validating results...\n");
	#pragma omp parallel for reduction(+:failed)
	for (uint32_t j = 0; j < num_pairs; j++)
	{
		// for (uint32_t i = 0; i<= cigar_length; i++){
		// 	if(cigar[j*cigar_length+i] == '\0') printf("0");
		// 	printf("%c", cigar[j*cigar_length+i]);
		// }
		// printf("\n");
		fprintf(output,"%d\t",dpu_distances[j]);
		cigar_debug = cigar_pprint2(&cigar[j*cigar_length], output);
		cigar_d = cigar_distance(&cigar[j*cigar_length]);
		if(cigar_d != dpu_distances[j])
		{
			#if DEBUG
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED"Mismatch result at index %d: host_distance=%d dpu_distance=%d cigar_distance=%d\n"ANSI_COLOR_RESET, j,host_distances[j],dpu_distances[j], cigar_d);
			printf("Printing wrowng CIGAR: %s\n",cigar_debug);
			#endif
			distance_failed++;
		}
		invalid = check_cigar(pattern_lengths[j], text_lengths[j],
							&patterns[j*longest_seq],
							&texts[j*longest_seq],
							&cigar[j*cigar_length]);
		if(invalid){
			#if DEBUG
			printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] "ANSI_COLOR_RED"CIGAR IS NOT VALID at index %d: host_distance=%d dpu_distance=%d cigar_distance=%d\n"ANSI_COLOR_RESET, j,host_distances[j],dpu_distances[j], cigar_d);
			#endif
			cigar_failed++;
		}
	}

	failed = cigar_failed;
	if(failed < distance_failed) failed = distance_failed;

	if (!failed) {
		printf("[" ANSI_COLOR_GREEN "-OK-" ANSI_COLOR_RESET "] Results are equal.\n");
	} else {
		printf("[" ANSI_COLOR_RED "ERR-" ANSI_COLOR_RESET "] Results differ! %d/%ld failed...\n", failed, num_pairs);

	}

	DPU_ASSERT(dpu_free(dpu_set));

	free(patterns);
	free(texts);
  	free(rev_patterns);
	free(rev_texts);
  return 0;
}
	double cc = 0;
    double cc_min = 0;
