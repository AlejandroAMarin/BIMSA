/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2019 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignments Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/*
 * Translate k and offset to coordinates h,v
 */
#define EWAVEFRONT_V(k,offset) ((offset)-(k))
#define EWAVEFRONT_H(k,offset) (offset)

#define EWAVEFRONT_DIAGONAL(h,v) ((h)-(v))
#define EWAVEFRONT_OFFSET(h,v)   (h)

#define MAX(a,b) (((a)>=(b))?(a):(b))
#define ABS(a) (((a)>=0)?(a):-(a))

#define SWAP(a,b) do {__typeof__(a) aux = a; a = b; b = aux;} while (0)

/*
 * Wavefront
 */
typedef int16_t ewf_offset_t;  // Edit Wavefront Offset
typedef struct {
  int lo;                      // Effective lowest diagonal (inclusive)
  int hi;                      // Effective highest diagonal (inclusive)
  ewf_offset_t* offsets;       // Offsets
} edit_wavefront_t;
/*
 * Edit Wavefronts
 */
typedef struct {
  // Dimensions
  int pattern_length;
  int text_length;
  int max_distance;
  // Waves Offsets
  edit_wavefront_t* wavefront;
  edit_wavefront_t* next_wavefront;
} edit_wavefronts_t;
void edit_wavefronts_init(
    edit_wavefronts_t* const wavefronts,
    const int pattern_length,
    const int text_length) {
  // Dimensions
  wavefronts->pattern_length = pattern_length;
  wavefronts->text_length = text_length;
  wavefronts->max_distance = pattern_length + text_length;
  // Allocate wavefronts
  wavefronts->wavefront = malloc(sizeof(edit_wavefront_t));
  wavefronts->wavefront->offsets = (ewf_offset_t*)calloc(2*wavefronts->max_distance,sizeof(ewf_offset_t)) + wavefronts->max_distance;
  wavefronts->next_wavefront = malloc(sizeof(edit_wavefront_t));
  wavefronts->next_wavefront->offsets = (ewf_offset_t*)calloc(2*wavefronts->max_distance,sizeof(ewf_offset_t)) + wavefronts->max_distance;
}
/*
 * Extend Wavefront
 */
void edit_wavefronts_extend_wavefront(
    edit_wavefront_t* const wavefront,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int distance) {
  // Parameters
  ewf_offset_t* const offsets = wavefront->offsets;
  const int k_min = wavefront->lo;
  const int k_max = wavefront->hi;
  // Extend diagonally each wavefront point
  int k;
  for (k=k_min;k<=k_max;++k) {
    int v = EWAVEFRONT_V(k,offsets[k]);
    int h = EWAVEFRONT_H(k,offsets[k]);
    while (v<pattern_length && h<text_length && pattern[v++]==text[h++]) {
      ++(offsets[k]);
    }
  }
}
/*
 * Edit Wavefront Compute
 */
void edit_wavefronts_compute_wavefront(
    edit_wavefront_t* const wavefront,
    edit_wavefront_t* const next_wavefront,
    const int pattern_length,
    const int text_length,
    const int distance) {
  // Fetch wavefronts
  const int hi = wavefront->hi;
  const int lo = wavefront->lo;
  next_wavefront->hi = hi+1;
  next_wavefront->lo = lo-1;
  // Fetch offsets
  ewf_offset_t* const offsets = wavefront->offsets;
  ewf_offset_t* const next_offsets = next_wavefront->offsets;

  // Loop peeling (k=lo-1)
  next_offsets[lo-1] = offsets[lo];
  // Loop peeling (k=lo)
  const ewf_offset_t bottom_upper_del = ((lo+1) <= hi) ? offsets[lo+1] : -1;
  next_offsets[lo] = MAX(offsets[lo]+1,bottom_upper_del);
  // Compute next wavefront starting point
  int k;
  #pragma GCC ivdep
  for (k=lo+1;k<=hi-1;++k) {

      const int del = offsets[k+1]; // Upper
      const int sub = offsets[k] + 1; // Mid
      const int ins = offsets[k-1] + 1; // Lower
      next_offsets[k] = MAX(sub,MAX(ins,del)); // MAX

    //const ewf_offset_t max_ins_sub = MAX(offsets[k],offsets[k-1]) + 1;
    //next_offsets[k] = MAX(max_ins_sub,offsets[k+1]);
  }
  // Loop peeling (k=hi)
  const ewf_offset_t top_lower_ins = (lo <= (hi-1)) ? offsets[hi-1] : -1;
  next_offsets[hi] = MAX(offsets[hi],top_lower_ins) + 1;
  // Loop peeling (k=hi+1)
  next_offsets[hi+1] = offsets[hi] + 1;
}
/*
 * Edit distance alignment using wavefronts
 */
int edit_wavefronts_distance(
    edit_wavefronts_t* const wavefronts,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length) {
  // Parameters
  const int max_distance = pattern_length + text_length;
  const int target_k = EWAVEFRONT_DIAGONAL(text_length,pattern_length);
  const int target_k_abs = ABS(target_k);
  const ewf_offset_t target_offset = EWAVEFRONT_OFFSET(text_length,pattern_length);
  // Init wavefronts
  int distance;
  wavefronts->wavefront->lo = 0;
  wavefronts->wavefront->hi = 0;
  wavefronts->wavefront->offsets[0] = 0;
  // Compute wavefronts for increasing distance
  for (distance=0;distance<max_distance;++distance) {
    // Extend diagonally each wavefront point
    edit_wavefronts_extend_wavefront(
        wavefronts->wavefront,pattern,pattern_length,
        text,text_length,distance);
    // Exit condition
    if (target_k_abs <= distance &&
        wavefronts->wavefront->offsets[target_k] == target_offset) break;
    // Compute next wavefront starting point
    edit_wavefronts_compute_wavefront(
        wavefronts->wavefront,wavefronts->next_wavefront,
        pattern_length,text_length,distance+1);
    // Swap
    SWAP(wavefronts->wavefront,wavefronts->next_wavefront);
  }
  // Return distance
  return distance;
}
int main(int argc,char* argv[]) {
  // Buffers
  char* pattern_mem =
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
      "TCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGTTCTTTACTCGCGCGTTGGAGAAATACAATAGT"
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY";
  char* text_mem    =
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"
      "TCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGTTCTATACTGCGCGTTTGGAGAAATAAAATAGT"
      "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY";
  // Pattern & Text
  char* pattern = pattern_mem + 64;
  char* text = text_mem + 64;
  const int pattern_length = strlen(pattern_mem)-2*64;
  const int text_length = strlen(text_mem)-2*64;
  const int reps = 10000000;
  // Init Wavefronts
  edit_wavefronts_t wavefronts;
  edit_wavefronts_init(&wavefronts,pattern_length,text_length);
  int i, dummy = 0 ;
  for (i=0;i<reps;++i) {
    dummy += edit_wavefronts_distance(&wavefronts,pattern,pattern_length,text,text_length);
  }
  printf("Result: %d\n", dummy);
}
