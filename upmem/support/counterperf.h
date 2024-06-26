#include <perfcounter.h>

// Timer
typedef struct perfcounter_cycles{
    perfcounter_t start;
    perfcounter_t end;
    perfcounter_t end2;

}perfcounter_cycles;

void timer_start(perfcounter_cycles *cycles){
    // if(me() == 0){
    // #if PERF_CYCLES
    // perfcounter_config(COUNT_CYCLES, true);
    // #endif
    // #if PERF_INSTRUCTIONS
    // perfcounter_config(COUNT_INSTRUCTIONS, true);
    // #endif
    // }
    cycles->start = perfcounter_get(); // START TIMER
}

uint64_t timer_stop(perfcounter_cycles *cycles){
    //if(me() == 0){
    cycles->end = perfcounter_get(); // STOP TIMER
    cycles->end2 = perfcounter_get();
    if(cycles->end < cycles->start){
        printf("[DPU ]COUNTER OVERFLOW\n");
        return 0;
    } else{
        return cycles->end - cycles->start;
        //return(((uint64_t)((uint32_t)(((cycles->end >> 4) - (cycles->start >> 4)) - ((cycles->end2 >> 4) - (cycles->end >> 4))))) << 4);
    }
    //cycles->end2 = perfcounter_get(); // STOP TIMER
    //return(((uint64_t)((uint32_t)(((cycles->end >> 4) - (cycles->start >> 4)) - ((cycles->end2 >> 4) - (cycles->end >> 4))))) << 4);
    //}
    return 0;
}