#include "timer.h"
#include <stddef.h>
#include <stdio.h>

void timer_init(Timer *timer, int i){
    timer->time[i] = 0.0;
}

void start(Timer *timer, int i) {
    gettimeofday(&timer->startTime[i], NULL);
}

void stop(Timer *timer, int i) {
    gettimeofday(&timer->stopTime[i], NULL);
    timer->time[i] += (timer->stopTime[i].tv_sec - timer->startTime[i].tv_sec) * 1000000.0 +
                      (timer->stopTime[i].tv_usec - timer->startTime[i].tv_usec);
}

void print(Timer *timer, int i, int REP) { printf("%f\t", timer->time[i] / (1000 * REP)); }