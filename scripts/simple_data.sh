#!/bin/bash

NR_DPUS=2500
NR_TASKLETS=12
EXP_FILE="/home/upmem0046/aalonso/experiments/scalability/BIMSA-$1-tl.out"

echo "\n" >  $EXP_FILE
for tl in 1 4 8 12 16 14
do
for i in 1 2 3
do
    cd	/home/upmem0046/aalonso/ulsapim-github/ULSAPIM/upmem/
        make clean
        make NR_TASKLETS=$tl NR_DPUS=$NR_DPUS PERF=1 DYNAMIC=1 MAX_DISTANCE_THRESHOLD=100000
    echo "\n Running experiment for NR_TASKLETS=$tl NR_DPUS=$NR_DPUS file=$1 repetition=$i\n" >> $EXP_FILE
        ./bin/ulsapim_host -i /home/upmem0046/aalonso/inputs/$1 -s "$2" >> $EXP_FILE
done
done

echo "\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!!! EXECUTION FINISHED $1 !!!!!!!!!!!!!!!!!!!!!! \n\n\n\n\n" 
