#!/bin/bash

NR_DPUS=2500
NR_TASKLETS=12
EXP_FILE="/home/upmem0046/aalonso/experiments/scalability/BIMSA-$1-tl-part3.out"

echo "\n" >  $EXP_FILE
for tl in 4 8 12
do
for i in 1 2 3
do
    cd	/home/upmem0046/aalonso/ulsapim-github/ULSAPIM/upmem/
        make clean
        make NR_TASKLETS=$tl NR_DPUS=$NR_DPUS DYNAMIC=0 MAX_DISTANCE_THRESHOLD=50000
    echo "\n Running experiment for NR_TASKLETS=$tl NR_DPUS=$NR_DPUS file=$1 repetition=$i\n" >> $EXP_FILE
        ./bin/bimsa_host -i /home/upmem0046/aalonso/inputs/$1 -s "$2" >> $EXP_FILE
done
done
for tl in 14 16 18
do
for i in 1 2 3
do
    cd	/home/upmem0046/aalonso/ulsapim-github/ULSAPIM/upmem/
        make clean
        make NR_TASKLETS=$tl NR_DPUS=$NR_DPUS DYNAMIC=0 MAX_DISTANCE_THRESHOLD=50000 CIGART=7 TASKT=7
    echo "\n Running experiment for NR_TASKLETS=$tl NR_DPUS=$NR_DPUS file=$1 repetition=$i\n" >> $EXP_FILE
        ./bin/bimsa_host -i /home/upmem0046/aalonso/inputs/$1 -s "$2" >> $EXP_FILE
done
done

echo "\n\n\n\n\n !!!!!!!!!!!!!!!!!!!!!! EXECUTION FINISHED $1 !!!!!!!!!!!!!!!!!!!!!! \n\n\n\n\n" 
