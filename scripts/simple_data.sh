#!/bin/bash

NR_DPUS=1000
BL=8
BLI=5
SZ=1024
NR_SEQ=5000000
NR_TASKLETS=24
EXP_FILE="/home/upmem0046/aalonso/experiments/scalability/biwfa-iterative.out"

echo "\n" >  $EXP_FILE
for seqlen in 150 1000
do
    for error in 2 5 10
    do
    cd	/home/upmem0046/aalonso/wfa-upmem/upmem/bi_wfa
        make clean
        make NR_TASKLETS="$NR_TASKLETS" NR_DPUS="$NR_DPUS" BL="$BL" BLI="$BLI" STACK_SIZE_DEFAULT="$SZ"
    echo "\n Running experiment for NR_TASKLETS=$NR_TASKLETS NR_DPUS=$NR_DPUS BL=$BL BLI=$BLI STACK_SIZE_DEFAULT=$SZ seqlen=$seqlen error=$error \n" >> $EXP_FILE
        ./bin/biwfa_host -i /home/upmem0046/aalonso/inputs/n"$NR_SEQ"_l"$seqlen"_e"$error".seq >> $EXP_FILE
    #cd /home/upmem0046/aalonso/WFA2-lib
        #./bin/align_benchmark -a edit-wfa --num-threads 32 --wfa-memory-mode ultralow -i /home/upmem0046/aalonso/inputs/n"$NR_SEQ"_l"$seqlen"_e"$error".seq -o cpu_output 2>>  /home/upmem0046/aalonso/simple_test.out
    done
done
