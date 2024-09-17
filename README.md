# BIMSA: Bidirectional In Memory Sequence Aligner
Implementation of a distance-based Bidirectional Wavefront Algorithm tailored for the UPMEM architecture.

## Setting up the UPMEM local simulator
Python is required if not installed.
```
sudo apt install python3
```

Download the UPMEM SDK from the following source:

https://sdk.upmem.com/

1. Untar the binary package, for instance into `$HOME/upmem-sdk` directory.
2. Source the script `$HOME/upmem-sdk/upmem_env.sh` to set environment variables to appropriate values.
```
mkdir $HOME/upmem-sdk
tar --strip-components=1 -zxvf $HOME/upmem-2023.2.0-Linux-x86_64.tar.gz -C $HOME/upmem-sdk
source $HOME/upmem-sdk/upmem_env.sh simulator
```

## Getting started with BIMSA
Clone this repository and init the submodule:
```
git clone git@github.com:AlejandroAMarin/BIMSA.git
cd BIMSA
git submodule init
git submodule update
```

Run a unitest:
```
python3 run_bimsa.py -d 2 -t 4 -f $HOME/BIMSA/inputs/wfa.utest.seq -s 20000
```

## Generating new synthetic inputs
To generate new synthetic inputs, first compile the tool, second run the script indicating your file sizes and error.
```
cd $HOME/BIMSA/tools/generate_dataset/
make
cd $HOME/BIMSA/scripts/
python3 gen_input.py -n 10 -l 100 -e 0.02
```
Run BIMSA on the new created inputs:
```
cd $HOME/BIMSA
python3 run_bimsa.py -d 2 -t 2 -f $HOME/projects/BIMSA/inputs/n10_l100_e2.seq -s 150
```
Be aware that the file generator creates files with sequence lengths aproximate to the indicated length, so the size `(-s)` on the BIMSA arguments must be at least 5% beyond the dataset sequence size. If the size is inferior to a sequence, this sequence will be discarded.

Generated files can be modified in number of pairs, sequence length and error rate.
```
usage: gen_input.py [-h] -n NR_PAIRS -l SEQ_LENGTH -e ERROR

Generate datasets with specified parameters

arguments:
  -h, --help            show this help message and exit
  -n NR_PAIRS, --nr_pairs NR_PAIRS
                        Number of pairs (NR_PAIRS)
  -l SEQ_LENGTH, --seq_length SEQ_LENGTH
                        Sequence length (SEQ_LENGTH)
  -e ERROR, --error ERROR
                        Error rate (ERROR)
```
## UPMEM simulator limitations
The UPMEM simulator is only able to run on 8 DPUs max and 24 tasklets.
The times probided by the simulator are not equivalent to a real hardware execution. The results optained from the alignment are the same as an UPMEM server.

## Running BIMSA on a UPMEM server
On a UPMEM server all libraries should be installed, so BIMSA should run straighforward by cloning the repositorie. If there are any SDK problems, contact the UPMEM team for troubleshooting.
```
git clone git@github.com:AlejandroAMarin/BIMSA.git
cd BIMSA
git submodule init
git submodule update
```
Create a large enough file to feed all the DPUs
```
cd $HOME/BIMSA/scripts/
python3 gen_input.py -n 1000000 -l 200 -e 0.10
```

Run the file indicating more DPUs and more tasklets using -d and -t respectively.
```
cd $HOME/BIMSA
python3 run_bimsa.py -d 2500 -t 12 -f $HOME/projects/BIMSA/inputs/n1000000_l200_e10.seq -s 250
```

## Executing BIMSA-Hybrid for real datasets
To execute real datasets with heterogeneus alignment sizes optimally. It is recomended to use BIMSA-Hybrid by tunning the parameters `-m` `-b` and using `--dynamic`.
- `-m` Tunes the nominal distance threshold where the alingments are send to the CPU recovery. Unfortunately, the number of alignments which are recovered in CPU can not be indicated, but tunning must be done to the threshold in order to achieve the desired percentage of alignments recovered.
- `-b` Tunes the number of alignments that fit in a batch so the number of batches is determined by the total number of alignments/batch_size. If the batches are used with the CPU recovery, the CPU and DPU kernel executions will overlap for each batch.
- `--dynamic` Activates the dynamic asignment of alignments between tasklets, which improves performance for heterogeneus datasets.

## Advanced BIMSA configurations
The user can configure BIMSA's WRAM sctructure sizes by using the arguments `--wf_trans`, `--seq_trans`, `--cigar_trans` indicating the size in a number from 3 to 10 which refers to a power of 2 for the size in MB.
- `--wf_trans` Controls the wavefront WRAM size.
- `--seq_trans` Controls the number of characters read from the sequence for the WRAM at a given time.
- `--cigar_trans` Controls the number of CIGAR characters stored in WRAM before they are written in MRAM.

If BIMSA arises an error comunicating that the WRAM space is surpassed, lowering `--wf_trans`, `--seq_trans`, `--cigar_trans` or the number of tasklets `-t` will trade performance for WRAM space utillization, which will make the file executable.

If BIMSA arises an error comunicating that the MRAM space is surpassed, only decreassing `-t` will trade performance for MRAM utilization. Additionally, if the number of sequences is enough to feed more DPUs, increasing `-d` will lower MRAM utilization. If MRAM utilization is superior to 100% with `-t 1` and using the maximum numbeber of DPUs, it means that the file is not supported for BIMSA execution.

The arguments `-p`, `-db`, `-bd`, `-ad` are only for debugging and developing pourposes.

## All BIMSA script arguments
```
usage: run_bimsa.py [-h] -f FILE [-t TASKLETS] [-d DPUS] [-wt WF_TRANS]
                    [-st SEQ_TRANS] [-ct CIGAR_TRANS] -s SIZE
                    [-m MAX_DISTANCE] [-b BATCH_SIZE] [-dn] [-p] [-db] [-bd]
                    [-ad]

Python program to run make and execute a C program with an input file.

arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Path to the input file
  -t TASKLETS, --tasklets TASKLETS
                        Number of tasklets (default 12)
  -d DPUS, --dpus DPUS  Number of DPUs (dafault 2500)
  -wt WF_TRANS, --wf_trans WF_TRANS
                        Wavefront transfer and cache size in power of 2
  -st SEQ_TRANS, --seq_trans SEQ_TRANS
                        Sequence and CIGAR transfer and cache size in power of
                        2
  -ct CIGAR_TRANS, --cigar_trans CIGAR_TRANS
                        Sequence and CIGAR transfer and cache size in power of
                        2
  -s SIZE, --size SIZE  Define maximum sequence size
  -m MAX_DISTANCE, --max_distance MAX_DISTANCE
                        Maximum wavefront distance (default 5000)
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Number of pairs to be executed in batches (default 0)
  -dn, --dynamic        Asign pairs dynamically to threads
  -p, --print           Print extra dpu execution information
  -db, --debug          Print extra validation information
  -bd, --banded         Use banded heuristic
  -ad, --adaptive       Use adaptive heuristic
```

## For developers
