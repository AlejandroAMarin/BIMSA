# BIMSA: Bidirectional In Memory Sequence Aligner
Implementation of a distance-based Bidirectional Wavefront Algorithm tailored for the UPMEM architecture.

## Features

* Calculates **edit distance (Levenshtein distance)** and the **optimal alignment path (CIGAR)**.
* It scales to **long and noisy sequences** as long as the allocated memory does not surpass the **MRAM limit** (64MB).
* Includes a synthetic file generator from [WFA2lib](https://github.com/smarco/WFA2-lib).
* Provides optional **batching**, optional **dynamic thread asignment** and optional **CPU recovery**.
* The user can configure the size for the different memory structures independently, as well as the number of comput units and threads.

Contents
--------

* [Features](#features)
* [Contents](#contents)
* [Setting up the UPMEM local simulator](#setting-up-the-upmem-local-simulator)
* [Getting started with BIMSA](#Getting-started-with-BIMSA)
* [Generating new synthetic inputs](#Generating-new-synthetic-inputs)
* [Running BIMSA on a UPMEM server](#Running-BIMSA-on-a-UPMEM-server)
* [Executing BIMSA-Hybrid for real datasets](#Executing-BIMSA-Hybrid-for-real-datasets)
  * [Reproducing the BIMSA-Hybrid configuration for the paper](#Reproducing-the-BIMSA-Hybrid-configuration-for-the-paper)
* [Advanced BIMSA configurations](#Advanced-BIMSA-configurations)
* [All BIMSA script arguments](#All-BIMSA-script-arguments)
* [For developers](#For-developers)
  * [UPMEM code file layout](#UPMEM-code-file-layout)
* [Access to a UPMEM server](#Access-to-a-UPMEM-server)

## Setting up the UPMEM local simulator
>[!NOTE]
> Python is required if not installed.
```
sudo apt install python3
```

1. Download the UPMEM SDK from the following source: [UPMEM SDK](https://sdk.upmem.com/)

2. Unpack the binary package, for instance into `$HOME/upmem-sdk` directory.
3. Source the script `$HOME/upmem-sdk/upmem_env.sh` to set environment variables to appropriate values.
```
mkdir $HOME/upmem-sdk
tar --strip-components=1 -zxvf $HOME/upmem-2023.2.0-Linux-x86_64.tar.gz -C $HOME/upmem-sdk
source $HOME/upmem-sdk/upmem_env.sh simulator
```

## Getting started with BIMSA
Clone this repository and initialize the submodule:
```
git clone git@github.com:AlejandroAMarin/BIMSA.git
cd BIMSA
git submodule init
git submodule update
```

Run a unit test:
```
python3 run_bimsa.py -d 2 -t 4 -f $PWD/inputs/wfa.utest.seq -s 20000
```

## Generating new synthetic inputs
To generate new synthetic inputs, first compile the tool, second run the script indicating your file sizes and error.
```
cd tools/generate_dataset/
make
cd ../../scripts/
python3 gen_input.py -n 10 -l 100 -e 0.02
```
Run BIMSA on the new created inputs:
```
cd ../
python3 run_bimsa.py -d 2 -t 2 -f $PWD/inputs/n10_l100_e2.seq -s 150
```
> [!IMPORTANT]
> Be aware that the file generator creates files with sequence lengths approximate to the indicated length, so the size `(-s)` on the BIMSA arguments must be at least 5% beyond the dataset sequence size. If the size is inferior to a sequence, this sequence will be discarded.

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
> [!WARNING]
> The UPMEM simulator is only able to run on 8 DPUs max and 24 tasklets.
> The times provided by the simulator are not equivalent to a real hardware execution. The results obtained from the alignment are the same as an UPMEM server.

## Running BIMSA on a UPMEM server
> [!NOTE]
> BIMSA was developed and tested in UPMEM UPMEM-v1A chips using SDK 2024.1.0.
> To check your sdk version run the `dpu-diag` command.

On a UPMEM server all libraries should be installed, so BIMSA should run straightforwardly by cloning the repository. If there are any SDK problems, contact the UPMEM team for troubleshooting.
```
git clone git@github.com:AlejandroAMarin/BIMSA.git
cd BIMSA
git submodule init
git submodule update
```
Create a large enough file to feed all the DPUs
```
cd tools/generate_dataset/
make
cd ../../scripts/
python3 gen_input.py -n 1000000 -l 200 -e 0.10
```

Run the file indicating more DPUs and more tasklets using -d and -t respectively.
```
cd ../
python3 run_bimsa.py -d 2500 -t 12 -f $PWD/inputs/n1000000_l200_e10.seq -s 250
```

## Executing BIMSA-Hybrid for real datasets
To execute real datasets with heterogeneous alignment sizes optimally. It is recommended to use BIMSA-Hybrid by tunning the parameters `-m` `-b` and using `--dynamic`.
- `-m` Tunes the nominal distance threshold where the alingments are sent to the CPU recovery. Unfortunately, the number of alignments which are recovered in CPU can not be indicated, but tunning must be done to the threshold in order to achieve the desired percentage of alignments recovered.
- `-b` Tunes the number of alignments that fit in a batch so the number of batches is determined by the total number of alignments/batch_size. If the batches are used with the CPU recovery, the CPU and DPU kernel executions will overlap for each batch.
- `--dynamic` Activates the dynamic assignment of alignments between tasklets, which improves performance for heterogeneus datasets.

### Reproducing the BIMSA-Hybrid configuration for the paper
>[!WARNING]
> The Illumina, PacBio and Nanopore files are not included in the repositorie due to file size. The datasets and similar ones can be found at [Genome in a bottle](https://github.com/genome-in-a-bottle/giab_data_indexes)

```
python3 run_bimsa.py -s 21709 -d 2500 -t 12 -f $PWD/inputs/Nanopore.bowden.1M.seq  -dn -m 600 -b 200000
```
```
python3 run_bimsa.py -s 18318 -d 2500 -t 12 -f $PWD/inputs/PacBio.CSS.1M.seq -m 50 -dn
```
## Advanced BIMSA configurations
The user can configure BIMSA's WRAM structure sizes by using the arguments `--wf_trans`, `--seq_trans`, `--cigar_trans` indicating the size in a number from 3 to 10 which refers to a power of 2 for the size in MB.
- `--wf_trans` Controls the wavefront WRAM size.
- `--seq_trans` Controls the number of characters read from the sequence for the WRAM at a given time.
- `--cigar_trans` Controls the number of CIGAR characters stored in WRAM before they are written in MRAM.

> [!IMPORTANT]
> If BIMSA raises an error communicating that the WRAM space is surpassed, lowering `--wf_trans`, `--seq_trans`, `--cigar_trans` or the number of tasklets `-t` will trade performance for WRAM space utilization, which will make the file executable.

> [!IMPORTANT]
> If BIMSA raises an error communicating that the MRAM space is surpassed, only decreasing `-t` will trade performance for MRAM utilization. Additionally, if the number of sequences is enough to feed more DPUs, increasing `-d` will lower MRAM utilization. If MRAM utilization is superior to 100% with `-t 1` and using the maximum numbeber of DPUs, it means that the file is not supported for BIMSA execution.

> [!WARNING]
> The arguments `-p`, `-db`, `-bd`, `-ad` are only for debugging and developing purposes. The heuristics code is not maintained.

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
This section is aimed at users who are curious about the structure and basic functionality of a UPMEM application.

UPMEM applications are composed mainly of a host file and a DPU kernel file. Both files are written in regular C using UPMEM library directives.

The DPU kernel file is compiled using the UPMEM SDK command.

The host file is compiled using the regular GCC, but using the UPMEM libraries includes inside the code.

Once both files are compiled, the host file will run and manage the DPU kernel code by using the UPMEM library directives.

The UPMEM programming paradigm is very similar to the CUDA programming paradigm regarding directives and code organization.

### UPMEM code file layout
```text
BIMSA/
├─ upmem/
│  ├─ dpu/
│  │  ├─ bimsa.c  <- Main DPU kernel
│  │  └─ bimsa.h  <- DPU kernel biwfa library functions
│  ├─ host/
│  │  ├─ bimsa_host.c  <- Main host file
│  │  ├─ timer.c
│  │  ├─ wfa_cpu.c  <- CPU recovery
│  │  └─ wfa_cpu.h
│  ├─ support/
│  │  ├─ common.h
│  │  ├─ counterperf.h
│  │  ├─ params.h
│  │  └─ timer.h
│  └─ Makefile
```

## Access to a UPMEM server
There are three options for accessing UPMEM infrastructure currently:
1. Renting cloud access to UPMEM’s servers.
2. Buying UPMEM DIMM modules. For now, UPMEM’s DIMMs are only compatible with the Intel® Server Board S2600WF Product Family.
3. Buying a pre-configured UPMEM server.

For more information contact [UPMEM](https://www.upmem.com/)