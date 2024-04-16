# BIMSA: Bidirectional In Memory Sequence Aligner
Implementation of a distance-based Bidirectional Wavefront Algorithm tailored for the UPMEM architecture.

## Prerequisites
You need to have the UPMEM SDK installed or access to a UPMEM, please follow the instructions here:

https://sdk.upmem.com/

And remember to source the enviromental variables:
```
source path-to-upmem-sdk/upmem_env.sh
```

## Getting Started
Clone this repository:
```
git clone https://gitlab.bsc.es/ifernand/wfa-upmem.git
```

## Generating some example inputs

```
cd tools/generate_dataset/
make
cd ../../scripts/
python3 generate_inputs.py
```
## Running BIMSA
BIMSA provides a python script to run with the different configurations.

Have in mind that some configurations might saturate the resources, this script does not limit parameters if they surpass the simulator or server capabilities.

```
usage: run_ulsapim.py [-h] -f FILE [-t TASKLETS] [-d DPUS] [-wt WF_TRANS]
                      [-st SEQ_TRANS] [-ct CIGAR_TRANS] -s SIZE
                      [-m MAX_DISTANCE] [-b BATCH_SIZE] [-dn] [-p] [-db] [-bd]
                      [-ad]

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