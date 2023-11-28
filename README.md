# WFA-UPMEM
Implementation of a distance-based WFA tailored for the UPMEM architecture.

<p style="text-align:center;"><img src="images/upmem.jpg" alt="upmem" class="center" title="UPMEM" width="250" height="180"/></p>

## Prerequisites
You need to have the UPMEM SDK installed, please follow the instructions here:

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
cd wfa-upmem/tools/generate_dataset/
make
cd ../../scripts/
python3 generate_inputs.py

```
## Running some experiments locally

```
cd wfa-upmem/scripts/
python3 run_upmem_short.py

```
## Checking the results
```
cd wfa-upmem/upmem/profile
cat *
```
You should see the outputs of the different executions.

## Running experiments on the cluster
Larger inputs are available under "/home/upmem0046/bscivanfv/wfa-upmem/inputs".