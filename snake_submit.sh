#!/bin/bash

#SBATCH --job-name=synleth

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jlm.nguyen@mail.utoronto.ca

source /cluster/home/julian/miniconda3/etc/profile.d/conda.sh
conda activate synleth

snakemake \
  -s workflow/Snakefile \
  --latency-wait 60 \
  --keep-going \
  --cluster "sbatch \
      -p {resources.partition} \
      --mem={resources.mem} \
      -t {resources.time} \
      -c {resources.cpus} \
      -o logs/{rule}.{wildcards}.out \
      -J {rule}" \
  -j 40 \
  -F
