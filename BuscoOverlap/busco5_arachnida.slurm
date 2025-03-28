#!/bin/bash

#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="busco5_arac"
#SBATCH --mail-user=david.luecke@usda.gov   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output="busco5_arachnida-%j-%N.stdout" # job standard output file (%j replaced by job id)
#SBATCH --error="busco5_arachnida-%j-%N.stderr" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load busco5/5.7.1

assembly=$1

busco -m genome -c 48 \
    -i ${assembly} \
    -o ${assembly%.*}_arachnida \
    -l /project/vpgru/databases/busco/busco_downloads/lineages/arachnida_odb10
