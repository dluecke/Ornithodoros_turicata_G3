#!/bin/bash
#SBATCH --job-name="blast-asm"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p short
#SBATCH -o "blast.stdout.%j.%N"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "blast.stderr.%j.%N"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails
#SBATCH -t 4:00:00

module load blast+/2.15.0

QUERY=$1
TARGET=$2

#makeblastdb -dbtype nucl -in $TARGET
blastn -num_threads 32 -db $TARGET -query $QUERY -evalue 1 -outfmt 6 -out ${QUERY}_on_${TARGET}-blast.out6

awk '$NF > 99' ${QUERY}_on_${TARGET}-blast.out6 > ${QUERY}_on_${TARGET}-blast_tophits.out6
