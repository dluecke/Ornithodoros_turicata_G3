#!/bin/bash
#SBATCH -J abyss-F1
#SBATCH -o "abyss.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "abyss.%j.%N.stderr"
#SBATCH -c 8
#SBATCH -p mem
#SBATCH -N 1

date

module load abyss

abyss-pe k=96 B=100G name=OturF_1 in='O_turi_F1.unmapped.fastq'

date
