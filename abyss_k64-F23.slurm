#!/bin/bash
#SBATCH -J abyss-F23
#SBATCH -o "abyss.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "abyss.%j.%N.stderr"
#SBATCH -c 12
#SBATCH -p mem
#SBATCH -N 1
#SBATCH -t 1-00:00:00

date

module load abyss

abyss-pe k=64 B=100G name=OturF_23-k64 lib='f2 f3' \
	f2='O_turi_F2.unmapped.fastq' \
	f3='O_turi_F3.unmapped.fastq'

date
