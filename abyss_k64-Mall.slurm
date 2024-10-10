#!/bin/bash
#SBATCH -J abyss-Mall
#SBATCH -o "abyss.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "abyss.%j.%N.stderr"
#SBATCH -c 8
#SBATCH -p mem
#SBATCH -N 1
#SBATCH -t 2-00:00:00
date

module load abyss/2.3.5

abyss-pe k=64 B=100G name=OturM_all-k64 lib='m1a m1b m2' \
	m1a='O_turi_M1a.unmapped.fastq' \
	m1b='O_turi_M1b.unmapped.fastq' \
	m2='O_turi_M2.unmapped.fastq'

date
