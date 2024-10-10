#!/bin/bash
#SBATCH --job-name="subseq"   #name of this job
#SBATCH -p short              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 4                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "seqtk-stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "seqtk-stderr.%j.%N"     #optional, prints our standard error
#SBATCH -t 1:00:00

date

module load seqtk/1.3

ASM=$1
LIST=$2

seqtk subseq $ASM $LIST > $LIST.fa

