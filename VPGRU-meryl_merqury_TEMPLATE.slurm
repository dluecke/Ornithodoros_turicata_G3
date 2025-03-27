#!/bin/bash
#SBATCH --job-name="merqury"   #name of this job
#SBATCH -p short              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                 #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "stderr.%j.%N"     #optional, prints our standard error

date

module load merqury
which meryl
which merqury.sh

# meryl directory names are seq file names with extension replaced by .meryl
# matches default meryl behavior 
READS_FN=$(basename $READS_FASTQ)
READS_MERYL="${READS_FN%.f*}.meryl"
ASM_FN=$(basename $ASM_FASTA)
ASM_MERYL="${ASM_FN%.f*a}.meryl"

echo "Meryl and Merqury jobs ran on -p short with -n ${THREADS}"
echo

echo -e "Running meryl count on reads:\n ${READS_FASTQ}"
echo -e "with kmer length $K_LEN\n"
meryl k=${K_LEN} count ${READS_FASTQ} output ${READS_MERYL}

echo -e "Running merqury comparing reads meryl output:\n ${READS_MERYL}"
echo -e "to assembly file:\n ${ASM_FASTA}"
echo -e "with output prefix ${MERQURY_OUT}\n"
merqury.sh ${READS_MERYL} ${ASM_FASTA} ${MERQURY_OUT}

if [[ -n $READS_ONLY ]]; then
    # ReadsOnly outputs names for fastq file name
    READSONLY_MERYL="${READS_FN%.f*}.ReadsOnly.meryl"
    READSONLY_FASTQ="${READS_FN%.f*}.ReadsOnly.fastq.gz"
    READSONLY_NAMES="${READS_FN%.f*}.ReadsOnly.names.txt"

    echo -e "Running meryl difference to find kmers observed only in hifi reads:\n ${READSONLY_MERYL}\n"
    meryl difference ${READS_MERYL} ${ASM_MERYL} output ${READSONLY_MERYL}

    echo -e "Running meryl-lookup to extract fastq.gz file of reads with kmers not observed in assembly:\n ${READSONLY_FASTQ}\n"
    meryl-lookup -include -sequence ${READS_FASTQ} -mers ${READSONLY_MERYL} | \
        pigz -c > ${READSONLY_FASTQ}

    echo -e "Extracting read names from ${READSONLY_FASTQ}:\n ${READSONLY_NAMES}\n"
    zcat ${READSONLY_FASTQ} | grep -e "^@" > ${READSONLY_NAMES}
fi

if [[ -n $TRIO ]]; then
    READS_MAT="${READS_FN%.f*}_mat.meryl"
    READS_PAT="${READS_FN%.f*}_pat.meryl"
    echo -e "\nTrio mode. Running meryl count on parental reads with kmer length ${K_LEN}."
    echo -e "Maternal:"
    ls $MATERNAL
    echo -e "output prefix ${READS_MAT}\n"
    meryl k=${K_LEN} count ${MATERNAL} output "${READS_MAT}" 
    echo -e "Paternal:"
    ls $PATERNAL
    echo -e "output prefix ${READS_PAT}\n"
    meryl k=${K_LEN} count ${PATERNAL} output "${READS_PAT}" 

    echo -e "Running merqury hapmers.sh on maternal, paternal, and F1 meryl dbs\n"
    $MERQURY/trio/hapmers.sh ${READS_MAT} ${READS_PAT} ${READS_MERYL}
fi

date
