#!/bin/bash

# RunMerylMerqury-VPGRU.sh runs meryl kmer analysis and merqury reads vs assembly tool

usage() { 
    echo "USAGE: $0 [-h|-o|-l|-t<-m-p>|-k|-c|-g] READS.FASTQ ASSEMBLY.FASTA"
    echo "  -o STRING FCS output directory, default ASSEMBLYNAME"
    echo "  -l FLAG write ReadsOnly files for leftovers"
    echo "  -t FLAG do trio binning hapmer steps, REQUIRES:"
    echo "    -m MATERNAL_READS.FASTA (quote globs)"
    echo "    -p PATERNAL_READS.FASTA (quote globs)"
    echo "  -k INT kmer length, default 21"
    echo "  -c INT cores, default 48"
    echo "  -g PATH to FlyComparativeGenomics git repo, default ~/FlyComparativeGenomics"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage
[[ "$*" == *" -h"* || $1 == "-h" ]] && usage


# last 2 args the reference and query seq files
READS_FASTQ="${@: -2:1}"
ASM_FASTA="${@: -1}"
ASM_FN="$(basename $ASM_FASTA)"
OUT_PREFIX="${ASM_FN%%.f*a}"
# default run parameters
READS_ONLY="" # empty string won't trigger ReadsOnly section
TRIO=""
K_LEN=21
N_THREAD=48
FCG_PATH=~/FlyComparativeGenomics

# get options, including call usage if -h flag
while getopts ":ho:ltm:p:k:c:g:" arg; do
    case $arg in
        o) # name for RunID and output directory, default assembly filename
            OUT_PREFIX="${OPTARG}"
            ;;
        l) # flag for writing leftover ReadsOnly files
            # any string will do
            READS_ONLY="l" 
            ;;
        t) # flag for trio binning hapmer analysis
            TRIO="t"
            ;;
        m) # maternal read file used in trio binning
            MATERNAL="${OPTARG}"
            ;;
        p) # paternal read file used in trio binning
            PATERNAL="${OPTARG}"
            ;;
        k) # kmer length for reads meryl
            K_LEN=${OPTARG}
            ;;
        c) # number of threads for SLURM submission
            N_THREAD=${OPTARG}
            ;;
        g) # path to annotation_tools/ git repo
            FCG_PATH=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

[[ -f $READS_FASTQ ]] || { echo "Can't find reads file ${READS_FASTQ}"; usage; }
[[ -f $ASM_FASTA ]] || { echo "Can't find assembly file ${ASM_FASTA}"; usage; }
[[ -d $FCG_PATH ]] || { echo "Can't find FlyComparativeGenomics repo at $FCG_PATH"; usage; }

if [[ -n $TRIO ]]; then
    find $MATERNAL -quit &> /dev/null || { echo "Can't find reads file ${MATERNAL}"; usage; }
    find $PATERNAL -quit &> /dev/null || { echo "Can't find reads file ${PATERNAL}"; usage; }
fi

echo -e "Running meryl and merqury with k=${K_LEN} on reads:\n ${READS_FASTQ}"
echo -e "and assembly:\n ${ASM_FASTA}\nRun in :\n $PWD"
echo -e "with merqury output prefix ${OUT_PREFIX}"

if [[ -n $READS_ONLY ]]; then
    echo -e "\nUsing meryl difference to generate ReadsOnly files of reads with kmers not in assembly"
fi

if [[ -n $TRIO ]]; then
    echo -e "\nTrio mode. Running separate meryl count on parental reads"
    echo -e "Maternal:"
    ls $MATERNAL
    echo -e "Paternal:"
    ls $PATERNAL
    echo -e "\nThen will run hapmers.sh to generate trio binned merqury plots"
fi

echo -e "\nSubmitting to the short partition with $N_THREAD tasks"
echo "with job name Meryl-${OUT_PREFIX}"

# launch slurm template with proper variables
sbatch --job-name="Meryl-${OUT_PREFIX}" \
    --mail-user="${USER}@usda.gov" \
    -n ${N_THREAD} \
    -o "Meryl-${OUT_PREFIX}.stdout.%j.%N" \
    -e "Meryl-${OUT_PREFIX}.stderr.%j.%N" \
    --export=ALL,READS_FASTQ=${READS_FASTQ},ASM_FASTA=${ASM_FASTA},K_LEN=${K_LEN},\
MERQURY_OUT=${OUT_PREFIX},READS_ONLY=${READS_ONLY},THREADS=${N_THREAD},FCG_REPO=${FCG_PATH},\
TRIO=${TRIO},MATERNAL="${MATERNAL}",PATERNAL="${PATERNAL}" \
    ${FCG_PATH}/VPGRU-meryl_merqury_TEMPLATE.slurm

