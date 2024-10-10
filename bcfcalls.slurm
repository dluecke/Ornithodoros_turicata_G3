#!/bin/bash
#SBATCH -J bcfcall_Otur
#SBATCH -o "bcf.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bcf.%j.%N.stderr"
#SBATCH -c 9
#SBATCH -p mem
#SBATCH -N 1
#SBATCH -t 2-00:00:00

module load samtools/1.17
module load bcftools/1.20

BAM=$1
ASSEMBLY="../GCF_037126465.1_ASM3712646v1_genomic.fna"

samtools index -M $BAM

bcftools mpileup --threads 8 -Ou -f $ASSEMBLY $BAM | bcftools call --threads 8 -cv -Ov -o $BAM.vcf

# write csv counting variant lines per scaffold
grep ^N $BAM.vcf | cut -f1 | sort | uniq -c | sed 's/ N/,N/' | tr -d ' ' > $BAM.nseg_by_scaffold.csv

# add lines to csv with NC_088201.1 variant count split at 90Mbp boundary
SPLIT1=$(grep ^NC_088201.1 $BAM.vcf | awk '$2 <= 90000000' | wc -l)
SPLIT2=$(grep ^NC_088201.1 $BAM.vcf | awk '$2 > 90000000' | wc -l)
echo "${SPLIT1},NC_088201.1:1-90M" >> $BAM.nseg_by_scaffold.csv
echo "${SPLIT2},NC_088201.1:90M-end" >> $BAM.nseg_by_scaffold.csv
