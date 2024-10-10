#!/bin/bash
#SBATCH -J bwa_Otur
#SBATCH -o "bwa.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bwa.%j.%N.stderr"
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -N 1

module load bwa_mem2/2.2.1
module load samtools/1.17

ASSEMBLY=../GCF_037126465.1_ASM3712646v1_genomic.fna
INDEX="Otur"
THREADS=32

bwa-mem2 index -p $INDEX $ASSEMBLY

# function for mapping steps to assembly, call with R1, R2, and OUT_PREFIX
Map_GetUnmapped_GetDepths() {
  R1=$1
  R2=$2
  OUT=$3
  bwa-mem2 mem -t $THREADS $INDEX $R1 $R2 | samtools sort -o $OUT.bam
  # flag -f12 for read and mate both unmapped
  samtools fastq -f12 $OUT.bam > $OUT.unmapped.fastq
  samtools depth -aa $OUT.bam > $OUT.depth.tsv
  awk '{ 
          sum[$1] += $3
          count[$1]++
      } 
      END {
          n=asorti(count, indices); 
          for (i=1; i<=n; i++) {
              print indices[i], "\t", sum[indices[i]]/count[indices[i]]
          }
      }' $OUT.depth.tsv > $OUT.depth_by_scaffold.tsv
}

# female read pairs
Map_GetUnmapped_GetDepths ../Illumina_gDNA/female/O_turi_F1_P15_S1_R1_001.fastq.gz ../Illumina_gDNA/female/O_turi_F1_P15_S1_R2_001.fastq.gz O_turi_F1
Map_GetUnmapped_GetDepths ../Illumina_gDNA/female/O_turi_F2_P17_S2_R1_001.fastq.gz ../Illumina_gDNA/female/O_turi_F2_P17_S2_R2_001.fastq.gz O_turi_F2
Map_GetUnmapped_GetDepths ../Illumina_gDNA/female/O_turi_F3_P24_S3_R1_001.fastq.gz ../Illumina_gDNA/female/O_turi_F3_P24_S3_R2_001.fastq.gz O_turi_F3

# male read pairs
Map_GetUnmapped_GetDepths ../Illumina_gDNA/male/O_turi_M1_P26_S4_R1_001.fastq.gz ../Illumina_gDNA/male/O_turi_M1_P26_S4_R2_001.fastq.gz O_turi_M1a
Map_GetUnmapped_GetDepths ../Illumina_gDNA/male/O_turi_M1_P29_S6_R1_001.fastq.gz ../Illumina_gDNA/male/O_turi_M1_P29_S6_R2_001.fastq.gz O_turi_M1b
Map_GetUnmapped_GetDepths ../Illumina_gDNA/male/O_turi_M2_P28_S5_R1_001.fastq.gz ../Illumina_gDNA/male/O_turi_M2_P28_S5_R2_001.fastq.gz O_turi_M2


