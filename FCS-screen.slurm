#!/bin/bash
#SBATCH --job-name="FCS-OturMall"   #name of this job
#SBATCH -p mem              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "fcs-stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "fcs-stderr.%j.%N"     #optional, prints our standard error
#SBATCH -t 8:00:00

date

module load fcs-genome
module load seqtk/1.3
module load seqkit/2.4.0

TAXID=$1
ASM=$2

#fcs.py screen genome --fasta $ASM --out-dir $ASM-FCSout --tax-id $TAXID

cat $ASM | fcs.py clean genome --action-report $ASM-FCSout/*fcs_gx_report.txt --output $ASM-clean.fa --contam-fasta-out $ASM-contam.fa

seqkit stats -a $ASM > $ASM.stats
seqkit stats -a $ASM-clean.fa > $ASM-clean.stats
seqkit stats -a $ASM-contam.fa > $ASM-contam.stats

# Main targets for web BLAST followup
# get insect-tagged (top hit) contigs larger than 1kb
cut -f1-8 $ASM-FCSout/*.taxonomy.rpt | grep anml:insect | awk '$2 > 999 && $2 != "taxonomy"' | cut -f1 > $ASM-insect.txt
seqtk subseq $ASM $ASM-insect.txt > $ASM-insect.fa

# additional candidate targets for bulk blast on Atlas
# get animal-tagged (any) contigs larger than 5kb, except established Gallus contamination in F1
grep anml $ASM-FCSout/*.taxonomy.rpt | awk '$2 > 4999 && $2 != "taxonomy" && $6 != "Gallus"' | cut -f1 > $ASM-anml_contig_IDs.txt
# get animal-tagged (top hit, not insect) contigs larger than 1kb, except insect or Gallus
cut -f1-8 $ASM-FCSout/*.taxonomy.rpt | grep anml | grep -v insect | grep -v Gallus | awk '$2 > 999 && $2 != "taxonomy"' | cut -f1 >> $ASM-anml_contig_IDs.txt
seqtk subseq $ASM <(sort -u $ASM-anml_contig_IDs.txt) > $ASM-anml_contigs.fa

