#!/bin/bash

# FindBuscoOverlaps.sh counts the shared BUSCO genes per chromosome pair between assemblies

# takes BUSCO full_table.tsv for each assembly (run for same lineage) 
FullTableOtur='full_table-Otur.tsv'
FullTableIsca='full_table-Isca.tsv'
FullTableRmic='full_table-Rmic.tsv'
# and file listing the sequence IDs (matching BUSCO results column 3) for each assembly
ChrIDsOtur='ChrIDs-Otur.txt'
ChrIDsIsca='ChrIDs-Isca.txt'
ChrIDsRmic='ChrIDs-Rmic.txt'
# outputs a CSV with overlapping BUSCOs per chromosome pair (Otur vs Isca/Rmic)
CSV_OUT='BuscoOverlap-Otur_IscaRmic.csv'

# Split full_table outputs by chromosome

# for focal species (Otur) only take "Complete" (not Duplicated)
# Otur
while read ChrID; do 
    awk -v chr="$ChrID" \
        '$2 == "Complete" && $3 == chr' \
    $FullTableOtur > FT-Otur-${ChrID}.tsv
done < $ChrIDsOtur

# for reference species take any hit
# Isca
while read ChrID; do 
    awk -v chr="$ChrID" \
        '$3 == chr' \
    $FullTableIsca > FT-Isca-${ChrID}.tsv
done < $ChrIDsIsca
#Rmic
while read ChrID; do 
    awk -v chr="$ChrID" \
        '$3 == chr' \
    $FullTableRmic > FT-Rmic-${ChrID}.tsv
done < $ChrIDsRmic

# Count overlapping BUSCO IDs for each chromosome pair

# CSV to hold overlap counts
echo "OturChr, RefSpp, RefChr, N_overlap" > $CSV_OUT

# pairwise comparisons of BUSCO ID column
while read OturChr; do 

    while read IscaChr; do

        N_OVERLAP=$(comm --total -123 \
            <(cut -f1 FT-Otur-${OturChr}.tsv | sort) \
            <(cut -f1 FT-Isca-${IscaChr}.tsv | sort -u) \
            | cut -f3)

        echo "${OturChr}, Isca, ${IscaChr}, ${N_OVERLAP}" >> $CSV_OUT

    done < $ChrIDsIsca

    while read RmicChr; do

        N_OVERLAP=$(comm --total -123 \
            <(cut -f1 FT-Otur-${OturChr}.tsv | sort) \
            <(cut -f1 FT-Rmic-${RmicChr}.tsv | sort -u) \
            | cut -f3)

        echo "${OturChr}, Rmic, ${RmicChr}, ${N_OVERLAP}" >> $CSV_OUT

    done < $ChrIDsRmic

done < $ChrIDsOtur
