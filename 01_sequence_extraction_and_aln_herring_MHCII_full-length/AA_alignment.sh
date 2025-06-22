#!/bin/bash

## Align the AA sequences

module load bioinfo-tools MAFFT

fasta_list=("DAA" "DAB" "DBA" "DBB")

for prefix in "${fasta_list[@]}"; do
    mafft --reorder "AA_${prefix}.fa" > "AA_${prefix}_aln.fa"
done
