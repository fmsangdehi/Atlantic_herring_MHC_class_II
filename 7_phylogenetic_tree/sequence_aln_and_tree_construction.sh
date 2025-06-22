#!/bin/bash

## Make an amino-acid-based ML phylogenetic tree with IQ-TREE

module load bioinfo-tools MAFFT iqtree

##---------------------------
## Perform sequence alignment for alpha and beta chains separately

mafft --reorder MHCII_herring_and_comparative_set_alpha.fa > MHCII_herring_and_comparative_set_alpha_aln.fa

mafft --reorder MHCII_herring_and_comparative_set_beta.fa > MHCII_herring_and_comparative_set_beta_aln.fa


##----------------------------
## Construct tree for alpha

AA_aln="MHCII_herring_and_comparative_set_alpha_aln.fa"

iqtree2 -s ${AA_aln} --seqtype AA --prefix ML_tree_alpha -T AUTO --threads-max 8 -m MFP -B 2000


##----------------------------
## Construct tree for beta

AA_aln="MHCII_herring_and_comparative_set_beta_aln.fa"

iqtree2 -s ${AA_aln} --seqtype AA --prefix ML_tree_beta -T AUTO --threads-max 8 -m MFP -B 2000
