#!/bin/bash
#SBATCH -J mafft
#SBATCH -p hour
#SBATCH -t 01:00:00
#SBATCH --mem=1000
#SBATCH -c 8

module load mafft/7.271

#OTUs
mafft \
 --thread 4 \
--auto ../OTUs/OTUssequences.fasta \
> OTUs.msa.fa
