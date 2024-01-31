#!/bin/bash
#SBATCH -J mafft
#SBATCH -p hour
#SBATCH -t 01:00:00
#SBATCH --mem=1000
#SBATCH -c 8

module load mafft/7.271

mafft --thread 8 --auto CRotuclusters_dna-sequences.fasta > CRotuclusters_dna-sequences_msa.fa
