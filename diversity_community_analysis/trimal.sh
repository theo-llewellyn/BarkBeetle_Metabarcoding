#!/bin/bash
#SBATCH -J trimal
#SBATCH -p hour
#SBATCH -t 01:00:00
#SBATCH --mem=1000
#SBATCH -c 1

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate trimAl-env

cd ~/guidance2/CRotuclusters_Guidance2

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Current working directory is `pwd`"

trimal -in MCRotuclusters_dna-sequences_msa.fa -out MSA.MAFFT.aln.With_Names_tAl_gappyout.fa -fasta -gappyout
trimal -in MCRotuclusters_dna-sequences_msa.fa -out MSA.MAFFT.aln.With_Names_tAl_gt0.2.fa -fasta -gt 0.2
