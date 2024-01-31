#!/bin/bash
#SBATCH -J iqtree
#SBATCH -p week
#SBATCH -t 7-00:00:00
#SBATCH --mem=100
#SBATCH -c 1

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate iqtree2-env

cd ~/iqtree

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Current working directory is `pwd`"

iqtree -s MSA.MAFFT.aln.With_Names_tAl_gappyout.fa \
 -m GTR+G \
 -B 1000 -bnni \
 --prefix CRotus_3820T_trimmed_guidance \
 -T 1
