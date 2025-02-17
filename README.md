# BarkBeetle_Metabarcoding
Bioinformatic scripts/code for 'Taxonomic and Geographic Diversity of Fungal Symbiont Communities Associated with Tropical Bark and Ambrosia Beetles'

All bash scripts were run on HPC servers using SLURM queueing system, therefore core/RAM/runtimes in .sh scripts are specified in this format. All R scripts were run locally in R v4.3.0 in RStudio v2023.06.0+421.

Order of analyses: 

[1. Sequence Denoising](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#1-sequence-denoising)

[2. Taxonomic identification](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#2-taxonomic-identification)

[3. Fungal Diversity and Community Detection Analysis](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#3-fungal-diversity-and-community-detection-analysis)

## 1. Sequence Denoising
`cd denoising`

The following scripts demultiplex and denoise reads, merge forward and reverse reads, detect consensus chimaeras and cluster reads in operational taxonomic units.
1. `./demultiplex.sh`
2. `./DADA2_denoising_chimaera.sh`
5. `./OTU_clustering.sh`

## 2. Taxonomic identification
`cd taxID`

These scripts assign taxonomy to fungal OTUs using a consensus approach of three methods: 1) the Bayesian classifier algorithm RDP trained on UNITE, 2) local BLASTn searches against UNITE, and 3) phylogenetic placement using TBAS
1. `./RDP.sh`
2. BLASTn and phylogenetic placement were run from the T-BAS GUI webpage (https://tbas.cifr.ncsu.edu/tbas2_3/pages/tbas.php) using the Fungi v3 reference tree and RAxML-EPA algorithm for phylogenetic placement.
3. `Rscript OTU_taxonomy_plot_TBAS_UNITE_ALLCR.R` This script finds consensus between three methods and produces summary plots of the taxonomy
4. `Rscript OTU_taxonomy_plot_TBAS_UNITE_ALLCR_scolytinae_vs_platypodinae.R` As above but also summarising differences between host subfamilies and locality

## 3. Fungal Diversity and Community Detection Analysis
`cd diversity_community_analysis`
The following scripts reconstruct a phylogeny of fungal OTUs, calculate alpha and beta diversity metrics, compare diversity between host subfamilies and localities, identify indicator species, and conduct species co-occurrence analysis

### 3.1 Phylogenetic reconstruction
1. `./fasttree.sh`
2. `sbatch mafft_BeetleOTU.sh`
3. `sbatch trimal.sh`
4. `sbatch iqtree.sh` `sbatch iqtree_trimmed.sh` `sbatch iqtree_trimmed_gt0.2.sh`
5. `Rscript OTU_ALLCR_treeplots.R` compares monophyly of fungal phyla across phylogenetic reconstruction methods. Requires the four tree files (FastTree and 3 IQTrees) and consensus taxonomy files from Section 2

### 3.2 Diversity Analysis
The following scripts were repeated for the four trees to test whether tree reconstruction method impacted results
1. `./calculate_diversity_metrics.sh`
2. `Rscript Multivariate_analysis_ALLCR.R` Beta diversity analysis
3. `Rscript MPD_MNTD.R` Calculate phylogenetic alpha diversity metrics
4. `Rscript Multivariate_analysis_ALLCR_alpha.R` alpha diversity analysis
5. `./calculate_principal_unifrac.sh` calculate unifrac distances for traps instead of samples
6. `Rscript Multivariate_analysis_ALLCR_trap_effect.R` test the effect of traps on beta diversity results

### 3.3 Indicator and Community Detection
1. `Rscript significant_OTUs_Borneo_FG_envfit.R` pulls OTUs that significantly correlate with the PCoA clusters obtained from step 3.2.2
2. `Rscript cooccurrence_analysis.R` requires a modified version of the cooccur package plot.cooccur function to allow for two-tail testing
