# BarkBeetle_Metabarcoding
Bioinformatic scripts/code for 'Taxonomic and Geographic Diversity of Fungal Symbiont Communities Associated with Tropical Bark and Ambrosia Beetles'

All bash scripts were run on HPC servers using SLURM queueing system, therefore core/RAM/runtimes in .sh scripts are specified in this format. All R scripts were run locally in R v4.3.0 in RStudio v2023.06.0+421.

Order of analyses: 

[1. Sequence Denoising](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#1-sequence-denoising)

[2. Taxonomic identification](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#2-taxonomic-identification)

[3. Fungal Diversity and Community Detection Analysis](https://github.com/theo-llewellyn/BarkBeetle_Metabarcoding#3-fungal-diversity-and-community-detection-analysis)

## 1. Sequence Denoising
`cd denoising`

The following scripts demultiplex and denoise reads, merge forward and reverse reads, detect consensus chimaeras and cluster reads into ASVs.
1. `./demultiplex.sh`
2. `./vsearch.sh`

## 2. Taxonomic identification
`cd taxID`

These scripts assign taxonomy to fungal ASVs using taxon-specific thresholds of ITS2 clustering. ASVs are then clustered into OTUs dynamically.
1. `./dnabarcoder.sh`
2. `./dynamicclustering.R`

## 3. Fungal Diversity and Community Detection Analysis
`cd diversity_community_analysis`
The following scripts calculates alpha and beta diversity metrics, compare diversity between host subfamilies and localities, identify indicator species, and conduct species co-occurrence analysis.

### 3.1 Phylogenetic reconstruction
1. `sbatch mafft_BeetleOTU.sh`
2. `./fasttree.sh`

### 3.2 Diversity Analysis
1. `Rscript alpha_div_2025.R` alpha diversity calculation and analysis
2. `Rscript MPD_MNTD.R` Calculate phylogenetic alpha diversity metrics
3. `Rscript Multivariate_analysis.R` Beta diversity analysis
4. `Rscript trap_effect.R` test the effect of traps on beta diversity results

### 3.3 Indicator and Community Detection
1. `Rscript significant_OTUs_Borneo_FG_envfit.R` pulls OTUs that significantly correlate with the PCoA clusters obtained from step 3.2.2
2. `Rscript cooccurrence_analysis.R` requires a modified version of the cooccur package plot.cooccur function to allow for two-tail testing
3. `Rscript FUNGUILD_cooccurrences.R` assesses and visualised functional ecology of OTUs from co-occurrence analysis
4. `Rscript PIME.R` uses prevalence filtering and Random Forests to detect core OTUs for countries
