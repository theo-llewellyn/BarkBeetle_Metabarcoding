#TESTING TRAP EFFECT
library(tidyverse)
library(vegan)
library(cowplot)

setwd("~/OneDrive - Natural History Museum/01_PUBLICATIONS/ITS_Fungi_Metabarcoding/RESULTS/")

#read in dissimilarity matrices
OTU_table <- read_tsv('OTUs/OTUsotu_table.txt') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID")

#remove OTUs with less than 10 reads in total
OTU_table <- OTU_table[rowSums(OTU_table) >= 10, ]
#remove OTUs in less than 3 samples
OTU_table <- OTU_table[rowSums(OTU_table > 0) >= 3, ]

#transpose
OTU_table <- OTU_table %>%
  t() %>%
  as.data.frame()

#hellinger transformation
OTU_table <- decostand(OTU_table, method = "hellinger")

#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  rename(index = fungal_metabarcode_ID) %>%
  replace_na(list(subfamily = 'Platypodinae')) %>%
  filter(index %in% rownames(OTU_table)) %>%
  arrange(index)

#Danum Valley to remove
Danum_Valley <- c('FG_Sco_P2_A09','FG_Sco_P2_A10','FG_Sco_P2_B09','FG_Sco_P2_B10','FG_Sco_P2_C09','FG_Sco_P2_C10','FG_Sco_P2_D09','FG_Sco_P2_D10','FG_Sco_P2_E09','FG_Sco_P2_E10','FG_Sco_P2_F09','FG_Sco_P2_F10','FG_Sco_P2_G09','FG_Sco_P2_H08','FG_Sco_P2_H09')
sample_data <- sample_data[ ! sample_data$index %in% Danum_Valley, ] %>%
  arrange(index)

#find out which samples dont have metadata for and remove them
missing_metadata <- setdiff(rownames(OTU_table), sample_data$index)
OTU_table_reduced <- OTU_table[ ! rownames(OTU_table) %in% missing_metadata, ]

#Calculate UniFrac indices
#read tree
iqtree <- read.tree('TREES/FastTree_OTUs.tree')
iqtree$tip.label <- gsub('_',':',iqtree$tip.label)
Physeq <- phyloseq(otu_table(OTU_table_reduced, taxa_are_rows=FALSE),phy_tree(iqtree))
UnwUniFrac_diss <- UniFrac(Physeq, weighted = FALSE, parallel = TRUE)
WeUniFrac_diss <- UniFrac(Physeq, weighted = TRUE, parallel = TRUE)

###############################################################
#test the effect of subfamily on Borneo and FG separately after accounting for any variation due to principal_ID
#separate the data
borneo_samples <- subset(sample_data, sample_data$country == 'Malaysia')
#subset dissimilarity matrices by country
borneo_subset <- OTU_table_reduced[borneo_samples$index,]
#reorder rows of metadata to match dissimilarity matrices
borneo_samples <-borneo_samples[ order(match(borneo_samples$index, rownames(borneo_subset))), ]

#same for FG
FG_samples <- subset(sample_data, sample_data$country == 'French Guiana')
#subset dissimilarity matrices by country
FG_subset <- OTU_table_reduced[FG_samples$index,]
#reorder rows to match dissimilarity matrices
FG_samples <-FG_samples[ order(match(FG_samples$index, rownames(OTU_table_reduced))), ]

UnwUniFrac_diss_borneo <- as.matrix(UnwUniFrac_diss)[borneo_samples$index,borneo_samples$index]
WeUniFrac_diss_borneo <- as.matrix(WeUniFrac_diss)[borneo_samples$index,borneo_samples$index]
UnwUniFrac_diss_FG <- as.matrix(UnwUniFrac_diss)[FG_samples$index,FG_samples$index]
WeUniFrac_diss_FG <- as.matrix(WeUniFrac_diss)[FG_samples$index,FG_samples$index]

####
#separating by country and see if effect significant with family*trap
Borneo_model_BrayCurtis <- adonis2(borneo_subset ~ principal_id + subfamily,
                                   strata = borneo_samples$principal_id,
                                   method = "bray", data = borneo_samples)
Borneo_model_Jaccard <- adonis2(borneo_subset ~ principal_id + subfamily,
                                strata = borneo_samples$principal_id,
                                method = "jaccard", data = borneo_samples)

Borneo_model_UnwUniFrac <- adonis2(UnwUniFrac_diss_borneo ~ principal_id+FAMILY,
                                   strata = borneo_samples$principal_id, 
                                   data = borneo_samples)
Borneo_model_WeUniFrac <- adonis2(WeUniFrac_diss_borneo ~ principal_id+FAMILY,
                                  strata = borneo_samples$principal_id, 
                                  data = borneo_samples)

FG_model_BrayCurtis <- adonis2(FG_subset ~ principal_id + subfamily,
                                   strata = FG_samples$principal_id,
                                   method = "bray", data = FG_samples)
FG_model_Jaccard <- adonis2(FG_subset ~ principal_id + subfamily,
                                strata = FG_samples$principal_id,
                                method = "jaccard", data = FG_samples)

FG_model_UnwUniFrac <- adonis2(UnwUniFrac_diss_FG ~ principal_id+FAMILY,
                               strata = FG_samples$principal_id, 
                               data = FG_samples)
FG_model_WeUniFrac <- adonis2(WeUniFrac_diss_FG ~ principal_id+FAMILY,
                              strata = FG_samples$principal_id, 
                              data = FG_samples)


###############################################################
#Collapse principals to test effect of country removing any principal effect
#merge samples within each trap total read counts across samples within trap and repeat stats
OTU_table_reduced %>%
  rownames_to_column() %>%
  left_join(sample_data[,c(19,28)], by = c('rowname' = 'index')) %>%
  group_by(principal_id) %>%
  summarise(across(-1, ~ sum(., na.rm = TRUE))) %>%
  filter(principal_id != '#N/A') -> trap_collapsed_data

#add country metadata
trap_collapsed_metadata <-  trap_collapsed_data %>%
  #add to country metadata
  left_join(sample_data[,c(6,28)]) %>%
  dplyr::select(principal_id,Country) %>%
  dplyr::distinct()

write_tsv(trap_collapsed_pivoted,'dynamicOTUs/OTUsotu_table_trap_collapse.tsv')

#collapse Unifrac matrices
Physeq_trap <- phyloseq(otu_table(trap_collapsed_data %>% column_to_rownames('principal_id'), taxa_are_rows=FALSE),phy_tree(iqtree))
UnwUniFrac_diss_trap <- UniFrac(Physeq_trap, weighted = FALSE, parallel = TRUE)
WeUniFrac_diss_trap <- UniFrac(Physeq_trap, weighted = TRUE, parallel = TRUE)


trap_collapsed_model_bray <- adonis2(formula = trap_collapsed_data[,-1] ~ Country,
                                     method = 'bray',
                                     data = trap_collapsed_metadata)

trap_collapsed_model_jaccard <- adonis2(formula = trap_collapsed_data[,-1] ~ Country,
                                        method = 'jaccard',
                                        data = trap_collapsed_metadata)

trap_collapsed_model_ununifrac <- adonis2(UnwUniFrac_diss_trap ~ Country,
                                          data = trap_collapsed_metadata)

trap_collapsed_model_weunifrac <- adonis2(WeUniFrac_diss_trap ~ Country,
                                          data = trap_collapsed_metadata)
