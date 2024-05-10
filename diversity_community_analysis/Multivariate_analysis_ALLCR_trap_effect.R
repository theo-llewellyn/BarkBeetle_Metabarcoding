#ADONIS testing the effect of traps on dissimilarity patterns
library(tidyverse)
library(vegan)
library(cowplot)

#read in dissimilarity matrices
BrayCurtis_diss <- read_tsv('CRcoremetrics/bray_curtis_distance_matrix/data/distance-matrix.tsv') %>%
  column_to_rownames("...1")
Jaccard_diss <- read_tsv('CRcoremetrics/jaccard_distance_matrix/data/distance-matrix.tsv') %>%
  column_to_rownames("...1")
UnwUniFrac_diss <- read_tsv('IQtreeCRcoremetrics/unweighted_unifrac/unweighted_unifrac_distance_matrix/distance-matrix.tsv') %>%
  column_to_rownames("...1")
WeUniFrac_diss <- read_tsv('IQtreeCRcoremetrics/weighted_unifrac/weighted_unifrac_distance_matrix/distance-matrix.tsv') %>%
  column_to_rownames("...1")

#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  #remove samples not in dissimilarity matrix
  filter(fungal_metabarcode_ID %in% colnames(BrayCurtis_diss))

#reorder rows to match dissimilarity matrices
sample_data <-sample_data[ order(match(sample_data$fungal_metabarcode_ID, colnames(BrayCurtis_diss))), ]

#check Danum valley samples removed
danum_subset <- subset(sample_data, sample_data$principal_id == 'Danum_Valley_1')
for(i in danum_subset$`Image no`){
  print(paste('Sample',i, sep = ' '))
  print(grep(i,colnames(BrayCurtis_diss)))
}

#original model with no trap ID batch effect
country_fam_model <- adonis2(as.dist(BrayCurtis_diss) ~ Country*FAMILY,
                             data = sample_data)

################################################################################
# TEST 1
# test the effect of subfamily on Borneo and French Guiana separately after 
# accounting for any variation due to principal_ID separate the data
borneo_samples <- subset(sample_data, sample_data$Country == 'Malaysia')

#subset dissimilarity matrices by country
BrayCurtis_diss_borneo <- subset(BrayCurtis_diss[,borneo_samples$fungal_metabarcode_ID], rownames(BrayCurtis_diss) %in% borneo_samples$fungal_metabarcode_ID)
Jaccard_diss_borneo <- subset(Jaccard_diss[,borneo_samples$fungal_metabarcode_ID], rownames(Jaccard_diss) %in% borneo_samples$fungal_metabarcode_ID)
UnwUniFrac_diss_borneo <- subset(UnwUniFrac_diss[,borneo_samples$fungal_metabarcode_ID], rownames(UnwUniFrac_diss) %in% borneo_samples$fungal_metabarcode_ID)
WeUniFrac_diss_borneo <- subset(WeUniFrac_diss[,borneo_samples$fungal_metabarcode_ID], rownames(WeUniFrac_diss) %in% borneo_samples$fungal_metabarcode_ID)

#reorder rows of metadata to match dissimilarity matrices
sample_data_borneo <-borneo_samples[ order(match(borneo_samples$fungal_metabarcode_ID, colnames(BrayCurtis_diss_borneo))), ]

#same for French Guiana
FG_samples <- subset(sample_data, sample_data$Country == 'French Guiana')

#subset dissimilarity matrices by country
BrayCurtis_diss_FG <- subset(BrayCurtis_diss[,FG_samples$fungal_metabarcode_ID], rownames(BrayCurtis_diss) %in% FG_samples$fungal_metabarcode_ID)
Jaccard_diss_FG <- subset(Jaccard_diss[,FG_samples$fungal_metabarcode_ID], rownames(Jaccard_diss) %in% FG_samples$fungal_metabarcode_ID)
UnwUniFrac_diss_FG <- subset(UnwUniFrac_diss[,FG_samples$fungal_metabarcode_ID], rownames(UnwUniFrac_diss) %in% FG_samples$fungal_metabarcode_ID)
WeUniFrac_diss_FG <- subset(WeUniFrac_diss[,FG_samples$fungal_metabarcode_ID], rownames(WeUniFrac_diss) %in% FG_samples$fungal_metabarcode_ID)

#reorder rows to match dissimilarity matrices
sample_data_FG <-FG_samples[ order(match(FG_samples$fungal_metabarcode_ID, colnames(BrayCurtis_diss_FG))), ]

####
#separating by country and see if effect of family significant after accounting for trap
Borneo_model_BrayCurtis <- adonis2(BrayCurtis_diss_borneo ~ principal_id+FAMILY,
                                   strata = sample_data_borneo$principal_id, 
                                   data = sample_data_borneo)
Borneo_model_Jaccard <- adonis2(Jaccard_diss_borneo ~ principal_id+FAMILY,
                                strata = sample_data_borneo$principal_id, 
                                data = sample_data_borneo)
Borneo_model_UnwUniFrac <- adonis2(UnwUniFrac_diss_borneo ~ principal_id+FAMILY,
                                   strata = sample_data_borneo$principal_id, 
                                   data = sample_data_borneo)
Borneo_model_WeUniFrac <- adonis2(WeUniFrac_diss_borneo ~ principal_id+FAMILY,
                                  strata = sample_data_borneo$principal_id, 
                                  data = sample_data_borneo)

FG_model_BrayCurtis <- adonis2(BrayCurtis_diss_FG ~ principal_id + FAMILY,
                               strata = sample_data_FG$principal_id, 
                               data = sample_data_FG)
FG_model_Jaccard <- adonis2(Jaccard_diss_FG ~ principal_id+FAMILY,
                            strata = sample_data_FG$principal_id, 
                            data = sample_data_FG)
FG_model_UnwUniFrac <- adonis2(UnwUniFrac_diss_FG ~ principal_id+FAMILY,
                               strata = sample_data_FG$principal_id, 
                               data = sample_data_FG)
FG_model_WeUniFrac <- adonis2(WeUniFrac_diss_FG ~ principal_id+FAMILY,
                              strata = sample_data_FG$principal_id, 
                              data = sample_data_FG)


################################################################################
# TEST 2
# Collapse principals to test effect of country removing any principal effect

#read in OTU reads counts
taxa <- read_tsv('ALLCRfiltered_feature-table.tsv',skip = 1) %>%
  #transpose
  pivot_longer(-`#OTU ID`) %>% 
  pivot_wider(names_from=`#OTU ID`, values_from=value) %>%
  #remove beetles with zero reads after denoising in qiime
  subset(name %in% sample_data$fungal_metabarcode_ID) %>%
  column_to_rownames(var = 'name')

# merge samples within each trap total read counts across samples within trap 
# and repeat stats
taxa %>%
  rownames_to_column() %>%
  left_join(sample_data[,c(19,28)], by = c('rowname' = 'fungal_metabarcode_ID')) %>%
  group_by(principal_id) %>%
  summarise(across(-1, ~ sum(., na.rm = TRUE))) %>%
  filter(principal_id != '#N/A') -> trap_collapsed_data

#add country metadata
trap_collapsed_metadata <-  trap_collapsed_data %>%
  #add to country metadata
  left_join(sample_data[,c(6,28)]) %>%
  dplyr::select(principal_id,Country) %>%
  dplyr::distinct()

#save a pivoted version to calculate unifrac distances for traps
trap_collapsed_data[,-3821] %>%
  pivot_longer(-principal_id) %>% 
  pivot_wider(names_from=principal_id, values_from=value) -> trap_collapsed_pivoted

write_tsv(trap_collapsed_pivoted,'ALLCRfiltered_feature-table_trap_collapse.tsv')

#run bash script calculate_principal_unifrac.sh

#read in unifrac dist matrices from qiime
UnwUniFrac_diss_trap <- read_tsv('IQtreeCRcoremetrics_trap_collapse/unweighted_unifrac_distance_matrix/distance-matrix.tsv') %>%
  column_to_rownames("...1")
WeUniFrac_diss_trap <- read_tsv('IQtreeCRcoremetrics_trap_collapse/weighted_unifrac_distance_matrix/distance-matrix.tsv') %>%
  column_to_rownames("...1")


trap_collapsed_model_bray <- adonis2(formula = trap_collapsed_data[,-c(1,3821)] ~ Country,
                                     method = 'bray',
                                     data = trap_collapsed_metadata)

trap_collapsed_model_jaccard <- adonis2(formula = trap_collapsed_data[,-c(1,3821)] ~ Country,
                                        method = 'jaccard',
                                        data = trap_collapsed_metadata)

trap_collapsed_model_ununifrac <- adonis2(UnwUniFrac_diss_trap ~ Country,
                                          data = trap_collapsed_metadata)

trap_collapsed_model_weunifrac <- adonis2(WeUniFrac_diss_trap ~ Country,
                                          data = trap_collapsed_metadata)
