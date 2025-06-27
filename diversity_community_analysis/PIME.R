#PIME
library(tidyverse)
library(vegan)
library(phyloseq)
library(pime)


setwd("~/OneDrive - Natural History Museum/01_PUBLICATIONS/ITS_Fungi_Metabarcoding/RESULTS/")

#read in dissimilarity matrices
OTU_table <- read_tsv('OTUs/OTUsotu_table.txt') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID")

#transpose
OTU_table <- OTU_table %>%
  t() %>%
  as.data.frame()

#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  dplyr::select(c(fungal_metabarcode_ID,subfamily,country)) %>% 
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
OTU_table_reduced <- OTU_table[ ! rownames(OTU_table) %in% missing_metadata, ] %>%
  mutate_if(is.numeric, ~1 * (. > 0))


#link data into phyloseq object
taxonomy <- read_tsv('OTUs/OTUstaxonomy.txt') %>%
  select(-c(ReferenceID, rank, score, cutoff, abundance)) %>%
  column_to_rownames('OTU_ID')
ps <- phyloseq(otu_table(OTU_table_reduced, taxa_are_rows=FALSE), 
               tax_table(as.matrix(taxonomy)))
sample_metadata <- phyloseq::sample_data(sample_data %>% column_to_rownames('index'))
physeq <- merge_phyloseq(ps,sample_metadata)



###########################
# Run PIME
#estimate OOB error rate
pime.oob.error(physeq, "country")
# 0.02281369

#split tables by country
per_variable_obj <- pime.split.by.variable(physeq, "country")
per_variable_obj

#estimate highest prevalence, calculates prevalence for taxa
prevalences <- pime.prevalence(per_variable_obj)
prevalences

#get best prevalence value
set.seed(8)
best.prev <- pime.best.prevalence(prevalences, "country")


#likelihood of introducing bias
#randomises sample labels into arbitrary groups and calculates OOB. Used to check whether diff. in original groups due to chance
randomized <- pime.error.prediction(physeq, "country", bootstrap = 100, parallel = FALSE, max.prev = 10)
randomized$Plot
randomized$`Results table`

#replicate without randomisation to show how consistent results are
replicated.oob.error <- pime.oob.replicate(prevalences, "country", bootstrap = 100, parallel = TRUE)
replicated.oob.error$Plot
replicated.oob.error$`Results table`

#show random and real together
png('FIGURES/dynamicOTUs/PIME_prediction_error.png', res = 400, height = 2000, width = 2000)
ggplot(data = randomized$`Results table`, 
       mapping = aes(x='5', y = `Prevalence5%`)) +
  geom_boxplot(aes(col = "Randomised")) +
  geom_boxplot(data = replicated.oob.error$`Results table`, aes(col = 'Original')) +
  ylim(0,1) +
  xlab('Prevalence %') +
  ylab('OOB error') +
  theme_bw() +
  coord_fixed() +
  scale_color_manual(name = "Dataset",
                     values = c("Randomised" = "black", "Original" = "red"))
dev.off()

imp5 <- best.prev$`Importance`

#add column showing whether that ASV contributes to defining Borneo or FG
imp5_table <- imp5$`Prevalence 5`
imp5_table$indicator_country <- NULL

for(i in 1:nrow(imp5_table)){
  if(imp5_table[i,'French.Guiana']>imp5_table[i,'Malaysia']){
    imp5_table[i,'indicator_country'] = 'Borneo'
  } else{
    imp5_table[i,'indicator_country'] = 'French Guiana'
  }
}

write_csv(imp5_table, 'dynamicOTUs/PIME_Results_country.csv')
head(imp5)

