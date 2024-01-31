#------------------------------------------------------------------------------------#
# Indicator species analysis

library(cooccur)
library(tidyverse)
library(reshape2)
library(indicspecies)

#Read in UNITE blastn results
taxonomy_data <- read_csv("UNITE_TBAS_CR.csv") %>%
  select(c(query, consensus_three))

#read in sample data
sample_data <- read_csv('ALLCR_level-5.csv') %>%
  select(c(index,subfamily,locality,country,species,place)) %>%
  filter(country %in% c('Borneo','French Guiana')) %>%
  arrange(index)

#read in data with OTUs and their presence in each sample
OTU_2_Sample <- read_tsv("ALLCRfiltered_feature-table.tsv", skip = 1) %>% 
  #convert to presence absence
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  #merge with taxonomy info
  left_join(.,taxonomy_data, by = c(`#OTU ID` = 'query')) %>%
  #remove OTU_ID_col
  select(-"#OTU ID") %>%
  #make taxa names unique
  mutate(consensus_three = paste(row_number(),consensus_three,sep = "_")) %>%
  #set OTU ID as rownames
  column_to_rownames(var = "consensus_three")

#replace NAs with Unclassified in taxa names and then cut to the lowest identified taxon rank
rownames(OTU_2_Sample) <- rownames(OTU_2_Sample) %>% 
  gsub('NA',"Unclassified",.) %>% 
  gsub('__Unclassified.*','_sp.',.)  %>% 
  str_remove(.,'Fungi__[:alpha:].*__')

#tranpose the table so that beetle samples are rows and fungi are columns and remove Jerrys unmatched samples
transposed_OTU_2_Sample <- OTU_2_Sample %>% 
  rownames_to_column() %>%
  as_tibble() %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "rowname", values_from = "value") %>%
  filter(name %in% sample_data$index) %>%
  arrange(name) %>%
  column_to_rownames(var = 'name')

#check rows are in same order
rownames(transposed_OTU_2_Sample) == sample_data$index

#indicator species analysis
#for country
#testing each species individually
indval <- multipatt(transposed_OTU_2_Sample, sample_data$country, 
                   control = how(nperm=999))
#show only indicators with high specificity and sensitivity
summary(indval, indvalcomp = TRUE, At=0.5, Bt=0.2)

#testing for combinations of indicators for Borneo
Borneo_indic <- indicators(X=transposed_OTU_2_Sample, cluster= sample_data$country, group='Borneo', 
                 max.order = 3, verbose=TRUE, 
                 At=0.5, Bt=0.2,
                 func = 'IndVal')
print(Borneo_indic, sqrtIVt = 0.6)

#calculate the proportion of sites of the target site group where one or another indicator is found
coverage(Borneo_indic, At=0.8, alpha =0.05)

#prune indicators
Borneo_indic2 <- pruneindicators(Borneo_indic, At=0.8, Bt=0.2, verbose=TRUE)
print(Borneo_indic2)
write.table(print(Borneo_indic2), sep = ",", 'Borneo_indicators.csv', col.names=NA)

#testing for combinations of indicators for FG
FG_indic <- indicators(X=transposed_OTU_2_Sample, cluster= sample_data$country, 
                       group='French Guiana',
                       max.order = 3, verbose=TRUE,
                       At=0.5, Bt=0.2)
print(FG_indic, sqrtIVt = 0.6)

#calculate the proportion of sites of the target site group where one or another indicator is found
coverage(FG_indic, At=0.8, alpha =0.05)

#prune indicators
FG_indic2 <- pruneindicators(FG_indic, At=0.8, Bt=0.2, verbose=TRUE)
print(FG_indic2)
write.table(print(FG_indic2), sep = ",", 'FG_indicators.csv', col.names=NA)

#####################
#for subfamily
indval_subfamily <- multipatt(transposed_OTU_2_Sample, sample_data$subfamily, 
                    control = how(nperm=999)) 
summary(indval_subfamily, indvalcomp = TRUE, At=0.5, Bt=0.2)

#testing for combinations of indicators for Scolytinae
Scolytinae_indic <- indicators(X=transposed_OTU_2_Sample, cluster= sample_data$subfamily, group='Scolytinae', 
                           max.order = 3, verbose=TRUE, 
                           At=0.5, Bt=0.2,
                           func = 'IndVal.g')
print(Scolytinae_indic, sqrtIVt = 0.6)

#calculate the proportion of sites of the target site group where one or another indicator is found
coverage(Scolytinae_indic, At=0.8, alpha =0.05)

#prune indicators
Scolytinae_indic2 <- pruneindicators(Scolytinae_indic, At=0.8, Bt=0.2, verbose=TRUE)
print(Scolytinae_indic2)
write.table(print(Scolytinae_indic2), sep = ",", 'Scolytinae_indicators.csv', col.names=NA)

#testing for combinations of indicators for Platypodinae
Platypodinae_indic <- indicators(X=transposed_OTU_2_Sample, cluster= sample_data$subfamily, group='Platypodinae', 
                               max.order = 3, verbose=TRUE, 
                               At=0.5, Bt=0.2,
                               func = 'IndVal.g')
print(Platypodinae_indic, sqrtIVt = 0.6)

#calculate the proportion of sites of the target site group where one or another indicator is found
coverage(Platypodinae_indic, At=0.8, alpha =0.05)

#prune indicators
Platypodinae_indic2 <- pruneindicators(Platypodinae_indic, At=0.8, Bt=0.2, verbose=TRUE)
print(Platypodinae_indic2)
