#Alpha diversity multivariate analysis
library(tidyverse)
library(vegan)
library(cowplot)
library(nlme)
library(lme4)
library(flexplot)

#read in data with OTUs and their presence in each sample
OTU_2_Sample <- read_tsv("ALLCRfiltered_feature-table.tsv", skip = 1) %>% 
  #convert to presence absence
  #transpose
  pivot_longer(Borneo_Sco_P1_A01:FG_Sco_P3_H08) %>% 
  pivot_wider(names_from = `#OTU ID`, values_from = value)

#visualise sequencing depth
OTU_2_Sample %>%
  column_to_rownames(var = "name") %>%
  mutate(total = rowSums(.,)) %>%
ggplot(aes(x = total)) + 
  geom_histogram(binwidth = 500) +
  scale_y_continuous(labels = scales::comma)

#see where bulk of data sits
OTU_2_Sample %>%
  column_to_rownames(var = "name") %>%
  mutate(total = rowSums(.,)) %>%
ggplot(aes(x=1,y=total)) +
  geom_jitter() +
  scale_y_continuous(labels = scales::comma)

#sort by the total number of sequences
OTU_2_Sample %>%
  column_to_rownames(var = "name") %>%
  mutate(total = rowSums(.,)) %>%
  arrange(total) %>%
  select(total) %>%
  head(50)

#rarefy to 500 seqs
rarefied_richness <- OTU_2_Sample %>%
  column_to_rownames(var = "name") %>%
  rarefy(., 1000) %>%
  as_tibble(rownames ='name') %>%
  select(name, rarefied_richness = value)

#non-rarefied richness
richness <- OTU_2_Sample %>%
  column_to_rownames(var = "name") %>%
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  mutate(total = rowSums(.,)) %>%
  select(total) %>%
  as_tibble(rownames ='name') %>%
  select(name, richness = total) %>%
  inner_join(rarefied_richness)


#read phylogenetic diversity indices
MPD_MNTD <- read_csv('MPD_MNTD_results.csv')
#merge the richness and phylogenetic diversity
alpha_stats <- left_join(richness, MPD_MNTD, by = c('name' = 'rowname'))

#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  #remove samples not in dissimilarity matrix
  filter(fungal_metabarcode_ID %in% alpha_stats$name)

#join sample and alpha data
alpha_stats_metadata <- inner_join(alpha_stats, sample_data, by = c('name' = 'fungal_metabarcode_ID'))

#GLM
#full model allowing random effect to vary both slope and intercept. Negative binomial distribution
richness_glm <- glmer.nb(rarefied_richness ~ FAMILY.x*Country + (1+FAMILY.x|principal_id), data = alpha_stats_metadata)

#check if significant
summary(richness_glm)

visualize(richness_glm, plot = 'model', sample = 30) +
  theme(legend.position = "none")

################################################################################
#same with MPD
#full model allowing random effect to vary both slope and intercept. Negative binomial distribution
MPD_glm <- glmer.nb(mpd.obs ~ FAMILY.x*Country + (1+FAMILY.x|principal_id), data = alpha_stats_metadata)

#check if significant
summary(MPD_glm)

################################################################################
#same with MNTD
#full model allowing random effect to vary both slope and intercept. Negative binomial distribution
MNTD_glm <- glmer.nb(mntd.obs ~ FAMILY.x*Country + (1+FAMILY.x|principal_id), data = alpha_stats_metadata)

#check if significant
summary(MNTD_glm)

visualize(MNTD_glm, plot = 'model', sample = 30) +
  theme(legend.position = "none")
