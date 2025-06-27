#Alpha diversity multivariate analysis
library(tidyverse)
library(vegan)
library(cowplot)
library(nlme)
library(lme4)
library(flexplot)


setwd("~/OneDrive - Natural History Museum/01_PUBLICATIONS/ITS_Fungi_Metabarcoding/RESULTS/")

#read in data with OTUs and their presence in each sample
OTU_table <- read_tsv('OTUs/OTUsotu_table.txt') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID") %>% 
  t() %>%
  as.data.frame()

#visualise sequencing depth
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x = total)) + 
  geom_histogram(binwidth = 500) +
  scale_y_continuous(labels = scales::comma)

#see where bulk of data sits
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  ggplot(aes(x=1,y=total)) +
  geom_jitter() +
  scale_y_log10()

#sort by the total number of sequences
OTU_table %>%
  mutate(total = rowSums(.,)) %>%
  arrange(total) %>%
  select(total) %>%
  head(50)

#plot rarefaction
rarecurve(OTU_table, step=1000, label=FALSE)
#zoom in
rarecurve(OTU_table, step=1000, label=FALSE, xlim = c(0,10000))

#rarefy to 1000 seqs
rarefied_richness <- OTU_table %>%
  rarefy(., 1000) %>%
  as_tibble(rownames ='name') %>%
  select(name, rarefied_richness = value)

#non-rarefied richness
alpha_stats <- tibble(name = rownames(OTU_table), richness = rowSums(OTU_table > 0)) %>%
  inner_join(rarefied_richness)


#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  #remove samples not in dissimilarity matrix
  filter(fungal_metabarcode_ID %in% alpha_stats$name)

#join sample and alpha data
alpha_stats_metadata <- inner_join(alpha_stats, sample_data, by = c('name' = 'fungal_metabarcode_ID'))

#GLM
#full model allowing random effect to vary both slope and intercept. Negative binomial distribution
richness_glm <- glmer.nb(rarefied_richness ~ FAMILY*Country + (1+FAMILY|principal_id), data = alpha_stats_metadata)

#check if significant
summary(richness_glm)
#Nakagawa & Schielzeth's (2013) r squared
library(r2glmm)
r2beta(richness_glm, method = 'sgv', partial = TRUE)

visualize(richness_glm, plot = 'model', sample = 30) +
  theme(legend.position = "none")

################################################################################
#read phylogenetic diversity indices
MPD_MNTD <- read_csv('MPD_MNTD_results.csv')
#merge the richness and phylogenetic diversity
alpha_stats <- left_join(alpha_stats, MPD_MNTD, by = c('name' = 'rowname'))

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
