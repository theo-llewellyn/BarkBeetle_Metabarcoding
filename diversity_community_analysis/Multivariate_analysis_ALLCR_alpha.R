#Alpha diversity multivariate analysis
library(tidyverse)
library(vegan)
library(cowplot)

#read in alpha diversity stats
faiths <- read_tsv('IQtreeCRcoremetrics/faith_pd_vector/alpha-diversity.tsv')
shannon <- read_tsv('IQtreeCRcoremetrics/shannon_vector/alpha-diversity.tsv')
#merge the two
alpha_stats <- left_join(faiths, shannon, by = c('#SampleID' = '...1'))

#read in sample data
sample_data <- read_csv('Sco_Pla_FG_Borneo_metadata.csv') %>%
  #remove samples not in dissimilarity matrix
  filter(fungal_metabarcode_ID %in% alpha_stats$`#SampleID`)

#join sample and alpha data
alpha_stats_metadata <- inner_join(alpha_stats, sample_data, by = c('#SampleID' = 'fungal_metabarcode_ID'))

#check Danum valley samples removed
subset(alpha_stats_metadata, alpha_stats_metadata$principal_id == 'Danum_Valley_1')

#visualise differences between principals
ggplot(alpha_stats_metadata, aes(x = principal_id, y = faith_pd, colour = Country)) + 
  geom_boxplot()+ 
  theme_classic()

#GLM
#full model allowing random effect to vary both slope and intercept
#Faiths PD
faiths_glm_full <- glmer(faith_pd ~ FAMILY*Country + (1|principal_id),
                         family = Gamma,
                         data = alpha_stats_metadata)
summary(faiths_glm_full)

#shannon
shannon_glm_full <- glmer(shannon_entropy ~ FAMILY*Country + (1|principal_id),
                          family = Gamma,
                          data = alpha_stats_metadata)
summary(shannon_glm_full)

# no significant terms
