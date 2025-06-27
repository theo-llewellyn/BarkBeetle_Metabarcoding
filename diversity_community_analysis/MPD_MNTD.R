#SES MPD MNTD

setwd("~/OneDrive - Natural History Museum/01_PUBLICATIONS/ITS_Fungi_Metabarcoding/RESULTS/")

library(tidyverse)
library(reshape2)
library(ape)
library(picante)
library(cowplot)

#read in dissimilarity matrices
OTU_table <- read_tsv('OTUs/OTUsotu_table.txt') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID") %>%
  t() %>%
  as.data.frame()

#read in sample data
sample_data <- read_csv('../NCBI_submission/Sco_Pla_FG_Borneo_metadata.csv') %>%
  dplyr::select(c(fungal_metabarcode_ID,subfamily,country)) %>% 
  rename(index = fungal_metabarcode_ID) %>%
  replace_na(list(subfamily = 'Platypodinae')) %>%
  filter(index %in% rownames(OTU_table)) %>%
  arrange(index)


#read tree
iqtree <- read.tree('TREES/FastTree_OTUs.tree')
iqtree$tip.label <- gsub('_',':',iqtree$tip.label)
#convert tree to phylogenetic distance matrix
phydist <- cophenetic(iqtree)

#calculate MPD and MNTD. Compare phylogenetic diversity to null model
ses.mpd.result <- ses.mpd(OTU_table, phydist, null.model="taxa.labels",
                          abundance.weighted=FALSE, runs=999)
ses.mntd.result <- ses.mntd(OTU_table, phydist, null.model="taxa.labels",
                            abundance.weighted=FALSE, runs=999)
ses.mpd.result
ses.mntd.result

#merge phylogenetic diversity with sample info
PD_sample_metadata <- ses.mpd.result %>%
  rownames_to_column() %>%
  #merge with metadata
  left_join(.,sample_data, by = c('rowname' = 'fungal_metabarcode_ID')) %>% 
  subset(., !is.na(country))

MNTD_sample_metadata <- ses.mntd.result %>%
  rownames_to_column() %>%
  #merge with metadata
  left_join(.,sample_data, by = c('rowname' = 'fungal_metabarcode_ID')) %>% 
  subset(., !is.na(country))

PD_country <- ggplot(PD_sample_metadata, aes(x = country, y = mpd.obs.z, colour = country)) +
  scale_x_discrete(labels=c("Malaysia" = "Borneo")) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width=0.1,show.legend = FALSE)

PD_subfam <- ggplot(PD_sample_metadata, aes(x = FAMILY, y = mpd.obs.z, colour = FAMILY)) +
  scale_x_discrete(labels=c("Curculionidae" = "Scolytinae", 'Platypodidae' = 'Platypodinae')) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width=0.1,show.legend = FALSE)


MNTD_country <- ggplot(MNTD_sample_metadata, aes(x = country, y = mntd.obs.z, colour = country)) +
  scale_x_discrete(labels=c("Malaysia" = "Borneo")) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width=0.1,show.legend = FALSE)

MNTD_subfam <- ggplot(MNTD_sample_metadata, aes(x = FAMILY, y = mntd.obs.z, colour = FAMILY)) +
  scale_x_discrete(labels=c("Curculionidae" = "Scolytinae", 'Platypodidae' = 'Platypodinae')) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width=0.1,show.legend = FALSE)

png("FIGURES/MPD_MNTD.png", res = 400, height = 2480, width = 2480)
plot_grid(PD_country, MNTD_country, 
          PD_subfam,MNTD_subfam,
          ncol = 2)
dev.off()

#classify into clustered, dispersed and random
#join mpd and mntd
MNTD_PD_sample_metadata <- left_join(PD_sample_metadata, MNTD_sample_metadata, by = "rowname", suffix = c("", ".annoying_duplicate_column")) %>%
  select(-ends_with(".annoying_duplicate_column")) %>%
  drop_na(mpd.obs.p)

#save the results
write_csv(MNTD_PD_sample_metadata, 'MPD_MNTD_results.csv')
#assign distribution
#Positive SES.MPD/SES.MNTD values and high quantiles (P>0.95) indicate phylogenetic dispersion, while negative and low quantiles (P<0.05) indicate phylogenetic clustering. Non-significant values indicate phylogenetic randomization.
#PD
for(i in 1:nrow(MNTD_PD_sample_metadata)){
  #PD random = if pval between 0.05 and 0.95
  if(MNTD_PD_sample_metadata[i,"mpd.obs.p"]>0.05 & MNTD_PD_sample_metadata[i,'mpd.obs.p']<0.95){
    MNTD_PD_sample_metadata[i,19] <- "Random"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(MNTD_PD_sample_metadata[i,'mpd.obs.z']<0 & MNTD_PD_sample_metadata[i,'mpd.obs.p']<=0.05){
    MNTD_PD_sample_metadata[i,19] <- "Clustered"
  }
  else if(MNTD_PD_sample_metadata[i,'mpd.obs.z']>0 & MNTD_PD_sample_metadata[i,'mpd.obs.p']>=0.95){
    MNTD_PD_sample_metadata[i,19] <- "Dispersed"
  }
}
#MNTD
for(i in 1:nrow(MNTD_PD_sample_metadata)){
  #PD random = if pval between 0.05 and 0.95
  if(MNTD_PD_sample_metadata[i,'mntd.obs.p']>0.05 & MNTD_PD_sample_metadata[i,'mntd.obs.p']<0.95){
    MNTD_PD_sample_metadata[i,20] <- "Random"
  }
  #Random = if D is greater than 0, pval is NOT dif to 1 and pval is dif to 0
  else if(MNTD_PD_sample_metadata[i,'mntd.obs.z']<0 & MNTD_PD_sample_metadata[i,'mntd.obs.p']<=0.05){
    MNTD_PD_sample_metadata[i,20] <- "Clustered"
  }
  else if(MNTD_PD_sample_metadata[i,'mntd.obs.z']>0 & MNTD_PD_sample_metadata[i,'mntd.obs.p']>=0.95){
    MNTD_PD_sample_metadata[i,20] <- "Dispersed"
  }
}


PD_dotplot_country <- ggplot(MNTD_PD_sample_metadata, aes(y = mpd.obs.z, colour = country)) +
  geom_boxplot( aes(x = country), width = 0.1, outlier.shape = NA) +
  geom_jitter(aes(shape = V19, x = country), alpha = 0.5, width = 0.4, size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_x_discrete(labels=c("Malaysia" = "Borneo")) +
  scale_colour_manual(values = c("#71B379","#B25690"), guide = "none")
MNTD_dotplot_country <- ggplot(MNTD_PD_sample_metadata, aes(y = mntd.obs.z, colour = country)) +
  geom_boxplot( aes(x = country), width = 0.1, outlier.shape = NA) +
  geom_jitter(aes(shape = V20, x = country), alpha = 0.5, width = 0.4, size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_x_discrete(labels=c("Malaysia" = "Borneo")) +
  scale_colour_manual(values = c("#71B379","#B25690"), guide = "none")
PD_dotplot_subfam <- ggplot(MNTD_PD_sample_metadata, aes(y = mpd.obs.z, colour = FAMILY)) +
  geom_boxplot( aes(x = FAMILY), width = 0.1, outlier.shape = NA) +
  geom_jitter(aes(shape = V19, x = FAMILY), alpha = 0.5, width = 0.4, size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_x_discrete(labels=c("Curculionidae" = "Scolytinae", 'Platypodidae' = 'Platypodinae')) +
  scale_colour_manual(values = c("#E69F00","#56B4E9"), guide = "none")
MNTD_dotplot_subfam <- ggplot(MNTD_PD_sample_metadata, aes(y = mntd.obs.z, colour = FAMILY)) +
  geom_boxplot( aes(x = FAMILY), width = 0.1, outlier.shape = NA) +
  geom_jitter(aes(shape = V20, x = FAMILY), alpha = 0.5, width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  scale_x_discrete(labels=c("Curculionidae" = "Scolytinae", 'Platypodidae' = 'Platypodinae')) +
  scale_colour_manual(values = c("#E69F00","#56B4E9"), guide = "none")

png("FIGURES/MPD_MNTD_dots.png", res = 400, height = 2480, width = 2480)
plot_grid(PD_dotplot_country, MNTD_dotplot_country, 
          PD_dotplot_subfam,MNTD_dotplot_subfam,
          ncol = 2)
dev.off()


ggplot(MNTD_PD_sample_metadata,aes(x = 1, fill = ...19)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Clustered","Dispersed","Random")) +
  ylab("proportion") 

ggplot(MNTD_PD_sample_metadata,aes(x = country, fill = ...19)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Clustered","Dispersed","Random")) +
  ylab("MPD proportion") -> country_plot

ggplot(MNTD_PD_sample_metadata,aes(x = FAMILY, fill = ...19)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Clustered","Dispersed","Random")) +
  ylab("MPD proportion") -> family_plot

ggplot(MNTD_PD_sample_metadata,aes(x = country, fill = ...20)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Clustered","Dispersed","Random")) +
  ylab("MNTD proportion") -> country_plot_MNTD

ggplot(MNTD_PD_sample_metadata,aes(x = FAMILY, fill = ...20)) + 
  geom_bar(position = "fill") +
  geom_text(aes(label = ..count..), stat = "count", position = position_fill(.5)) +
  scale_fill_manual(values=c('#9f7441',"red", "#327ec1"), name = "Distribution", labels = c("Clustered","Dispersed","Random")) +
  ylab("MNTD proportion") -> family_plot_MNTD

#get legend
legend_subfam <- get_legend(
  # create some space to the left of the legend
  family_plot + theme(legend.box.margin = margin(0, 0, 0, 12), legend.direction = "horizontal")
)

barplots <- plot_grid(country_plot + theme_minimal() + theme(legend.position="none",
                                                             axis.text.x = element_blank(),
                                                             axis.title.x = element_blank()), 
                      family_plot + theme_minimal() + theme(legend.position="none",
                                                            axis.text.x = element_blank(),
                                                            axis.title.x = element_blank(),
                                                            axis.text.y = element_blank(),
                                                            axis.title.y = element_blank()),
                      country_plot_MNTD + theme_minimal() + theme(legend.position="none"), 
                      family_plot_MNTD + theme_minimal() + theme(legend.position="none",
                                                                 axis.text.y = element_blank(),
                                                                 axis.title.y = element_blank()) +
                        xlab('family')
)

png("FIGURES/MPD_MNTD_barcharts.png", res = 400, height = 2480*0.811, width = 2480)
plot_grid(legend_subfam, barplots, nrow = 2, rel_heights = c(0.1,1))
dev.off()


