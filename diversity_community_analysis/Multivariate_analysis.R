#ADONIS
library(tidyverse)
library(vegan)
library(cowplot)
library(phyloseq)
library(ape)

setwd("~/OneDrive - Natural History Museum/01_PUBLICATIONS/ITS_Fungi_Metabarcoding/RESULTS/")

#read in dissimilarity matrices
OTU_table <- read_tsv('OTUs/OTUsotu_table.txt') %>% 
  slice(-1) %>%
  column_to_rownames(var = "OTU_ID")

#mean beetle occupancy
summary(rowSums(OTU_table > 0))
#mean no. OTUs per beetle
summary(colSums(OTU_table > 0))

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
OTU_table_reduced <- OTU_table[ ! rownames(OTU_table) %in% missing_metadata, ]


#Calculate UniFrac indices
#read tree
iqtree <- read.tree('TREES/FastTree_OTUs.tree')
iqtree$tip.label <- gsub('_',':',iqtree$tip.label)
Physeq <- phyloseq(otu_table(OTU_table_reduced, taxa_are_rows=FALSE),phy_tree(iqtree))
UnwUniFrac_diss <- UniFrac(Physeq, weighted = FALSE, parallel = TRUE)
WeUniFrac_diss <- UniFrac(Physeq, weighted = TRUE, parallel = TRUE)

#PCoA of distances
pco_BrayCurtis <- wcmdscale(vegdist(OTU_table_reduced, method = 'bray'), eig = TRUE)
pco_Jaccard <- wcmdscale(vegdist(OTU_table_reduced, method = 'jaccard'), eig = TRUE)
pco_UnwUniFrac <- wcmdscale(UnwUniFrac_diss, eig = TRUE)
pco_WeUniFrac <- wcmdscale(WeUniFrac_diss, eig = TRUE)


#add metdatat
metadat_PcoA <- inner_join(rownames_to_column(as.data.frame(pco_BrayCurtis$points)),
                           sample_data, by = c('rowname' = 'index'))
metadat_PcoA_Jaccard <- inner_join(rownames_to_column(as.data.frame(pco_Jaccard$points)),
                                   sample_data, by = c('rowname' = 'index'))
metadat_PcoA_UnwUniFrac <- inner_join(rownames_to_column(as.data.frame(pco_UnwUniFrac$points)),
                           sample_data, by = c('rowname' = 'index'))
metadat_PcoA_WeUniFrac <- inner_join(rownames_to_column(as.data.frame(pco_WeUniFrac$points)),
                           sample_data, by = c('rowname' = 'index'))

#plot PcoA
metadat_PcoA %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_BrayCurtis_subfam
#by country
metadat_PcoA %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_BrayCurtis_country

#plot PcoA
metadat_PcoA_Jaccard %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point() +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_Jaccard_subfam
#by country
metadat_PcoA_Jaccard %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_point() +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_Jaccard_country

#plot PcoA
metadat_PcoA_UnwUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_UnwUniFrac_subfam
#by country
metadat_PcoA_UnwUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_UnwUniFrac_country

#plot PcoA
metadat_PcoA_WeUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_WeUniFrac_subfam
#by country
metadat_PcoA_WeUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_WeUniFrac_country


#get legend
legend_subfam <- get_legend(
  # create some space to the left of the legend
  PCoA_Jaccard_subfam + theme(legend.box.margin = margin(0, 0, 0, 12), legend.direction = "horizontal")
)
legend_country <- get_legend(
  # create some space to the left of the legend
  PCoA_Jaccard_country + theme(legend.box.margin = margin(0, 0, 0, 12), legend.direction = "horizontal")
)

png("FIGURES/dynamicOTUs/dynamicOTUS_PCOAs_all.png", res = 400, height = 3508, width = 2480)
plot_grid(legend_subfam, legend_country,
          PCoA_BrayCurtis_subfam, PCoA_BrayCurtis_country,
          PCoA_Jaccard_subfam + theme(legend.position="none"),  PCoA_Jaccard_country + theme(legend.position="none"),
          PCoA_UnwUniFrac_subfam + theme(legend.position="none"),  PCoA_UnwUniFrac_country + theme(legend.position="none"),
          PCoA_WeUniFrac_subfam + theme(legend.position="none"),  PCoA_WeUniFrac_country + theme(legend.position="none"),
          labels=c('','','a','b','c','d','e','f'), 
          ncol = 2, rel_heights = c(0.25,1,1,1,1)
)
dev.off()


#PERMANOVA testing effect of 
adonis2(OTU_table_reduced ~ country * subfamily, method = "jaccard", data = sample_data)
adonis2(OTU_table_reduced ~ country * subfamily, method = "bray", data = sample_data)
adonis2(UnwUniFrac_diss ~ country * subfamily, data = sample_data)
adonis2(WeUniFrac_diss ~ country * subfamily, data = sample_data)

#PERMDISP
#Jaccard
jaccard_subfam <- betadisper(vegdist(OTU_table_reduced, method = 'jaccard'),sample_data$subfamily)
permutest(jaccard_subfam)
jaccard_country <- betadisper(vegdist(OTU_table_reduced, method = 'jaccard'),sample_data$country)
permutest(jaccard_country)
#BrayCurtis
bd_subfam <- betadisper(vegdist(OTU_table_reduced, method = 'bray'), sample_data$subfamily)
permutest(bd_subfam)
bd_country <- betadisper(vegdist(OTU_table_reduced, method = 'bray'),sample_data$country)
permutest(bd_country)
#UnwUniFrac
UnwUniFrac_subfam <- betadisper(UnwUniFrac_diss, sample_data$subfamily)
permutest(UnwUniFrac_subfam)
UnwUniFrac_country <- betadisper(UnwUniFrac_diss,sample_data$country)
permutest(UnwUniFrac_country)
#BrayCurtis
WeUniFrac_subfam <- betadisper(WeUniFrac_diss, sample_data$subfamily)
permutest(WeUniFrac_subfam)
WeUniFrac_country <- betadisper(WeUniFrac_diss,sample_data$country)
permutest(WeUniFrac_country)

