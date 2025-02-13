#ADONIS
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
sample_data <- read_csv('Sco_Pla_FG_Borneo_metadata.csv') %>%
  dplyr::select(c(fungal_metabarcode_ID,subfamily,country)) %>% 
  rename(index = fungal_metabarcode_ID) %>%
  replace_na(list(subfamily = 'Platypodinae')) %>%
  filter(index %in% colnames(BrayCurtis_diss))

#Danum Valley to remove
Danum_Valley <- c('FG_Sco_P2_A09','FG_Sco_P2_A10','FG_Sco_P2_B09','FG_Sco_P2_B10','FG_Sco_P2_C09','FG_Sco_P2_C10','FG_Sco_P2_D09','FG_Sco_P2_D10','FG_Sco_P2_E09','FG_Sco_P2_E10','FG_Sco_P2_F09','FG_Sco_P2_F10','FG_Sco_P2_G09','FG_Sco_P2_H08','FG_Sco_P2_H09')
#remove Danum, NA index and P2_A12 as its not in filtered OTU data
sample_data <- sample_data[ ! sample_data$index %in% Danum_Valley, ] %>%
  arrange(index)

#PERMANOVA testing effect of 
adonis2(as.dist(BrayCurtis_diss) ~ country * subfamily, data = sample_data)
adonis2(as.dist(Jaccard_diss) ~ country * subfamily, data = sample_data)
adonis2(as.dist(UnwUniFrac_diss) ~ country * subfamily, data = sample_data)
adonis2(as.dist(WeUniFrac_diss) ~ country * subfamily, data = sample_data)

#PERMDISP
#BrayCurtis
bd_subfam <- betadisper(as.dist(BrayCurtis_diss),sample_data$subfamily)
permutest(bd_subfam)
bd_country <- betadisper(as.dist(BrayCurtis_diss),sample_data$country)
permutest(bd_country)
#Jaccard
jaccard_subfam <- betadisper(as.dist(Jaccard_diss),sample_data$subfamily)
permutest(jaccard_subfam)
jaccard_country <- betadisper(as.dist(Jaccard_diss),sample_data$country)
permutest(jaccard_country)
#Unw UniFrac
UnwUniFrac_subfam <- betadisper(as.dist(UnwUniFrac_diss),sample_data$subfamily)
permutest(UnwUniFrac_subfam)
UnwUniFrac_country <- betadisper(as.dist(UnwUniFrac_diss),sample_data$country)
permutest(UnwUniFrac_country)
#We UniFrac
WeUniFrac_subfam <- betadisper(as.dist(WeUniFrac_diss),sample_data$subfamily)
permutest(WeUniFrac_subfam)
WeUniFrac_country <- betadisper(as.dist(WeUniFrac_diss),sample_data$country)
permutest(WeUniFrac_country)



#PCoA of distances
pco_BrayCurtis <- wcmdscale(as.dist(BrayCurtis_diss), eig = TRUE)
pco_Jaccard <- wcmdscale(as.dist(Jaccard_diss), eig = TRUE)
pco_UnwUniFrac <- wcmdscale(as.dist(UnwUniFrac_diss), eig = TRUE)
pco_WeUniFrac <- wcmdscale(as.dist(WeUniFrac_diss), eig = TRUE)

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
  scale_colour_manual(values = c("#B25690","#71B379")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_BrayCurtis_country

#plot PcoA
metadat_PcoA_Jaccard %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_Jaccard_subfam

#by country
metadat_PcoA_Jaccard %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#B25690","#71B379")) + 
  geom_point(show.legend=FALSE) +
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
  scale_colour_manual(values = c("#B25690","#71B379")) + 
  geom_point(show.legend=FALSE) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_UnwUniFrac_country

#plot PcoA
metadat_PcoA_WeUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = subfamily)) +
  scale_colour_manual(values = c("#E69F00","#56B4E9")) + 
  geom_point() +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_WeUniFrac_subfam

#by country
metadat_PcoA_WeUniFrac %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#B25690","#71B379")) + 
  geom_point() +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') -> PCoA_WeUniFrac_country

#get legend
legend_subfam <- get_legend(
  # create some space to the left of the legend
  PCoA_WeUniFrac_subfam + theme(legend.box.margin = margin(0, 0, 0, 12), legend.direction = "horizontal")
)
legend_country <- get_legend(
  # create some space to the left of the legend
  PCoA_WeUniFrac_country + theme(legend.box.margin = margin(0, 0, 0, 12), legend.direction = "horizontal")
)

png("FIGURES/ALLCR_PCoAs_IQTree.png", res = 400, height = 3508, width = 2480)
plot_grid(legend_subfam, legend_country,
          PCoA_BrayCurtis_subfam, PCoA_BrayCurtis_country,
          PCoA_Jaccard_subfam, PCoA_Jaccard_country, 
          PCoA_UnwUniFrac_subfam, PCoA_UnwUniFrac_country, 
          PCoA_WeUniFrac_subfam + theme(legend.position="none"), PCoA_WeUniFrac_country + theme(legend.position="none"),
          labels=c('','','a','b','c','d','e','f','g','h'), 
          ncol = 2, rel_heights = c(0.25,1,1,1,1)
          )
dev.off()
