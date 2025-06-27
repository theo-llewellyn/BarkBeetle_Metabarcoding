#ADONIS
library(tidyverse)
library(vegan)
library(cowplot)
library(FUNGuildR)

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

#PCoA of distances
pco_BrayCurtis <- wcmdscale(vegdist(OTU_table_reduced, method = 'bray'), eig = TRUE)
pco_Jaccard <- wcmdscale(vegdist(OTU_table_reduced, method = 'jaccard'), eig = TRUE)
pco_UnwUniFrac <- wcmdscale(UnwUniFrac_diss, eig = TRUE)


#add metdatat
metadat_PcoA <- inner_join(rownames_to_column(as.data.frame(pco_BrayCurtis$points)),
                           sample_data, by = c('rowname' = 'index'))
metadat_PcoA_Jaccard <- inner_join(rownames_to_column(as.data.frame(pco_Jaccard$points)),
                                   sample_data, by = c('rowname' = 'index'))
metadat_PcoA_UnwUniFrac <- inner_join(rownames_to_column(as.data.frame(pco_UnwUniFrac$points)),
                                      sample_data, by = c('rowname' = 'index'))


#as we cant access loadings in PCoA due to information loss when converting to dissimilarity, we can fit the OTUs as environmental variables for the biplot to see which ones correlate with the dispersal of points from PCoA
#bray curtis
efit <- envfit(pco_BrayCurtis, OTU_table_reduced)

#Jaccard
efit_Jaccard <- envfit(pco_Jaccard, OTU_table_reduced)


#unweighted unifrac
efit_unifrac <- envfit(pco_UnwUniFrac, OTU_table_reduced)


#extract a table with those ones for each subfamily
spp.scrs <- as.data.frame(vegan::scores(efit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, OTU = rownames(spp.scrs), Pvalues = efit$vectors$pvals, R_squared = efit$vectors$r, ArrowLength = abs(spp.scrs$Dim1 - spp.scrs$Dim2))

spp.scrs_Jaccard <- as.data.frame(vegan::scores(efit_Jaccard, display = "vectors"))
spp.scrs_Jaccard <- cbind(spp.scrs_Jaccard, OTU = rownames(spp.scrs_Jaccard), Pvalues = efit_Jaccard$vectors$pvals, R_squared = efit_Jaccard$vectors$r, ArrowLength = abs(spp.scrs_Jaccard$Dim1 - spp.scrs_Jaccard$Dim2))

spp.scrs_unifrac <- as.data.frame(vegan::scores(efit_unifrac, display = "vectors"))
spp.scrs_unifrac <- cbind(spp.scrs_unifrac, OTU = rownames(spp.scrs_unifrac), Pvalues = efit_unifrac$vectors$pvals, R_squared = efit_unifrac$vectors$r, ArrowLength = abs(spp.scrs_unifrac$Dim1 - spp.scrs_unifrac$Dim2))

# select significant P-values and then subset those that correlate with the countries based on their position in the biplot
# Extract scores (arrows), R² and p-values
arrows <- scores(efit, display = "vectors")
r2_values <- efit$vectors$r
pvals <- efit$vectors$pvals
# Combine into one data frame
efit_df <- data.frame(arrows, R2 = r2_values, P = pvals)
# Filter for high R² and significant p-value
efit_filtered <- efit_df[efit_df$R2 >= 0.2 & efit_df$P <= 0.01, ]

#show plot
metadat_PcoA %>%
  ggplot(aes(x= Dim1, y = Dim2, colour = country)) +
  scale_colour_manual(values = c("#71B379","#B25690")) + 
  geom_point() +
  #geom_text(aes(label = rowname)) +
  stat_ellipse(type = "norm", show.legend = FALSE) +
  theme_minimal() +
  xlab('PCoA axis 1') +
  ylab('PCoA axis 2') +
  geom_segment(data = efit_filtered,
               aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red") +
  geom_text(data = efit_filtered,
            aes(x = Dim1, y = Dim2, label = rownames(efit_filtered)),
            hjust = 0, nudge_x = 0.05, color = "red") +
  xlim(-0.5,1)



spp.scrs <- subset(spp.scrs, Pvalues < 0.01 & R_squared > 0.2)
#to select only top 10 add the slice_max
spp.scrs.Borneo <- subset(spp.scrs, Dim1 > 0 & Dim2 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs.FG <- subset(spp.scrs, Dim1 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs.Borneo, "dynamicOTUs/significant_OTUs_Borneo.csv")
write_csv(spp.scrs.FG, "dynamicOTUs/significant_OTUs_FG.csv")


#jaccard
spp.scrs_Jaccard <- subset(spp.scrs_Jaccard, Pvalues < 0.01 & R_squared > 0.2)
spp.scrs_Jaccard.Borneo <- subset(spp.scrs_Jaccard, Dim1 > 0 & Dim2 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs_Jaccard.FG <- subset(spp.scrs_Jaccard, Dim1 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs_Jaccard.Borneo, "dynamicOTUs/significant_OTUs_Jaccard_Borneo.csv")
write_csv(spp.scrs_Jaccard.FG, "dynamicOTUs/significant_OTUs_Jaccard_FG.csv")

#unifrac
spp.scrs_unifrac <- subset(spp.scrs_unifrac, Pvalues < 0.01 & R_squared > 0.2)
spp.scrs_unifrac.Borneo <- subset(spp.scrs_unifrac, Dim1 < 0 & Dim2 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs_unifrac.FG <- subset(spp.scrs_unifrac, Dim1 < 0 & Dim2 > 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs_unifrac.Borneo, "dynamicOTUs/significant_OTUs_unifrac_Borneo.csv")
write_csv(spp.scrs_unifrac.FG, "dynamicOTUs/significant_OTUs_unifrac_FG.csv")
