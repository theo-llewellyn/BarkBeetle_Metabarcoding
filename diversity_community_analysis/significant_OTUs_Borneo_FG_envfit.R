#pull out which OTUs significantly correlate with the clustering of points on PCoA plots
library(tidyverse)
library(vegan)
library(cowplot)
library(FUNGuildR)

#read in dissimilarity matrices
BrayCurtis_diss <- read_tsv('CRcoremetrics/bray_curtis_distance_matrix/data/distance-matrix.tsv') %>%
  column_to_rownames("...1")
UnwUniFrac_diss <- read_tsv('IQtreeCRcoremetrics/unweighted_unifrac/unweighted_unifrac_distance_matrix/distance-matrix.tsv') %>%
  column_to_rownames("...1")
Jaccard_diss <- read_tsv('CRcoremetrics/jaccard_distance_matrix/data/distance-matrix.tsv') %>%
  column_to_rownames("...1")
#read in sample data
sample_data <- read_csv('Sco_Pla_FG_Borneo_metadata.csv') %>%
  dplyr::select(c(fungal_metabarcode_ID,subfamily,country)) %>% 
  rename(index = fungal_metabarcode_ID) %>%
  replace_na(list(subfamily = 'Platypodinae')) %>%
  filter(index %in% colnames(BrayCurtis_diss))

#PCoA of distances
pco_BrayCurtis <- wcmdscale(as.dist(BrayCurtis_diss), eig = TRUE)
pco_UnwUniFrac <- wcmdscale(as.dist(UnwUniFrac_diss), eig = TRUE)
pco_Jaccard <- wcmdscale(as.dist(Jaccard_diss), eig = TRUE)

#add metdatat
metadat_PcoA <- inner_join(rownames_to_column(as.data.frame(pco_BrayCurtis$points)),
                           sample_data, by = c('rowname' = 'index'))
metadat_PcoA_UnwUniFrac <- inner_join(rownames_to_column(as.data.frame(pco_UnwUniFrac$points)),
                                      sample_data, by = c('rowname' = 'index'))
metadat_PcoA_Jaccard <- inner_join(rownames_to_column(as.data.frame(pco_Jaccard$points)),
                                   sample_data, by = c('rowname' = 'index'))


#as we cant access loadings in PCoA due to information loss when converting to dissimilarity, we can fit the OTUs as environmental variables for the biplot to see which ones correlate with the dispersal of points from PCoA
#bray curtis
envfit_data <- read_csv('ALLCR_level-5.csv') %>%
  filter(index %in% colnames(BrayCurtis_diss))
envfit_data_OTUs <- envfit_data[,-c(421:426)] %>% column_to_rownames('index')
efit <- envfit(pco_BrayCurtis, envfit_data)

#unweighted unifrac
envfit_data_unifrac <- read_csv('ALLCR_level-5.csv') %>%
  filter(index %in% colnames(UnwUniFrac_diss))
envfit_data_unifrac_OTUs <- envfit_data_unifrac[,-c(421:426)] %>% column_to_rownames('index')
efit_unifrac <- envfit(pco_UnwUniFrac, envfit_data_unifrac)

#Jaccard
envfit_data_Jaccard <- read_csv('ALLCR_level-5.csv') %>%
  filter(index %in% colnames(Jaccard_diss))
envfit_data_Jaccard_OTUs <- envfit_data_Jaccard[,-c(421:426)] %>% column_to_rownames('index')
efit_Jaccard <- envfit(pco_Jaccard, envfit_data_Jaccard)


#function to extract variable and make a list of colours
taxon_colour_function <- function(colours,taxa_list, variable){
  taxon_vector <- c()
  #for each taxon on the list of taxa
  for(taxon in taxa_list){
    #find which group its in at a particular taxonomic rank
    group <- subset(sample_data, index == taxon)[,variable][[1]]
    taxon_vector <- c(taxon_vector,group)
  }
  taxon_vector <- as.factor(taxon_vector)
  #make a vector of colours using the colours given and the taxa list
  colour_vector <- colours[taxon_vector]
  #make a vector given the levels of the taxa
  level_vector <- levels(taxon_vector)
  #combine into data frame to be able to access
  output <- data.frame(taxon_vector,colour_vector)
}

colours_country <- taxon_colour_function(colours = c("#B25690","#71B379"),
                                                taxa_list = sample_data$index,
                                                variable = "country")


#we can then plot the significant ones on top of the biplot
png("FIGURES/BrayCurtis_PCoA_envfittOTUs.png",  res = 300, width = 2000, height = 2000)
plot(pco_BrayCurtis,type = "n")
points(pco_BrayCurtis$points, display = "sites",
       scaling = "symmetric",
       col = colours_country$colour_vector)
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_BrayCurtis$points, groups = na.omit(colours_country$taxon_vector),
            kind = "ehull", col = c("#B25690","#71B379"),
            scaling = "symmetric", lwd = 2, lty = 2)
legend("topright",
       legend = unique(colours_country$taxon_vector),
       col = unique(colours_country$colour_vector),
       pch = 19)
plot(efit, col = "red", cex = 0.5, p.max=0.05)
dev.off()

png("FIGURES/UnwUnifrac_PCoA_envfitOTUs.png",  res = 300, width = 2000, height = 2000)
plot(pco_UnwUniFrac,type = "n")
points(pco_UnwUniFrac$points, display = "sites",
       scaling = "symmetric",
       col = colours_country$colour_vector)
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_UnwUniFrac$points, groups = na.omit(colours_country$taxon_vector),
            kind = "ehull", col = c("#B25690","#71B379"),
            scaling = "symmetric", lwd = 2, lty = 2)
legend("topright",
       legend = unique(colours_country$taxon_vector),
       col = unique(colours_country$colour_vector),
       pch = 19)
plot(efit_unifrac, col = "red", cex = 0.5, p.max=0.05)
dev.off()

png("FIGURES/Jaccard_PCoA_envfittOTUs.png",  res = 300, width = 2000, height = 2000)
plot(pco_Jaccard,type = "n")
points(pco_Jaccard$points, display = "sites",
       scaling = "symmetric",
       col = colours_country$colour_vector)
abline(h = 0, v = 0, lty = 2)
ordiellipse(pco_Jaccard$points, groups = na.omit(colours_country$taxon_vector),
            kind = "ehull", col = c("#B25690","#71B379"),
            scaling = "symmetric", lwd = 2, lty = 2)
legend("topright",
       legend = unique(colours_country$taxon_vector),
       col = unique(colours_country$colour_vector),
       pch = 19)
plot(efit_Jaccard, col = "red", cex = 0.5, p.max=0.05)
dev.off()

#extract a table with those ones for each subfamily
spp.scrs <- as.data.frame(vegan::scores(efit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, OTU = rownames(spp.scrs), Pvalues = efit$vectors$pvals, R_squared = efit$vectors$r, ArrowLength = abs(spp.scrs$Dim1 - spp.scrs$Dim2))

spp.scrs_unifrac <- as.data.frame(vegan::scores(efit_unifrac, display = "vectors"))
spp.scrs_unifrac <- cbind(spp.scrs_unifrac, OTU = rownames(spp.scrs_unifrac), Pvalues = efit_unifrac$vectors$pvals, R_squared = efit_unifrac$vectors$r, ArrowLength = abs(spp.scrs_unifrac$Dim1 - spp.scrs_unifrac$Dim2))

spp.scrs_Jaccard <- as.data.frame(vegan::scores(efit_Jaccard, display = "vectors"))
spp.scrs_Jaccard <- cbind(spp.scrs_Jaccard, OTU = rownames(spp.scrs_Jaccard), Pvalues = efit_Jaccard$vectors$pvals, R_squared = efit_Jaccard$vectors$r, ArrowLength = abs(spp.scrs_Jaccard$Dim1 - spp.scrs_Jaccard$Dim2))

# select significant P-values and then subset those that correlate with the countries based on their position in the biplot
spp.scrs <- subset(spp.scrs, Pvalues < 0.01)
#to select only top 10 add the slice_max
spp.scrs.Borneo <- subset(spp.scrs, Dim1 < 0 ) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs.FG <- subset(spp.scrs, Dim1 > 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs.Borneo, "significant_OTUs_Borneo.csv")
write_csv(spp.scrs.FG, "significant_OTUs_FG.csv")

#unifrac
spp.scrs_unifrac <- subset(spp.scrs_unifrac, Pvalues < 0.01)
spp.scrs_unifrac.Borneo <- subset(spp.scrs_unifrac, Dim1 < 0 & Dim2 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs_unifrac.FG <- subset(spp.scrs_unifrac, Dim1 < 0 & Dim2 > 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs_unifrac.Borneo, "significant_OTUs_unifrac_Borneo.csv")
write_csv(spp.scrs_unifrac.FG, "significant_OTUs_unifrac_FG.csv")

#jaccard
spp.scrs_Jaccard <- subset(spp.scrs_Jaccard, Pvalues < 0.01)
spp.scrs_Jaccard.FG <- subset(spp.scrs_Jaccard, Dim1 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
spp.scrs_Jaccard.Borneo <- subset(spp.scrs_Jaccard, Dim1 > 0 & Dim2 < 0) #%>% slice_max(order_by = ArrowLength, n = 10)
write_csv(spp.scrs_Jaccard.Borneo, "significant_OTUs_Jaccard_Borneo.csv")
write_csv(spp.scrs_Jaccard.FG, "significant_OTUs_Jaccard_FG.csv")
