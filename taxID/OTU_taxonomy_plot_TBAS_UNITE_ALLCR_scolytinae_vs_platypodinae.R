########################
## OTU CLASSIFICATION ##
########################
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(stringi)

#Read in UNITE blastn results
unite <- read_csv("ALLCR_combined_report426CU6C3.csv")

#Remove underscores from species
unite$blast_species <- gsub("_", " ", unite$blast_species)
unite$`Taxon assignment` <- gsub("_", " ", unite$`Taxon assignment`)

#Add column for whether UNITE and T-BAS agree
unite[,c('phylum.match','class.match','order.match','family.match','genus.match','species.match', 'phylum.consensus','class.consensus','order.consensus','family.consensus','genus.consensus','species.consensus','consensus')] <- NA


for (i in 1:nrow(unite)) {
  #if there is no UNITE match then skip everything
  if(!is.na(unite$blast_match[i])){
    #if blast species matches TBAS species then add Y(es) to match column else add N(o)
    if(grepl(unite$blast_species[i], unite$`Taxon assignment`[i], fixed = TRUE) == TRUE){
      unite$species.match[i] <- 'Y'
    }else{
      unite$species.match[i] <- 'N'
    }
    if(grepl(unite$blast_genus[i], unite$`Genus-level assignment`[i], fixed = TRUE) == TRUE){
      unite$genus.match[i] <- 'Y'
    }else{
      unite$genus.match[i] <- 'N'
    }
    if(grepl(unite$blast_family[i], unite$`Family-level assignment`[i], fixed = TRUE) == TRUE){
      unite$family.match[i] <- 'Y'
    }else{
      unite$family.match[i] <- 'N'
    }
    if(grepl(unite$blast_order[i], unite$`Order-level assignment`[i], fixed = TRUE) == TRUE){
      unite$order.match[i] <- 'Y'
    }else{
      unite$order.match[i] <- 'N'
    }
    if(grepl(unite$blast_class[i], unite$`Class-level assignment`[i], fixed = TRUE) == TRUE){
      unite$class.match[i] <- 'Y'
    }else{
      unite$class.match[i] <- 'N'
    }
    if(grepl(unite$blast_phylum[i], unite$`Phylum-level assignment`[i], fixed = TRUE) == TRUE){
      unite$phylum.match[i] <- 'Y'
    }else{
      unite$phylum.match[i] <- 'N'
    }
    #calculate the minimum taxonomic rank at which they agree
    min_taxon_match <- max(which(unite[i,33:38]=='Y'))
    #if they dont agree on any make consensus fungi
    if(min_taxon_match != -Inf){
      #add consensus taxa into consensus columns
      unite[i,c(39:(38+min_taxon_match))] <- unite[i,c(8:(7+min_taxon_match))]
    }
    # if TBAS and UNITE agree on genus but not on species and unite match is >99% then use UNITE species for consensus
    if(unite$genus.match[i] == 'Y' & unite$species.match[i] == 'N' & unite$blast_percent[i] >= 99 & !grepl(" sp", unite$blast_species[i])){
      unite$species.consensus[i] <- unite$blast_species[i]
    }
    unite$consensus[i] <- paste('Fungi', paste(unite[i,c(39:44)],collapse = '__'),sep = '__')
  }
}

#save consensus in new table and remove any sequences without matches
unite_tbas_consensus <- unite[,c(1,2,4,39:44)] %>% 
  filter(!is.na(blast_match)) %>%
  dplyr::select(-blast_match) %>% 
  mutate(across(-query, stri_replace_first_fixed, '_',' '))

### QIIME data
#read in the QIIME taxonomy
CR_taxonomy <- read_tsv("ALLCR_taxonomy.tsv")
#split the taxonomy string into a column for each taxon rank and reformat to remove p__ etc
CR_taxonomy_split <- data.frame(str_split_fixed(CR_taxonomy$Taxon,';',7)) %>% mutate(across(, str_replace, '^.*__','')) %>% mutate(across(, stri_replace_first_fixed, '_',' '))
CR_taxonomy_split <- cbind(query = CR_taxonomy$`Feature ID`, CR_taxonomy$Confidence, CR_taxonomy_split)





#join unite tbas and qiime results
UNITE_TBAS_CR <- left_join(unite_tbas_consensus, CR_taxonomy_split, by = 'query')
UNITE_TBAS_CR[UNITE_TBAS_CR == ''] <- 'Unclassified'
UNITE_TBAS_CR[is.na(UNITE_TBAS_CR)] <- 'Unclassified'

#Add column for whether UNITE and T-BAS agree
UNITE_TBAS_CR[,c('phylum.match','class.match','order.match','family.match','genus.match','species.match', 'phylum.consensus_3','class.consensus_3','order.consensus_3','family.consensus_3','genus.consensus_3','species.consensus_3','consensus_three')] <- NA

for (i in 1:nrow(UNITE_TBAS_CR)) {
  #if blast species matches TBAS species then add Y(es) to match column else add N(o)
  if(grepl(UNITE_TBAS_CR$species.consensus[i], UNITE_TBAS_CR$X7[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$species.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$species.match[i] <- 'N'
  }
  if(grepl(UNITE_TBAS_CR$genus.consensus[i], UNITE_TBAS_CR$X6[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$genus.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$genus.match[i] <- 'N'
  }
  if(grepl(UNITE_TBAS_CR$family.consensus[i], UNITE_TBAS_CR$X5[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$family.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$family.match[i] <- 'N'
  }
  if(grepl(UNITE_TBAS_CR$order.consensus[i], UNITE_TBAS_CR$X4[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$order.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$order.match[i] <- 'N'
  }
  if(grepl(UNITE_TBAS_CR$class.consensus[i], UNITE_TBAS_CR$X3[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$class.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$class.match[i] <- 'N'
  }
  if(grepl(UNITE_TBAS_CR$phylum.consensus[i], UNITE_TBAS_CR$X2[i], fixed = TRUE) == TRUE){
    UNITE_TBAS_CR$phylum.match[i] <- 'Y'
  }else{
    UNITE_TBAS_CR$phylum.match[i] <- 'N'
  }
  #calculate the minimum taxonomic rank at which they agree
  min_taxon_match <- max(which(UNITE_TBAS_CR[i,17:22]=='Y'))
  #if they dont agree on any make consensus fungi
  if(min_taxon_match != -Inf){
    #add consensus taxa into consensus columns
    UNITE_TBAS_CR[i,c(23:(22+min_taxon_match))] <- UNITE_TBAS_CR[i,c(3:(2+min_taxon_match))]
  }
  # if TBAS and UNITE agree on genus but not on species and unite match is >99% then use UNITE species for consensus
  if(UNITE_TBAS_CR$genus.match[i] == 'Y' & UNITE_TBAS_CR$species.match[i] == 'N' & UNITE_TBAS_CR$blast_percent[i] >= 99){
    UNITE_TBAS_CR$species.consensus_3[i] <- UNITE_TBAS_CR$species.consensus[i]
  }
  UNITE_TBAS_CR$consensus_three[i] <- paste('Fungi', paste(UNITE_TBAS_CR[i,c(23:28)],collapse = '__'),sep = '__')
}

##################
# Combine the above data with Sample information so we can subset the data by beetle subfamily

#read in data with OTUs and their presence in each sample
OTU_2_Sample <- read_tsv("ALLCRfiltered_feature-table.tsv", skip = 1)

#make a table of OTU id and consensus classification
data.frame(UNITE_TBAS_CR$query, str_split_fixed(UNITE_TBAS_CR$consensus_three,'__',7)) -> UNITE_TBAS_QIIME_CR_consensus 
#replace empty values as unclassified
UNITE_TBAS_QIIME_CR_consensus[UNITE_TBAS_QIIME_CR_consensus == 'NA'] <- 'Unclassified'

#histogram of number of Classified OTUs at each rank
png("ALLCR_TBAS_UNITE_Classified_OTUs.png",  res = 400, width = 2480, height = 2480)
UNITE_TBAS_QIIME_CR_consensus[,-1] %>%
  pivot_longer(cols = starts_with('X')) %>%
  count(name, value) %>%
  pivot_wider(names_from = name, values_from = n) %>%
  filter(value == 'Unclassified') %>%
  pivot_longer(cols = starts_with('X'), names_to = "name", names_repair = "unique") %>%
  ggplot(aes(x = name, y = 100-(value...3/nrow(UNITE_TBAS_QIIME_CR_consensus))*100)) +
  geom_col() +
  geom_text(aes(label = nrow(UNITE_TBAS_QIIME_CR_consensus)-value...3), vjust = -0.5) +
  scale_x_discrete(labels=c("X1" = "Kingdom","X2" = "Phylum","X3" = "Class","X4" = "Order","X5" = "Family","X6" = "Genus","X7" = "Species"), limits = c("X2","X3","X4","X5","X6","X7")) +
  scale_y_continuous(limits = c(0,100)) +
  ylab("OTUs Classified (%)") +
  theme_minimal() +
  theme(axis.title.x = element_blank())
dev.off()

#join the taxonomy with the sample info tables
UNITE_TBAS_QIIME_CR_consensus_tax <- right_join(OTU_2_Sample,UNITE_TBAS_QIIME_CR_consensus, by = c(`#OTU ID` = 'UNITE_TBAS_CR.query'))

#sample and beetle subfamily
OTU_taxa <- read_csv("ALLCR_level-5.csv")
subfamily_info <- OTU_taxa[,c(1,422)]


UNITE_TBAS_QIIME_CR_consensus_tax_long <- UNITE_TBAS_QIIME_CR_consensus_tax %>%
  pivot_longer(
    cols = starts_with(c("Borneo","FG"))
  )

#remove absent taxa
UNITE_TBAS_QIIME_CR_consensus_tax_long <- UNITE_TBAS_QIIME_CR_consensus_tax_long[UNITE_TBAS_QIIME_CR_consensus_tax_long$value > 0, ]

left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) -> UNITE_TBAS_QIIME_CR_consensus_tax_long_csv

write_csv(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, 'UNITE_TBAS_QIIME_CR_consensus_tax_long.csv')

#subset by subfamily
Scolytinae_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Scolytinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

Platypodinae_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Platypodinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

#plot Venn diagram of shared OTUs
library(VennDiagram)
venn_data <- list(Scolytinae = Scolytinae_OTUs$`#OTU ID`, Platypodinae = Platypodinae_OTUs$`#OTU ID`)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(venn_data, fill = c("#E69F00","#56B4E9"), cat.pos = c(0,0))
venn.diagram(venn_data, filename = "Sco_Plat_CR_venn.png", fill = c("#E69F00","#56B4E9"), cat.pos = c(0,10), category.names = c('Fungal OTUs\n Scolytinae','Fungal OTUs\n Platypodinae'))

#Venn diagram of Locality
#subset by subfamily
Borneo_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('Borneo', name)) %>% 
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

FG_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('FG', name)) %>% 
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

#plot Venn diagram of shared OTUs
venn_data_locality <- list(Borneo = Borneo_OTUs$`#OTU ID`, FG = FG_OTUs$`#OTU ID`)
display_venn(venn_data_locality, fill = c("#B25690","#71B379"), cat.pos = c(0,0))
venn.diagram(venn_data_locality, filename = "Borneo_FG_CR_venn.png", fill = c("#B25690","#71B379"), cat.pos = c(0,0), category.names = c('Fungal OTUs\n Borneo','Fungal OTUs\n French Guiana'))

#split the Sco and Platy into two sites as well
Scolytinae_Borneo_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Scolytinae') %>%
  filter(grepl('Borneo', name)) %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

Scolytinae_FG_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Scolytinae') %>%
  filter(grepl('FG', name)) %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

venn_data_locality_Sco <- list(Borneo = Scolytinae_Borneo_OTUs$`#OTU ID`, FG = Scolytinae_FG_OTUs$`#OTU ID`)
display_venn(venn_data_locality_Sco, fill = c("#B25690","#71B379"), cat.pos = c(0,0))
venn.diagram(venn_data_locality_Sco, filename = "Scolytinae_Borneo_FG_CR_venn.png", fill = c("#B25690","#71B379"), cat.pos = c(0,0), category.names = c('Fungal OTUs\n Borneo','Fungal OTUs\n French Guiana'))


Platypodinae_Borneo_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Platypodinae') %>%
  filter(grepl('Borneo', name)) %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

Platypodinae_FG_OTUs <- left_join(UNITE_TBAS_QIIME_CR_consensus_tax_long, subfamily_info, by = c('name' = 'index')) %>%
  filter(subfamily == 'Platypodinae') %>%
  filter(grepl('FG', name)) %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

venn_data_locality_Pla <- list(Borneo = Platypodinae_Borneo_OTUs$`#OTU ID`, FG = Platypodinae_FG_OTUs$`#OTU ID`)
display_venn(venn_data_locality_Pla, fill = c("#B25690","#71B379"), cat.pos = c(0,0), inverted = TRUE)
venn.diagram(venn_data_locality_Pla, filename = "Platypodinae_Borneo_FG_CR_venn.png", fill = c("#B25690","#71B379"), cat.pos = c(0,0), category.names = c('Fungal OTUs\n Borneo','Fungal OTUs\nFrench Guiana'), inverted = TRUE)

#split localities into subfamilies
Borneo_Platy_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('Borneo', name)) %>%
  filter(subfamily == 'Platypodinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

Borneo_Scoly_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('Borneo', name)) %>%
  filter(subfamily == 'Scolytinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

#plot Venn diagram of shared OTUs
venn_data_Borneo_subfam <- list(Platypodinae = Borneo_Platy_OTUs$`#OTU ID`, Scolytinae = Borneo_Scoly_OTUs$`#OTU ID`)
display_venn(venn_data_Borneo_subfam, fill = c("#56B4E9","#E69F00"), cat.pos = c(0,0))
venn.diagram(venn_data_Borneo_subfam, filename = "Borneo_Sco_Plat_CR_venn.png", fill = c("#56B4E9","#E69F00"), cat.pos = c(0,0), category.names = c('Fungal OTUs\n Platypodinae','Fungal OTUs\n Scolytinae'))

FG_Platy_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('FG', name)) %>%
  filter(subfamily == 'Platypodinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

FG_Scoly_OTUs <- filter(UNITE_TBAS_QIIME_CR_consensus_tax_long_csv, grepl('FG', name)) %>%
  filter(subfamily == 'Scolytinae') %>%
  dplyr::select(-c('name','value','subfamily')) %>%
  distinct()

#plot Venn diagram of shared OTUs
venn_data_FG_subfam <- list(Platypodinae = FG_Platy_OTUs$`#OTU ID`, Scolytinae = FG_Scoly_OTUs$`#OTU ID`)
display_venn(venn_data_FG_subfam, fill = c("#56B4E9","#E69F00"), cat.pos = c(0,0))
venn.diagram(venn_data_FG_subfam, filename = "FG_Sco_Plat_CR_venn.png", fill = c("#56B4E9","#E69F00"), cat.pos = c(0,0), category.names = c('Fungal OTUs\n Platypodinae','Fungal OTUs\n Scolytinae'))


Scolytinae_OTUs[,-c(1,2,7,8)] %>%
  add_count(X2, name = "X2_count") %>%
  add_count(X3, name = "X3_count") %>%
  add_count(X4, name = "X4_count") %>%
  add_count(X5, name = "X5_count") %>%
  arrange(desc(X2_count),desc(X3_count),desc(X4_count),desc(X5_count)) %>%
  select(-c(X2_count,X3_count,X4_count,X5_count)) %>%
  mutate(phylum = X2) -> Scolytinae_OTUs_sorted
#if it says unclassified change it to the name of the upper taxonomic rank at which it was classified, that way unclassified taxa in different higher taxa arent merged
for(col in 2:ncol(Scolytinae_OTUs_sorted)){
  for(i in 1:nrow(Scolytinae_OTUs_sorted)){
    if(Scolytinae_OTUs_sorted[i,col]=='Unclassified'){
      Scolytinae_OTUs_sorted[i,col] = paste(Scolytinae_OTUs_sorted[i,col-1],'_unclass.',sep='')
    }
  }
}

#we will assign colours to the names of the phyla manually so they are the same colours across all plots and datasets
unique_phyla <- c('Ascomycota','Basidiomycota','Chytridiomycota','Entomophthoromycota','Fungi_unclass.','Kickxellomycota','Mortierellomycota','Mucoromycota')

palette1_named = setNames(object = brewer.pal(8, "Set2"), nm = unique_phyla)

#donut plot
inset <- Scolytinae_OTUs_sorted %>%
  pivot_longer(1:4) %>%
  group_by(name, value) %>%
  mutate(width = n()) %>%
  unique() %>%
  arrange(phylum,name) %>%
  group_by(name) %>%
  mutate(ymid = as.numeric(sub("\\D+", "", name)),
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = c(0, head(cumsum(width), -1)),
         xmax = cumsum(width),
         xmid = (xmax + xmin) / 2) %>%
  ggplot(aes(xmid, ymid, fill = phylum == "Unclassified_unclass._unclass._unclass._unclass.")) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = phylum  == "Unclassified_unclass._unclass._unclass._unclass.")) +
  scale_fill_manual(values = c('grey','grey20'), labels = c("Classified","Unclassified")) +
  scale_colour_manual(values = c('grey','grey20'), guide = "none") +
  coord_polar() +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.key.size = unit(.25, 'cm'),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))

main <- Scolytinae_OTUs_sorted %>%
  filter(X2 != 'Unclassified') %>%
  pivot_longer(1:4) %>%
  group_by(name, value) %>%
  mutate(width = n()) %>%
  unique() %>%
  arrange(phylum,name) %>%
  group_by(name) %>%
  mutate(ymid = as.numeric(sub("\\D+", "", name)),
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = c(0, head(cumsum(width), -1)),
         xmax = cumsum(width),
         xmid = (xmax + xmin) / 2) %>%
  ggplot(aes(xmid, ymid, fill = phylum)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                alpha = name, color = phylum )) +
  scale_alpha_manual(values = seq(1,0.3,length.out = 5), guide = "none") +
  scale_fill_manual(values = palette1_named) +
  scale_colour_manual(values = rep("white",7), guide = "none") +
  geom_text_repel(data = . %>% filter(width/4510 > 0.005 & !grepl('_unclass.', value) & ymid<5),
                  aes(label = value, group = value, y = ymid), 
                  size = 2,
                  max.overlaps = 30,
                  show.legend = FALSE,
                  min.segment.length = .1) +
  geom_label_repel(data = . %>% filter(width/11292 > 0.001 & !grepl('_unclass.', value) & ymid == 5),
                   aes(label = value, group = value, y = ymax),
                   nudge_y = 0.5,
                   size = 2,
                   max.overlaps = Inf,
                   show.legend = FALSE,
                   min.segment.length = 0) +
  coord_polar() + 
  theme_void() +
  theme(legend.title = element_blank(),
        legend.key.size = unit(.25, 'cm'),
        legend.text = element_text(size=5),
        #plot.margin=grid::unit(c(-20,0,-20,0), "mm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))

library(cowplot)
plot.with.inset <-
  ggdraw() +
  draw_plot(main) +
  draw_plot(inset, x = .7, y = .76, width = .3, height = .3)

png("ALLCR_OTU_taxonomy_TBAS_UNITE_Scolytinae_inset.png",  res = 400, width = 2480, height = 2480)
plot.with.inset
dev.off()

#can then repeat previous section replacing Scolytinae for Platypodinae
