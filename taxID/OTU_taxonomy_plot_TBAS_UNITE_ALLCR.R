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
  select(-blast_match) %>% 
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

write_csv(UNITE_TBAS_CR, 'UNITE_TBAS_CR.csv')

# prepare data for pie chart
data.frame(str_split_fixed(UNITE_TBAS_CR$consensus_three,'__',7)) -> UNITE_TBAS_CR_consensus 
UNITE_TBAS_CR_consensus[UNITE_TBAS_CR_consensus == 'NA'] <- 'Unclassified'

UNITE_TBAS_CR_consensus %>%
  add_count(X2, name = "X2_count") %>%
  add_count(X3, name = "X3_count") %>%
  add_count(X4, name = "X4_count") %>%
  add_count(X5, name = "X5_count") %>%
  arrange(desc(X2_count),desc(X3_count),desc(X4_count),desc(X5_count)) %>%
  select(-c(X2_count,X3_count,X4_count,X5_count)) %>%
  mutate(phylum = X2) -> UNITE_TBAS_CR_consensus_sorted

#if it says unclassified change it to the name of the upper taxonomic rank at which it was classified, that way unclassified taxa in different higher taxa arent merged
for(col in 2:ncol(UNITE_TBAS_CR_consensus_sorted)){
  for(i in 1:nrow(UNITE_TBAS_CR_consensus_sorted)){
    if(UNITE_TBAS_CR_consensus_sorted[i,col]=='Unclassified'){
      UNITE_TBAS_CR_consensus_sorted[i,col] = paste(UNITE_TBAS_CR_consensus_sorted[i,col-1],'_unclass.',sep='')
    }
  }
}

#pie chart
png("ALLCR_OTU_taxonomy_TBAS_UNITE_all.png",  res = 400, width = 2480, height = 2480)
UNITE_TBAS_CR_consensus_sorted[,c(2:5,8)] %>%
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
                alpha = name, color = phylum)) +
  scale_alpha_manual(values = seq(1,0.3,length.out = 5), guide = "none") +
  scale_fill_manual(values = brewer.pal(8, "Set2"), labels = sort(unique(UNITE_TBAS_CR_consensus_sorted$X2))) +
  scale_colour_manual(values = rep("white",13), guide = "none") +
  geom_text_repel(data = . %>% filter(width/3718 > 0.005 & !grepl('_unclass.', value) & ymid<5),
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
dev.off()
