#funguild of communities
library(FUNGuildR)

### FUNGUILD OF SIGNIFICANT OTUs
#link to funguild
fung <- read_csv('/Users/tbl19/OneDrive - Natural History Museum/01_PUBLICATIONS/Ophiostomatoid_Tree/00_DATA/Ophiostomatoid_seqs/MICROASCALES/FUNGuild/FUNGuild_database_2024.csv')

OTUtaxonomy <- read_tsv('OTUs/OTUstaxonomy.txt')

#separate taxonomy into ranks as family level seems to be blocking species level searches
Borneo_fungi <- read_csv('cooccurrence/dynamicOTUs/pos_cooccurrence_modules_Borneo.csv') %>%
  left_join(OTUtaxonomy, by=c('label'='OTU_ID'))

#use species to search
Borneo_results <- funguild_assign(Borneo_fungi, db = fung, tax_col = 'species')

#if result is NA go up a taxon level
for(i in 1:nrow(Borneo_results)){
  print(Borneo_results$species[i])
  #if its not finding the species label try going up a taxon rank
  if(is.na(Borneo_results[i,17])){
    # use the order
    print('trying genus')
    Borneo_results[i,] <- funguild_assign(Borneo_fungi[i,], db = fung, tax_col = 'genus')
    # if its not still finding the taxa
    if(is.na(Borneo_results[i,17])){
      print('trying family')
      # use class
      Borneo_results[i,] <- funguild_assign(Borneo_fungi[i,], db = fung, tax_col = 'family')
      # if its not still finding the taxa
      if(is.na(Borneo_results[i,17])){
        print('trying class')
        # use class
        Borneo_results[i,] <- funguild_assign(Borneo_fungi[i,], db = fung, tax_col = 'class')
        #if its still not finding the class use phylum
        if(is.na(Borneo_results[i,17])){
          print('trying phylum')
          #use the family
          Borneo_results[i,] <- funguild_assign(Borneo_fungi[i,], db = fung, tax_col = 'phylum')
      }
    }
  }
    }else{
  }
}

#plot the proportions of different FUNGUILD trophic modes/lifestyles and facet by community module
df_prop <- subset(Borneo_results, Borneo_results$community != 283 & !is.na(trophicMode)) %>%
  mutate(community = as.factor(community)) %>%
  mutate(community = fct_recode(community, core = "276", accessory = '999')) %>%
  group_by(community, trophicMode) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(community) %>%
  mutate(proportion = count / sum(count))

trophic_colors <- c(
  "Saprotroph" = "#66c2a5",   # teal green
  "Pathotroph-Saprotroph" = "#fc8d62",  # soft orange
  "Pathotroph-Saprotroph-Symbiotroph" = "#8da0cb",  # soft blue
  "Pathotroph" = "#e78ac3",   # pink-mauve
  "Saprotroph-Symbiotroph" = "#a6d854",  # light green
  "Symbiotroph" = "#ffd92f",  # yellow
  "Pathotroph-Symbiotroph" = "#e5c494",  # beige-tan
  "Pathotroph-Pathotroph-Saprotroph" = "#b3b3b3"  # medium gray
)

# Reorder trophicMode within each community
df_prop <- df_prop %>%
  group_by(community) %>%
  arrange(desc(proportion)) %>%
  mutate(
    trophic_fill = paste(community, trophicMode, sep = "_"),  # unique per group
    trophic_fill = factor(trophic_fill, levels = trophic_fill)  # ordered per group
  ) %>%
  ungroup() %>%
  arrange(community)

Borneo <- ggplot(df_prop, aes(y = community, x = proportion, fill = trophic_fill)) +
  geom_col(position = "fill") +
  scale_fill_manual(
    values = setNames(trophic_colors[df_prop$trophicMode], df_prop$trophic_fill),
    labels = df_prop$trophicMode[match(levels(df_prop$trophic_fill), df_prop$trophic_fill)]) +
  theme(legend.position = 'none')



#separate taxonomy into ranks as family level seems to be blocking species level searches
FG_fungi <- read_csv('cooccurrence/dynamicOTUs/pos_cooccurrence_modules_FG.csv') %>%
  left_join(OTUtaxonomy, by=c('label'='OTU_ID'))

#use species to search
FG_results <- funguild_assign(FG_fungi, db = fung, tax_col = 'species')

#if result is NA go up a taxon level
for(i in 1:nrow(FG_results)){
  print(FG_results$species[i])
  #if its not finding the species label try going up a taxon rank
  if(is.na(FG_results[i,17])){
    # use the order
    print('trying genus')
    FG_results[i,] <- funguild_assign(FG_fungi[i,], db = fung, tax_col = 'genus')
    # if its not still finding the taxa
    if(is.na(FG_results[i,17])){
      print('trying family')
      # use class
      FG_results[i,] <- funguild_assign(FG_fungi[i,], db = fung, tax_col = 'family')
      # if its not still finding the taxa
      if(is.na(FG_results[i,17])){
        print('trying class')
        # use class
        FG_results[i,] <- funguild_assign(FG_fungi[i,], db = fung, tax_col = 'class')
        #if its still not finding the class use phylum
        if(is.na(FG_results[i,17])){
          print('trying phylum')
          #use the family
          FG_results[i,] <- funguild_assign(FG_fungi[i,], db = fung, tax_col = 'phylum')
        }
      }
    }
  }else{
  }
}

#plot the proportions of different FUNGUILD trophic modes/lifestyles and facet by community module
df_prop_FG <- FG_results %>%
  drop_na(trophicMode) %>%
  mutate(community = as.factor(community)) %>%
  mutate(community = fct_recode(community, core = "93", accessory = '999')) %>%
  group_by(community, trophicMode) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(community) %>%
  mutate(proportion = count / sum(count))

# Reorder trophicMode within each community
df_prop_FG <- df_prop_FG %>%
  group_by(community) %>%
  arrange(desc(proportion)) %>%
  mutate(
    trophic_fill = paste(community, trophicMode, sep = "_"),  # unique per group
    trophic_fill = factor(trophic_fill, levels = trophic_fill)  # ordered per group
  ) %>%
  ungroup() %>%
  arrange(community)

FG <- ggplot(df_prop_FG, aes(y = community, x = proportion, fill = trophic_fill)) +
  geom_col(position = "fill") +
  scale_fill_manual(
    values = setNames(trophic_colors[df_prop_FG$trophicMode], df_prop_FG$trophic_fill),
    labels = df_prop_FG$trophicMode[match(levels(df_prop_FG$trophic_fill), df_prop_FG$trophic_fill)]) +
  theme(legend.position = 'none')

legend_plot <- get_legend(ggplot(df_prop_FG %>% distinct(trophicMode), aes(x = trophicMode, fill = trophicMode)) +
                            geom_bar() +
                            scale_fill_manual(values = trophic_colors) +
                            theme(legend.direction = "horizontal") +
                            guides(fill = guide_legend(title = "Trophic Mode", ncol = 2)))

png('FIGURES/dynamicOTUs/Borneo_FG_networks_posneg_FUNGUILD.png', res = 400, height = 2400, width = 3000)
plot_grid(Borneo, FG, legend_plot, ncol = 1, rel_heights = c(2,2,1), labels = c('a: Borneo',' b: FG'), hjust = -.1, label_size = 11)
dev.off()

##### INSERT HERE

write_tsv(Borneo_results, 'dynamicOTUs/Borneo_poscooccurrence_FUNGuild_results.txt')
write_tsv(FG_results, 'dynamicOTUs/FG_poscooccurrence_FUNGuild_results.txt')
