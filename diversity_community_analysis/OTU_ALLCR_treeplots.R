########################
## TREE PLOTTING ##
########################
library(tidyverse)
library(RColorBrewer)
library(ggtree)
library(ape)
library(aplot)
library(ggnewscale)
library(ggtreeExtra)
library(picante)
library(cowplot)

#Read in UNITE blastn results
taxonomy_data <- read_csv("UNITE_TBAS_CR.csv")
#read in tree data
#FastTree
tree_data <- read.tree('CRunrooted-tree.nwk')
#iqtree
iqtree <- read.tree('CRotus_3820T_untrimmed_guidance.treefile')
iqtree_trimmed <- read.tree('CRotus_3820T_trimmed_guidance.treefile')
iqtree_trimmed_0.2 <- read.tree('CRotus_3820T_trimmed_gt0.2_guidance.treefile')

#read in sample data
sample_data <- read_csv('UNITE_TBAS_QIIME_CR_consensus_tax_long.csv')

ggplot(funguild_data, aes(x="", fill=trophicMode)) +
  geom_bar() +
  coord_polar("y", start=0) +
  theme_void()


#sort taxonomy amd funguild table by tip labels
taxonomy_data <- taxonomy_data[match(gsub('\'','',tree_data$tip.label), taxonomy_data$query),]
taxonomy_data_iqtree <- taxonomy_data[match(gsub('\'','',iqtree$tip.label), taxonomy_data$query),]
taxonomy_data_iqtree_trimmed <- taxonomy_data[match(gsub('\'','',iqtree_trimmed$tip.label), taxonomy_data$query),]
taxonomy_data_iqtree_trimmed_0.2 <- taxonomy_data[match(gsub('\'','',iqtree_trimmed_0.2$tip.label), taxonomy_data$query),]
  
#make a column showing if OTU is one or both subfamilies and localilities
sample_data[,c(1,11)] %>%
  distinct() %>%
  mutate(present = subfamily) %>%
  pivot_wider(values_from = present, names_from = subfamily) %>%
  mutate_if(is.character, ~replace_na(.,"")) %>%
  mutate(subfamily = paste(Scolytinae,Platypodinae,sep = "")) %>%
  dplyr::select(`#OTU ID`,subfamily) -> sample_subfamilies



#sort sample data the same way
sample_subfamilies <- sample_subfamilies[match(taxonomy_data$query, sample_subfamilies$`#OTU ID`),]

#change tip labels to the taxonomy
#tree_data$tip.label <- taxonomy_data$consensus_three
tree_data$tip.label <- gsub('\'','',tree_data$tip.label)
iqtree$tip.label <- gsub('\'','',iqtree$tip.label)
iqtree_trimmed$tip.label <- gsub('\'','',iqtree_trimmed$tip.label)
iqtree_trimmed_0.2$tip.label <- gsub('\'','',iqtree_trimmed_0.2$tip.label)

#find which OTUs are in Chytridiomycota
outgroup <- subset(taxonomy_data, taxonomy_data$phylum.consensus_3 == 'Chytridiomycota')
tree_data <- root(tree_data, outgroup = outgroup$query[2])
iqtree <- root(iqtree, outgroup = outgroup$query)
iqtree_trimmed <- root(iqtree_trimmed, outgroup = outgroup$query)
iqtree_trimmed_0.2 <- root(iqtree_trimmed_0.2, outgroup = outgroup$query)

tree1 <- ggtree(tree_data, layout = 'circular',branch.length="none")
tree2 <- ggtree(iqtree, layout = 'circular',branch.length="none")
tree3 <- ggtree(iqtree_trimmed, layout = 'circular',branch.length="none")
tree4 <- ggtree(iqtree_trimmed_0.2, layout = 'circular',branch.length="none")

#compare phyla monophyly of four trees
tree_A <- tree1 +
  geom_fruit(data=taxonomy_data, geom=geom_tile,
             mapping=aes(y=query, x='phylum', fill=phylum.consensus_3),
             pwidth = 10) +
  ggtitle(subtitle = 'FastTree untrimmed') +
  theme(legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        plot.title = element_text(hjust = 0.5))
tree_B <- tree2 +
  geom_fruit(data=taxonomy_data, geom=geom_tile,
             mapping=aes(y=query, x='phylum', fill=phylum.consensus_3),
             pwidth = 10) +
  ggtitle(subtitle = 'IQTree untrimmed') +
  theme(legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        plot.title = element_text(hjust = 0.5))
tree_C <- tree3 +
  geom_fruit(data=taxonomy_data, geom=geom_tile,
             mapping=aes(y=query, x='phylum', fill=phylum.consensus_3),
             pwidth = 10) +
  ggtitle(subtitle = 'IQTree trimmed') +
  theme(legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        plot.title = element_text(hjust = 0.5))
tree_D <- tree4 +
  geom_fruit(data=taxonomy_data, geom=geom_tile,
             mapping=aes(y=query, x='phylum', fill=phylum.consensus_3),
             pwidth = 10) +
  ggtitle(subtitle = 'IQTree trimmed gt0.2') +
  theme(legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        plot.title = element_text(hjust = 0.5))
#extract legend
legend <- get_legend(tree4 +
                       geom_fruit(data=taxonomy_data, geom=geom_tile,
                                  mapping=aes(y=query, x='phylum', fill=phylum.consensus_3),
                                  pwidth = 10,
                                  offset = 0.075) +
                       theme(legend.position = 'bottom',
                             legend.direction = "horizontal",
                             legend.box = "vertical",
                             legend.key.size = unit(.25, 'cm'),
                             legend.text = element_text(size=6),
                             legend.title = element_text(size=6),
                             legend.box.spacing = unit(0, "pt"),
                             legend.margin=margin(0,0,0,0),
                             legend.box.margin=margin(0,0,0,0)
                       ))

#plot four trees and monophyly of phyla
top_plot <- plot_grid(tree_A, tree_B, tree_D, tree_C, labels="auto")

png("Tree_comparison_plot.png",  res = 400, width = 2480, height = 2480)
plot_grid(top_plot, legend, ncol = 1, rel_heights = c(1,0.1))
dev.off()
