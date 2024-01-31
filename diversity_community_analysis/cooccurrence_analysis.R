#------------------------------------------------------------------------------------#
# Community Co-occurrence analysis

library(cooccur)
library(tidyverse)
library(reshape2)
source(file = 'plot.coocur.R')

#Read in UNITE blastn results
taxonomy_data <- read_csv("UNITE_TBAS_CR.csv") %>%
  select(c(query, consensus_three))

#read in data with OTUs and their presence in each sample
OTU_2_Sample <- read_tsv("ALLCRfiltered_feature-table.tsv", skip = 1) %>% 
  #convert to presence absence
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  #merge with taxonomy info
  left_join(.,taxonomy_data, by = c(`#OTU ID` = 'query')) %>%
  #remove OTU_ID_col
  select(-"#OTU ID") %>%
  #make taxa names unique
  mutate(consensus_three = paste(row_number(),consensus_three,sep = "_")) %>%
  #set OTU ID as rownames
  column_to_rownames(var = "consensus_three")

#replace NAs with Unclassified in taxa names and then cut to the lowest identified taxon rank
rownames(OTU_2_Sample) <- rownames(OTU_2_Sample) %>% 
  gsub('NA',"Unclassified",.) %>% 
  gsub('__Unclassified.*','_sp.',.)  %>% 
  str_remove(.,'Fungi__[:alpha:].*__')

#plot histogram of number of beetles each OTU was found in
ggplot(tibble(rowSums(tibble(OTU_2_Sample)),rownames(OTU_2_Sample)), aes(x = `rowSums(tibble(OTU_2_Sample))`)) +
  geom_histogram()
#get mean, sd, range
summary(rowSums(tibble(OTU_2_Sample)))

#plot histogram of number of OTUs each beetle had
ggplot(tibble(colSums(tibble(OTU_2_Sample)),colnames(OTU_2_Sample)), aes(x = `colSums(tibble(OTU_2_Sample))`)) +
  geom_histogram()
#get mean, sd, range
summary(colSums(tibble(OTU_2_Sample)))

#remove singleton OTUs
OTU_2_Sample <- OTU_2_Sample[-which(rowSums(OTU_2_Sample[sapply(OTU_2_Sample, is.numeric)]) == 1),]
#remove OTUs in less than 5% of samples
OTU_2_Sample <- OTU_2_Sample[-which(rowSums(OTU_2_Sample[sapply(OTU_2_Sample, is.numeric)]) < 27),]
#remove OTUs in less than 15% of samples
OTU_2_Sample_15perc <- OTU_2_Sample[-which(rowSums(OTU_2_Sample[sapply(OTU_2_Sample, is.numeric)]) < 82),]

#Borneo
#if we want the same taxa in the full analysis
OTU_2_Sample_Borneo <- OTU_2_Sample_15perc %>% select(Borneo_Sco_P1_A01:Borneo_Sco_P3_H12)
#if we want the 15% taxa across Borneo sites
OTU_2_Sample_Borneo_15perc <- OTU_2_Sample_Borneo[-which(rowSums(OTU_2_Sample_Borneo[sapply(OTU_2_Sample_Borneo, is.numeric)]) < 43),]

#French Guiana
OTU_2_Sample_FG <- OTU_2_Sample_15perc %>% select(FG_Sco_P1_A01:FG_Sco_P3_H08)
#if we want the 15% taxa across Borneo sites
OTU_2_Sample_FG_15perc <- OTU_2_Sample_FG[-which(rowSums(OTU_2_Sample_FG[sapply(OTU_2_Sample_FG, is.numeric)]) < 38),]

set.seed(071123)
cooccur.OTUs <- cooccur(mat = OTU_2_Sample, type = "spp_site", thresh = TRUE, spp_names = TRUE)
cooccur.OTUs_15perc <- cooccur(mat = OTU_2_Sample_15perc, type = "spp_site", thresh = TRUE, spp_names = TRUE)
cooccur.OTUs_Borneo <- cooccur(mat = OTU_2_Sample_Borneo_15perc, type = "spp_site", thresh = TRUE, spp_names = TRUE)
cooccur.OTUs_FG <- cooccur(mat = OTU_2_Sample_FG_15perc, type = "spp_site", thresh = TRUE, spp_names = TRUE)


#summary number of significant cooccurrences
summary(cooccur.OTUs_15perc)
#summary at 0.025 level
nrow(subset(prob.table(cooccur.OTUs_15perc), p_lt <0.025))
nrow(subset(prob.table(cooccur.OTUs_15perc), p_gt <0.025))
nrow(prob.table(cooccur.OTUs))

#plot matrix
tiff('ALLCR_cooccur_0.025_5perc.tiff', res = 400, height = 1748, width = 1748)
OTU_heatmap <- plot(cooccur.OTUs)
(OTU_heatmap <- OTU_heatmap + theme(
  legend.position = c(0.875,0.5),
  legend.background = element_blank(),
  legend.key.size = unit(.25, 'cm'),
  legend.text = element_text(size=5),
  plot.title = element_blank(),
  plot.margin=margin(-20,-20,-60,0)
))
dev.off()

png('ALLCR_cooccur_0.025_15perc.png', res = 400, height = 1748, width = 1748)
OTU_heatmap_15perc <- plot(cooccur.OTUs_15perc)
(OTU_heatmap_15perc <- OTU_heatmap_15perc + theme(
                                    legend.position = c(0.875,0.5),
                                    legend.background = element_blank(),
                                    legend.key.size = unit(.25, 'cm'),
                                    legend.text = element_text(size=5),
                                    plot.title = element_blank(),
                                    plot.margin=margin(-20,-20,-60,0)
                                    ))
dev.off()


#make effect size matrix with no threshold
eff_sizes <- effect.sizes(cooccur.OTUs, standardized = TRUE, matrix = TRUE)
prob.table(cooccur.OTUs)
eff_sizes_15perc <- effect.sizes(cooccur.OTUs_15perc, standardized = TRUE, matrix = TRUE)
eff_sizes_Borneo <- effect.sizes(cooccur.OTUs_Borneo, standardized = TRUE, matrix = TRUE)
eff_sizes_FG <- effect.sizes(cooccur.OTUs_FG, standardized = TRUE, matrix = TRUE)

#show matrix as network
library(visNetwork)
#make a dataframe of the nodes i.e. species names
nodes <- data.frame(id = 1:nrow(OTU_2_Sample),
                    label = rownames(OTU_2_Sample),
                    color = "grey") 
nodes_15perc <- data.frame(id = 1:nrow(OTU_2_Sample_15perc),
                    label = rownames(OTU_2_Sample_15perc),
                    color = "grey") 

#make a dataframe of the edges i.e. cooccurrences
edges <- data.frame(from = cooccur.OTUs$results$sp1, to = cooccur.OTUs$results$sp2,
                    color = NA)
edges_15perc <- data.frame(from = cooccur.OTUs_15perc$results$sp1, to = cooccur.OTUs_15perc$results$sp2,
                    color = NA)

for(i in 1:nrow(edges)){
  if(cooccur.OTUs$results$p_lt[i] <= 0.025){
    edges$color[i] = '#ffcc66'
  }
  if(cooccur.OTUs$results$p_gt[i] <= 0.025){
    edges$color[i] = '#acd8e6'
  }
  if(cooccur.OTUs$results$p_lt[i] > 0.025 & cooccur.OTUs$results$p_gt[i] > 0.025){
    edges$color[i] = 'white'
  }
}

for(i in 1:nrow(edges_15perc)){
  if(cooccur.OTUs_15perc$results$p_lt[i] <= 0.025){
    edges_15perc$color[i] = '#ffcc66'
  }
  if(cooccur.OTUs_15perc$results$p_gt[i] <= 0.025){
    edges_15perc$color[i] = '#acd8e6'
  }
  if(cooccur.OTUs_15perc$results$p_lt[i] > 0.025 & cooccur.OTUs_15perc$results$p_gt[i] > 0.025){
    edges_15perc$color[i] = 'white'
  }
}
#dont plot random associations
edges <- subset(edges, edges$color != 'white')
edges_15perc <- subset(edges_15perc, edges_15perc$color != 'white')

visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout() %>%
  visSave('test.html')

visNetwork(nodes = nodes_15perc, edges = edges_15perc) %>%
  visIgraphLayout()

write_csv(edges, "ALLCR_cooccurrence_network_edges_5perc.csv")
write_csv(nodes, "ALLCR_cooccurrence_network_nodes_5perc.csv")
write_csv(edges_15perc, "ALLCR_cooccurrence_network_edges_15perc.csv")
write_csv(nodes_15perc, "ALLCR_cooccurrence_network_nodes_15perc.csv")

library(vegan)
library(fpc)
library(cluster)
library(factoextra)
library(cowplot)

#function to calculate the calinski harabasz index to select the best k value for PAM
fviz_ch <- function(data,kvals) {
  ch <- c()
  for (i in kvals) {
    pam <- pam(data, k = i, diss = TRUE) # PAM clustering
    ch[i] <- calinhara(data, # data
                       pam$clustering, # cluster assignments
                       cn=max(pam$clustering) # total cluster number
    )
  }
  ch.k.df <- tibble(kvalues = kvals, chvals = ch)
  return(ch.k.df)
}
#run function with kvalues (num. clusters) being 2 to number of taxa - 1
ch.k.df <- fviz_ch(1 - eff_sizes, c(2:nrow(OTU_2_Sample)-1))
ch.k.df_15perc <- fviz_ch(1 - eff_sizes_15perc, c(2:nrow(OTU_2_Sample_15perc)-1))
ch.k.df_Borneo <- fviz_ch(1 - eff_sizes_Borneo, c(1:29))
ch.k.df_FG <- fviz_ch(1 - eff_sizes_FG, c(1:20))


png("ALLCR_cooccur_0.025_5perc_kvals.png",  res = 400, width = 2480, height = 2480)
(ch.k.df_5perc_plot <- ggplot(ch.k.df, aes(x = kvalues, y = chvals)) +
  geom_col(fill = 'grey') +
  theme_minimal() +
  labs(x="k Clusters", y = "Caliński-Harabasz index"))
dev.off()

png("ALLCR_cooccur_0.025_15perc_kvals.png",  res = 400, width = 2480, height = 2480)
(ch.k.df_15perc_plot <- ggplot(ch.k.df_15perc, aes(x = kvalues, y = chvals)) +
  geom_col(fill = 'grey') +
  theme_minimal() +
    theme(panel.background = element_rect(colour = "grey", fill=NA)) +
  labs(x="k Clusters", y = "Caliński-Harabasz index"))
dev.off()

png("ALLCR_cooccur_0.025_15perc_kvals_Borneo.png",  res = 400, width = 2480, height = 2480)
(ch.k.df_15perc_plot_Borneo <- ggplot(ch.k.df_Borneo, aes(x = kvalues, y = chvals)) +
    geom_col(fill = 'grey') +
    theme_minimal() +
    labs(x="k Clusters", y = "Caliński-Harabasz index"))
dev.off()

png("ALLCR_cooccur_0.025_15perc_kvals_FG.png",  res = 400, width = 2480, height = 2480)
(ch.k.df_15perc_plot_FG <- ggplot(ch.k.df_FG, aes(x = kvalues, y = chvals)) +
    geom_col(fill = 'grey') +
    theme_minimal() +
    labs(x="k Clusters", y = "Caliński-Harabasz index"))
dev.off()

#with the best K visualise the PAM clustering
cluster_output <- pam(1 - eff_sizes, k= ch.k.df[which.max(ch.k.df$chvals),1], diss = TRUE)
cluster_output$data <- 1 - eff_sizes
OTU_cluster <- fviz_cluster(cluster_output,
                            show.legend=FALSE, 
                            labelsize = 0) + 
  scale_colour_manual(values=c("#CC79A7","#009E73")) +
  scale_fill_manual(values=alpha(c("#CC79A7","#009E73"),0.4)) +
  theme(plot.title = element_blank())

blank_gg <- ggplot() + theme_void()

png("ALLCR_cooccur_5perc_multiplot.png",  res = 400, width = 3508, height = 2480)
plot_grid(OTU_heatmap, blank_gg, ch.k.df_5perc_plot, OTU_cluster, labels="auto")
dev.off()

cluster_output_15perc <- pam(1 - eff_sizes_15perc, k= ch.k.df_15perc[which.max(ch.k.df_15perc$chvals),1], diss = TRUE)
cluster_output_15perc$data <- 1 - eff_sizes_15perc
#plot clustering analysis
(OTU_cluster_15perc <- fviz_cluster(cluster_output_15perc, 
                                   repel = TRUE, 
                                   show.legend=FALSE, 
                                   labelsize = 5) + 
  scale_colour_manual(values=c("#FC8D62","#8DA0CB")) +
  scale_fill_manual(values=alpha(c("#FC8D62","#8DA0CB"),0.4)) +
  theme_minimal() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(colour = "grey", fill=NA))
)


png("ALLCR_cooccur_15perc_multiplot.png",  res = 400, width = 3508, height = 2480)
plot_grid(OTU_heatmap_15perc, blank_gg, ch.k.df_15perc_plot, OTU_cluster_15perc, labels="auto")
dev.off()


#Borneo
cluster_output_Borneo <- pam(1 - eff_sizes_Borneo, k= ch.k.df_Borneo[which.max(ch.k.df_Borneo$chvals),1], diss = TRUE)
cluster_output_Borneo$data <- 1 - eff_sizes_Borneo
#plot clustering analysis
OTU_cluster_Borneo <- fviz_cluster(cluster_output_Borneo, 
                                   repel = TRUE, 
                                   show.legend=FALSE, 
                                   labelsize = 5) + 
  scale_colour_manual(values=c("#009E73","#CC79A7")) +
  scale_fill_manual(values=alpha(c("#009E73","#CC79A7"),0.4)) +
  theme(plot.title = element_blank())

#FG
cluster_output_FG <- pam(1 - eff_sizes_FG, k= ch.k.df_FG[which.max(ch.k.df_FG$chvals),1], diss = TRUE)
cluster_output_FG$data <- 1 - eff_sizes_FG
#plot clustering analysis
OTU_cluster_FG <- fviz_cluster(cluster_output_FG, 
                                   repel = TRUE, 
                                   show.legend=FALSE, 
                                   labelsize = 5) + 
  scale_colour_manual(values=c("#CC79A7","#009E73")) +
  scale_fill_manual(values=alpha(c("#CC79A7","#009E73"),0.4)) +
  theme(plot.title = element_blank())


png("ALLCR_PAM_clusters_all_Borneo_FG.png",  res = 400, height = 3508, width = 2480)
plot_grid(OTU_cluster_15perc, OTU_cluster_Borneo, OTU_cluster_FG, labels=c("all","Borneo","FG"), ncol = 1)
dev.off()
