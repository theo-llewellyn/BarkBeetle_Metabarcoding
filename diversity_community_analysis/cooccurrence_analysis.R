#------------------------------------------------------------------------------------#
# Community Co-occurrence analysis

library(cooccur)
library(tidyverse)
library(reshape2)

source(file = '../R_SCRIPTS/plot.coocur.R')

###################################################################################################################
#read and format data

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

#get sample metadata
sample_data <- read_csv('Sco_Pla_FG_Borneo_metadata.csv')


##Borneo
#if we want the same taxa in the full analysis
OTU_2_Sample_Borneo <- OTU_2_Sample[, names(OTU_2_Sample) %in% 
                                      subset(sample_data, sample_data$Country == 'Malaysia')$fungal_metabarcode_ID]
#remove OTUs in less than 10 samples in Borneo
OTU_2_Sample_Borneo <- OTU_2_Sample_Borneo[-which(rowSums(OTU_2_Sample_Borneo[sapply(OTU_2_Sample_Borneo, is.numeric)]) < 10),]

#French Guiana
OTU_2_Sample_FG <- OTU_2_Sample[, names(OTU_2_Sample) %in% 
                                  subset(sample_data, sample_data$Country == 'French Guiana')$fungal_metabarcode_ID]
#remove OTUs in less than 10 samples in FG
OTU_2_Sample_FG <- OTU_2_Sample_FG[-which(rowSums(OTU_2_Sample_FG[sapply(OTU_2_Sample_FG, is.numeric)]) < 10),]


###################################################################################################################
set.seed(071123)
cooccur.OTUs_Borneo <- cooccur(mat = OTU_2_Sample_Borneo, type = "spp_site", thresh = TRUE, spp_names = TRUE)
cooccur.OTUs_FG <- cooccur(mat = OTU_2_Sample_FG, type = "spp_site", thresh = TRUE, spp_names = TRUE)


#summary number of significant cooccurrences Borneo
summary(cooccur.OTUs_Borneo)
#summary at 0.025 level
nrow(subset(prob.table(cooccur.OTUs_Borneo), p_lt <0.025))
nrow(subset(prob.table(cooccur.OTUs_Borneo), p_gt <0.025))
nrow(prob.table(cooccur.OTUs_Borneo))

#summary number of significant cooccurrences Borneo
summary(cooccur.OTUs_FG)
#summary at 0.025 level
nrow(subset(prob.table(cooccur.OTUs_FG), p_lt <0.025))
nrow(subset(prob.table(cooccur.OTUs_FG), p_gt <0.025))
nrow(prob.table(cooccur.OTUs_FG))

###################################################################################################################
#make effect size matrix with no threshold
eff_sizes_Borneo <- effect.sizes(cooccur.OTUs_Borneo, standardized = TRUE, matrix = TRUE)
eff_sizes_FG <- effect.sizes(cooccur.OTUs_FG, standardized = TRUE, matrix = TRUE)

#show matrix as network
library(visNetwork)
#make a dataframe of the nodes i.e. species names
nodes_Borneo <- data.frame(id = 1:nrow(OTU_2_Sample_Borneo),
                           label = rownames(OTU_2_Sample_Borneo),
                           color = "grey") 
nodes_FG <- data.frame(id = 1:nrow(OTU_2_Sample_FG),
                       label = rownames(OTU_2_Sample_FG),
                       color = "grey") 
#make a dataframe of the edges i.e. cooccurrences
edges_Borneo <- data.frame(from = cooccur.OTUs_Borneo$results$sp1, 
                           to = cooccur.OTUs_Borneo$results$sp2,
                           color = NA)
edges_FG <- data.frame(from = cooccur.OTUs_FG$results$sp1, 
                       to = cooccur.OTUs_FG$results$sp2,
                       color = NA)

for(i in 1:nrow(edges_Borneo)){
  if(cooccur.OTUs_Borneo$results$p_lt[i] <= 0.025){
    edges_Borneo$color[i] = '#ffcc66'
  }
  if(cooccur.OTUs_Borneo$results$p_gt[i] <= 0.025){
    edges_Borneo$color[i] = '#acd8e6'
  }
  if(cooccur.OTUs_Borneo$results$p_lt[i] > 0.025 & cooccur.OTUs_Borneo$results$p_gt[i] > 0.025){
    edges_Borneo$color[i] = 'white'
  }
}

for(i in 1:nrow(edges_FG)){
  if(cooccur.OTUs_FG$results$p_lt[i] <= 0.025){
    edges_FG$color[i] = '#ffcc66'
  }
  if(cooccur.OTUs_FG$results$p_gt[i] <= 0.025){
    edges_FG$color[i] = '#acd8e6'
  }
  if(cooccur.OTUs_FG$results$p_lt[i] > 0.025 & cooccur.OTUs_FG$results$p_gt[i] > 0.025){
    edges_FG$color[i] = 'white'
  }
}
#dont plot random associations
edges_Borneo <- subset(edges_Borneo, edges_Borneo$color != 'white')
edges_FG <- subset(edges_FG, edges_FG$color != 'white')

#save node info as table
write_csv(nodes_Borneo, "cooccurrence/ALLCR_cooccurrence_network_nodes_Borneo.csv")
write_csv(nodes_FG, "cooccurrence/ALLCR_cooccurrence_network_nodes_FG.csv")

#save the effect sizes as list of pairs
eff_size_pairs_Borneo <- as.data.frame(as.table(as.matrix(eff_sizes_Borneo)))
eff_size_pairs_FG <- as.data.frame(as.table(as.matrix(eff_sizes_FG)))
#change variable names to the OTU codes used in the nodes and edges tables
eff_size_pairs_Borneo$Var1 <- rep(c(1:nrow(nodes_Borneo)),nrow(nodes_Borneo))
eff_size_pairs_Borneo$Var2 <- rep(1:nrow(nodes_Borneo),each=nrow(nodes_Borneo))
eff_size_pairs_FG$Var1 <- rep(c(1:nrow(nodes_FG)),nrow(nodes_FG))
eff_size_pairs_FG$Var2 <- rep(1:nrow(nodes_FG),each=nrow(nodes_FG))
#merge with edges table to have a column with the effect size for each edge

edges_Borneo_ES <- left_join(edges_Borneo,eff_size_pairs_Borneo, by = c('from'='Var1', 'to'='Var2'))
edges_FG_ES <- left_join(edges_FG,eff_size_pairs_FG, by = c('from'='Var1', 'to'='Var2'))

#save edges tables
write_csv(edges_Borneo_ES, "cooccurrence/ALLCR_cooccurrence_network_edges_Borneo.csv")
write_csv(edges_FG_ES, "cooccurrence/ALLCR_cooccurrence_network_edges_FG.csv")

################################################################################
# Network plotting and analysis

library(vegan)
library(fpc)
library(cluster)
library(factoextra)
library(cowplot)
library(igraph)

#make network with igraph
graph <- graph_from_data_frame(subset(edges_Borneo_ES, edges_Borneo_ES$color == '#acd8e6'))
#calculate communities with edge betweenness method
community_modules <- cluster_edge_betweenness(graph, weights = E(graph)$Freq)
modularityBorneo <- modularity(graph, membership = community_modules$membership)
transitivityBorneo <- transitivity(graph,type = 'weighted',weights = E(graph)$Freq)
betweennessBorneo <- betweenness(graph, weights = E(graph)$Freq)
degreeBorneo <- degree(graph)

#make FG network with igraph
graph_FG <- graph_from_data_frame(subset(edges_FG_ES, edges_FG_ES$color == '#acd8e6'))
#calculate communities with edge betweenness method
community_modules_FG <- cluster_edge_betweenness(graph_FG, weights = E(graph_FG)$Freq)
modularityFG <- modularity(graph_FG, membership = community_modules_FG$membership)
transitivityFG<-transitivity(graph_FG,type = 'weighted',weights = E(graph_FG)$Freq)
betweennessFG<-betweenness(graph_FG, weights = E(graph_FG)$Freq)
degreeFG<-degree(graph_FG)

#visualise and calculate average network metrics for both countries
transFG<-data.frame(transitivityFG,rep('FG',length(transitivityFG)))
colnames(transFG)<-c('transitivity','locality')
transBorneo<-data.frame(transitivityBorneo,rep('Borneo',length(transitivityBorneo)))
colnames(transBorneo)<-c('transitivity','locality')
trans<-rbind(transFG,transBorneo)
transplot <- ggplot(trans,aes(x=locality,y=transitivity)) + 
  geom_boxplot()
betweenFG<-data.frame(betweennessFG,rep('FG',length(betweennessFG)))
colnames(betweenFG)<-c('betweenness','locality')
betweenBorneo<-data.frame(betweennessBorneo,rep('Borneo',length(betweennessBorneo)))
colnames(betweenBorneo)<-c('betweenness','locality')
betweenness<-rbind(betweenFG,betweenBorneo)
betweennessplot <- ggplot(betweenness,aes(x=locality,y=betweenness)) + 
  geom_boxplot()
degrFG<-data.frame(degreeFG,rep('FG',length(degreeFG)))
colnames(degrFG)<-c('degree','locality')
degrBorneo<-data.frame(degreeBorneo,rep('Borneo',length(degreeBorneo)))
colnames(degrBorneo)<-c('degree','locality')
degree<-rbind(degrFG,degrBorneo)
degreeplot <- ggplot(degree,aes(x=locality,y=degree)) + 
  geom_boxplot()

plot_grid(transplot,betweennessplot,degreeplot,nrow=3, ncol = 1,  align = 'v', align.axis = TRUE, labels="auto")

networkstats <- list(transitivityFG,transitivityBorneo,
                     betweennessFG,betweennessBorneo,
                     degreeFG,degreeBorneo)
names(networkstats) <- c('transitivityFG','transitivityBorneo','betweennessFG','betweennessBorneo','degreeFG','degreeBorneo')
for(i in 1:length(networkstats)){
  print(paste('mean',names(networkstats)[i],'=',mean(na.omit(networkstats[[i]])), sep = ' '))
}

#test if localities differ
wilcox.test(x = transitivityFG, y = transitivityBorneo, paired = FALSE, alternative = "two.sided")
wilcox.test(x = betweennessFG, y = betweennessBorneo, paired = FALSE, alternative = "two.sided")
wilcox.test(x = degreeFG, y = degreeBorneo, paired = FALSE, alternative = "two.sided")

###############################################################################
#plot negative cooccurrences
graph_neg <- graph_from_data_frame(subset(edges_Borneo_ES, edges_Borneo_ES$color == '#ffcc66'))
community_modules_neg <- cluster_edge_betweenness(graph_neg, weights = -(E(graph_neg)$Freq))
#plot negative cooccurrences FG
graph_neg_FG <- graph_from_data_frame(subset(edges_FG_ES, edges_FG_ES$color == '#ffcc66'))
community_modules_neg_FG <- cluster_edge_betweenness(graph_neg_FG, weights = -(E(graph_neg_FG)$Freq))


##FUNCTION TO PLOT NETWORK LAYOUT WITH MORE DISTANCE BETWEEN THAN WITHIN COMMUNITIES
weight.community=function(row,membership,weigth.within,weight.between){
  if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
    weight=weigth.within
  }else{
    weight=weight.between
  }
  return(weight)
}
###############################################################################

png('FIGURES/Borneo_FG_networks_posneg.png', res = 400, height = 2400, width = 2400)
par(mar=c(0,0,0,0), mai = c(0, 0, 0, 0))
layout(matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2), 
       heights = c(2, 1),
       widths = c(2, 2)) 

# remove singletons
singletons <- which(table(community_modules$membership) < 2)
community_modules1<-community_modules
community_modules1$membership[community_modules1$membership %in% singletons] <- 999
#layout with weighting by communities
E(graph)$weight=apply(get.edgelist(graph),1,weight.community,membership(community_modules),10,1)
graph$layout=layout.fruchterman.reingold(graph,weights=E(graph)$weight)
#make all singletons as transparent so the polygon doesnt show
colvec <- c('#F5B84C','#E7A7B5','#6D6D6D')[as.numeric(as.factor(community_modules1$membership))]

#colour the edges by communities. if they are within communities then blue if not then grey
test_table <- data.frame(community_modules1$names,community_modules1$membership)
colnames(test_table) <- c('V1','community')
test_table2<-left_join(as.data.frame(get.edgelist(graph)),test_table)
test_table2<- left_join(test_table2,test_table, by=c('V2'='V1'))
for(i in 1:nrow(test_table2)){
  test_table2$colour[i] <- ifelse(test_table2$community.x[i]==999|test_table2$community.y[i]==999, 
                                  '#6D6D6D66', '#acd8e666')
}

plot(community_modules1, graph, 
     col=colvec,
     vertex.size=2,
     vertex.frame.color='#00000000',
     edge.arrow.size=0, arrow.mode = 0, edge.curved=.1,
     vertex.label=NA, vertex.size = 1,
     edge.color=test_table2$colour,
     mark.col = '#00000000', 
     mark.border = '#00000000')
title("Borneo", line = -1.5, font.main = 1)
title("a", line = -1, adj = 0)

#######################
#FG
# remove singletons
singletons <- which(table(community_modules_FG$membership) < 2)
community_modules1<-community_modules_FG
community_modules1$membership[community_modules1$membership %in% singletons] <- 999
#layout with weighting by communities
E(graph_FG)$weight=apply(get.edgelist(graph_FG),1,weight.community,membership(community_modules_FG),5,1)
graph_FG$layout=layout.fruchterman.reingold(graph_FG,weights=E(graph_FG)$weight)
#make all singletons as transparent so the polygon doesnt show
colvec <- c('#F5B84C','#E7A7B5','#6D6D6D')[as.numeric(as.factor(community_modules1$membership))]

#colour the edges by communities. if they are within communities then blue if not then grey
test_table <- data.frame(community_modules1$names,community_modules1$membership)
colnames(test_table) <- c('V1','community')
test_table2<-left_join(as.data.frame(get.edgelist(graph_FG)),test_table)
test_table2<- left_join(test_table2,test_table, by=c('V2'='V1'))
for(i in 1:nrow(test_table2)){
  test_table2$colour[i] <- ifelse(test_table2$community.x[i]==999|test_table2$community.y[i]==999, 
                                  '#6D6D6D66', '#acd8e666')
}

plot(community_modules1, graph_FG, 
     col=colvec,
     vertex.size=2,
     vertex.frame.color='#00000000',
     edge.arrow.size=0, arrow.mode = 0, edge.curved=.1,
     vertex.label=NA, vertex.size = 1,
     edge.color=test_table2$colour,
     mark.col = '#00000000',
     mark.border = '#00000000')
title("French Guiana", line = -1.5, font.main=1)
title("b", line = -1, adj = 0)

##############################################################################################################
#negative plots
#remove singleton
singletons <- which(table(community_modules_neg$membership) < 2)
community_modules1<-community_modules_neg
community_modules1$membership[community_modules1$membership %in% singletons] <- 999
#layout with weighting by communities
E(graph_neg)$weight=apply(get.edgelist(graph_neg),1,weight.community,membership(community_modules_neg),10,1)
graph_neg$layout=layout.fruchterman.reingold(graph_neg,weights=E(graph_neg)$weight)
#list of colours for communities
colvec <- c(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'#6D6D6D')[as.factor(community_modules1$membership)]

#colour the edges by communities. if they are within communities then blue if not then grey
test_table <- data.frame(community_modules1$names,community_modules1$membership)
colnames(test_table) <- c('V1','community')
test_table2<-left_join(as.data.frame(get.edgelist(graph_neg)),test_table)
test_table2<- left_join(test_table2,test_table, by=c('V2'='V1'))
for(i in 1:nrow(test_table2)){
  test_table2$colour[i] <- ifelse(test_table2$community.x[i]==test_table2$community.y[i], 
                                  '#ffcc6666','#6D6D6D66')
  if(test_table2$community.x[i]==999|test_table2$community.y[i]==999){
    test_table2$colour[i] <- '#6D6D6D66'}
}

plot(community_modules1, graph_neg, 
     col=colvec,
     vertex.size=3,
     vertex.frame.color='#00000000',
     edge.arrow.size=0, arrow.mode = 0, edge.curved=.1,
     vertex.label=NA, vertex.size = 1,
     edge.color=test_table2$colour,
     mark.col = c(paste(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'66',sep=''),'#00000000'),
     mark.border = c(paste(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'66',sep=''),'#00000000'))
title("c", line = -1, adj = 0)

#FG negative plot
# remove singletons
singletons <- which(table(community_modules_neg_FG$membership) < 2)
community_modules1<-community_modules_neg_FG
community_modules1$membership[community_modules1$membership %in% singletons] <- 999
#layout with weighting by communities
E(graph_neg_FG)$weight=apply(get.edgelist(graph_neg_FG),1,weight.community,membership(community_modules_neg_FG),2,1)
graph_neg_FG$layout=layout.fruchterman.reingold(graph_neg_FG,weights=E(graph_neg_FG)$weight)
#list of colours for communities
colvec <- c(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'#6D6D6D')[as.factor(community_modules1$membership)]

#colour the edges by communities. if they are within communities then blue if not then grey
test_table <- data.frame(community_modules1$names,community_modules1$membership)
colnames(test_table) <- c('V1','community')
test_table2<-left_join(as.data.frame(get.edgelist(graph_neg_FG)),test_table)
test_table2<- left_join(test_table2,test_table, by=c('V2'='V1'))
for(i in 1:nrow(test_table2)){
  test_table2$colour[i] <- ifelse(test_table2$community.x[i]==test_table2$community.y[i], 
                                  '#ffcc6666','#6D6D6D66')
  if(test_table2$community.x[i]==999|test_table2$community.y[i]==999){
    test_table2$colour[i] <- '#6D6D6D66'}
}

plot(community_modules1, graph_neg_FG, 
     col=colvec,
     vertex.size=3,
     vertex.frame.color='#00000000',
     edge.arrow.size=0, arrow.mode = 0, edge.curved=.1,
     vertex.label=NA, vertex.size = 1,
     edge.color=test_table2$colour,
     mark.col = c(paste(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'66',sep=''),'#00000000'),
     mark.border = c(paste(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(community_modules1$membership))-1),'66',sep=''),'#00000000'))
title("d", line = -1, adj = 0)

dev.off()





#################################################################################################################################
