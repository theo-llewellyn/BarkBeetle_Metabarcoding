#!/bin/bash
###IQ tree####

conda activate qiime2-amplicon-2023.9

#import IQTree phylogeny into qiime data format
qiime tools import \
  --input-path CRotus_3820T_untrimmed_guidance.treefile \
  --output-path iqtree.qza \
  --type 'Phylogeny[Rooted]'


#untrimmed IQ stats#

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny iqtree.qza \
  --i-table CRfilteredfeaturetable.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file nomalaysiametadata.tsv \
  --output-dir IQtreeCRcoremetrics

  #high cut off##
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny iqtree.qza \
  --i-table CRfilteredfeaturetable.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file nomalaysiametadata.tsv \
  --output-dir IQtreehighcutoffCRcoremetrics

#phylogenetic alpha diversity
qiime diversity alpha-phylogenetic \
  --i-table CRfilteredfeaturetable.qza \
  --i-phylogeny iqtree.qza \
  --p-metric faith_pd \
  --o-alpha-diversity IQtreeCRfaith_pd_vector.qza 

##visualization (Kruskal Wallis)
qiime diversity alpha-group-significance \
--i-alpha-diversity IQtreeCRfaith_pd_vector.qza \
--m-metadata-file ALLsequencesmetadata.tsv.txt \
--o-visualization IQtreeCRfaith-pd-group-significance.qzv
