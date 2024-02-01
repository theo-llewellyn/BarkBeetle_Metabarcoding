###fasttree script###

###generate a tree for phylogenetic diversity tests (https://cryptick-lab.github.io/NGS-Analysis/_site/QIIME2-DiversityPhylogeny.html)

###### CR OTU #########
qiime alignment mafft \
--i-sequences CRotuclusters.qza \
--o-alignment CRaligned-sequences.qza

qiime alignment mask \
--i-alignment CRaligned-sequences.qza \
--o-masked-alignment CRmasked-aligned-sequences.qza

qiime phylogeny fasttree \
--i-alignment CRmasked-aligned-sequences.qza \
--o-tree CRunrooted-tree.qza

qiime phylogeny midpoint-root \
--i-tree CRunrooted-tree.qza \
--o-rooted-tree CRrooted-tree.qza

#to export tree
qiime tools export \
  --input-path CRrooted-tree.qza \
  --output-path CRexported-rooted-tree/
