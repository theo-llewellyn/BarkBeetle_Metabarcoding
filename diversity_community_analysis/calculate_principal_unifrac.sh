
#convert feature table to biom format
biom convert -i ALLCRfiltered_feature-table_trap_collapse.tsv -o ALLCRfiltered_feature-table_trap_collapse.biom --table-type "Table" --to-hdf5

#import table to qiime
qiime tools import --input-path ALLCRfiltered_feature-table_trap_collapse.biom --type FeatureTable[Frequency] --output-path ALLCRfiltered_feature-table_trap_collapse.qza

#calculate unifrac distances for principals using iqtree phylogeny
qiime diversity core-metrics-phylogenetic --i-phylogeny iqtree.qza --i-table ALLCRfiltered_feature-table_trap_collapse.qza --p-sampling-depth 1000 --output-dir IQtreeCRcoremetrics_trap_collapse --m-metadata-file trap_collapse_metadata.tsv
