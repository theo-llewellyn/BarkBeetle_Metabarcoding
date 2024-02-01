#Taxonomic Assignment
#download fungal database https://doi.plutof.ut.ee/doi/10.15156/BIO/2483915 
#unzip downloaded reference sequence database
tar -xvf UNITEreferencesequences.tgz
tar -xf UNITEreferencesequences.tgz -C .

#import database
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path sh_refs_qiime_ver9_dynamic_29.11.2022.fasta \
--output-path UNITEfasta.qza

#Taxonomic classification, using the results from fit-classifier-naive-bayes
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path sh_taxonomy_qiime_ver9_dynamic_29.11.2022.txt \
--output-path UNITEtaxonomy.qza

#classify-sklearn

#Training the Bayesian classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads UNITEfasta.qza \
--i-reference-taxonomy UNITEtaxonomy.qza \
--o-classifier bayesian_classifier.qza

#Applying the trained algorithm to my data
####closed reference OTU####
qiime feature-classifier classify-sklearn \
--i-reads CRotuclusters.qza \
--i-classifier bayesian_classifier.qza \
--o-classification ALLCRtaxonomic_classifications_sklearn.qza

###Filter out unassigned sequences from table
qiime taxa filter-table \
--i-table CRotutable.qza \
--i-taxonomy ALLCRtaxonomic_classifications_sklearn.qza \
--p-exclude Unassigned \
--p-mode contains \
--o-filtered-table ALLCRfilteredfeaturetable.qza \

#view filtered table
qiime metadata tabulate \
--m-input-file ALLCRfilteredfeaturetable.qza \
--o-visualization ALLCRfilteredfeaturetable.qzv

#filter out unmatched (no picture ID) samples with metadata sheet
qiime feature-table filter-samples \
    --i-table ALLCRfilteredfeaturetable.qza \
    --m-metadata-file ALLsequencesmetadata.tsv.txt \
    --o-filtered-table ALLCRfinalfilteredtable.qza

###Filtered visualization
qiime taxa barplot \
--i-table ALLCRfinalfilteredtable.qza \
--i-taxonomy ALLCRtaxonomic_classifications_sklearn.qza \
--m-metadata-file ALLsequencesmetadata.tsv.txt \
--o-visualization ALLCRbarplot_sklearn.qzv
