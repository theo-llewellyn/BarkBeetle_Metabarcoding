####OTU Clustering####

###closed reference OTU clustering
qiime vsearch cluster-features-closed-reference \
  --i-table ALLdada2-featuretable.qza \
  --i-sequences ALLdada2-repseq.qza \
  --i-reference-sequences UNITEfasta.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table CRotutable.qza \
  --o-clustered-sequences CRotuclusters.qza \
  --o-unmatched-sequences CRotuunmatched.qza

#view CR table and sequences
qiime feature-table tabulate-seqs \
--i-data CRotuclusters.qza \
--o-visualization CRclustersequences.qzv

qiime feature-table tabulate-seqs \
--i-data CRotuunmatched.qza \
--o-visualization CRunmatched.qzv

#filter out unmatched (no picture ID) samples with metadata sheet
qiime feature-table filter-samples \
    --i-table CRotutable.qza \
    --m-metadata-file ALLsequencesmetadata.tsv.txt \
    --o-filtered-table filteredCRotutable.qza

qiime feature-table summarize \
--i-table filteredCRotutable.qza \
--m-sample-metadata-file ALLsequencesmetadata.tsv.txt \
--o-visualization CROTUsummary.qzv
