###DADA2 Denoising (and Chimera detection)###

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ALLimportedsequences.qza \
--p-trim-left-f 60 \
--p-trunc-len-f 0 \
--p-trim-left-r 60 \
--p-trunc-len-r 0 \
--o-representative-sequences ALLdada2-repseq.qza \
--o-table ALLdada2-featuretable.qza \
--o-denoising-stats ALLdada2-stats.qza --verbose

#Checking denoising statistics
qiime metadata tabulate \
--m-input-file ALLdada2-stats.qza \
--o-visualization ALLdada2-stats.qzv

qiime feature-table tabulate-seqs \
  --i-data ALLdada2-repseq.qza \
  --o-visualization ALLpostdadalength-stats.qzv

##### view in viewqiime2.org
qiime metadata tabulate \
--m-input-file ALLdada2-featuretable.qza \
--o-visualization ALLdada2featuretable.qzv

###### view in viewqiime2.org
qiime feature-table tabulate-seqs \
--i-data ALLdada2-repseq.qza \
--o-visualization ALLdada2sequences.qzv

#filter out unmatched samples (no picture ID in metadata) with metadata sheet
qiime feature-table filter-samples \
    --i-table ALLdada2-featuretable.qza \
    --m-metadata-file ALLsequencesmetadata.tsv.txt \
    --o-filtered-table filteredALLdata2-featuretable.qza

qiime feature-table summarize \
--i-table filteredALLdata2-featuretable.qza \
--m-sample-metadata-file ALLsequencesmetadata.tsv.txt \
--o-visualization ALLsequencesummary.qzv
