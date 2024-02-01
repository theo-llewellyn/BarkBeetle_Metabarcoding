###Demultiplex data###

#move all R1 & R2 of each sample to one folder 
#rename ****_1.fq.gz to ****_R1.fq.gz and ****_2.fq.gz to ****_R2.fq.gz (for var in *1.fq.gz; do mv "$var" "${var%1.fq.gz}R1.fq.gz"; done); 
#generate a demux table 
#copy muldemux.py to your folder; 
#run 'python3 muldemux.py -i 0_rawdata/*.fq.gz -d EUQN_235705_02_demux.tsv -o 0_demux/ -s 0_demux_outcounts.tsv' (0_rawdata is  #your rawdata folder, EUQN_235705_02_demux.tsv is your demux table, 0_demux is the output folder).

python3 demultiplex_reads.py -i *.fq.gz -d UKdemux.tsv -o 0_demuxUK/ -s 0_demux_outcountsUK.tsv

###input data
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path combinedmanifest_new.csv \
--output-path ALLimportedsequences.qza \
--input-format PairedEndFastqManifestPhred33

###view results of demux
qiime demux summarize \
--i-data ALLimportedsequences.qza \
--o-visualizationALLdemux.qzv
