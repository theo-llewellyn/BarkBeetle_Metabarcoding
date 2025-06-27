#!/bin/sh                                                                       

# Process samples                                                               
conda activate cutadapt

for f in *_R1*.fastq.gz; do

    s=$(cut -d_ -f1,2,3,4 <<< "$f")
    r=${s}_*R2*.fastq.gz

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================
    
    # Commands to demultiplex and remove tags and primers                       
    # using e.g. cutadapt may be added here.         
    echo Removing ITS86F from forward reads   
    cutadapt -g GTGAATCATCGAATCTTTGAA -o $s.R1.fastq.gz $f
    echo Removing ITS4 from reverse reads
    cutadapt -g TCCTCCGCTTATTGATATGC -o $s.R2.fastq.gz $r
    
    echo
    echo Merging FWD and REV

    ~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch --threads 15 \
        --fastq_mergepairs $s.R1.fastq.gz \
        --reverse $s.R2.fastq.gz \
        --fastqout $s.merged.fastq \
        --fastq_eeout                       

    echo
    echo Quality filtering

    ~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch --threads 1 \
        --fastq_filter $s.merged.fastq \
        --fastq_maxee 0.5 \
        --fastq_maxns 0 \
        --sample $s \
        --fastqout $s.filtered.fastq \
        --fasta_width 0

done

rm *.R2.fastq.gz
rm *.R1.fastq.gz
rm *.merged.fastq
                     
echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

cat *.filtered.fastq > all.fq

echo
echo Dereplicate across samples and remove singletons

~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch --threads 1 \
 --fastx_uniques all.fq \
 --sizein \
 --sizeout \
 --fastqout all_derep.fq

echo Unique non-singleton sequences: $(grep -c "^@" all_derep.fq) >> log.txt

echo
echo remove seqs with less than 4 reads to remove very low abundance
echo cluster/denoise into 0pc ASVs or zOTUs
echo 

~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch \
 --cluster_unoise all_derep.fq \
 --sizein --sizeout \
 --minsize 4 \
 --centroids all_denoised.fa

echo Unique sequences after preclustering: $(grep -c "^>" all_denoised.fa) >> log.txt

echo
echo De novo chimera detection

~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch \
 --uchime3_denovo all_denoised.fa \
 --sizein \
 --sizeout \
 --fasta_width 0 \
 --nonchimeras all_nonchimeras.fa

echo Unique sequences after chimera detection: $(grep -c "^>" all_nonchimeras.fa) >> log.txt

echo
echo Generate ASV table
echo convert fastq to fasta
~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch \
 --fastq_filter all.fq \
 -fastaout all.fa
 
echo replace space in fasta header
sed -i 's/ /_/g' all.fa

echo
echo Generate table
~/software/vsearch-2.30.0-linux-x86_64/bin/vsearch \
 --usearch_global all.fa \
 --db all_nonchimeras.fa \
 --id 1.0 \
 --sizein \
 --sizeout \
 --threads 32 \
 --otutabout ASV_frequency_table.tsv
