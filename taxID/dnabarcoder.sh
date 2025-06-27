cd RESULTS/ITSx 

~/bin/ITSx_1.1.3/ITSx \
 -i ../ALLimportedsequences/all_nonchimeras.fa \
 --saveregions ITS2 \
 -o all_nochimeras_formatted \
 --cpu 8 \
 --detailed_results T \
 --preserve T
 
#get list of only Fungal ITSx sequences
cat all_nochimeras_formatted.extraction.results | awk '{if ($3 == "F") print $1;}' > Fungal_ITS_headers.txt

#subset fasta file for just Fungal ITS2
cat all_nochimeras_formatted.ITS2.fasta | seqkit grep -f Fungal_ITS_headers.txt > all_nochimeras_formatted.Fungi.ITS2.fasta

cd ../dnabarcoder
source activate dnabarcoder-env

#search for best matches against UNITE ITS2 db
~/miniconda3/envs/dnabarcoder-env/bin/dnabarcoder/dnabarcoder.py search \
 -i all_nochimeras_formatted.Fungi.ITS2.fasta \
 -r unite2024ITS2.fasta -ml 50

#assign sequences to diff taxonomic groups
~/miniconda3/envs/dnabarcoder-env/bin/dnabarcoder/dnabarcoder.py classify \
 -i dnabarcoder/all_nochimeras_formatted.Fungi.ITS2.unite2024ITS2_BLAST.bestmatch \
 -c unite2024ITS2.classification \
 -cutoffs unite2024ITS2.unique.cutoffs.best.json

#Need to edit a the dnabarcoder .classification file and ASV_frequency_table to fit the format. See the .dynamic.format files.
#also need to convert the unite json file into a csv with the headers TaxonomicName,TaxonRank,Cut-off


#remove size info from dynamic format files to match the ASV table
sed 's/;size.*//g' ../ITSx/all_nochimeras_formatted.Fungi.ITS2.fasta > ../ITSx/all_nochimeras_formatted.Fungi.ITS2.nosize.fasta
sed 's/;size=[0-9A-Za-z]*//g' dnabarcoder/all_nochimeras_formatted.Fungi.ITS2.unite2024ITS2_BLAST.classification.dynamic.format > dnabarcoder/all_nochimeras_formatted.Fungi.ITS2.unite2024ITS2_BLAST.classification.dynamic.format.nosize

#dynamic clustering based on dnabarcoer cutoffs
#run RStudio from conda environment
conda activate dynamic_clustering
/Applications/RStudio.app/Contents/MacOS/RStudio
