###fasttree script###

###generate a tree for phylogenetic diversity tests

sed 's/:/_/g' ../MAFFT/OTUs.msa.fa > ../MAFFT/OTUs.msa.FastTree.fa
~/bin/FastTree -gtr -nt -gamma ../MAFFT/OTUs.msa.FastTree.fa > FastTree_OTUs.tree
