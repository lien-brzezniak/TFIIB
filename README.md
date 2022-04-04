#### This repository contains the code used to reproduce the figures.
#### The initial data required to run the scripts can be downloaded from NCBI Gene Expression Omnibus (GEO; https://www.ncbi.nlm.nih.gov/geo/), accession number: GSE199058

#### functions.R contain the functions called in the scripts.
#### sample_info_3RNAseq.tsv contains metadata describing samples analysed in 3'RNA-seq experiment
#### sample_GEO_HiC.tsv contains metadata describing samples analysed in Hi-C experiment
#### reproduce_*.Rmd are markdown files with R code chunks that generate figures presented in the manuscript
#### *.html contain reports from the execution of .Rmd files
#### protein_coding_genes.bed contains annotation for protein coding genes of A. thaliana according to tair10
#### PCSD_allData.RDS contain lists of genes with various chromatin state enrichments, described in Liu Y, et al., 2018. PCSD: a plant chromatin state database. Nucleic Acids Res 46: D1157–D1167.
#### drought_genes.tsv - table with trainability index for genes induced by drought, calculated based on 3'RNAseq data
#### calculate_interaction_freq_100bp - a description on how interaction frequencies were calculated for the distance normalization 
####  
