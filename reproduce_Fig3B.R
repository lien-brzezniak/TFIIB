# set path to repository directory
setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
source('functions.R')

#specify the directory, where GMatrix.RDS files are located:
dir='/home/lien/data/TFIIB/GEO_upload/HiC/'

# provide the list of gene AGI numbers
genes=c('AT1G05490')

#data are converted from upper-triangle sparse format to a symmetric matrix
gene_matrixes <- singleGene_map(sample='1',genes=genes, format='raw')
gene_matrixes[[1]][1:10,1:10]

# it is possible to save the results as pdf:
# barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
# output_path = paste0(barcode,'/')
# system(paste0('mkdir ',output_path),intern=F)

 # pdf_file=paste0(output_path,paste0(names(gene_matrix),collapse='_'),'.pdf')
 # pdf(pdf_file)

 lapply(genes,function(GeneID){
    f_draw_matrix(raw_array=gene_matrixes[[GeneID]],title=GeneID,threshold=0.995)
  })

# dev.off()
