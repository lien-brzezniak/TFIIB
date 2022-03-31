setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
library("gplots")
library(gridExtra)
library(grid)
source('functions.R')



#output_path = '/home/lien/data/TFIIB/Hi_C/TagDirs/matrices/raw_relative/GeneMatrix_allSizes/calc_matrix/'
sample_info <- read_tsv('/home/lien/data/skrypty/TFIIB/samples_GEO_HiC.tsv',col_types='ffffff')
group_sets=list(
  'other'=c('All'))
  #'geneBody'=c(3,12,25,26))
down_thresh=5000
up_thresh=50000
 
#grupy=c('CDS_states','R58C_TrainabilityLoss')
#grupy=c('top_genes','deltaMem','non_trainable')
#grupy=c('single')

samples=c('1','chang_rep1')

Groups <- create_Groups(group_sets,down_thresh,up_thresh,complement_rest=F,exclusive=F)

####********MAIN
####*
####*

output_path='/home/lien/data/skrypty/TFIIB/gMpIw1648136307_FigSupp3AB/'
matrices <-  lapply(names(Groups), function(group){
 matrices <-lapply(samples,function(sample){
        norm_matrix <- read_tsv(paste0(output_path,'Meta_norm_',sample,'_',group,'_',down_thresh/1000,'to',up_thresh/1000,'kb_nonscaled.txt'), col_names=F) %>%
          as.matrix()
         scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
     
 })
 names(matrices) <- samples
 return(matrices)
})
 names(matrices)=names(Groups) 
 raw_matrices <-  lapply(names(Groups), function(group){
 matrices <-lapply(samples,function(sample){
        norm_matrix <- read_tsv(paste0(output_path,'Meta_raw_',sample,'_',group,'_',down_thresh/1000,'to',up_thresh/1000,'kb_nonscaled.txt'), col_names=F) %>%
          as.matrix()
         
     
 })
 names(matrices) <- samples
 return(matrices)
})
 names(raw_matrices)=names(Groups) 
  threshold=0.995
max_norm = quantile(abs((unlist(matrices))),probs=threshold,na.rm=T)
max_raw = quantile(abs((unlist(raw_matrices))),probs=threshold,na.rm=T)

  
  pdf_file=paste0(output_path,paste0(names(Groups),collapse='_'),'singleSamples.pdf')
  pdf(pdf_file)
lapply(names(Groups),function(group){ 

  lapply(samples,function(sample){
    norm_map <- matrices[[group]][[sample]]
     raw_map <- raw_matrices[[group]][[sample]]
    title=paste0(group,' \nThresh: ',threshold,'\nSample ',sample)
    
    f_draw_matrix(norm_array=norm_map,raw_array=raw_map, threshold=1,title=title,do_scaling=F,max_scale_norm=max_norm,max_scale_raw=NULL)
  })
})

  dev.off()

  