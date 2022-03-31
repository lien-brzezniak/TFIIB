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
 
 
#grupy=c('CDS_states','R58C_TrainabilityLoss')
#grupy=c('top_genes','deltaMem','non_trainable')
#grupy=c('single')

samples=c(1:12)



####********MAIN
####*
####*

barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'/')
system(paste0('mkdir ',output_path),intern=F)

res = 100
reach = 1500
down_thresh=3000
up_thresh=50000
write_matrix_file=T
Groups <- create_Groups(group_sets,down_thresh,up_thresh,complement_rest=F,exclusive=F)

metamatrix_data <- lapply(1:length(Groups),function(G){
 
  lapply(samples,function(s){
    
   genes=Groups[[G]]
  result <- calculate_metamatrix(genes,s,reach,res)
  result$sample=paste0('Sample',s)
  result$name <- unite(sample_info,c(2:4),col='name',sep='_' )[which(sample_info$sample==result$sample),] %>% pull(name)
  result$group=names(Groups)[G]
 
  if(write_matrix_file){
    title=paste0(result$sample,'_',result$group,'_',down_thresh/1000,' kb to ',up_thresh/1000,' kb.matrix')
     write.table(result$mat_raw, file=paste0(output_path,'Meta_raw_',title), sep='\t', row.names = F, col.names =  F)
     write.table(result$mat_norm, file=paste0(output_path,'Meta_norm_',title), sep='\t', row.names = F, col.names =  F)

  }
   return(result)
  })
}) %>% reduce(c)

 pdf_file=paste0(output_path,paste0(names(Groups),collapse='_'),'.pdf')
 pdf(pdf_file)
lapply(metamatrix_data, function(x){

  title=paste0(x$sample,'   ',x$name,'\n',x$group,'; ',down_thresh/1000,' kb to ',up_thresh/1000,' kb\n',x$nr_of_genes,' genes')
     f_draw_matrix(x$mat_raw, x$mat_norm,do_scale=T,title,threshold=0.995)
 
})
dev.off()

######draw deltas from saved files
output_path='/home/lien/data/skrypty/TFIIB/aaAib1647967542Fig4C/'
matrices <-  lapply(names(Groups), function(group){
 lapply(samples,function(sample){

        norm_matrix <- read_tsv(paste0(output_path,'Meta_norm_Sample',sample,'_',group,'_',down_thresh/1000,' kb to ',up_thresh/1000,' kb.matrix'), col_names=F) %>%
          as.matrix()
         scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
     
 })
})
 names(matrices)=names(Groups) 
  threshold=0.995
 
max_norm = quantile(abs((unlist(matrices))),probs=threshold,na.rm=T)

  
  pdf_file=paste0(output_path,paste0(names(Groups),collapse='_'),'singleSamples.pdf')
  pdf(pdf_file)
lapply(names(Groups),function(group){ 

  lapply(samples,function(sample){
    norm_map <- matrices[[group]][[sample]]
    title=paste0(group,' \nThresh: ',threshold,'\nSample ',sample)
    
    f_draw_matrix(norm_array=norm_map, title=title, do_scaling=F,max_scale_norm=max_norm)
  })
})

  dev.off()
