setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
# library("gplots")
# library(gridExtra)
# library(grid)
source('functions.R')


write_matrix_file=T
#output_path = '/home/lien/data/TFIIB/Hi_C/TagDirs/matrices/raw_relative/GeneMatrix_allSizes/calc_matrix/'
sample_info <- read_tsv('/home/lien/data/skrypty/TFIIB/samples_GEO_HiC.tsv',col_types='ffffff')
group_sets=list('geneBody'=c(3,12,25,26))
#group_sets=list('other'=c('R58C_TrainabilityLoss'))

s='1'

dir='/home/lien/data/TFIIB/GEO_upload/HiC/'
barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'/')
system(paste0('mkdir ',output_path),intern=F)

down_thresh=3000
up_thresh=50000
res = 100
reach = down_thresh/2

Groups <- create_Groups(group_sets,down_thresh,up_thresh,complement_rest=F,exclusive=T)
metamatrix_data <- lapply(1:length(Groups),function(G){
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


threshold=0.99
max_norm = quantile(abs(scale(unlist(lapply(metamatrix_data,function(x){x$'mat_norm'})))),probs=threshold,na.rm=T)

 pdf_file=paste0(output_path,paste0(unlist(group_sets),collapse='_'),'.pdf')
 pdf(pdf_file)
lapply(metamatrix_data, function(x){
  title=paste0(x$sample,'   ',x$name,'\n',x$group,'; ',down_thresh/1000,' kb to ',up_thresh/1000,' kb\n',x$nr_of_genes,' genes')
     f_draw_matrix(raw_array=x$mat_raw, norm_array=x$mat_norm,title=title,threshold=threshold,max_scale=c('norm'=max_norm))
 
})
dev.off()
