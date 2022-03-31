setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
source('functions.R')




sample_info <- read_tsv('samples_GEO_HiC.tsv',col_types='ffffff')

# read functions.R::create_Groups() to  see other defined gene-groups. A custom list of genes may also be provided.
group_sets=list('other'=c('All_PCG'))

#specify the directory, where GMatrix.RDS files are located:
dir='/home/lien/data/TFIIB/GEO_upload/HiC/'
down_thresh=5000
up_thresh=50000
Groups <- create_Groups(group_sets,down_thresh,up_thresh)
s='1'
metamatrix_data <- lapply(1:length(Groups),function(G){
  genes=Groups[[G]]
  result <- calculate_metamatrix(genes,s,reach=2500)
  result$sample=paste0('Sample',s)
  result$name <- unite(sample_info,c(2:4),col='name',sep='_' )[which(sample_info$sample==result$sample),] %>% pull(name)
  result$group=names(Groups)[G]
  return(result)
})

# the results will be saved as pdf report and/or matrix .txt files
barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'/')
system(paste0('mkdir ',output_path),intern=F)
pdf_file=paste0(output_path,paste0(unlist(group_sets),collapse='_'),'.pdf')
pdf(pdf_file)
lapply(metamatrix_data, function(x){
  title=paste0(x$sample,'   ',x$name,'\n',x$group,'; ',down_thresh/1000,' kb to ',up_thresh/1000,' kb\n',x$nr_of_genes,' genes')
  f_draw_matrix(x$mat_raw, x$mat_norm,do_scaling=T,title,threshold=0.995)
  write.table(mat_raw, file=paste0(output_path,'Meta_raw_Sample',s,'_',names(Groups)[G],'_',descr1,'.txt'), sep='\t', row.names = F, col.names =  F)
  write.table(mat_norm, file=paste0(output_path,'Meta_norm_Sample',s,'_',names(Groups)[G],'_',descr1,'.txt'), sep='\t', row.names = F, col.names =  F)

})
dev.off()
