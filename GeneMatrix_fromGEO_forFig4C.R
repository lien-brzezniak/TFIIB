setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
library("gplots")
library(gridExtra)
library(grid)
source('functions.R')



#output_path = '/home/lien/data/TFIIB/Hi_C/TagDirs/matrices/raw_relative/GeneMatrix_allSizes/calc_matrix/'
sample_info <- read_tsv('/home/lien/data/skrypty/TFIIB/samples_GEO_HiC.tsv',col_types='ffffff')
group_sets=list(
  'other'=c('All'),
  'geneBody'='acetylated')
  #'other'='R58C_TrainabilityLoss')
 
 
#grupy=c('CDS_states','R58C_TrainabilityLoss')
#grupy=c('top_genes','deltaMem','non_trainable')
#grupy=c('single')

samples=c(1:12)



####********MAIN
####*
####*

barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'Fig4C/')
system(paste0('mkdir ',output_path),intern=F)

res = 100
reach = 1500
down_thresh=3000
up_thresh=50000
write_matrix_file=F
Groups <- create_Groups(group_sets,down_thresh,up_thresh,complement_rest=T,exclusive=F)

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

threshold=0.995
max_norm = quantile(abs(scale(unlist(lapply(metamatrix_data,function(x){x$'mat_norm'})))),probs=threshold,na.rm=T)


 pdf_file=paste0(output_path,paste0(names(Groups),collapse='_'),'.pdf')
 pdf(pdf_file)
lapply(metamatrix_data, function(x){

  title=paste0(x$sample,'   ',x$name,'\n',x$group,'; ',down_thresh/1000,' kb to ',up_thresh/1000,' kb\n',x$nr_of_genes,' genes')
     f_draw_matrix(x$mat_raw, x$mat_norm,do_scaling=T,title=title,max_scale_norm=max_norm,threshold=threshold)
 return(title)
})
dev.off()

##combine replicates

output_path='/home/lien/data/skrypty/TFIIB/aaAib1647967542Fig4C/'
Samples=list(WT=list(water=c(1,3,9,11),drought=c(2,4,10)),R58C=list(water=c(5,7),drought=c(6,8,12)))
matrices <- lapply(names(Groups),function(group){
  # group=names(Groups)[3]
   lapply(Samples,function(genotype){
    lapply(genotype,function(condition){
      lapply(condition,function(sample){
        norm_matrix <- read_tsv(paste0(output_path,'Meta_norm_Sample',sample,'_',group,'_',down_thresh/1000,' kb to ',up_thresh/1000,' kb.matrix'), col_names=F) %>%
          as.matrix()
         scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
      }) %>% Reduce('+',.)/length(condition)
    })
  })
})
names(matrices)<-names(Groups)
  
  threshold=0.9999
  
conditions=c('drought','water')
deltas <- lapply(names(Groups),function(group){lapply(conditions,function(condition){
 d_array = matrices[[group]][['R58C']][[condition]]-matrices[[group]][['WT']][[condition]]
 title=paste0(group,'\n',condition,'\nR58C-WT\nDifferential matrix')
 return(list(delta_array=d_array, group=group, condition=condition))
 })})%>% reduce(c)


max_norm = quantile(abs((unlist(matrices))),probs=threshold,na.rm=T)
#max_delta = quantile(abs((unlist(deltas))),na.rm=T,probs=threshold)
max_delta = quantile(abs((unlist(lapply(deltas,function(x){x$delta_array})))),probs=threshold,na.rm=T)

  
  pdf_file=paste0(output_path,paste0(names(Groups),collapse='_'),'.pdf')
  pdf(pdf_file)
lapply(names(Groups),function(group){ 
lapply(names(Samples),function(genotype){
  lapply(conditions,function(condition){
    norm_map <- matrices[[group]][[genotype]][[condition]]
    title=paste0(group,' \nThresh: ',threshold,'\n',genotype,' ',condition)
    
    f_draw_matrix(norm_array=norm_map, title=title, do_scaling=F,max_scale_norm=max_norm)
  })
})
})
lapply(deltas, function(d){
  title=paste0(d$group,'\n',d$condition,'\nR58C-WT')
  
  f_draw_matrix(delta_array=d$delta_array, title=title, do_scaling=F, max_scale_delta=max_delta)
})
  dev.off()
  
