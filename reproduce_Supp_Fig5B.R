# The script takes list of genes and draws normalized or raw (norm=c('norm','raw')) heatmaps from one sample, local scaled.
setwd('/home/lien/data/skrypty/TFIIB/')
library("tidyverse")
library("gplots")
library(gridExtra)
library(grid)
source('functions.R')

output_path = ''

group_sets=list('custom_gene_list'=c('AT1G03080','AT3G12580','AT5G20900','AT4G39980','AT2G39800',
                                     'AT1G09970','AT1G62660','AT4G16760',
                                     'AT4G18950','AT4G23050','AT5G10650','AT5G51070'))

Groups <- create_Groups(group_sets)

dir='/home/lien/data/TFIIB/GEO_upload/HiC/'


sample_sets=list('Wt'=as.character(c(2,4,10)),'R58C'=as.character(c(6,8,12)))

# extract a contact matrix for each gene, each sample
data <- lapply(names(Groups),function(group){
  genes = Groups[[group]]
  d <- lapply(names(sample_sets),function(genotype){
    d <-lapply(sample_sets[[genotype]],function(sample){
      # browser()
      gene_matrix <- singleGene_map(sample,genes,format='norm')
      
    })
    names(d) <- sample_sets[[genotype]]
    return(d)
    
  })
  names(d) <- names(sample_sets)
  return(d)
})
names(data) <- names(Groups)

# set the threshold for heatmap scale
threshold=0.995

# average tha matrixes for each gene across the sample sets 
combined_repl_data <- lapply(names(Groups),function(group){
  d <-lapply(Groups[[group]],function(gene){
    d <-lapply(names(sample_sets), function(set){
      lapply(sample_sets[[set]], function(sample){
        norm_matrix <- data[[group]][[set]][[sample]][[gene]]
        scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
      })%>% Reduce('+',.)/length(sample_sets[[set]])
    })
    names(d) <- names(sample_sets)
    d$max_norm = quantile(abs((unlist(d))),probs=threshold,na.rm=T)
    return(d)
  })
  names(d) <- Groups[[group]]
  return(d)
})
names(combined_repl_data) <- names(Groups)
# plot the first gene:
group=names(Groups)[1]
gene=Groups[[1]][1]
lapply(names(sample_sets), function(set){
  matrix <- combined_repl_data[[group]][[gene]][[set]]
  thresh=combined_repl_data[[group]][[gene]][['max_norm']]
  title=paste0(gene,' from group: ',group,'\n',set,' (',paste0(sample_sets[[set]],collapse=', '),')')
  f_draw_matrix(norm_array=matrix,title=title, max_scale_norm=thresh)
})

# plot all heatmaps in a pdf report, 1 file for 1 group of genes defined in group_sets
barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'_SupplementalFig_5B/')
system(paste0('mkdir ',output_path),intern=F)
lapply(1:length(Groups), function(g){
  group=names(Groups)[g]
  genes=Groups[[g]]
  pdf_file=paste0(output_path,group,'_G',g,'_single_metagenMaps.pdf')
  pdf(pdf_file)
  
  lapply(genes,function(gene){
    lapply(names(sample_sets), function(set){
      matrix <- combined_repl_data[[group]][[gene]][[set]]
      thresh=combined_repl_data[[group]][[gene]][['max_norm']]
      title=paste0(gene,' from group: ',group,'\n',set,' (',paste0(sample_sets[[set]],collapse=', '),')')
      f_draw_matrix(norm_array=matrix,title=title, max_scale_norm=thresh)
    })
  })
  dev.off() 
})