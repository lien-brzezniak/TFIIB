setwd('/home/lien/data/skrypty/TFIIB')
library(gplots)
library(tidyverse)
source('functions.R')



data_dir='/home/lien/data/TFIIB/figury/Fig3/qfuKz1645798144/'
out_dir=paste0(data_dir,'draw/')
system(paste0('mkdir ',out_dir),intern=F)

groups=c('All')

Samples=list(WT=list(Water=c(1,9),S1=c(2,10),Recovery=c(3,11),S2=4),R58C=list(Water=5,S1=c(6,14),Recovery=7, S2=8))

matrices <-  lapply(Samples,function(genotype){
  lapply(genotype,function(condition){
    lapply(condition,function(sample){
      norm_matrix <- read_tsv(paste0(data_dir,'Meta_norm_',sample,'_',groups,'_3to50kb_nonscaled.txt'), col_names=F) %>% 
        as.matrix()
      scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
    }) %>% Reduce('+',.)/length(condition)
  })
})

threshold=0.998

conditions=c('Water','S1','Recovery','S2')
deltas <- lapply(conditions,function(condition){
 matrices[['R58C']][[condition]]-matrices[['WT']][[condition]]
 })
names(deltas)<-conditions

max_norm = quantile(abs((unlist(matrices))),probs=threshold,na.rm=T)
max_delta = quantile(abs((unlist(deltas))),na.rm=T,probs=threshold)


pdf_title=paste0(out_dir,groups, '_maxthresh_',threshold,'.pdf')
pdf(pdf_title)

lapply(names(matrices),function(genotype){
  lapply(names(matrices[[genotype]]),function(condition){
    norm_map <- matrices[[genotype]][[condition]]
    title=paste0(groups,' \nThresh: ',threshold,'\n',genotype,' ',condition)
    
    f_draw_matrix(norm_array=norm_map, title=title, do_scaling=F,max_scale_norm=max_norm)
  })
})
lapply(names(deltas), function(d){
  title=paste0(groups,'\nR58C-WT\n',d)
  
  f_draw_matrix(delta_array=deltas[[d]], title=title, do_scaling=F, max_scale_delta=max_delta)
})
dev.off()



