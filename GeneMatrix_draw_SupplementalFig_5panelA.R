setwd('/home/lien/data/skrypty/TFIIB')
library(gplots)
library(tidyverse)
source('functions.R')



data_dir='/home/lien/data/skrypty/TFIIB/GhvTX1648049065_Fig5/'
out_dir=paste0(data_dir,'draw/')
system(paste0('mkdir ',out_dir),intern=F)

groups=c('Trainability_Loss','Trainability_Gain')

Samples=list(Wt=list(Water=c(1,9),S1=c(2,10),Recovery=c(3,11),S2=4),R58C=list(Water=5,S1=c(6,12),Recovery=7, S2=8))

matrices <-  lapply(groups,function(group){
  lapply(Samples,function(genotype){
  lapply(genotype,function(condition){
    lapply(condition,function(sample){
      # browser()
      norm_matrix <- read_tsv(paste0(data_dir,'Meta_norm_Sample',sample,'_',group,'_3 kb to 50 kb.matrix'), col_names=F) %>% 
        as.matrix()
      scaled_norm_matrix <- (norm_matrix-mean(norm_matrix,na.rm=T))/sd(norm_matrix,na.rm=T)
    }) %>% Reduce('+',.)/length(condition)
  })
})
})
names(matrices)=groups


conditions=c('Water','S1','Recovery','S2')

deltas <- lapply(groups, function(group){
  deltas <-lapply(conditions,function(condition){
 matrices[[group]][['R58C']][[condition]]-matrices[[group]][['Wt']][[condition]]
 })
names(deltas)<-conditions
return(deltas)
})
names(deltas)=groups

threshold=0.99
max_norm = quantile(abs((unlist(matrices))),probs=threshold,na.rm=T)
max_delta = quantile(abs((unlist(deltas))),na.rm=T,probs=threshold)


pdf_title=paste0(out_dir,groups, '_maxthresh_',threshold,'Supp_Fig5A.pdf')
pdf(pdf_title)

lapply(groups, function(group){
lapply(c('Wt','R58C'),function(genotype){
  lapply(conditions,function(condition){
    
    norm_map <- matrices[[group]][[genotype]][[condition]]
    title=paste0(group,' \nThresh: ',threshold,'\n',genotype,' ',condition)
    
    f_draw_matrix(norm_array=norm_map, title=title, do_scaling=F,max_scale_norm=max_norm)
  })
})
})
lapply(groups, function(group){
lapply(conditions, function(condition){
  title=paste0(group,'\nR58C-Wt\n',condition)
  
  f_draw_matrix(delta_array=deltas[[group]][[condition]], title=title, do_scaling=F, max_scale_delta=max_delta)
})
})
dev.off()



