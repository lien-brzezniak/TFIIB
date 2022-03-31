library(HiTC)
library(tidyverse)
library(scico)
source('/home/lien/data/skrypty/TFIIB/functions.R')

lists_samples=list(Wt_NT=as.character(c(1,9)),
                   Wt_S1=as.character(c(2,10)),
                   Wt_Rec=as.character(c(3,11)),
                   Wt_S2=as.character(c(4)),
                   R58C_NT=as.character(c(5)),
                   R58C_S1=as.character(c(6,14)),
                   R58C_Rec=as.character(c(7)),
                   R58C_S2=as.character(c(8)))
lapply(names(lists_samples),f_convertAndSaveCdata_combined)


ICEnormed <-  lapply(names(lists_samples), function(sample){
  data=readRDS(paste0('/home/lien/data/TFIIB/Hi_C/HiTC/',sample,'_100kbWholeGenomeMatrix'))
  normICE(data, max_iter=200,sparse.filter=0.02)
})
names(ICEnormed)<- names(lists_samples)
saveRDS(ICEnormed,paste0('/home/lien/data/TFIIB/Hi_C/HiTC/ICEnormed_comb_replicates.RDS'))

all_data <-   readRDS(paste0('/home/lien/data/TFIIB/Hi_C/HiTC/ICEnormed_allData.RDS'))
ICEnormed <-  readRDS(paste0('/home/lien/data/TFIIB/Hi_C/HiTC/ICEnormed_WTandR58C_Data.RDS'))
ICE_comb_repl <-  readRDS(paste0('/home/lien/data/TFIIB/Hi_C/HiTC/ICEnormed_comb_replicates.RDS'))

kolory = c('white',scico(palette='lajolla',n=31))
# mapC(ICEnormed[['WT']],minrange=0,maxrange=12,log.data=T,show.zero=F,show.na=F,col.pos=kolory)
# mapC(ICEnormed[['R58C']],minrange=0,maxrange=10,log.data=T,show.zero=F,show.na=F,col.pos=c("white", "yellow","orange", "red", "black"))
#mapC(ICEnormed[['WT']],ICEnormed[['R58C']],minrange=0,maxrange=10,log.data=T,show.zero=F,show.na=F,col.pos=kolory)
mapC(ICE_comb_repl[['Wt_Rec']],ICEnormed[['Wt_S2']],minrange=0,maxrange=10,log.data=T,show.zero=F,show.na=F,col.pos=kolory)


