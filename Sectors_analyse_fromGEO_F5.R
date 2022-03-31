#*This script analyses sector calculations created by GeneMatrix_calculateSectors.R

setwd('/home/lien/data/skrypty/TFIIB/')
library(tidyverse)
# library(stats)
 library(ggpubr)
library(ggplot2)
# library("boot")
# library(rstatix)
library(gridExtra)
library(grid)
# library(scico)
source('functions.R')

#*******************GeneGroups
#*
# all_genes <- read_tsv('/home/lien/data/data_mining/TAIR10/proten_coding_genes.bed') %>%
#   transmute(GeneID, gLength=floor((end-start))) 
group_sets=list(
 'other'=c('Trainability_Loss','Trainability_Gain'))

Groups <- create_Groups(group_sets,down_thresh=3000,up_thresh=50000,complement_rest=F,exclusive=F)



sample_info <- read_tsv('/home/lien/data/skrypty/TFIIB/samples_GEO_HiC.tsv',col_types='ffffff')

data <- lapply(1:12,function(i){
  read_tsv(paste0('/home/lien/data/TFIIB/GEO_upload/HiC/Sample',i,'_sectors.tsv')) %>%
   filter(complete.cases(.))
}) %>% reduce(bind_rows) 

summary <- make_contact_summary(data,Groups)
plots <- report_contact_summary(summary,Groups)
plots <- reduce(plots,c)
plots$AsymmetryRatio_Trainability_Loss <- plots$AsymmetryRatio_Trainability_Loss+
 scale_y_continuous(limits=c(0,1.25),breaks=seq(-1,4,0.25))
  plots$AsymmetryRatio_Trainability_Gain <- plots$AsymmetryRatio_Trainability_Gain+
 scale_y_continuous(limits=c(0,1.25),breaks=seq(-1,4,0.25))
#plots$boot_plot <- bootANDplot(data,Groups,merge_samples=c(6,8,14))
#print(plots[5])

barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
output_path = paste0(barcode,'_plotsFig5/')
system(paste0('mkdir ',output_path),intern=F)
pdf_name=paste0(output_path, 'AsRatio_Fig5.pdf')
pdf(pdf_name,width = 4, height = 4)
print(plots)
lapply(names(Groups),function(x){
  table <- summary %>% select(group, sample,mean_AsRatio) %>%
    subset(group==x)%>%
    tableGrob(theme=ttheme_default(base_size=6))
   grid.newpage()
  grid.draw(table)
})

dev.off()

lapply(1:length(plots),function(x){
  svg(paste0(output_path,x,'.svg'),width = 4, height = 4)
  print(plots[[x]])
  dev.off()
})

