#*This script analyses sector calculations created by GeneMatrix_calculateSectors.R
library(tidyverse)
# library(stats)
 library(ggpubr)
library(ggplot2)
# library("boot")
# library(rstatix)
library(gridExtra)
library(grid)
# library(scico)
source('/home/lien/data/skrypty/TFIIB/functions.R')
setwd('/home/lien/data/skrypty/TFIIB/')
#*******************GeneGroups
#*
all_genes <- read_tsv('/home/lien/data/data_mining/TAIR10/proten_coding_genes.bed') %>%
  transmute(GeneID, gLength=floor((end-start))) 
group_sets=list(
  'other'=c('All'),
  'geneBody'=c(3,12,25,26))
#'other'=c('R58C_TrainabilityLoss','R58C_TrainabilityGain'))
down_thresh=3000
up_thresh=50000
exclusive=T
complement_rest=T
Groups <- create_Groups(group_sets,down_thresh,up_thresh,complement_rest,exclusive)

#out_dir='/home/lien/data/TFIIB/figury/Fig4/more_tests/'
out_dir=''
sample_info <- read_tsv('/home/lien/data/skrypty/TFIIB/samples_GEO.tsv',col_types='ffffff')

data <- lapply(1:12,function(i){
  read_tsv(paste0('/home/lien/data/TFIIB/Hi_C/TagDirs/matrices/raw_relative/GeneMatrix_allSizes/GEO_submission/Sample',i,'_sectors.tsv'))
}) %>% reduce(bind_rows) 

summary <- make_contact_summary(data,Groups)
plots <- report_contact_summary(summary,Groups)
boot_plot <- bootANDplot(data,Groups,merge_samples=c(1:4,9:11))

pdf_name=paste0(out_dir, 'AsRatio_Fig4B.pdf')
pdf(pdf_name,width = 4, height = 4)
print(plots)
print(boot_plot)
lapply(names(Groups),function(x){
  table <- summary %>% select(group, sample,mean_AsRatio) %>%
    subset(group==x)%>%
    tableGrob(theme=ttheme_default(base_size=6))
   grid.newpage()
  grid.draw(table)
})

dev.off()
