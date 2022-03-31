
setlibrary("tidyverse")
library(ggpubr)
library("scales")
setwd('/home/lien/data/skrypty/TFIIB/')
source('functions.R')

#*******************Groups
#*

group_sets=list(
  'geneBody'=c(3,12,25,26),
 'other'=c('All'))
#'other'=c('R58C_TrainabilityLoss','R58C_TrainabilityGain'))

Groups <- create_Groups(group_sets,down_thresh=3000,up_thresh=50000,complement_rest=F,exclusive=T)
#****************************************************
kolory=c("TFIIB_WT"="#4261bb","TFIIB_R58C"="#d15992","TFIIB_mut"="#d15992","WT"="#4261bb","R58C"="#d15992")

nor <- read_tsv('readCounts_table.tsv')


sample_info <- read_tsv('sample_info_3RNAseq.tsv', col_types ='cffdffff') %>%
  mutate(comparison=paste0(stress,'_',timepoint),
         genotype=factor(genotype, levels=c('Wt','R58C')))

calculate_mean_rpm <- function(nor, filter_nuclear=T, filter_protCoding=T
                               
WT_NT_mean_rpm <- mean_rpm[c('GeneID','Wt_S1_0')]

expression_data <- lapply(names(Groups),function(g){
WT_NT_mean_rpm %>%
  subset(GeneID %in% Groups[[g]])%>%
  mutate(group=g) 
}) %>% purrr::reduce(bind_rows)
# kolors <- c('#F8766D','#7CAE00','#454545','#00BFC4','#C77CFF')
expression_data$group <- factor(expression_data$group, levels=c("State 26","State 12", "All","State 25", "State 3" ))

expressionPlot <- expression_data %>%
  ggviolin(x="group",y="Wt_S1_0",alpha=0.3,color='group',fill='group') +
  scale_color_manual(values=kolory, aesthetics = c("color","fill"))+
   scale_y_log10(name='Expression level [rpm]',breaks=c(1,10,100,1000), labels=c('1','10','100','1000'))+
  scale_x_discrete(limits=c("State 26","State 12", "All","State 25", "State 3" ))+
   geom_boxplot(aes(x=group,y=Wt_S1_0,color=group),width=0.1)+ 
  newtheme+
  theme( aspect.ratio=1,
         #strip.text=element_blank(),
         axis.title=element_text(size = 10),
         axis.text=element_text(size = 10),
         axis.text.x=element_text(angle = 75,hjust=1))

svg('expression_singlePlot.svg', width=3)
print(expressionPlot)
dev.off()


