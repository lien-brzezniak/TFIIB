library(DESeq2)
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(scales)
setwd('/home/lien/data/skrypty/TFIIB/')
source('functions.R')


nor <- read_tsv("readCounts_table.tsv")

sample_info_3RNAseq <- read_tsv('sample_info_3RNAseq.tsv', col_types ='cffdffff') %>%
  mutate(comparison=paste0(stress,'_',timepoint),
         genotype=factor(genotype, levels=c('Wt','R58C')))

nor_nuclear <- nor[-grep('AT[CM]',nor$GeneID),]
nor_nuclearPCG <- nor_nuclear[which(nor_nuclear$GeneID %in% all_genes$GeneID),] 
cts <- as.matrix(nor_nuclearPCG[, sample_info_3RNAseq$sample])
rownames(cts) <- nor_nuclearPCG$GeneID





## multifactorial analysis - PCA plot
 ddsMF <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sample_info_3RNAseq ,
                              design = ~ genotype + comparison)

ddsMF <- estimateSizeFactors(ddsMF)
vsd <- vst(ddsMF, blind=FALSE)
plotPCA(vsd, intgroup=c("comparison", "genotype"))
pcaData <- plotPCA(vsd, intgroup=c("comparison", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) +
  scale_color_manual(values=kolory)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  newtheme

# identify genes differentially expressed in TFIIB1-R58C compared to TFIIB1-WT
comparisons <- sample_info$comparison %>% unique()
cts_list <- as.vector(comparisons, mode='list')
names(cts_list)=comparisons
 
cts_list <- lapply(cts_list, function(x){
 # x=comparisons[2]
  samples=sample_info$sample[which(sample_info$comparison==x)]
  coldata <- sample_info[which(sample_info$sample %in% samples),]
  cts <- cts[,coldata$sample]
  return(list('coldata'=coldata,'cts'=cts))
})

dds_list <- lapply(cts_list, function(x) {
  dds <- DESeqDataSetFromMatrix(countData = x$cts,
                              colData = x$coldata,
                              design = ~ genotype)
  keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$comparison <- factor(dds$comparison, levels=c('Wt','R58C'))
dds$comparison <- relevel(dds$comparison, ref = "Wt")
dds <- DESeq(dds)
})

res_list <- lapply(dds_list,function(x){
res <- results(x)
resOrdered <- res[order(res$pvalue),]
})
deseq_tables <- lapply(res_list,function(x){
    frame <- as.data.frame(x)
    #%>%
   # mutate(comparison=x)
    frame$GeneID <- rownames(frame)
    rownames(frame)=NULL
    colnames(frame) = c('baseMean','log2fc','lfcSE','stat','pvalue','padj','GeneID')
  return(as_tibble(frame))
})

DEGs_tables <- lapply(deseq_tables, function(x){
  x %>% subset(padj<0.05 & abs(log2fc)>=0)
})
DEGs_tables
DEGs_No <- lapply(names(DEGs_tables), function(c){
  DEGs_tables[[c]] %>% mutate(comparison=c,count=1)
}) %>% reduce(bind_rows) %>%
  right_join(distinct(sample_info[,c('comparison','condition')])) %>%
  select(GeneID,condition,count) %>%
  distinct() %>%
  group_by(condition) %>%
  summarize(n_of_DEGSs=sum(count, na.rm=T))
DEGs_No

All_DEGs_list <-reduce(DEGs_tables,bind_rows) %>%
  select(GeneID) %>%
  distinct() 
#write_tsv(All_DEGs_list,'DEGs.txt')

All_DEGs_list


# calculate trainability Index

mean_rpm <- calculate_mean_rpm(nor, filter_nuclear=T, filter_protCoding=T) %>%
  pivot_longer(2:13, values_to='mean_rpm') %>%
  subset(mean_rpm>5) %>%
  pivot_wider(values_from='mean_rpm') %>%
  subset(complete.cases(.)) 
# %>%
#   subset(rowSums(.[2:13])>250)

normalized_counts <-calculate_mean_normCTS(nor, filter_nuclear=T, filter_protCoding=T)%>%
   pivot_longer(2:13, values_to='mean_rpm') %>%
  subset(mean_rpm>5) %>%
  pivot_wider(values_from='mean_rpm') %>%
  subset(complete.cases(.))%>%
  subset(rowSums(.[2:13])>500)

# calculated <- mean_rpm %>% 
#   mutate(fc_Wt_S1=(Wt_S1_90/Wt_S1_0+Wt_S1_180/Wt_S1_0)/2) %>%
#   mutate(fc_Wt_S4=(Wt_S4_90/Wt_S4_0+Wt_S4_180/Wt_S4_0)/2) %>%
#   mutate(fc_R58C_S1=(R58C_S1_90/R58C_S1_0+R58C_S1_180/R58C_S1_0)/2) %>%
#   mutate(fc_R58C_S4=(R58C_S4_90/R58C_S4_0+R58C_S4_180/R58C_S4_0)/2) %>%
#   mutate(Trainability_Wt = ((Wt_S4_180-Wt_S4_0+2*Wt_S4_90-2*Wt_S4_0)/2-(Wt_S1_180-Wt_S1_0+2*Wt_S1_90-2*Wt_S1_0)/2) )  %>%
#   mutate(Trainability_R58C =((R58C_S4_180-R58C_S4_0+2*R58C_S4_90-2*R58C_S4_0)/2-(R58C_S1_180-R58C_S1_0+2*R58C_S1_90-2*R58C_S1_0)/2) )

calculated <- mean_rpm %>% 
   mutate(fc_Wt_1=(Wt_S1_90/Wt_S1_0+Wt_S1_180/Wt_S1_0)/2) %>%
  mutate(fc_Wt_4=(Wt_S4_90/Wt_S4_0+Wt_S4_180/Wt_S4_0)/2) %>%
  mutate(fc_R58C_1=(R58C_S1_90/R58C_S1_0+R58C_S1_180/R58C_S1_0)/2) %>%
  mutate(fc_R58C_4=(R58C_S4_90/R58C_S4_0+R58C_S4_180/R58C_S4_0)/2) %>%
    mutate(Trainability_Wt = ((Wt_S4_180-Wt_S4_0+2*Wt_S4_90-2*Wt_S4_0)/2-(Wt_S1_180-Wt_S1_0+2*Wt_S1_90-2*Wt_S1_0)/2) )  %>%
  mutate(Trainability_R58C =((R58C_S4_180-R58C_S4_0+2*R58C_S4_90-2*R58C_S4_0)/2-(R58C_S1_180-R58C_S1_0+2*R58C_S1_90-2*R58C_S1_0)/2) )

calculated

#filter only genes that are induced in any of measured treatments
filtered_up <- calculated %>% 
  ungroup() %>%
  filter(rowAny(across(contains('fc_'),~.x>3)))
  
filtered_up

arranged <- filtered_up %>%
  arrange(desc(Trainability_Wt)) %>%

  transmute(GeneID, Wt=Trainability_Wt,R58C=Trainability_R58C) %>%
  pivot_longer(2:3,names_to="Genotype") %>%
  mutate(value=scale(value)[,1]) %>%
  pivot_wider(names_from="Genotype",values_from="value")
  

arranged$Index <- 1:nrow(arranged)

arranged
x <- arranged %>%
  select(GeneID,Index,Wt,R58C)%>%
  mutate(delabel=ifelse(GeneID %in% c('AT5G66400','AT5G52300','AT5G52310','AT2G42540'),GeneID,NA))

arranged <- arranged %>%
  select(GeneID,Index,Wt,R58C) %>%
  mutate(delta=R58C-Wt,
         Type=ntile(Wt,n=4),
         Type_L=case_when(Type==1~'<25%',
                        (Type>1 & Type<4) ~ '25-75%',
                        Type==4 ~ '>75%'),
         delabel=ifelse(GeneID %in% c('AT5G66400','AT5G52300','AT5G52310','AT2G42540'),GeneID,NA))
drought_upregulated_table <- arranged %>%
  transmute(GeneID,
            Tr_order_Wt=Index,
            Tr_Index_Wt = Wt,
            Tr_Index_R58C = R58C,
            'Tr_delta(R58C-Wt)'=delta) %>%
  left_join(filtered_up[,c('GeneID',paste0(unique(sample_info_3RNAseq$red_sample)))])

# write_tsv(drought_upregulated_table,'/home/lien/data/TFIIB/figury/Trainability_table.tsv')
# write_tsv(arranged,'drought_genes.tsv')

ggplot(arranged, aes(x=Index, y=Wt,color=Type_L)) + 
  ggtitle("Trainability. Icebox-2") + 
 # geom_vline(xintercept=c(-1,699/4,699/4*2,699/4*3,700),color='grey')+
  geom_hline(yintercept=0,size=1)+
  geom_point(size=2)+
  geom_point(data=x[which(!is.na(x$delabel)),],color='white',size=1)+
  scale_color_manual(values=kolory2)+
  scale_y_continuous(breaks=seq(-12.5,12.5,2.5))+
  xlab("Drought-induced genes ordered by Trainability Index in WT") + 
  geom_text_repel(aes(label=delabel))+
  ylab("Trainability Index") +
  newtheme



stacked_effect <- arranged %>%
  mutate(effect=ifelse(delta<0,'loss','gain')) %>%
  arrange(abs(delta))

xxxc<- stacked_effect %>%
  group_by(Type, effect) %>%
  summarize(n=n())
xxxc
stacked_effect$Type <- factor(stacked_effect$Type, levels=c('1','2','3','4'))

stacked_effect %>% ggplot(aes(x=Type, y=delta,group=effect,colour=Type_L))+
   geom_col(width=0.8,size=1 , fill=NA)+
  scale_x_discrete( limits=as.character(c(4:1)), labels=c(1:6))+
  scale_color_manual(aesthetics=c('color'),values=kolory2)+
  geom_hline(yintercept = 0, colour = "black",linetype=2) + 
  newtheme

expression_data <-mean_rpm %>%
  subset(GeneID %in% filtered_up$GeneID)%>%
  pivot_longer(2:13, names_to='red_sample')%>%
  left_join(distinct(sample_info_3RNAseq[,c('red_sample','comparison','genotype')]))


expression_data$comparison <- factor(expression_data$comparison,levels=c('S1_0','S1_90','S1_180','S4_0','S4_90','S4_180'))

stats <- expression_data %>%
  ggviolin(x="genotype",y="value",fill="genotype",color="genotype",alpha=0.3) +
  scale_color_manual(values=kolory, aesthetics = c("color","fill"))+
    scale_y_log10()+
    stat_compare_means(label = "p.format",size=3, method='wilcox.test', paired=F)+
  geom_boxplot(width=0.5,size=0.5,aes(color=genotype))+ 
 facet_wrap(~comparison,nrow=1, scales="fixed")+
 ggtitle('Expression of drought-induced genes')+
    newtheme
stats
 # svg('/home/lien/data/TFIIB/figury/SuppFig1/expression_droughtGenes.svg', width=6,height=3)
 # print(stats)
 # dev.off()
  
  
  
rpm <- nor %>%
  mutate_if(is.numeric, rpm_fun)






single_plots <- lapply(c(1:30),function(x){
  
gene=arranged$GeneID[x]
TrI=arranged$Index[x]
plot_gene <- rpm %>%
  filter(GeneID %in% gene) %>%
    pivot_longer(2:37, names_to = "sample")%>%
  left_join(sample_info,by = "sample")%>%
  group_by(genotype,stress,timepoint) %>%
  summarize(GeneID,mean_rpm=mean(value),stdev=sd(value))%>%
  distinct() %>%
  mutate(stress=if_else(stress=='S1','Primary Stress','Subsequent Stress'))
plot_gene$timepoint <- factor(plot_gene$timepoint,levels=c('0','90','180'))
plot <-ggplot(plot_gene, aes(timepoint,mean_rpm,colour=genotype,group=genotype))+
  geom_point() +
  geom_line() +
  ggtitle(paste0(gene,", ",TrI)) +
  scale_color_manual(values=kolory)+
  scale_y_continuous(name=NULL)+
  scale_x_discrete(name=NULL)+
  geom_errorbar(aes(ymin=mean_rpm-stdev, ymax=mean_rpm+stdev), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~stress)+
  newtheme+
  theme(axis.line.x=element_line(size=0.7),
        axis.line.y=element_blank(),
        text = element_text(size = 8))

})

single_plots[1]
# pdf('/home/lien/data/TFIIB/figury/icebox/single_gene_rpm.pdf',width=3.5,height=3)
# print(single_plots)
# dev.off()
# 
# svg_dir = '/home/lien/data/TFIIB/figury/icebox/single_gene_rpm_svgs/'
# system(paste0('mkdir ',svg_dir))
# lapply(1:length(single_plots), function(x){
# svg(paste0(svg_dir,x,'.svg'),width=2.2,height=2)
#  print(plot)
# dev.off()
# })
