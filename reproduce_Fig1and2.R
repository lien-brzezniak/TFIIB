#This code reproduces the analysis of 3'RNAseq data
library(DESeq2)
library(tidyverse)
library(ggpubr)
library(ggrepel)
setwd('/home/lien/data/skrypty/TFIIB/')
source('functions.R')

#readCounts_table.tsv can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199056
nor <- read_tsv("readCounts_table.tsv")

sample_info_3RNAseq <- read_tsv('sample_info_3RNAseq.tsv', col_types ='cffdffff') %>%
  mutate(comparison=paste0(stress,'_',timepoint),
         genotype=factor(genotype, levels=c('Wt','R58C')))

#We used deseq2 package to make PCA analysis and identify DEGs. We took only nuclear protein coding genes:
nor_nuclear <- nor[-grep('AT[CM]',nor$GeneID),]
nor_nuclearPCG <- nor_nuclear[which(nor_nuclear$GeneID %in% all_genes$GeneID),] 
cts <- as.matrix(nor_nuclearPCG[, sample_info_3RNAseq$sample])
rownames(cts) <- nor_nuclearPCG$GeneID

## multifactorial analysis - PCA plot (Fig. 1C)
ddsMF <- DESeqDataSetFromMatrix(countData = cts,
                                colData = sample_info_3RNAseq ,
                                design = ~ genotype + comparison)

ddsMF <- estimateSizeFactors(ddsMF)
vsd <- vst(ddsMF, blind=FALSE)
# change scales to percentage and add custom colors
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
comparisons <- sample_info_3RNAseq$comparison %>% unique()
comparisons
cts_list <- as.vector(comparisons, mode='list')
names(cts_list)=comparisons

cts_list <- lapply(cts_list, function(x){
  # x=comparisons[2]
  samples=sample_info_3RNAseq$sample[which(sample_info_3RNAseq$comparison==x)]
  coldata <- sample_info_3RNAseq[which(sample_info_3RNAseq$sample %in% samples),]
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
  right_join(distinct(sample_info_3RNAseq[,c('comparison','condition')])) %>%
  select(GeneID,condition,count) %>%
  distinct() %>%
  group_by(condition) %>%
  summarize(n_of_DEGSs=sum(count, na.rm=T))
DEGs_No

All_DEGs_list <-reduce(DEGs_tables,bind_rows) %>%
  select(GeneID) %>%
  distinct() 
#write_tsv(All_DEGs_list,'DEGs.txt')
#this list was used as a query for GO term enrichment analysis
All_DEGs_list


# calculate trainability Index (Fig. 2)
# calculate mean rpm values for all conditions and genotypes (average from 3 biological repeats);
# filter out lowly expressed genes - reject genes that had mean rpm less than 5 in any experimental condition 
mean_rpm <- calculate_mean_rpm(nor, filter_nuclear=T, filter_protCoding=T) %>%
  pivot_longer(2:13, values_to='mean_rpm') %>%
  subset(mean_rpm>5) %>%
  pivot_wider(values_from='mean_rpm') %>%
  subset(complete.cases(.)) 

#calculate average foldchanges during primary and subsequent stress (for filtering drought-induced genes), and the trainability index
calculated <- mean_rpm %>% 
  mutate(fc_Wt_1=(Wt_S1_90/Wt_S1_0+Wt_S1_180/Wt_S1_0)/2) %>%
  mutate(fc_Wt_4=(Wt_S4_90/Wt_S4_0+Wt_S4_180/Wt_S4_0)/2) %>%
  mutate(fc_R58C_1=(R58C_S1_90/R58C_S1_0+R58C_S1_180/R58C_S1_0)/2) %>%
  mutate(fc_R58C_4=(R58C_S4_90/R58C_S4_0+R58C_S4_180/R58C_S4_0)/2) %>%
  mutate(Trainability_Wt = ((Wt_S4_180-Wt_S4_0+2*Wt_S4_90-2*Wt_S4_0)/2-(Wt_S1_180-Wt_S1_0+2*Wt_S1_90-2*Wt_S1_0)/2) )  %>%
  mutate(Trainability_R58C =((R58C_S4_180-R58C_S4_0+2*R58C_S4_90-2*R58C_S4_0)/2-(R58C_S1_180-R58C_S1_0+2*R58C_S1_90-2*R58C_S1_0)/2) )

calculated

#filter only genes that are induced in either primary or subsequent stress
filtered_up <- calculated %>% 
  ungroup() %>%
  filter(rowAny(across(contains('fc_'),~.x>3)))

filtered_up

#order by trainability in TFIIB1-WT
arranged <- filtered_up %>%
  arrange(desc(Trainability_Wt)) %>%
  transmute(GeneID, Wt=Trainability_Wt,R58C=Trainability_R58C) %>%
  pivot_longer(2:3,names_to="Genotype") %>%
  mutate(value=scale(value)[,1]) %>%
  pivot_wider(names_from="Genotype",values_from="value")

arranged$Index <- 1:nrow(arranged)
arranged

# prepare labels and color by trainability
arranged <- arranged %>%
  select(GeneID,Index,Wt,R58C) %>%
  mutate(delta=R58C-Wt,
         Type=ntile(Wt,n=4),
         Type_L=case_when(Type==1~'<25%',
                          (Type>1 & Type<4) ~ '25-75%',
                          Type==4 ~ '>75%'),
         delabel=ifelse(GeneID %in% c('AT5G66400','AT5G52300','AT5G52310','AT2G42540'),GeneID,NA))

write_tsv(arranged,'drought_genes.tsv')

#Fig. 2B upper panel
ggplot(arranged, aes(x=Index, y=Wt,color=Type_L)) + 
  ggtitle("Trainability. Icebox-2") + 
  # geom_vline(xintercept=c(-1,699/4,699/4*2,699/4*3,700),color='grey')+
  geom_hline(yintercept=0,size=1)+
  geom_point(size=2)+
  geom_point(data=arranged[which(!is.na(arranged$delabel)),],color='white',size=1)+
  scale_color_manual(values=kolory2)+
  scale_y_continuous(breaks=seq(-12.5,12.5,2.5))+
  xlab("Drought-induced genes ordered by Trainability Index in WT") + 
  geom_text_repel(aes(label=delabel))+
  ylab("Trainability Index") +
  newtheme

#prepare a Supplemental Table S8
drought_upregulated_table <- arranged %>%
  transmute(GeneID,
            Tr_order_Wt=Index,
            Tr_Index_Wt = Wt,
            Tr_Index_R58C = R58C,
            'Tr_delta(R58C-Wt)'=delta) %>%
  left_join(filtered_up[,c('GeneID',paste0(unique(sample_info_3RNAseq$red_sample)))])
drought_upregulated_table
# write_tsv(drought_upregulated_table,'/home/lien/data/TFIIB/figury/Trainability_table.tsv')

# prepare Fig 2B, lower panel
stacked_effect <- arranged %>%
  mutate(effect=ifelse(delta<0,'loss','gain')) %>%
  arrange(abs(delta))


stacked_effect$Type <- factor(stacked_effect$Type, levels=c('1','2','3','4'))

stacked_effect %>% ggplot(aes(x=Type, y=delta,group=effect,colour=Type_L))+
  geom_col(width=0.8,size=1 , fill=NA)+
  scale_x_discrete( limits=as.character(c(4:1)), labels=c(4:1))+
  scale_color_manual(aesthetics=c('color'),values=kolory2)+
  geom_hline(yintercept = 0, colour = "black",linetype=2) + 
  newtheme
# count genes changed in the same direction
numbers<- stacked_effect %>%
  group_by(Type, effect) %>%
  summarize(n=n())
numbers

# Supplemental Fig. S1
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

# Supplemental Fig. S2 - Single gene rpm graphs  
# transfrom counts table to rpm table  
rpm <- nor %>%
  mutate_if(is.numeric, rpm_fun)
# Top 30 most trainable genes:
single_plots <- lapply(c(1:30),function(x){
  gene=arranged$GeneID[x]
  TrI=arranged$Index[x]
  plot_gene <- rpm %>%
    filter(GeneID %in% gene) %>%
    pivot_longer(2:37, names_to = "sample")%>%
    left_join(sample_info_3RNAseq,by = "sample")%>%
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

