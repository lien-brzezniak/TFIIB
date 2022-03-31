library(tidyverse)
library(boot)
library(scico)
library(gplots)

library(grid)
library(gridExtra)
library(ggpubr)
library(ggrepel)

# load length data for all protein coding genes, according to TAIR10 annotation 
all_genes <- read_tsv('proten_coding_genes.bed') %>%
  mutate( gLength=floor((end-start))) 

#*create_Groups() prepares a list of gene-lists from pre-defined sets.
#* group_sets must be a list. Examples:
# group_sets=list(
#   'geneBody'=c(1:36),
#   'promoter'='acetylated',
#   'promoter'=c(1:8,15:21),
#   'custom_gene_list'=c('AT1G03080','AT3G12580','AT5G20900','AT4G39980','AT2G39800',
#                        'AT1G09970','AT1G62660','AT4G16760',
#                        'AT4G18950','AT4G23050','AT5G10650','AT5G51070'),
#   'other'=c('All_PCG','Trainability_Loss','Trainability_Gain','droughtgenes'))
# down_thresh=0     filter genes by min length
# up_thresh=50000   filter genes by max length
# complement_rest=F creates an additional list of genes, called 'rest', that contains other protein 
#   coding genes not listed in provided groups
# exclusive=F       leave only genes that are unique in each group. Affects all groups except for 'All_PCG'
# create_Groups() returns a list of character vectors containing AGI identifiers
create_Groups <- function(group_sets,down_thresh=0,up_thresh=50000,complement_rest=F,exclusive=F){
    Groups <- lapply(names(group_sets),function(set){
     grupy=group_sets[[set]]
      if (set %in% c('promoter','exon','intron','downstream','geneBody')) {
# gene lists from PCSD database (Liu Y, et al., 2018. Nucleic Acids Res 46: D1157–D1167.). State 1 to 36, enriched in specific gene region
      Groups <- readRDS('PCSD_allData.RDS')[[set]]
          if(grupy=='acetylated'){
        Groups=list('acetylated' = Groups[which(names(Groups) %in% paste0('At_genes_S',c(22:28)))] %>% unlist() %>% unique())
      } else{
        Groups <- Groups[which(names(Groups) %in% paste0('At_genes_S',grupy))]
      }
      
    }  else if (set == 'custom_gene_list') {
      Groups <- list('custom_gene_list'=grupy)
      
    }  else {
      drought_genes <- read_tsv('drought_genes.tsv')
      subSet = list()
      subSet$All_PCG <- pull(all_genes, GeneID)
      subSet$Trainability_Loss <-  drought_genes %>% subset(delta<0) %>%  pull(GeneID)
      subSet$Trainability_Gain <-  drought_genes %>% subset(delta>0) %>%  pull(GeneID)
       subSet$droughtGenes <-  drought_genes %>%  pull(GeneID)
     
      Groups = NULL
      for (G in 1:length(grupy)) {
        Groups[[G]] <-
          all_genes[which(all_genes$GeneID %in% subSet[[which(names(subSet) == grupy[G])]]), ] %>% pull(GeneID)
      }
      names(Groups) = grupy
     
    }
     return(Groups)
  }) %>% purrr::reduce(c)
  
  
  if(exclusive){
    
    Groups2 <- Groups[which(!names(Groups)=='All_PCG')]
    names<- paste0('Excl_',names(Groups2))
    Groups2 <- lapply(1:length(Groups2),function(i){
      
      x=Groups2[[i]]
      y= Groups2[which(!names(Groups2)==names(Groups2)[i])] %>% reduce(c)
      setdiff(x,y)
    })
    names(Groups2) <- names
    Groups <- c(Groups2, Groups[which(names(Groups)=='All')])
    
  }
  if (complement_rest){
    Groups$rest <- all_genes$GeneID[which(!all_genes$GeneID %in% unlist(Groups[which(!names(Groups)=='All_PCG')]))]
  }
  Groups <- lapply(Groups, function(n) {
    #browser()
    all_genes %>%  subset((gLength > down_thresh) &
                            (gLength < up_thresh) &
                            GeneID %in% n) %>% pull(GeneID)
  })
  
  new_names <- names(Groups)
  new_names <- sub('Excl_At_genes_S','Excl_State ',new_names)
  names(Groups) <- new_names
  return(Groups)
}

# create_meta_pos() builds a table necessary to assign positions on a meta-matrix
# reach=   set up the size of a meta-matrix. It should be not greater than half of min_thresh of analyzed gene length:
#   if genes longer than 3 kb are analysed, reach is set to 1500
create_meta_pos <- function(reach,res=100,gLength){
  tibble(label=c(as.character(seq(-4900,-100,100)),'TSS',as.character(seq(100,reach,100)),
                 rep('GeneBody',1),
                 as.character(seq(-1*(reach),-100,100)),'TTS',as.character(seq(100,4900,100))),
         pos=c(1:((reach/res*2)+101)),
         relTSS=c(seq(-4900,reach,100),NA,seq(gLength-reach,gLength+4900,100)))
  #,         relTTS=c(rep(NA,(51+(reach/res))), seq(-1*reach,4900,100)))
}

# calculate_metamatrix() calculate a meta-gene matrix, where each value is a Hi-C signal averaged across a group of genes.
# It takes Hi-C data for each gene, converts from upper-triangular sparse format to symmetric matrix,
# Next, genes are aligned at the TSS and TTS. Genes are not scaled. Instead, only fragments 
# from -5000 to <reach> with respect to TSS
# and from -<reach> to +5000 with respect to TTS are analysed.  
# genes=    a vector containing AGIs of genes from the analysed group. Usually the list is an element of <Groups>
# sample=    a number id of analysed sample, from 1 to 12, according to the samples_GEO_HiC.tsv
# reach=    should be less than half of the shortest gene analysed (usually reach=down_thresh/2)
# calculate_metamatrix() produces a list containing 2 meta-matrixes : raw Hi-C signal and norm - normalized to distance,
# and the number of genes that were actually analysed
calculate_metamatrix <- function(genes,sample,reach,res=100){
#   genes=c('AT1G03080','AT2G39800',
#                        'AT1G09970','AT1G62660','AT4G16760',
#                        'AT4G18950','AT4G23050','AT5G10650','AT5G51070')
# sample='1'
#   reach=1500
#   res=100
  # initialize empty meta_matrixes
  genes <- Groups[[1]]
  meta_pos_e=create_meta_pos(reach,res,gLength=(reach*2))
  matrix_raw_sum <- matrix(data = 0, nrow = nrow(meta_pos_e), ncol = nrow(meta_pos_e))
  matrix_norm_sum <- matrix(data = 0, nrow = nrow(meta_pos_e), ncol = nrow(meta_pos_e))
  # read the GMatrix raw data
  data <- readRDS(file=paste0(data_dir,'GMatrix_Sample',sample,'.RDS'))
  data <- data[which(names(data) %in% genes)]
  if(length(data)==0){next}
  genes_an <- names(data)
  TTSpos <- all_genes %>%  subset(GeneID %in% genes_an) %>%
    transmute(GeneID, TTSrel=floor((end-start)/res)*res)
 
  for (n in 1:nrow(TTSpos)){
 
      GeneID=TTSpos$GeneID[n]
    geneLength=TTSpos$TTSrel[n] 
    meta_pos= create_meta_pos(reach, res,geneLength)
    matrix_raw <- matrix(data = 0, nrow = nrow(meta_pos), ncol = nrow(meta_pos))
    matrix_norm <- matrix(data = 0, nrow = nrow(meta_pos), ncol = nrow(meta_pos))
    data_lowerTriangle <- data[[GeneID]] %>%
      subset(!pos1==pos2) %>%
      transmute(GeneID,relTSS1=pos2,
                relTSS2=pos1, norm_signal,raw_signal)
    
    data_s <- data[[GeneID]]%>%
      transmute(GeneID, relTSS1=pos1, relTSS2=pos2,norm_signal,raw_signal) %>%
      bind_rows(data_lowerTriangle) %>%
        right_join(meta_pos, by=c('relTSS1'='relTSS')) %>%
      transmute(meta1=pos, raw_signal, norm_signal, relTSS2) %>%
      right_join(meta_pos, by=c('relTSS2'='relTSS')) %>%
      transmute(meta2=pos, meta1,raw_signal, norm_signal)

  
    for (j in 1:nrow(data_s)) {
      matrix_raw[data_s$meta1[j],data_s$meta2[j]] <- data_s$raw_signal[j]
      matrix_norm[data_s$meta1[j],data_s$meta2[j]] <- data_s$norm_signal[j]
      
    }
    
    matrix_raw_sum <- matrix_raw_sum+matrix_raw
    matrix_norm_sum <- matrix_norm_sum+matrix_norm
    
  }
  
  matrix_raw_mean=matrix_raw_sum/length(genes_an)
  matrix_norm_mean=matrix_norm_sum/length(genes_an)
  matrix_norm_mean[which(meta_pos_e$label=='GeneBody'),] <- NA
  matrix_norm_mean[,which(meta_pos_e$label=='GeneBody')] <- NA
  matrix_raw_mean[which(meta_pos_e$label=='GeneBody'),] <- NA
  matrix_raw_mean[,which(meta_pos_e$label=='GeneBody')] <- NA
  
  rownames(matrix_raw_mean)=meta_pos_e$label
  colnames(matrix_raw_mean)=meta_pos_e$label
  rownames(matrix_norm_mean)=meta_pos_e$label
  colnames(matrix_norm_mean)=meta_pos_e$label
  
  list(mat_raw=matrix_raw_mean, mat_norm = matrix_norm_mean, nr_of_genes = length(genes_an))
  
}

singleGene_map <- function(sample,genes,res=100,format){
 
  TTSpos <- all_genes %>%  transmute(GeneID, TTSrel=floor((end-start)/res)*res, bins_nr=100+TTSrel/res)
  data <- readRDS(file=paste0(data_dir,'GMatrix_Sample',sample,'.RDS'))
  data <- data[which(names(data) %in% genes)]
  
  arrays <- lapply(names(data), function(GeneID){
   
    data_lowerTriangle <- data[[GeneID]] %>%
      subset(!pos1==pos2) %>%
      transmute(GeneID,relTSS1=pos2,
                relTSS2=pos1, norm_signal,raw_signal)
    
    data_s <- data[[GeneID]]%>%
      transmute(GeneID, relTSS1=pos1, relTSS2=pos2,norm_signal,raw_signal) %>%
      bind_rows(data_lowerTriangle)
    nr_bins <- TTSpos$bins_nr[which(TTSpos$GeneID==GeneID)]
    array <- matrix(data=0,
                    ncol = nr_bins,
                    nrow = nr_bins)
    rownames(array)=seq(-4900,(nr_bins-50)*100,100)
    colnames(array)=seq(-4900,(nr_bins-50)*100,100)
    if (format=='norm'){
      for(i in 1:nrow(data_s)){
        array[which(rownames(array)==data_s$relTSS1[i]),which(colnames(array)==data_s$relTSS2[i])] <- data_s$norm_signal[i]
      }
      
    } else if (format=='raw') {
      
      for(i in 1:nrow(data_s)){
        array[which(rownames(array)==data_s$relTSS1[i]),which(colnames(array)==data_s$relTSS2[i])] <- data_s$raw_signal[i]
      }
    }
    return(array)
  })
  names(arrays) <- names(data)
  return(arrays)
}

f_draw_matrix <- function(raw_array=NULL, norm_array=NULL, delta_array=NULL,add_frame=T, do_scaling=T,title='Metagene Matrix', threshold=0.98, max_scale_raw=NULL, max_scale_norm=NULL, max_scale_delta=NULL){
 # browser()
    if(!is.null(raw_array)){
    if(is.null(max_scale_raw)){ max_scale=quantile(raw_array,probs=threshold,na.rm=T)} else max_scale=max_scale_raw
    col_breaks_raw = seq(0,max_scale, length=12)
    raw_palette <- c('white',scico(10,palette='lajolla',direction=1))
    #raw_palette <- c('white',scico(30,palette='lajolla',direction=1))
    if(add_frame){frame_pos=c(50,nrow(raw_array)-50)}
    par(cex.main=0.7)
    heatmap.2(raw_array, col=raw_palette, breaks=col_breaks_raw,na.color='grey', colsep=frame_pos,rowsep=frame_pos,sepcolor="black", cexRow=0.3, cexCol=0.3,dendrogram='none',density.info="none", Rowv=FALSE, Colv=FALSE,trace='none',main=paste0(title,'\nRaw'))
  }
  if(!is.null(norm_array)){
     if(do_scaling){
       totalMean=mean(norm_array, na.rm=T)
      totalSD=sd(norm_array, na.rm=T)
      norm_array <- (norm_array-totalMean)/totalSD
    }
    if(is.null(max_scale_norm)){ max_scale=quantile(norm_array,probs=threshold,na.rm=T)  } else max_scale=max_scale_norm
    
    #norm_palette <- colorRampPalette(c("black","darkblue","cyan","white","white","yellow","firebrick1","darkred"))(n=31)
    norm_palette <- scico(31,palette='roma',direction=-1)
    col_breaks_norm = seq(-1*max_scale,max_scale, length=32)
    if(add_frame){frame_pos=c(50,(nrow(norm_array)-50))}
    par(cex.main=0.7)
    heatmap.2(norm_array, col=norm_palette, na.color='white', colsep=frame_pos,rowsep=frame_pos,sepcolor="black", breaks=col_breaks_norm,cexRow=0.3, cexCol=0.3,density.info="none",dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',main=paste0(title,'\nNormalized to distance'))
  }
  if(!is.null(delta_array)){
     if(do_scaling){
       totalMean=mean(delta_array, na.rm=T)
      totalSD=sd(delta_array, na.rm=T)
      delta_array <- (delta_array-totalMean)/totalSD
     }
    if(is.null(max_scale_delta)){ max_scale=quantile(delta_array,probs=threshold,na.rm=T)} else max_scale=max_scale_delta
    col_breaks_delta = seq(-1*max_scale,max_scale, length=32)
    delta_palette <- scico(palette='vikO',n=31)
    if(add_frame){frame_pos=c(50,(nrow(delta_array)-50))}
    par(cex.main=0.7)
    heatmap.2(delta_array, col=delta_palette,  breaks=col_breaks_delta, na.color='white', colsep=frame_pos, rowsep=frame_pos, sepcolor="black", cexRow=0.3, cexCol=0.3,dendrogram='none',density.info="none", Rowv=FALSE, Colv=FALSE,trace='none',main=paste0(title,'\nDifferential Map'))
  }
  
}

bootANDplot <- function(data,Groups,merge_samples){

  boot <- lapply(1:length(Groups),function(G){
    dat <- data %>%
      subset(GeneID %in% Groups[[G]]) %>%
      subset(sample %in% paste0('Sample',merge_samples)) %>%
      group_by(GeneID) %>%
      summarize(TSSprox=mean(TSSprox,na.rm=T),
                TTSprox=mean(TTSprox,na.rm=T)) %>%
      ungroup()
    f_ratio <- function(dat,idx){
      mean_ratio = log2(mean(dat$TTSprox[idx],na.rm=T)/mean(dat$TSSprox[idx]))
      return(mean_ratio)
    }
    bootstrap <- boot(dat, f_ratio, R = 100)
    a <- boot.ci(boot.out = bootstrap, index=1, type=c('norm'))
    min <- a$normal[1,2]
    max <- a$normal[1,3]
    return(tibble(mean=bootstrap$t0,
                  min=a$normal[1,2],
                  max = a$normal[1,3],
                  group=names(Groups)[G]))
  }) %>% reduce(bind_rows)
  boot_order <- boot %>%
    arrange(mean) %>% pull(group)
  boot$group <- factor(boot$group, levels=boot_order)
  
  kolory<- c("#8C0172", "#97523A","#4e4e4e","#99A021", "#29bc7e")
  
  boot_plot <-boot %>%
    ggplot(aes(x=group,y=mean,color=group))+
    geom_point(size=2, position=position_dodge(width=0.5))+
    geom_errorbar(aes(ymin=min, ymax=max) ,position=position_dodge(width=0.5), width=0.2,size=1)+
    scale_y_continuous(limits=c(min(boot$min),max(boot$max)),breaks=seq(-1,4,0.25))+
    scale_color_manual(values=kolory)+
    ylab('Asymmetry Ratio')+
    newtheme+
    theme(aspect.ratio=4,
          axis.text.y=element_text(angle = 90,hjust=0.5),
          axis.text.x=element_text(angle = 90,hjust=1),
          strip.text.y=element_text(angle = 90))
  table <- boot %>%
    tableGrob(theme=ttheme_default(base_size=10))
  return(list(boot_plot,table))
}

make_contact_summary <- function(data,Groups){
  summary <- lapply(names(Groups), function(G){
    data %>% subset(GeneID %in% Groups[[G]]) %>%
      group_by(sample) %>%
      summarize(TSSprox=mean(TSSprox*100000),
                TTSprox=mean(TTSprox*100000),
                TSSTTS=mean(TSSTTS*100000), 
                mean_AsRatio=log2((TTSprox)/(TSSprox)),
                group=G,
                n=n())
  }) %>% reduce(bind_rows) %>% left_join(sample_info) %>% ungroup()
  summary$group <- factor(summary$group, levels=names(Groups))
  summary$stress <- factor(summary$stress, levels=c('drought','water'))
  return(summary) 
}

report_contact_summary <- function(contact_summary,Groups){
  sectors_plots <-lapply(names(Groups),function(G){
 
    ss <- summary %>%
      subset(group==G)%>%
      pivot_longer(c(2:4),names_to='type')
    ss %>%
      ggplot(aes(x=stress,y=value,color=genotype,fill=genotype))+
      geom_col(width=0.9,position=position_dodge2(padding=0.3),size=0.8 )+
      scale_color_manual(values=genotypes_colors,aesthetics = c("color","fill") )+
      #ylim(0,max(summary$mean)+0.01)+
      scale_shape_manual(values=c(21,1))+
      ylab('mean score')+
      facet_wrap(~type, scales='fixed')+
      newtheme+
      ggtitle(paste0(G))
  })
   names(sectors_plots) <- paste0('Sectors_',names(Groups))
  ratio_plots <- lapply(names(Groups),function(G){
    ss <- summary %>%
      subset(group==G)%>%
      ggplot(aes(x=genotype,y=mean_AsRatio,color=genotype,fill=genotype))+
      geom_jitter(aes(shape=plants),size=1.5,stroke=2,
                  position=position_jitterdodge(dodge.width=0.4,jitter.width=0.4))+
      #scale_y_continuous(limits=c(min(summary$mean_AsRatio)-0.01,max(summary$mean_AsRatio)+0.01),breaks=seq(-1,4,0.25))+
      scale_y_continuous(limits=c(-0.15,0.9),breaks=seq(-1,4,0.25))+
      scale_shape_manual(values=c(21,1))+
      scale_color_manual(values=genotypes_colors, aesthetics = c("color","fill"))+
      ylab('Asymmetry score')+
      facet_grid(group~stress, scales='fixed')+
      newtheme+
      theme(aspect.ratio=6,
            axis.text.y=element_text(angle = 90,hjust=0.5),
            axis.text.x=element_text(angle = 90,hjust=1),
            strip.text.y=element_text(angle = 90))
  })
   names(ratio_plots) <- paste0('AsymmetryRatio_',names(Groups))
  stat_plots <-lapply(names(Groups),function(g){
    summary %>%
      subset(group==g)%>%
      ggboxplot(x='genotype',y='mean_AsRatio',color='genotype')+
      geom_jitter(aes(shape=plants,color=genotype))+
      scale_color_manual(values=genotypes_colors, aesthetics = c("color","fill"))+
      stat_compare_means(label = "p.format",size=5,method='t.test',vjust=1)+
      #ylim(min(summary$mean)-0.01,max(summary$mean)+0.01)+
      scale_shape_manual(values=c(16,1))+
      ylab('Mean Asymmetry score')+
      facet_wrap(~stress, scales='fixed')+
      newtheme+
      ggtitle(paste0(g,' genes\nn=',min(summary$n[which(summary$group==g)]),'\n'))
  }) 
  names(stat_plots) <- paste0('Stats_',names(Groups))
  report <- list(sectors_plots,ratio_plots,stat_plots)
  names(report) <- c('sectors_plots','ratio_plots','stat_plots')
  return(report)
}

rpm_fun <- function(x, na.rm = TRUE) x / sum(x, na.rm) * 1000000

calculate_mean_rpm <- function(nor, filter_nuclear=T, filter_protCoding=T){
  nuclear <- nor$GeneID[-grep('AT[CM]',nor$GeneID)]
  if(filter_nuclear){nor <- subset(nor, GeneID %in% nuclear)}
if(filter_protCoding){nor <- subset(nor, GeneID %in% all_genes$GeneID)}


mean_rpm <- nor %>%
  mutate_if(is.numeric, rpm_fun) %>%
  pivot_longer(2:37, names_to='sample', values_to='rpm') %>%
  left_join(sample_info_3RNAseq[,c('sample','red_sample')]) %>%
  group_by(GeneID, red_sample) %>%
  summarize(mean_rpm=mean(rpm,na.rm=F)) %>%
  pivot_wider(values_from='mean_rpm', names_from='red_sample') %>%
  subset(complete.cases(.))


return(mean_rpm)
}

calculate_mean_normCTS <- function(nor, filter_nuclear=T, filter_protCoding=T){
library(DESeq2)
  nuclear <- nor$GeneID[-grep('AT[CM]',nor$GeneID)]
  if(filter_nuclear){nor <- subset(nor, GeneID %in% nuclear)}
if(filter_protCoding){nor <- subset(nor, GeneID %in% all_genes$GeneID)}
  
 
cts <- as.matrix(nor[, sample_info_3RNAseq$sample])
rownames(cts) <- nor$GeneID
ddsMF <- DESeqDataSetFromMatrix(countData = cts,
                              colData = sample_info_3RNAseq ,
                              design = ~ genotype + condition)
ddsMF <- estimateSizeFactors(ddsMF)
counts <- counts(ddsMF, normalized=TRUE)
normalized_counts <- as_tibble(counts) %>% mutate(GeneID=rownames(counts),.before=1) %>%
    pivot_longer(2:37, names_to='sample', values_to='rpm') %>%
  left_join(sample_info_3RNAseq[,c('sample','red_sample')]) %>%
  group_by(GeneID, red_sample) %>%
  summarize(mean_rpm=mean(rpm,na.rm=T)) %>%
  #filter(mean_rpm>5) %>%
  pivot_wider(values_from='mean_rpm', names_from='red_sample') %>%
  subset(complete.cases(.))

return(normalized_counts)
}

f_convertAndSaveCdata_combined <- function(list_of_samples){
  
  xgi.bed='/home/lien/data/TFIIB/Hi_C/HiTC/matrix_for_import/xgi_100kb.bed'
  #list_of_samples='R58C'
    data <-lapply(lists_samples[[list_of_samples]],function(sample){
   
 table <- read_tsv(paste0('/home/lien/data/TFIIB/Hi_C/HiTC/matrix_for_import/',sample,'_100kbWholeGenomeMatrix.tsv')) %>%
          select(-c(1))
 matrix<-table[,2:length(table)] %>% as.matrix()
 rownames(matrix)<-pull(table[1])
 return(matrix)
}) %>% Reduce('+',.)
  table <- as_tibble(data)
     sparse_matrix <-table %>%
       mutate(bin1=rownames(data)) %>%
  pivot_longer(1:(length(table)-1),names_to='bin2',values_to='score' )%>%
  transmute(bin1, bin2,score) %>%
  subset(score != 0)
write_tsv(sparse_matrix,'/home/lien/data/TFIIB/Hi_C/HiTC/tmp.tsv',col_names =FALSE)
con='/home/lien/data/TFIIB/Hi_C/HiTC/tmp.tsv'
Cdata=importC(con, xgi.bed, ygi.bed=NULL, allPairwise=T,rm.trans=FALSE, lazyload=FALSE)
saveRDS(Cdata,paste0('/home/lien/data/TFIIB/Hi_C/HiTC/',list_of_samples,'_100kbWholeGenomeMatrix'))
}

rowAny <- function(x) rowSums(x) > 0

genotypes_colors=c("Wt"="#3b3b3b","R58C"="#c87970")
kolory2=c("<25%"="#191d57","25-75%"="#949496",">75%"="#c1841e")

newtheme <- theme_minimal() + theme(text = element_text(size = 14),
                                    #aspect.ratio=1,
                                    #strip.text=element_blank(),
                                    axis.title=element_text(size = 10),
                                    #axis.text=element_text(size = 6),
                                    #axis.text.x=element_text(angle = 75,hjust=1),
                                    #axis.text.x=element_blank(),
                                    axis.line.x=element_line(size=1),
                                    axis.line.y=element_line(size=1),
                                    axis.ticks.y=element_line(size=1),
                                    legend.text=element_text(size = 10),
                                    legend.title=element_blank(),
                                    #legend.position = "top",
                                    legend.position = "none",
                                    legend.box = "vertical",
                                    panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank())
