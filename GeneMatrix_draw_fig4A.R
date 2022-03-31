library(gplots)
library(tidyverse)


 data_dir='/home/lien/data/TFIIB/figury/Fig4/EIipR1645810458/'
out_dir=paste0(data_dir,'draw/')
system(paste0('mkdir ',out_dir),intern=F)


#barcode=paste0(paste0(sample(c(letters,LETTERS),5,replace=T),collapse=""),as.integer(as.POSIXct(Sys.time())))
res=100
reach=1500
meta_pos <- tibble(label=
                     c(as.character(seq(-4900,-100,100)),'TSS',as.character(seq(100,reach,100)),
                       rep('GeneBody',1),
                       as.character(seq(-1*(reach),-100,100)),'TTS',as.character(seq(100,4900,100))),
                   pos=c(1:((reach/res*2)+101)),
                   relTSS=c(seq(-4900,reach,100), rep(NA,(51+(reach/res)))), 
                   relTTS=c(rep(NA,(51+(reach/res))), seq(-1*reach,4900,100)))
frame_pos=c(meta_pos$pos[which(meta_pos$label=='TSS')],meta_pos$pos[which(meta_pos$label=='TTS')+1])
groups=c('drought_genes')

# samples1=c('1','2','3','4','9','10','11','2','10','9','2','6')
# samples2=c('5','6','7','8','5','6','7','14','14','1','10','14')
samples=c(1:11,14)
read_matrix <- function(x,descr) {
  norm_matrix <- read_tsv(paste0(data_dir,'Meta_norm_',x,'_',descr,'_3to50kb_nonscaled.txt'), col_names=F) %>% 
       as.matrix() 
  aver=mean(norm_matrix,na.rm=T)
  dev=sd(norm_matrix,na.rm=T)
  norm_matrix <-( norm_matrix-aver)/dev
  # raw_matrix <- read_tsv(paste0(data_dir,'Meta_raw_',x,'_',descr,'_3to50kb_nonscaled.txt'), col_names=F) %>% 
  #      as.matrix() 
  return(norm_matrix)
}
norm_palette <- colorRampPalette(c("black","darkblue","cyan","white","white","yellow","firebrick1","darkred"))(n=31)

matrices <- sapply(groups,function(x){
 sapply(samples, read_matrix, simplify=F,x, USE.NAMES=T)
})

  
threshold=0.995
 max_norm=quantile(unlist(matrices),probs=threshold,na.rm=T)

pdf_title=paste0(out_dir, 'maxthresh_',threshold,'.pdf')
pdf(pdf_title)
lapply(c(1:length(matrices)),function(x){
 
  # aver=mean(matrices[[x]],na.rm=T)
  # dev=sd(matrices[[x]],na.rm=T)
  # norm_map <-( matrices[[x]]-aver)/dev
  # max_norm=quantile(norm_map,probs=threshold,na.rm=T)
  norm_map <- matrices[[x]]
  col_breaks_norm = seq(-1*max_norm,max_norm, length=32)
  heatmap.2(norm_map,col=norm_palette, breaks=col_breaks_norm,colsep=frame_pos,rowsep=frame_pos,sepcolor='black', na.color='grey',cexRow=0.2, cexCol=0.2,density.info="none",dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',main=paste0('Thresh: ',threshold,'\n',names(matrices[x])))
 
  # heatmap.2(filtered_map,colsep=c(50,82),rowsep=c(50,82),sepcolor='black',col=norm_palette, breaks=col_breaks_norm, na.color='grey',cexRow=0.2, cexCol=0.2,density.info="none",dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',main=paste0(descr,'\nThresh: ',threshold,'\n',names(matrices[x])))
  })
dev.off()



