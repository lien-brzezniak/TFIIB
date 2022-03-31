# dpn_data <- read_tsv('/home/lien/data/TFIIB/Hi_C/TagDirs/matrices/raw_relative/GeneMatrix_allSizes/slices/dpn2metaSlice.tsv')
# dpn <- lapply(genes,function(x){
#   # x=genes[1]
#   values <- dpn_data[which(dpn_data$GeneID==x),2:length(dpn_data)] %>% as.double()
#   values <- values[which(!is.na(values))]
# })
# names(dpn) <- genes
# dpnM <- matrix(0,ncol=max(lengths(dpn)),nrow=length(dpn))
# rownames(dpnM) <- genes
# for(x in genes){
#   # Cx=genes[1]
#   dpnM[x,] <- c(dpn[[x]],rep(NA,max(lengths(dpn))-length(dpn[[x]])))
# }
# dpn_palette <- colorRampPalette(c("white","black"))(n=2)
# heatmap.2(dpnM,sepcolor='black', colsep=0,rowsep=0,col=dpn_palette, breaks=c(0,0.5,1.5),na.color='blue',  cexRow=0.3, cexCol=0.3,key=T,density.info="none",dendrogram='none',Rowv=FALSE, Colv=FALSE,trace='none',main=paste0('DpnI sites'))
