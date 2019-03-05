setwd("~/Documents/MRC16")

require(XLConnect)
require(pheatmap)
require(ggplot2)
require(reshape)
require(StingRay)

reformat <- function (x,lab="\n", spt="_") { sapply(x, function (c) { paste(unlist(strsplit(as.character(c) , split=spt)),collapse=lab) }) }

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# BUILD INPUT DATA
wb <- loadWorkbook("data/tRNA abundance/original data/tRNA_steady_state_abundances_update.xlsx")
tRNAab <- list()
tRNAab$abundance <- readWorksheet(wb, sheet = getSheets(wb)[2],header=T, startRow=2, endRow=44)
tRNAab$abundance[,2:ncol(tRNAab$abundance)] <- sapply(tRNAab$abundance[,-1], as.numeric)
tRNAab$conditions <- readWorksheet(wb, sheet = getSheets(wb)[4],header=T)
tRNAab$conditions$label <- paste( reformat(tRNAab$conditions$treatment,lab=" ",spt="_"), " t=",tRNAab$conditions$time, sep="" )

# Transform into a matrix (for pheatmap)
    rawM <- as.matrix(tRNAab$abundance[,-1]) ; rownames(rawM) <- tRNAab$abundance$tRNA_anticodon
    M <- rawM[complete.cases(rawM),]
    colnames(M) <- tRNAab$conditions$label
    
    # log2 transformation ( since raw value = abundance(x)/abundance(wt) )
    log2M <- log2(M)
    # scaled (centred on 0 and scaled to 1)
    scaledM <- apply(M, 2, scale); rownames(scaledM) <- rownames(M)
    # scalation of the log2 data
    tM <- apply(log2M, 2, scale); rownames(tM) <- rownames(M)
    
    save(tM,file="results/tRNA.abundance/tM.Rda")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Gather information on anticodons
require(StingRay); data(TAB)
  
      # Table with aa and aa group for each anticodon
      d1 <- unique( merge(TAB$aminoacid.properties[,c("aa","group")], TAB$genetic.code[,c("aa","anticodon")], by="aa"))
      # Table with tRNA copy number for each anticodon
      d2 <- data.frame( table(TAB$genetic.code[,2]) ); colnames(d2) <- c("anticodon","tRNA_gene_copy_number")
      # Table with number of synonymous codons for each anticodon
      d3 <- data.frame(table(unique(TAB$genetic.code[,1:2])[,2])); colnames(d3) <- c("anticodon","n_isoaccepting_codons")
      # Table with mean tAI of codons accepting a given tRNA. Doesn't account for codon frequency nor expression, so not sure it's relevant to use
      TAB$tAI$codon <- tolower(TAB$tAI$codon)
      d4 <- ddply( arrange(merge(TAB$tAI[,1:2], unique(TAB$genetic.code[,1:2]), by="codon"),anticodon), .(anticodon), summarise, mean.tAI.syn.codons = mean(Scer, na.rm=T) )

  tRNA_info <- merge_recurse(list(d1,d2,d3,d4),by="anticodon")
  tRNA_info <- unique(tRNA_info[,-3]); rownames(tRNA_info) <- tRNA_info$anticodon; tRNA_info <- tRNA_info[,-1]

save(tRNA_info, file="data/Lister/tRNA_info.Rda")
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Select relevant time points
  t20  <- subset(tRNAab$conditions, time==20 & type == "average" )$label
  t60  <- subset(tRNAab$conditions, time==60 & type == "average")$label
  t120 <- subset(tRNAab$conditions, time==120 & type == "average")$label

  # Clustering and heatmap
pdf("results/tRNA.abundance/response.to.stress.foldchange.abundance.pdf", width=15, height=5 )
  pUNS <- pheatmap(mat= t(tM[,t120]), scale="none", cluster_cols=T, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])
  pheatmap(mat= t(tM[,t60])[,pUNS$tree_col$order], scale="none", cluster_cols=F, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])
  pheatmap(mat= t(tM[,t20])[,pUNS$tree_col$order], scale="none", cluster_cols=F, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])
dev.off()

  # Scaled by column to highlight differences across conditions per anticodon
  pSCA <- pheatmap(mat= t(scaledM[,t120]), scale="column", cluster_cols=T, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])
  pheatmap(mat= t(scaledM[,t60])[,pSCA$tree_col$order], scale="column", cluster_cols=F, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])
  pheatmap(mat= t(scaledM[,t20])[,pSCA$tree_col$order], scale="column", cluster_cols=F, cluster_rows=T, cellwidth=15, cellheight=15, annotation=tRNA_info[,-1])


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Transform into a long-format dataframe
generate.D <- function(X=log2M, trna=tRNAab, method="ward"){
  require(StingRay)
  data(TAB)
  require(seqinr)
  
  d.x <- dist2(X)
    h.x <- hclust(d.x,method=method)
      o.x <- h.x$labels[h.x$order]
  d.y <- dist2(t(X))
    h.y <- hclust(d.y,method=method)
      o.y <- h.y$labels[h.y$order]
  trna$conditions$reference  <-  factor(trna$conditions$reference, levels=o.y)
  
  D <- melt(X); colnames(D) <- c("anticodon","variable","value")
  D <- merge(D, trna$conditions, by.x="variable", by.y="reference")
  D$anticodon <- factor(D$anticodon, levels= o.x )
  D$variable <- factor(D$variable, levels= o.y )
  D$treatment <- reformat(D$treatment, lab=" ")
  D$treatment <- factor(D$treatment, levels=reformat(unique(trna$conditions$treatment[h.y$order]),lab=" ")  )
  
  # Heavy and non-elegant approach to construct a master table
  TAB$anticodon2codon$anticodon <- chartr("U","T", TAB$anticodon2codon$anticodon)
  Add  <- unique(merge(TAB$anticodon2codon[,c(1,3)], TAB$tRNA.copy.number, by="codon"))
  Add2 <- merge(Add,TAB$nTE.codon.optimality[,c("codon","Scer")], by="codon")
    colnames(Add2)[ncol(Add2)] <- "optimal.nTE"
  Add3 <- merge(Add2, TAB$nTE[,c("codon","Scer")], by="codon")
    colnames(Add3)[ncol(Add3)] <- "nTE"
  Add4 <- merge(Add3, TAB$CAI[,c("codon","Scer")], by="codon")
    colnames(Add4)[ncol(Add4)] <- "CAI"
  Add5 <- merge(Add4, TAB$CAI.codon.optimality[,c("codon","optimality")], by="codon")
    colnames(Add5)[ncol(Add5)] <- "optimal.CAI"
  Add6 <- merge(Add5, TAB$tAI[,c("codon","Scer")], by="codon")
    colnames(Add6)[ncol(Add6)] <- "tAI"
  E <- merge(D,Add6, by="anticodon")
  E$AA <- sapply(E$codon, function(x){ aaa(translate(s2c(x)) )})
  E$anticodon <- as.character(E$anticodon)
  return(E)
}



coolColors <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)

  D <- generate.D(X=tM)
  ggplot(data=generate.D(X=tM), aes(x=time, y=value) ) + geom_violin(alpha=0.5) + geom_point() + facet_grid( . ~ treatment) + theme(strip.text.y = element_text(angle=0) )

  ggplot(data=subset(D, variable %in% t60), aes(x=treatment, y=anticodon, fill=value) ) + 
  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors)

  # Supply of tRNAs
  ggplot(data=D, aes(x=treatment, y=anticodon, fill=value) ) + ylab("tRNA's anticodon") +
  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors) + facet_grid(. ~ time) + theme(axis.text.x=element_text(angle=-45, hjust=0), legend.position="top") + labs(fill="Relative tRNA Abundance")

  # Supply of tRNA per type of amino-acid
  ggplot(data=D, aes(x=treatment, y=AA, fill=value) ) + ylab("anticodon") +  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors) + facet_grid(. ~ time,space="free_y",drop=T) + theme(axis.text.x=element_text(angle=-45, hjust=0), legend.position="top") + labs(fill="Relative tRNA Abundance")

  ggplot(data=D, aes(x=treatment, y=anticodon, fill=value) ) + ylab("anticodon") +  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors) + facet_grid(AA ~ time,space="free", scales="free",drop=T) + theme(axis.text.x=element_text(angle=-45, hjust=0), legend.position="top") + labs(fill="Relative tRNA Abundance")
