# This script is deprecated #
setwd("~/Documents/MRC16")
require(XLConnect)
require(pheatmap)
require(ggplot2)
require(reshape)
require(StingRay)

# Data collected by Marc Torrent (Old)
wb <- loadWorkbook("data/tRNA abundance/original data/tRNA_steady_state_abundances_old.xlsx")
tRNAab <- readWorksheet(wb, sheet = getSheets(wb)[2:4],header=T)
tRNAab$abundance[,2:ncol(tRNAab$abundance)] <- sapply(tRNAab$abundance[,-1], as.numeric)

# Transform into a matrix (for pheatmap)
rawM <- as.matrix(tRNAab$abundance[,-1]) ; rownames(rawM) <- tRNAab$abundance$tRNA_anticodon
M <- rawM[complete.cases(rawM),]

# log2 transformation ( since raw value = abundance(x)/abundance(wt) )
log2M <- log2(M)
# scaled (centred on 0 and scaled to 1)
scaledM <- apply(M, 2, scale); rownames(scaledM) <- rownames(M)
# scalation of the log2 data
tM <- apply(log2M, 2, scale); rownames(tM) <- rownames(M)

  # Select relevant time points
  t20 <- subset(tRNAab$conditions, time=="20min" & type == "average" & treatment != "temperature_shift_(30-45)")$reference
  t60 <- subset(tRNAab$conditions, time=="60min" & type == "average")$reference
  t120 <- subset(tRNAab$conditions, time=="120min" & type == "average")$reference
  # Clustering and heatmap
  pheatmap(mat= tM[,t20], scale="none", cluster_cols=F, cellwidth=25, cellheight=15)
  
  pheatmap(mat= M[,t20], scale="row", cluster_cols=F, cellwidth=25, cellheight=15)
  pheatmap(mat= scaledM[,t60], scale="row", cluster_cols=F, cellwidth=25, cellheight=15)

# Transform into a long-format dataframe
generate.D <- function(X=log2M, trna=tRNAab, method="ward.D"){
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
  Add5$codon <- tolower(Add5$codon)
  Add6 <- merge(Add5, TAB$tAI[,c("codon","Scer")], by="codon")
    colnames(Add6)[ncol(Add6)] <- "tAI"
  E <- merge(D,Add6, by="anticodon")
  E$aa <- sapply(E$codon, function(x){ aaa(seqinr::translate(s2c(x)) )})
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
  ggplot(data=D, aes(x=treatment, y=aa, fill=value) ) + ylab("anticodon") +  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors) + facet_grid(. ~ time,space="free_y",drop=T) + theme(axis.text.x=element_text(angle=-45, hjust=0), legend.position="top") + labs(fill="Relative tRNA Abundance")

  ggplot(data=D, aes(x=treatment, y=anticodon, fill=value) ) + ylab("anticodon") +  geom_tile(alpha=0.5) + geom_tile() + scale_fill_gradientn(colours=coolColors) + facet_grid(AA ~ time,space="free", scales="free",drop=T) + theme(axis.text.x=element_text(angle=-45, hjust=0), legend.position="top") + labs(fill="Relative tRNA Abundance")
