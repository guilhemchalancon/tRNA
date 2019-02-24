
require(NeatMap)
make.circularmap((M.tRNA),metric="euclidean",cluster.method="complete.linkage", normalize.profiles=FALSE)
make.heatmap1(M.tRNA[,91:95],row.method="PCA",column.method="average.linkage",row.normalize=T,column.metric="euclidean")
#make.profileplot3d(Mg_20[[2]],row.method="PCA",column.method="average.linkage")

# First raw analysis
# pdf("results/gash/codon.usage.average_top10.pdf", title="Codon usage (averaged among top10 highly expressed ORFs/condition)")
# 
# boxplot(gash.all,axes=F,main="normalised expression levels, log transformed - all ORFs")
# box()
# axis(2)
# 
# show.binning(gash.q10,gash.all)
# heatmap(gash.matrix_05[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 5% ORFs/condition)")
# heatmap(gash.matrix_10[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 10% ORFs/condition)")
# heatmap(gash.matrix_15[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 15% ORFs/condition)")
# heatmap(gash.matrix_20[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 20% ORFs/condition)")
# 
# heatmap(gash.matrix_05[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 5% ORFs/condition)")
# heatmap(gash.matrix_10[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 10% ORFs/condition)")
# heatmap(gash.matrix_15[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 15% ORFs/condition)")
# heatmap(gash.matrix_20[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 20% ORFs/condition)")
# dev.off()

# Refined analysis
setwd("~/Documents/MRC16")
source("scripts/conditions.R")
conditions_to_remove <- setdiff(colnames(gash.all),conditions_to_keep)


# Generate matrices with selected conditions, Trp and Met removed  
Mg_10 <- refine.matrices(gash.top10,conditions_to_keep);
Mg_20 <- refine.matrices(gash.top20,conditions_to_keep);
Mg_80 <- refine.matrices(gash.bottom20,conditions_to_keep);
Mg_90 <- refine.matrices(gash.bottom10,conditions_to_keep);

#####################################################

# Normalization
require(som)
Mg_20.n <- refine.matrices(gash.top20); tmp <- colnames(Mg_20.n[[2]])
Mg_20.n <- lapply(Mg_20.n, function (x) { normalize(x,byrow=T) } ); 
colnames(Mg_20.n[[1]]) <- tmp
colnames(Mg_20.n[[2]]) <- tmp

M <- t(Mg_20[[2]]) ; hMg_20a <- heatmap.3(matrix=M,SCALE="col",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, color scaling)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=col$codons, RowSideColors=col$conditions)
M <- t(Mg_20.n[[2]][,conditions_to_keep]) ; hMg_20b <- heatmap.3(matrix=M,CENTRE=1,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=col$codons, RowSideColors=col$conditions,centered.break=T)
M <- t(Mg_20.n[[2]][,conditions_to_keep]) ; hMg_20c <- heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=col.codons, RowSideColors=col.conditions)
cola <- color.clusters(hMg_20a)
colb <- color.clusters(hMg_20b)
colc <- color.clusters(hMg_20c)

# Testing color scaling and centring
pdf("results/gash/color.scaling.pdf")
# no normalisation, but color scaling (per column, i.e. codons)
M <- t(Mg_20[[2]]) ; hMg_20 <- heatmap.3(matrix=M,SCALE="col",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, color scaling)",LABROW=rdc(colnames(M)),a=3,b=8,METHOD=cm[1],ColSideColors=cola$codons, RowSideColors=cola$conditions)
# matrix normalisation + color centring at RSCU=1
M <- t(Mg_20.n[[2]][,conditions_to_keep]) ; hMg_20 <- heatmap.3(matrix=M,CENTRE=1,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=8,METHOD=cm[1],ColSideColors=colb$codons, RowSideColors=colb$conditions,centered.break=T)
# matrix normalisation, no color scaling
M <- t(Mg_20.n[[2]][,conditions_to_keep]) ; hMg_20 <- heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=colc$codons, RowSideColors=colc$conditions)
# hMg_20$carpet[rownames(Mg_20.n[[2]]),conditions_to_keep] - Mg_20.n[[2]][,conditions_to_keep]
boxplot(test.CAI$CAI ~ test.CAI$cluster, ylab="Codon Adaptation Index", main="Clusters based on codon usage\n(Up-regulated genes, top 20%)", notch=F, col=names(table(col.codons)), boxwex=0.5)
dev.off()



plot(hMg_20$breaks[1:length(palette)], col=palette,pch=16)
hist(M)

#####################################################


# Plot hierarchical clustering of average codon usage in up-regulated genes upon different stress conditions
op <- par(mar = par("mar")*2,oma=c(2,2,2,2))
pdf("results/gash/codon.usage.selected.conditions_rscu.pdf", title="RSCU (averaged among top10 highly expressed ORFs/condition)")
  M <- t(Mg_10[[2]]); heatmap.3(M,SCALE="col",ROWV=T,COLV=T,SEPW=NA,main="RSCU (top 10%)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=col.codons, RowSideColors=col.conditions)
  M <- t(Mg_20.n[[2]][,conditions_to_keep]) ; hMg_20 <- heatmap.3(matrix=M,CENTRE=1,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="RSCU (top 20%)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1],ColSideColors=col.codons, RowSideColors=col.conditions)
  M <- t(Mg_80[[2]][ hMg_20$colInd, rev(hMg_20$rowInd) ]) ; heatmap.3(M,SCALE="col",DENDRO="none",ROWV=F,COLV=F,SEPW=NA,main="RSCU (bottom 20%), same order",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1], ColSideColors=col.codons[hMg_20$colInd], RowSideColors=col.conditions[rev(hMg_20$rowInd)])
  M <- t(Mg_80[[2]]); heatmap.3(M,SCALE="col",ROWV=T,COLV=T,SEPW=NA,main="RSCU (bottom 20%)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1], ColSideColors=col.codons, RowSideColors=col.conditions)
dev.off()
par(op)


# What is the Codon Adaptation Index of the codons in each cluster? 
CAI <- caitab$sc; names(CAI) <- rownames(caitab)
test.CAI <- merge(clusters.codons, CAI, by=0); rownames(test.CAI) <- test.CAI$Row.names; test.CAI <- test.CAI[,-1]; colnames(test.CAI) <- c("cluster", "CAI")
pdf("results/gash/codon.adaptation.selected_conditions.pdf")
  par(cex.main=0.8)
  boxplot(test.CAI$CAI ~ test.CAI$cluster, ylab="Codon Adaptation Index", main="Clusters based on codon usage\n(Up-regulated genes, top 20%)", notch=F, col=names(table(col.codons)), boxwex=0.5)
dev.off()

# Message Passing clustering
library(apcluster)
x1 <- (Mg_10[[2]])
s1 <- negDistMat(x1, r=2)
ranges <- round(preferenceRange(s1,exact=F),2)
apres1a <- apcluster(s1, q=0.5);plot(apres1a, s1)
apres2a <- apclusterK(s1, K=2); plot(apres2a, s1)
aggres1a <- aggExCluster(s1); plot(aggres1a, s1)
cl1a <- cutree(aggres1a, k = 2)

pdf("results/gash/codon.usage.average_apcluster_sim.pdf", title="Similarity Matrix (averaged among top10 highly expressed ORFs/condition)")
plot(apres1a, s1)
plot(apres2a, s1)
dev.off()


# What genes are in there?
# require(Matrix)
# genes.in.Mg_10 <- matrix(NA, nr=length(orf_genes_seq), nc=ncol(Mg_10[[2]]))
# rownames(genes.in.Mg_10)
# genes.in.Mg_10 <- sapply( gash.top10, function(x) { genes.in.Mg_10[x, ] <- 1})
# orf_genes_seq
