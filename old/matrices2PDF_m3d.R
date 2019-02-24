# First raw analysis
pdf("results/m3d/codon.usage.average_top10.pdf", title="Codon usage (averaged among top10 highly expressed ORFs/condition)")

boxplot(m3d.log,axes=F,main="normalised expression levels, log transformed - all ORFs")
box()
axis(2)

show.binning(m3d.q10,m3d.log)
heatmap(m3d.matrix_05[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 5% ORFs/condition)")
heatmap(m3d.matrix_10[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 10% ORFs/condition)")
heatmap(m3d.matrix_15[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 15% ORFs/condition)")
heatmap(m3d.matrix_20[[1]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="Codon frequency (top 20% ORFs/condition)")

heatmap(m3d.matrix_05[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 5% ORFs/condition)")
heatmap(m3d.matrix_10[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 10% ORFs/condition)")
heatmap(m3d.matrix_15[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 15% ORFs/condition)")
heatmap(m3d.matrix_20[[2]],hclust=function(c){hclust(c, method=cm[1])},cexCol=0.3,cexRow=0.6,main="RSCU (top 20% ORFs/condition)")
dev.off()



# Refined analysis
experiments_to_remove <- c("Rpd3 and H3 H4 deletions","lfh1","Demos", "Rap1 and Abf1 mutants", "TFIIH mutants with methyl methanesulfonate", "NSM mRNA decay","allele specific expression","PMI40 deletion" ,"metabolic cycle")

############# Generate dataset
m3d.top10.readable <- m3d.top10; names(m3d.top10.readable) <- experiments[names(m3d.top10), 2]
m3d.top20.readable <- m3d.top20; names(m3d.top20.readable) <- experiments[names(m3d.top20), 2]

M_10 <- generate.matrices(m3d.top10.readable); 
M_10 <- lapply( M_10, function(x) experiment.subset(x, experiments_to_remove,2) )
M_10 <- lapply( M_10, function(x) codon.subset(x, c("tgg -> Trp", "atg -> Met")) )

M_20 <- generate.matrices(m3d.top20.readable); 
M_20 <- lapply( M_20, function(x) experiment.subset(x, experiments_to_remove,2) )
M_20 <- lapply( M_20, function(x) codon.subset(x, c("tgg -> Trp", "atg -> Met")) )

pdf("results/codon.usage.average_hcluster_rscu.pdf", title="RSCU (averaged among top10 highly expressed ORFs/condition)")
heatmap.3(matrix=M_10[[2]],ROWV=T,COLV=T,SEPW=NA,main="RSCU (top 10% ORFs/condition)",a=7,b=5,METHOD=cm[1], KEY=F)
heatmap.2(M_10[[2]])
heatmap.3(matrix=M_20[[2]],ROWV=T,COLV=T,SEPW=NA,main="RSCU (top 20% ORFs/condition)",a=7,b=5,METHOD=cm[1])
dev.off()

library(apcluster)
x1 <- (M_10[[2]])
s1 <- negDistMat(x1, r=2)
ranges <- round(preferenceRange(s1,exact=F),2)

apres1a <- apcluster(s1, q=0.5);plot(apres1a, s1)
apres2a <- apclusterK(s1, K=2); plot(apres2a, s1)


pdf("results/m3d/codon.usage.average_apcluster_sim.pdf", title="Similarity Matrix (averaged among top10 highly expressed ORFs/condition)")
plot(apres1a, s1)
plot(apres2a, s1)
dev.off()