
Mat <- t(w.Mg_20.n[,c("ypd.exp.growth",conditions_to_keep)]) ; h.w.Mg_20 <- heatmap.3(matrix=Mat,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1])
colc.w <- color.clusters(h.w.Mg_20)
colc.w$conditions  <- replace(colc.w$conditions, which(colc.w$conditions=="green2"), "red3");
colc.w$conditions  <- replace(colc.w$conditions, which(colc.w$conditions=="red2"), "green2"); 

CODONS <- colnames(M[, h.w.Mg_20$colInd])

pdf(file="results/gash/checkings.23Feb.pdf")

# PAGE 1
heatmap.3(matrix=Mat,SCALE="none",ROWV=T,COLV=T,SEPW=NA,main="Up-regulated genes (top 20%)\ncodon usage (RSCU, normalized)",LABROW=rdc(colnames(M)),a=3,b=7,METHOD=cm[1], ,ColSideColors=colc.w$codons, RowSideColors=colc.w$conditions)

# PAGE 2
boxplot(test.CAI$CAI ~ test.CAI$cluster, ylab="Codon Adaptation Index", main="Clusters based on codon usage\n(Up-regulated genes, top 20%)", notch=F, col=names(table(col.codons)), boxwex=0.45)

# PAGE 3
plot(wt.matrix.rscu[ CODONS ], type="l", lwd=2, xlab="", ylab="", ylim=c(0,6), col="grey", axes=F )
par(new=T)
plot(Mg_20[[2]][ CODONS,"steady_state_25_dec_C_ct-2" ], type="l", lwd=2, xlab="", ylab="", ylim=c(0,6), col="skyblue", axes=F )
par(new=T)
plot(Mg_20[[2]][ CODONS ,"Nitrogen_Depletion_30_min." ], type="l", lwd=2, xlab="", ylab="codon usage (RSCU)", ylim=c(0,6), col="orange",  axes=F)
axis(1, at = 1:length(CODONS), labels = formatC(CODONS, format="fg"), side=1, pos=-0.35, las=2, cex.axis=0.75)
axis(2, at = 0:6, labels = formatC(0:6, format="fg"), pos=-0.35, las=2, cex.axis=0.75)

# PAGE 4
plot(w.Mg_20.n[ CODONS,"steady_state_25_dec_C_ct-2" ], type="l", lwd=2, xlab="", ylab="", ylim=c(-4,4), col="skyblue", axes=F )
#abline(h=0, col="#0000045",lwd=5)
par(new=T)
my.plot <- plot(w.Mg_20.n[ CODONS ,"Nitrogen_Depletion_30_min." ], type="l", lwd=2, xlab="", ylab="Normalised codon usage (RSCU)", ylim=c(-4,4), col="orange",  axes=F)
par(new=T)
plot(w.Mg_20.n[ CODONS, "ypd.exp.growth" ], type="l", lwd=2, xlab="", ylab="", ylim=c(-4,4), col="grey", axes=F )
axis(1, at = 1:length(CODONS), labels = formatC(CODONS, format="fg"), side=1, pos=-4.35, las=2, cex.axis=0.75)
axis(2, at = -4:4, labels = formatC(-4:4, format="fg"), pos=-0.35, las=2, cex.axis=0.75)


dev.off()



