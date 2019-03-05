# Monitor stresses
Mg.all_10 <- lapply( gash.matrix_10, function(x) codon.subset(x, c("tgg", "atg")) )
Mg.all_20 <- lapply( gash.matrix_20, function(x) codon.subset(x, c("tgg", "atg")) )

pdf("results/gash/codon.usage.during_stresses.pdf", title="Codon usage during stresses")
par(cex.main=0.8)
# AA starvation
M <- Mg_20.n[[2]][,91:95]; heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\naa starvation\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,91:95]; heatmap.3(matrix=M,SCALE="row",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\naa starvation\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,91:95], CENTRE=0, SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="aa starvation\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1],centered.break=T)

# Nitrogen depletion
M <- Mg_20.n[[2]][,96:105] ;heatmap.3(matrix=M,CENTRE=1, SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nNitrogen depletion\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,96:105] ;heatmap.3(matrix=M,CENTRE=1, SCALE="row",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nNitrogen depletion\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,96:105],CENTRE=0,SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="Nitrogen depletion\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1], centered.break=T)

# Diauxic shift
M <- Mg_20.n[[2]][,106:112] ;heatmap.3(matrix=M,SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nDiauxic shift\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,106:112] ;heatmap.3(matrix=M,SCALE="row", ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nDiauxic shift\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,106:112],CENTRE=0, SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="Diauxic shift\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1],centered.break=T)
#plot(M[hMg_20$colInd,5])

# YPD
M <- Mg_20.n[[2]][,113:122] ;heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nYPD\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,113:122] ;heatmap.3(matrix=M,SCALE="row",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nYPD\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,113:122],CENTRE=0, SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="YPD\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1],centered.break=T)

# Heat shock
M <- Mg_20.n[[2]][,1:8] ;heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nHeat shock\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,1:8] ;heatmap.3(matrix=M,SCALE="row",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nHeat shock\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,1:8],CENTRE=0, SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="Heat shock\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1],centered.break=T)

# Sorbitol
M <- Mg_20.n[[2]][,30:35] ;heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes\nSorbitol - top 20%\nno scaling",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
M <- Mg_20.n[[2]][,30:35] ;heatmap.3(matrix=M,SCALE="row",ROWV=T,COLV=F,SEPW=NA,main="Up-regulated genes - top 20%\nSorbitol\ncolor relative to codon",LABCOL=rdc(rownames(M)),a=7,b=5,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=M.tRNA[,30:35],CENTRE=0, SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="Sorbitol\n highlights\nexpression levels (log2)",a=7,b=5,METHOD=cm[1],centered.break=T)

dev.off()
