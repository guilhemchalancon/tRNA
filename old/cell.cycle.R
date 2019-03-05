setwd("~/Documents/MRC16/")

CC.alpha <- read.table(file="data/Cell cycle (high resolution)/marina_alpha_orfs.tsv",sep="\t", row.names=1,header=1,check.names=F)
CC.cdc28 <- read.table(file="data/Cell cycle (high resolution)/marina_cdc28_orfs.tsv",sep="\t", row.names=1,header=1,check.names=F)

CC.alpha.q20 <- quantilize(CC.alpha,name="CC.alpha",bin.number=20)
CC.cdc28.q20 <- quantilize(CC.cdc28,name="CC.cdc28",bin.number=20)

CC.alpha.top20 <- select.bins(CC.alpha.q20,c("0-5","5-10","10-15","15-20"))
CC.alpha.bottom20 <- select.bins(CC.alpha.q20,c("80-85","85-90","90-95","95-100"))

CC.cdc28.top20 <- select.bins(CC.cdc28.q20,c("0-5","5-10","10-15","15-20"))
CC.cdc28.bottom20 <- select.bins(CC.cdc28.q20,c("80-85","85-90","90-95","95-100"))

save(CC.alpha.top20,file="data/Rdata/CC.alpha.top20.up_regulated.ORFs.Rda")
save(CC.cdc28.top20,file="data/Rdata/CC.cdc28.top20.up_regulated.ORFs.Rda")

M.alpha <- refine.matrices(CC.alpha.top20);
M.cdc28 <- refine.matrices(CC.cdc28.top20);

require(som)
M.alpha.n <- normalize(M.alpha[[2]]); colnames(M.alpha.n) <- colnames(M.alpha[[2]])
M.cdc28.n <- normalize(M.cdc28[[2]]); colnames(M.cdc28.n) <- colnames(M.cdc28[[2]])

##########################################################################################
query <- c("LOS1","GCN4","MTR10","SOL1","SOL2","CEX1","UTP8","TEF1","TEF2","RSP5")
CC.alpha.highlights <- highlights(query,CC.alpha)
CC.cdc28.highlights <- highlights(query,CC.cdc28)

pdf("results/cell.cycle/cell.cycle_codon.usage.pdf")
par(cex.main=0.8)
M <- M.alpha.n ; hM.alpha <- heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Codon usage during cell cycle\n(alpha factor G1 arrest)\nUp-regulated genes (top 20%)",LABCOL=rdc(rownames(M)),a=3,b=8,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=CC.alpha.highlights,CENTRE=0,SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="expression levels (log2)\n(alpha factor G1 arrest)\n highlights",a=7,b=5,METHOD=cm[1], centered.break=T)

M <- M.cdc28.n ; hM.cdc28 <- heatmap.3(matrix=M,SCALE="none",ROWV=T,COLV=F,SEPW=NA,main="Codon usage during cell cycle\n(38ºC heat shock in cdc28-13)\nUp-regulated genes (top 20%)",LABCOL=rdc(rownames(M)),a=3,b=8,METHOD=cm[1],RowSideColors=col.codons)
heatmap.3(matrix=CC.cdc28.highlights,CENTRE=0,SCALE="none", ROWV=T,COLV=F,SEPW=NA,main="expression levels (log2)\n(38ºC heat shock in cdc28-13)\n highlights",a=7,b=5,METHOD=cm[1], centered.break=T)
dev.off()

