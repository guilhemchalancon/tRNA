
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- COMMANDS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
require(StingRay)
source("scripts/commands.R")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- DATA ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
cell.cycle.regulators <- read.table("data/Cell cycle (high resolution)/SGD_regulation_of_cell_cycle_genes.txt", header=T, sep="\t", stringsAsFactors = F)
data(SEQ) # for conding sequences
data(TAB) # for scores on tRNA adaptability


setdiff( cell.cycle.regulators$gene, names(SEQ$ORF_CDS) ) # there is a sequence for every cell cycle regulator gene

# tRNA adaptation index
TAB$tAI$codon <- tolower(TAB$tAI$codon)


# counts of codons in all protein-coding genes in yeast
all.codons <- data.frame( gene = names(SEQ$ORF_CDS), do.call( rbind, lapply( names(SEQ$ORF_CDS), function(x) { uco(  s2c(toString(SEQ$ORF_CDS[[x]])), index = "rscu"  ) } ) ), stringsAsFactors=F)
save(all.codons, file="results/codon.counts.allORFs.Rda")

# clusters of codons we identified as preferred in top-20% most up-regulated genes in response to stress
cluster.codons <- read.table("results/lister/clusters_codons.txt",header=T, stringsAsFactor=F)


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- MASTER TABLE ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

CDS.length <- data.frame( gene = all.codons$gene, length = rowSums( all.codons[,-1] ) )
  
M1 <- merge(CDS.length, all.codons, by="gene")  
M2 <- melt(M1, id.vars=c("gene","length") ) 
  colnames(M2)[3:4] <- c("codon","rcsu")
M3 <- merge(M2, TAB$tAI[,c("codon","Scer")], by = "codon", all.x=T)
  colnames(M3)[5] <- "tAI"
M4 <- merge(M3, unique( TAB$genetic.code[, c("codon","aa")] ), by = "codon", all.x =T) 
M5 <- merge(M4, cluster.codons, by = "codon")

# Formatting and refining
  #M5$freq  <-  round( M5$rscu / M5$length, 2)
  M5$name  <- fNAME(M5$gene)
  M5 <- M5[,c("name","gene","length","aa","codon","up_in","rcsu","tAI")]
  M5 <- subset( M5, !codon %in% c("taa","tag","tga") )
  ordered.codons <- arrange( merge( cluster.codons, subset(TAB$tAI[,c("codon","Scer")], !is.nan(Scer) ), by = "codon" ) , up_in, Scer )$codon
  M5$codon <- factor( M5$codon, levels = ordered.codons )
save(M5, file="results/lister/master_table.Rda")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- ANALYSIS  ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#


M5.cell_cycle <- subset( M5, gene %in% cell.cycle.regulators$gene )
M5.others     <- subset( M5, !gene %in% cell.cycle.regulators$gene )


matrix <- as.matrix( cast( data= M5.cell_cycle, formula = name ~ codon, value = "rcsu"  ) )
for.annotations <- unique(M5[,c("aa","up_in","tAI","codon")] )
annotations <- for.annotations[,-ncol(for.annotations)]; rownames(annotations) <- for.annotations$codon



require(pheatmap)

pheatmap(mat= matrix , scale="row", cluster_cols=F, cluster_rows=T, cellwidth=10, cellheight=10, annotation= annotations, clustering_method= "mcquitty" )
ggplot( M5.cell_cycle, aes(x=codon, y = name, fill=rcsu) ) + geom_tile() + facet_grid( ~ up_in, drop=T )


# Ok I let this here for now ( to reiterate later - 20 May 2014)

cell.cycle.genes <-  read.table("data/Cell cycle (high resolution)/SGD_all_cell_cycle_genes.txt", header = 1)
cell.cycle.genes$gene <- fORF(cell.cycle.genes$name)

data(RNA)
codon.usage <- subset(RNA$codon.usage, gene %in% NAMES$systematic.name)
rownames(codon.usage) <- fNAME(codon.usage$gene)
codon.usage.quantiles <- quantilize( codon.usage[,-1], bin.number= 10 ) 
codon.usage.quantiles$name <- rownames(codon.usage.quantiles)
colnames(codon.usage.quantiles)[1:4] <- paste0(colnames(codon.usage.quantiles)[1:4], ".q")

M6 <- merge( codon.usage, codon.usage.quantiles, by = 0 )
M6$cell.cycle <- ifelse( M6$gene %in% cell.cycle.genes$gene, T, F )
low.tAI.cellcycle.guys <- subset( M6, gene %in% cell.cycle.genes$gene & tAI.q %in% c("90-100") )

other.tAI.cellcycle.guys <- subset( M6, gene %in% cell.cycle.genes$gene & tAI.q %in% c("80-90") )
nrow(low.tAI.cellcycle.guys) / nrow(cell.cycle.genes)

ggplot( data = M6, aes(x=factor(cell.cycle), y=tAI, group=cell.cycle, colour=cell.cycle) ) +
  geom_boxplot(notch=T, width=0.5) + theme_bw()

guys <- ( arrange( subset(low.tAI.cellcycle.guys, select=c(name, tAI.q, codon.bias, tAI) ), tAI ) )
write.table(guys, file="results/lister/cell_cycle_regulators_with_low_tAI.txt",sep="\t", quote=F, row.names=F)
