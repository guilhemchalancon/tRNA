PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
setwd("~/Documents/MRC16/")
source("scripts/commands.R"); source("scripts/SMoT.R")

require(XLConnect)
require(Biostrings)
require(seqinr)
require(plyr)
require(tidyr)

wb <- loadWorkbook("~/Documents/MRC16/data/qPCR/PCR_efficiencies.xlsx")
qPCR.efficiency <- readWorksheet(wb, sheet = getSheets(wb)[1],header=T, autofitCol=T)
  qPCR.efficiency$GC_forward  <- gc.content(qPCR.efficiency$forward_primer)
  qPCR.efficiency$GC_reverse  <- gc.content(qPCR.efficiency$reverse_primer)
  qPCR.efficiency$GC_both     <- gc.content( paste(qPCR.efficiency$forward_primer, qPCR.efficiency$reverse_primer,sep="")) # could also just take the mean
  qPCR.efficiency$length_forward <- nchar(qPCR.efficiency$forward_primer)
  qPCR.efficiency$length_reverse <- nchar(qPCR.efficiency$reverse_primer)
head(qPCR.efficiency)

primer.ID <- rbind( data.frame( ID = paste(qPCR.efficiency$anticodon, "_forward",sep=""), seq = qPCR.efficiency$forward_primer ),
                    data.frame( ID = paste(qPCR.efficiency$anticodon, "_reverse",sep=""), seq = qPCR.efficiency$reverse_primer ))
primer.ID <- matrixify(primer.ID)



alignments.to.do <- data.frame(t(combn( rownames(primer.ID), 2 )))
colnames(alignments.to.do) <- c("id1","id2")
alignments.to.do$type  <-  paste(gsub(alignments.to.do$id1, pattern="^(.*)_(.*)$", replacement="\\2"), gsub(alignments.to.do$id2, pattern="^(.*)_(.*)$", replacement="\\2"),sep="_")
nrow(alignments.to.do)
table(alignments.to.do$type)

# PCR effiency (PrimerBLAST)

ggplot(data=qPCR.efficiency, aes(x=PCR_efficiency)) + geom_dotplot(binwidth=1) + xlim(0,105) + 
  labs(x="RT2 PCR efficiency")

require(scales)
g <- ggplot(data=qPCR.efficiency, aes(x=PCR_efficiency/100, y = GC_both)) + 
  geom_rect(data=qPCR.efficiency[1,],xmax=Inf, xmin=-Inf, ymin=0.35,ymax=0.65,fill="#96969690") +
  geom_point(size=2.5) + 
  geom_vline(xintercept=0.90, lty=3) + 
  scale_x_continuous(labels=percent, limits=c(0,1.2)) +
  scale_y_continuous(labels=percent, limits=c(0,1)) +
  labs(x=expression(paste(RT^2," PCR efficiency")), y = "GC-content")
ggsave(plot=g, filename = paste0(PATH,"part1/qPCR_primers.efficiency.pdf"), width=6.74, height=2, useDingbats=F )

low.efficiency.anticodons <- qPCR.efficiency$anticodon[ qPCR.efficiency$PCR_efficiency < 90 ]

QA <- marc.data
QA <- ddply(QA, .(experiment, time), mutate, z.av = zscore(average), z.std = zscore(std) )
QA$below.efficiency  <- ifelse(QA$anticodon %in% low.efficiency.anticodons, "<90",">90")
# ggplot(data = QA, aes(x = below.efficiency, y = z.av) ) + geom_boxplot(notch=T) +
#   labs(x="", y="mean fold change (Z-score)")
g <- ggplot(data = QA, aes(x = below.efficiency, y = z.std) ) + geom_boxplot(notch=F) +
  labs(x="", y="std error (Z-score)")
ggsave(plot=g, filename = paste0(PATH,"part1/qPCR_primers.std.pdf"), width = 2.51, height=2, useDingbats=F )
wilcox.test( QA$z.std ~ QA$below.efficiency, alternative = "greater")

# run all pairwise alignments of primers and collect information on them
primer.alignments <- do.call(rbind, apply( alignments.to.do, 1, function(x){
 # c(primer.ID[x[[1]]],x[[2]])
    aln <- pairwiseAlignment( pattern = DNAString(primer.ID[x[[1]],]), 
                              subject = DNAString(primer.ID[x[[2]],]) 
                            )
    output <- data.frame(
      id1 = x[[1]],
      id2 = x[[2]],
      identity.all = pid(aln, type="PID1"),
      identity.aln.positions = pid(aln, type="PID2"),
      score = aln@score,
      aln.pattern = as.character(aln@pattern), 
      aln.subject = as.character(aln@subject)
    )
}))



primer.alignments$type <-  paste(gsub(primer.alignments$id1, pattern="^(.*)_(.*)$", replacement="\\2"), 
                                 gsub(primer.alignments$id2, pattern="^(.*)_(.*)$", replacement="\\2"),sep="_")
primer.alignments$anticodon1 <-  gsub(primer.alignments$id1, pattern="^(.*)_(.*)$", replacement="\\1") 
primer.alignments$anticodon2 <-  gsub(primer.alignments$id2, pattern="^(.*)_(.*)$", replacement="\\1")
primer.alignments[,c("anticodon1","anticodon2")] <- t(apply(primer.alignments[,c("anticodon1","anticodon2")], 1, function(x) sort(x)))
write.table(primer.alignments, "~/Documents/MRC16/data/qPCR/primer.alignements.txt", row.names=F, sep="\t", quote=F)

primer.alignments <- read.table("~/Documents/MRC16/data/qPCR/primer.alignements.txt",header=1)




# Compare how similar the primers of every possible pair of tRNAs are.
forward_forward <- subset(primer.alignments, type=="forward_forward", select=c(id1,id2, identity.all))
forward_forward$id1 <- gsub(forward_forward$id1, pattern="^(.*)_(.*)$", replacement="\\1")
forward_forward$id2 <- gsub(forward_forward$id2, pattern="^(.*)_(.*)$", replacement="\\1")
forward_forward[,c("id1","id2")] <- t(apply(forward_forward[,c("id1","id2")], 1, function(x) sort(x)))
colnames(forward_forward)[3] <- "identity.forward"

reverse_reverse <- subset(primer.alignments, type=="reverse_reverse", select=c(id1,id2, identity.all))
reverse_reverse$id1 <- gsub(reverse_reverse$id1, pattern="^(.*)_(.*)$", replacement="\\1")
reverse_reverse$id2 <- gsub(reverse_reverse$id2, pattern="^(.*)_(.*)$", replacement="\\1")
reverse_reverse[,c("id1","id2")] <- t(apply(reverse_reverse[,c("id1","id2")], 1, function(x) sort(x)))
colnames(reverse_reverse)[3] <- "identity.reverse"

aa_anticodon <- unique(subset(rosetta.anticodons, select=c(anticodon,aa)))
comparison.tRNA <- merge(forward_forward,reverse_reverse, by=c("id1","id2"))
comparison.tRNA$mean <- (comparison.tRNA$identity.forward +comparison.tRNA$identity.reverse)/2
comparison.tRNA <- arrange(comparison.tRNA, plyr::desc(mean))
comparison.tRNA$id1 <- factor(comparison.tRNA$id1, 
                              levels= as.character(aa_anticodon$anticodon),
                              labels= paste(as.character(aa_anticodon$anticodon), as.character(aa_anticodon$aa),sep="-"))
comparison.tRNA$id2 <- factor(comparison.tRNA$id2, 
                              levels= as.character(aa_anticodon$anticodon),
                              labels= paste(as.character(aa_anticodon$anticodon), as.character(aa_anticodon$aa),sep="-"))
                              
                              
                              
identity.matrix <- get.sym.matrix(df = primer.alignments, value.var = "identity.all")
identity.matrices <- dlply( primer.alignments, .(type), 
                            get.sym.matrix, value.var = "identity.all") 
identity.matrices$by.tRNA <- get.sym.matrix(df = subset(comparison.tRNA, id1!="CAU2" & id2!="CAU2" ) , value.var = "mean")


order.anticodons <- as.character(arrange(annotations, aa, anticodon)$anticodon)

identity.matrices$by.tRNA[is.na(identity.matrices$by.tRNA)] <- 0

require(pheatmap)
heatmaps <- list()
heatmaps$forward_forward <- pheatmap(identity.matrices$forward_forward, treeheight_row = 10, treeheight_col=10)
heatmaps$reverse_reverse <- pheatmap(identity.matrices$reverse_reverse, treeheight_row = 10, treeheight_col=10)
heatmaps$forward_reverse <- pheatmap(identity.matrices$forward_reverse, treeheight_row = 10, treeheight_col=10)

pdf(paste0(PATH,"part1/qPCR_primers.seq_similarity.pdf"), width = 7, height = 6)
heatmaps$by.tRNA         <- pheatmap(identity.matrices$by.tRNA[order.anticodons,order.anticodons], 
                                     color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Greys")))(100))[1:80],
                                     cluster_cols = F, cluster_rows = F,
                                     treeheight_row = 15, treeheight_col = 15, 
                                     annotation = data.frame(annotations.2[,c("aa"),drop=F]) )
dev.off()


# check correlation between tRNA abundance measurmeents in function of PCR primer similarity
correlations.all.points <- melt(cor(t(M.tRNAab$all.together)))
colnames(correlations.all.points) <- c("id1","id2","Pearson")
correlations.all.points[,c("id1","id2")] <- t(apply(correlations.all.points[,c("id1","id2")], 1, function(x) sort(x)))
correlations.all.points <- unique(correlations.all.points)
require(tidyr)
test.correlations.primer.similarity <- arrange(merge( comparison.tRNA, correlations.all.points, by=c("id1","id2")), plyr::desc(mean))
test.correlations.primer.similarity$label1 <- test.correlations.primer.similarity$id1
test.correlations.primer.similarity$label2 <- test.correlations.primer.similarity$id2
test.correlations.primer.similarity <- separate(test.correlations.primer.similarity, col = "id1", into = c("codon1","aa1"), sep="-" )
test.correlations.primer.similarity <- separate(test.correlations.primer.similarity, col = "id2", into = c("codon2","aa2"), sep="-" )
test.correlations.primer.similarity$same.aa <- ifelse(test.correlations.primer.similarity$aa1 == test.correlations.primer.similarity$aa2, T,F)
test.correlations.primer.similarity$aa.set <- apply(test.correlations.primer.similarity[,c("aa1","aa2")], 1, function(x){ paste(unique(x),collapse="+")})


primer.correlations <- ggplot(data=test.correlations.primer.similarity, aes(x=mean,y=Pearson) ) + 
  geom_point( data = subset(test.correlations.primer.similarity, same.aa==F), size=3, color="gray70") + 
  geom_point( data = subset(test.correlations.primer.similarity, same.aa==T), shape=1, size=3, pch=21) +
  geom_hline( yintercept =-0.5, lty=3) +
  geom_hline( yintercept = 0.5, lty=3) +
  geom_vline( xintercept = 90, lty=3) + coord_flip() +
  geom_text( data = subset(test.correlations.primer.similarity, mean > 90 ), aes(label= paste(codon1,codon2,sep="\n")), size =2.4 ) +
  labs(x="sequence similarity of PCR primers", y="Pearson correlation of changes in tRNA abundance")
ggsave(plot = primer.correlations, filename=paste0(PATH,"part1/qPCR_primers.correlations.pdf"), width = 6, height = 6, useDingbats=F)

distribution.seq.sim <-  ggplot(test.correlations.primer.similarity, aes(x=mean))  + geom_histogram(binwidth=3, fill="gray40") + 
  labs(x="sequence similarity of PCR primers") + theme_minimal()
ggsave(plot = distribution.seq.sim, filename=paste0(PATH,"part1/qPCR_primers.correlations_dist.pdf"), width = 3, height = 3, useDingbats=F)



subset(primer.alignments, anticodon1 == "CUG" & anticodon2=="UUG")

subset(primer.alignments, anticodon1 == "CCC" & anticodon2=="GCC")
