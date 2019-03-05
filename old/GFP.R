library(GeneRfold)
library(GeneGA)
library(seqinr)
library(Biostrings)
setwd("~/Documents/MRC16/data/")
GFP <- read.fasta(file="GFP seq/Scer.GFP.fasta")
GFP <- readDNAStringSet(file="GFP seq/Scer.GFP.fasta")

data(caitab) # CAI Table
w <- as.matrix(caitab[,"sc",drop=F])
cai(GFP[[1]],w)
translate( GFP[[1]] )

GFP.codons <- lapply(GFP, get.codons)
GFP.cai <- cai.sequence(GFP.codons[[1]])

require(ggplot2)

p <- ggplot(GFP.cai, aes(x=index, y=CAI) )
q <- p + geom_violin()
p + geom_smooth()
p + geom_hex() + opts(title = "adaptation score of GFP codons")

###########################################
inputdatfile <- system.file("sequences/input.dat", package = "seqinr")
input <- read.fasta(file = inputdatfile) # read the FASTA file

scucofile <- system.file("sequences/scuco.txt", package = "seqinr")
scuco.res <- read.table(scucofile, header = TRUE) # read codonW result file

cai.res <- sapply(input, cai, w = w)
plot(cai.res, scuco.res$CAI,
     main = "Comparison of seqinR and codonW results",
     xlab = "CAI from seqinR",
     ylab = "CAI from codonW",
     las = 1)
abline(c(0,1))