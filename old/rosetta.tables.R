require(plyr)
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)

# correspondence between codon, anticodons and amino acids (no major difference with rosettta.stone, just that in this 
# version I added info about the 1st position of the anticodon and the 3rd position of the codon)

# Correspondance between codons and anticodons
#   rosetta.stone.stAI <- unique(TAB$genetic.code[,-3])
#   rosetta.stone.stAI$pos.1_anticodon <- sapply( as.character(rosetta.stone.stAI$anticodon), function(x) {s2c(x)[1]} )
#   rosetta.stone.stAI$pos.2.3_anticodon <- sapply( as.character(rosetta.stone.stAI$anticodon), function(x) {c2s(s2c(x)[2:3])} )
#   rosetta.stone.stAI$pos.3_codon <- sapply( as.character(rosetta.stone.stAI$codon), function(x) {s2c(x)[3]} )
#   write.table(rosetta.stone.stAI,"data/info/rosetta.stone2.txt",sep="\t",row.name=F,quote=F)


rosetta.stone2 <- read.table("data/info/rosetta.stone2.txt",header=1,sep="\t")


# bring data on codon optimality. but because things will be anticodon-centric, will need to distinguish 3 types of anticodon:
# those that only recognise optimal synonymous codons
# those that recognise both optimal and non-optimal synonymous codons
# those that only recognise non-optimal codons

data(TAB)
# CAI score
CAI <- TAB$CAI[,c("codon","Scer")]
CAI$codon <- tolower(CAI$codon)
CAI <- subset(CAI, !codon %in% c("tga","taa","tag") ) # ignore STOP codons
colnames(CAI)[2] <- "CAI"
# tAI score
tAI <- TAB$tAI[,c("codon","Scer")]
colnames(tAI)[2] <- "tAI"
# nTE
nTE <- TAB$nTE[,c("codon","Scer")]
nTE$codon <- tolower(nTE$codon)
colnames(nTE)[2] <- "nTE"
# nTE codon optimality
codon.optimality.nTE <- TAB$nTE.codon.optimality[,c("codon","Scer")]
codon.optimality.nTE$codon <- tolower(codon.optimality.nTE$codon)
colnames(codon.optimality.nTE)[2] <- "nTE.O"
# CAI codon optimality
codon.optimality.CAI <- TAB$CAI.codon.optimality[,c("codon","optimality")]
codon.optimality.CAI$codon <- tolower(codon.optimality.CAI$codon)
codon.optimality.CAI[which(codon.optimality.CAI$codon == "ggc"),"optimality"] <- "O"
codon.optimality.CAI[which(codon.optimality.CAI$codon == "tgt"),"optimality"] <- "N"
# manual corrections (should modify the source table and recompile StingRay)
codon.optimality.CAI$optimality[which(codon.optimality.CAI$codon=="atg")] <- "O"
codon.optimality.CAI$optimality[which(codon.optimality.CAI$codon=="ttg")] <- "O"
colnames(codon.optimality.CAI)[2] <- "CAI.O"
 

# "topology" of the genetic code: I defined a finite set of architectures of aa-tRNA-codon graphs # Follow the Newick tree format

genetic.code.topology <- melt( list(
  "aa(tRNA(c))"               = c("Met","Trp"), # single tRNA, single codon
  "aa(tRNA(c,c))"             = c("Cys","His","Phe","Asp","Tyr","Asn"), # single tRNA, multiple codons
  "aa(tRNA(c),tRNA(c)+)"      = c("Gln","Glu","Ile","Lys"), # multiple tRNAs, each of which has a unique codons
  "aa(tRNA(c)+,tRNA(c,c+)+)"  = c("Thr","Ser","Val","Leu","Gly","Arg"), # multiple tRNAs, some of which have multiple codons
  "aa(tRNA(c,c),tRNA(c,c))"   = c("Pro","Ala") # multiple tRNAs, each of which has multiple codons
  ) )
genetic.code.topology$L1 <- factor(genetic.code.topology$L1, levels=c("aa(tRNA(c))", "aa(tRNA(c,c))", "aa(tRNA(c),tRNA(c)+)", "aa(tRNA(c)+,tRNA(c,c+)+)", "aa(tRNA(c,c),tRNA(c,c))") )
colnames(genetic.code.topology) <- c("aa","topology")
write.table(genetic.code.topology, file = "data/info/genetic.code.topology.txt",sep="\t",quote=F, row.names=F)

#-#-#-#-#-#-#-#-#-#-#-#-#-#- CODON-CENTRIC TABLE #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
rosetta.codons <- Reduce( merge, list( CAI, tAI, nTE, codon.optimality.CAI, codon.optimality.nTE, rosetta.stone2) )
rosetta.codons <- merge(rosetta.codons, genetic.code.topology, by="aa")

# integrate information about mild/severe stress analysis performed on Gasch et al.2000 (Art. Wuster and myself)
upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)
rosetta.codons <- merge( rosetta.codons, upstress_codons, by = "codon", all.x=T)

load("data/Lister/orf_all.codon.usage_freq.Rda")
freq <- colMeans(matrixify(as.data.frame(uco.freq)))
rosetta.codons <- merge( rosetta.codons, data.frame(codon = names(freq), genomic.freq = freq), by = "codon")


write.table(rosetta.codons, file="data/info/rosetta.codons.txt",sep="\t",quote=F,row.names=F)


#-#-#-#-#-#-#-#-#-#-#-#- ANTICODON-CENTRIC TABLE #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#



# collapse the table to have unique anticodons, with all synonymous codons grouped together
rosetta.anticodons <- unique(plyr::arrange(ddply( rosetta.codons[,c("aa","anticodon","codon","CAI.O","nTE.O","pos.2.3_anticodon","pos.1_anticodon","pos.3_codon")], 
                        .(anticodon), function(x) { x$codon <- paste(x$codon,collapse="+"); 
                                                    x$pos.3_codon  <- paste(x$pos.3_codon,collapse="+")
                                                    x$CAI.O <- paste(sort(unique(x$CAI.O)),collapse=".")
                                                    x$nTE.O <- paste(sort(unique(x$nTE.O)),collapse=".")
                                                    return(x) } 
                                  ),
                          aa)
                      )



# compute tRNA Gene Copy Number (tGCN)
# tGCN <- data.frame( table(unique(TAB$anticodon2codon[,1:2])$anticodon) )
# ensures the tRNA gene copy number is identical to the one used in the simulation 
tGCN <- read.table("data/SMoPT/stochastic.simulation/anticodons.txt",header=1)[,-3]
colnames(tGCN) <- c("anticodon","tGCN")

rosetta.anticodons <- subset(merge(rosetta.anticodons, tGCN, by="anticodon"), anticodon != "UCA") # mitochondrial tRNA to ignore here.

# relative availability: how many genes a particular tRNA has in relation to all tRNAs coding for the same amino acid?
rosetta.anticodons <- ddply( rosetta.anticodons, .(aa), mutate, relative.availability.normal = round(tGCN / sum(tGCN),2)  )

# how many tRNA types (or anticodons) code for each amino acid? 
n.anticodon.per.aa <- data.frame( table(unique(rosetta.anticodons[,c("anticodon","aa")])$aa) )
colnames(n.anticodon.per.aa) <- c("aa","n.anticodons")
n.anticodon.per.aa <- arrange(n.anticodon.per.aa, n.anticodons)
aa.with.one.anticodon <- as.character(subset(n.anticodon.per.aa, n.anticodons == 1)$aa)

# is the anticodon unique to its amino acid, or are there synonymous anticodons for that amino acid?
rosetta.anticodons$aa.with.1.anticodon  <- ifelse( rosetta.anticodons$aa %in% aa.with.one.anticodon, T, F )

# integrate information about genetic code topology (see above)
rosetta.anticodons <- merge(rosetta.anticodons, genetic.code.topology, by="aa")
write.table(rosetta.anticodons, file="data/info/rosetta.anticodons.txt",sep="\t",quote=F,row.names=F)

# remove all non-essential data
rm(n.anticodon.per.aa, nTE, tGCN, aa.with.one.anticodon, rosetta.stone2, codon.optimality.CAI, codon.optimality.nTE, CAI, tAI, genetic.code.topology)
