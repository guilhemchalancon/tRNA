setwd("~/Documents/MRC16")
require(reshape)
require(reshape2)
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(pheatmap)

source("scripts/commands.R")
source("scripts/SMoT.R")



# # correspondence between codon, anticodons and amino acids
# rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",header=1,sep="\t")
# rownames(rosetta.anticodons) <- rosetta.anticodons$anticodon
# 
#       # Fold changes in tRNA abundance quantified by Marc, with stress names (e.g. diauxic_20, temp_120) [wide format]
#       tRNAab.foldchange <- read.table("data/tRNA abundance/tRNAab_quantified.foldchange.txt",sep="\t",header=1)
#       tRNAab.foldchange$normal_0 <- 1
#       
#       # same thing, long format
#       tRNAab.foldchange.long <- reshape::melt(tRNAab.foldchange)
#       colnames(tRNAab.foldchange.long)[2:3] <- c("stress","foldchange")
#       tRNAab.foldchange.long$log2.foldchange <- log2(tRNAab.foldchange.long$foldchange)
#       
# 
#       # build a table that contains the necessary information to study variation of tRNA abundance during stress for indidivual anticodons
#       anticodon.supply.foldchange.long <- merge( tRNAab.foldchange.long, rosetta.anticodons, by = "anticodon"  )
#       anticodon.supply.foldchange.long$time       <- sapply( strsplit( as.character(anticodon.supply.foldchange.long$stress), split = "_"), function(x) x[[2]] )
#       anticodon.supply.foldchange.long$experiment <- sapply( strsplit( as.character(anticodon.supply.foldchange.long$stress), split = "_"), function(x) x[[1]] )
#       anticodon.supply.foldchange.long <- ddply( anticodon.supply.foldchange.long, .(stress, time), mutate, total.tRNAab = 3.33*10^6 * (tGCN * foldchange / sum(tGCN) )  )
# 
# 
#       additions <- list( 
#         ox      = subset(anticodon.supply.foldchange.long, experiment=="normal"), 
#         diauxic = subset(anticodon.supply.foldchange.long, experiment=="normal"),
#         temp    = subset(anticodon.supply.foldchange.long, experiment=="normal"),
#         osm     = subset(anticodon.supply.foldchange.long, experiment=="normal")
#       )
# 
#       additions <- lapply( names(additions), 
#                       function(x) { y <- additions[[x]]; y$stress <- paste0(x,"_0"); y$experiment <- x; return(data.frame(y, row.names=NULL) )  } ) 
# 
#       anticodon.supply.foldchange.long <- rbind( anticodon.supply.foldchange.long, do.call(rbind, additions) )
# 
#         write.table(anticodon.supply.foldchange.long, file="data/tRNA abundance/tRNAab_quantified.foldchange_long.txt",sep="\t", row.names=F, quote=F)
#       

# tRNA_molecule.copy.number <- read.table("data/tRNA abundance/tRNA_molecule.copy.number.txt",sep="\t",header=1)


# anticodon.supply.foldchange.long <- merge( tRNA_molecule.copy.number, rosetta.anticodons, by = "anticodon" )
# 
#   write.table(anticodon.supply.foldchange.long, file="data/tRNA abundance/anticodon.supply.foldchange.long.txt",sep="\t", row.names=F)

# build a table that contains the necessary information to study variation of tRNA abundance during stress for indidivual codons
# codon.potential.demand <- merge( tRNA_molecule.copy.number[, c("anticodon","experiment","time","stress",
#                                                                              "tGCN","adjusted.tGCN",
#                                                                              "foldchange","log2.foldchange","free.tRNAab","total.tRNAab","total.tRNAab.marc",
#                                                                              "free.tRNAab.normal","total.tRNAab.normal","delta_free.tRNAab",
#                                                                              "ratio_free.tRNAab")  
#                                                                            ], 
#                                                  rosetta.codons, by = c("anticodon"), all.y=T  )
# 
# codon.potential.demand <- codon.potential.demand[,c("codon",setdiff( colnames(codon.potential.demand), "codon" ) )]
# write.table(codon.potential.demand, file="data/tRNA abundance/codon.potential.demand.txt",sep="\t", row.names=F, quote=F)
# 
#   
#   head(codon.potential.demand)
