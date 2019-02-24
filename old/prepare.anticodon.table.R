start_from_scratch = T

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- ANTICODON MASTER TABLE ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(plyr)
require(reshape)
require(reshape2)
source("scripts/commands.R")
source("scripts/SMoT.R")

# correspondence between codon, anticodons and amino acids
rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",header=1,sep="\t")
rownames(rosetta.anticodons) <- rosetta.anticodons$anticodon

# ----- follow the time evolution (0, 20, 60, 120) ----
if (start_from_scratch == T ){
  
    # Fold changes in tRNA abundance quantified by Marc, with stress names (e.g. diauxic_20, temp_120) [wide format]
    tRNAab.foldchange <- read.table("data/tRNA abundance/tRNAab_quantified.foldchange.txt",sep="\t",header=1)
    tRNAab.foldchange$normal_0 <- 1
    
    # same thing, long format
    tRNAab.foldchange.long <- reshape::melt(tRNAab.foldchange, id.vars = "anticodon")
    colnames(tRNAab.foldchange.long)[2:3] <- c("stress","foldchange")
    tRNAab.foldchange.long$log2.foldchange <- log2(tRNAab.foldchange.long$foldchange)
    
    
    # build a table that contains the necessary information to study variation of tRNA abundance during stress for indidivual anticodons
    tRNAab_quantified.foldchange_long <- merge( tRNAab.foldchange.long, rosetta.anticodons, by = "anticodon"  )
    tRNAab_quantified.foldchange_long$time       <- sapply( strsplit( as.character(tRNAab_quantified.foldchange_long$stress), split = "_"), function(x) x[[2]] )
    tRNAab_quantified.foldchange_long$experiment <- sapply( strsplit( as.character(tRNAab_quantified.foldchange_long$stress), split = "_"), function(x) x[[1]] )
    tRNAab_quantified.foldchange_long <- ddply( tRNAab_quantified.foldchange_long, .(stress, time), mutate, total.tRNAab = 3.33*10^6 * (tGCN * foldchange / sum(tGCN) )  )
    
    
    additions <- list( 
      ox      = subset(tRNAab_quantified.foldchange_long, experiment=="normal"), 
      diauxic = subset(tRNAab_quantified.foldchange_long, experiment=="normal"),
      temp    = subset(tRNAab_quantified.foldchange_long, experiment=="normal"),
      osm     = subset(tRNAab_quantified.foldchange_long, experiment=="normal")
    )
    
    additions <- lapply( names(additions), 
                         function(x) { y <- additions[[x]]; y$stress <- paste0(x,"_0"); y$experiment <- x; return(data.frame(y, row.names=NULL) )  } ) 
    
    tRNAab_quantified.foldchange_long <- rbind( tRNAab_quantified.foldchange_long, do.call(rbind, additions) )
    
    write.table(tRNAab_quantified.foldchange_long, file="data/tRNA abundance/tRNAab_quantified.foldchange_long.txt",sep="\t", row.names=F, quote=F)
    rm(tRNAab.foldchange, additions,tRNAab.foldchange.long )
 } else {
   tRNAab_quantified.foldchange_long <- read.table(file="data/tRNA abundance/tRNAab_quantified.foldchange_long.txt", sep="\t", header=1)
}

##------- Anticodon Supply -----
# Build a table giving the average number of freely available tRNA molecule of each species
tRNA_molecule.copy.number <- compile.tRNAabundance.master.table()
write.table(tRNA_molecule.copy.number, file="data/tRNA abundance/tRNA_molecule.copy.number.txt", sep="\t", quote=F, row.names=F)

#-------- Adjusted tGCN (for stress) ----

    # wide format table showing tGCN values for the stress conditions
    adjusted.tGCN <- ddply( subset(tRNA_molecule.copy.number,time!=0), .(time), function(x){ reshape2::dcast( data =  subset(x, select=c(anticodon, adjusted.tGCN, experiment)), 
                                                                                                              formula = anticodon ~ experiment, value.var = "adjusted.tGCN")  })
    adjusted.tGCN <- reshape2::dcast( data =  subset(tRNA_molecule.copy.number,time!=0, select=c(anticodon, adjusted.tGCN, experiment, time)), formula = anticodon ~ experiment + time, value.var = "adjusted.tGCN")
    adjusted.tGCN$normal_0 <- arrange(unique(tRNA_molecule.copy.number[,c("anticodon","tGCN")]),anticodon)$tGCN # add normal conditions too
    write.table(adjusted.tGCN, "results/stress tAI/adjusted.tGCN.txt", sep="\t",quote=F)


####

rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1,sep="\t")
expressed.codons         <- read.table ("results/SMoPT/expressed.codons.txt",sep="\t", header=1)
expressed.codons.up_reg  <- read.table ("results/SMoPT/expressed.codons.up_reg.txt",sep="\t", header=1)


anticodon.demand.codons         <- subset( merge(expressed.codons,        subset(rosetta.codons, select=c(codon, anticodon,aa)), by="codon", all.x=TRUE), !codon %in% c("taa","tag","tga") )
anticodon.demand.codons.up_reg  <- subset( merge(expressed.codons.up_reg, subset(rosetta.codons, select=c(codon, anticodon,aa)), by="codon", all.x=TRUE), !codon %in% c("taa","tag","tga") )


anticodon.demand.all    <- ddply( anticodon.demand.codons,  .(anticodon, experiment, time), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), 
                                                                                                                     demand.mRNAab = sum(x$demand.mRNAab), 
                                                                                                                     demand.events = sum(x$demand.events),
                                                                                                                     demand.mRNAab.normal = sum(x$demand.mRNAab.normal,na.rm=T)#, 
                                                                                                                     #demand.events.normal = sun(x$demand.events.normal, na.rm=T) 
                                                                                                                     ) } ) 

anticodon.demand.up_reg <- ddply( anticodon.demand.codons.up_reg,  .(anticodon, experiment, time), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), 
                                                                                                                            demand.mRNAab.up = sum(x$demand.mRNAab.up), 
                                                                                                                            demand.events.up = sum(x$demand.events.up),
                                                                                                                            demand.mRNAab.up.normal = sum(x$demand.mRNAab.up.normal)#, 
  #                                                                                                                          demand.events.up.normal = sun(x$demand.events.up.normal)
                                                                                                                            ) } ) 

# all.x = T because anticodon.demand.all also contains "normal_0", which is by definition absent from the up_reg subset
anticodon.demand <- merge(anticodon.demand.all, anticodon.demand.up_reg, by=c("anticodon","experiment","time", "recognised"), all.x=T)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# merge data on anticodon demand and supply
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

anticodon.master.table <- merge( anticodon.demand, tRNA_molecule.copy.number, by=c("anticodon","experiment","time") )
anticodon.master.table <- anticodon.master.table[, !colnames(anticodon.master.table) %in% c("codon","label","copy.number") ]
anticodon.master.table$GC.anti <- gc.content( anticodon.master.table$anticodon )
anticodon.master.table$label <- paste(anticodon.master.table$aa, anticodon.master.table$anticodon, sep="-")
write.table(anticodon.master.table, file="results/master tables/anticodon.master.table.txt",sep="\t", quote=F, row.names=F)

rm(anticodon.demand, anticodon.demand.all, anticodon.demand.up_reg, anticodon.demand.codons, anticodon.demand.codons.up_reg, expressed.codons, expressed.codons.up_reg)
