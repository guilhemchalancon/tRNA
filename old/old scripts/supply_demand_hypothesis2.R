setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(plyr)
require(reshape)
source("scripts/commands.R")
source("scripts/SMoT.R")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- LOAD DATA ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# ---- Codon-based data ----
data(TAB)
TAB$nTE # normalised TE
TAB$CAI # Codon Adaptation Index
TAB$tAI # tRNA Adaptation Index


# ---- Sequence data ----
data(SEQ)
SEQ$ORF_CDS

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# mRNA abundance data (also compiled in master.table so not sure it's still of any use)
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

    # 
    # # identify the genes in the SMoPT data
    # genes <- read.table("~/Documents/MRC16/data/SMoPT/example/input/S.cer.mRNA.abndc.ini.tsv", sep="\t", header=1,stringsAsFactors=F)
    # genes$Gene <- 0:(nrow(genes)-1)
    # 
    # # table giving the (rounded up) mRNA abundance for each mRNA (calculated by Marc using the rand_mRNA column in the output of the SMoPT simulation for normal calculations, and adjusted with the foldchange from Gasch2000 microarray data)
    # mRNA_abundance <- read.table("data/SMoPT/stochastic.simulation/mRNA_abundance.txt", header=1, stringsAsFactors=F)
    # mRNA_abundance$genome_normal_0 <- genes$rand_mRNA
    # log2.mRNA_abundance <- read.table("data/SMoPT/stochastic.simulation/mRNA_abundance_change.txt", header=1, stringsAsFactors=F)
    # log2.mRNA_abundance$genome_normal_0 <- 0
    # 
    # mRNA_molecule.copy.number <- melt( merge( genes[,c("Gene","ORF","IniProb")], mRNA_abundance, by="Gene" ), id.vars = c("Gene","ORF","IniProb") )[,-1]
    # colnames(mRNA_molecule.copy.number)[3:4] <- c("experiment","mRNAab")
    # #  mRNA_molecule.copy.number <- rbind( mRNA_molecule.copy.number, data.frame( ORF = genes$ORF, experiment = "genome_normal_conditions", IniProb = genes$IniProb, mRNA_abundance = genes$rand_mRNA) )
    # 
    # 
    # mRNA_log2.change <- melt( log2.mRNA_abundance, id.vars = "ORF" )
    # colnames(mRNA_log2.change)[2:3] <- c("experiment","log2.mRNAab")
    # 
    # mRNA_molecule.copy.number <-  merge( mRNA_molecule.copy.number, mRNA_log2.change, by=c("ORF","experiment") )
    # mRNA_molecule.copy.number$time       <- sapply( strsplit( as.character(mRNA_molecule.copy.number$experiment), split = "_"), function(x) x[[3]] )
    # mRNA_molecule.copy.number$experiment <- sapply( strsplit( as.character(mRNA_molecule.copy.number$experiment), split = "_"), function(x) x[[2]] )
    # mRNA_molecule.copy.number$stress  <- paste(mRNA_molecule.copy.number$experiment, mRNA_molecule.copy.number$time, sep="_")  
    # 
    # mRNA_molecule.copy.number <- subset(mRNA_molecule.copy.number, select=c(ORF,stress,experiment,time,IniProb,mRNAab,log2.mRNAab))
    # 
    # write.table(mRNA_molecule.copy.number, file="data/mRNA_molecule.copy.number.txt",sep="\t", row.names=F)
    # 


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- GENE MASTER TABLE ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

#   # reduced set of features to include in the table
#   additions <- c("Num_of_events","av.initiation_time","av.initiation_time.normal","ratio_initiation","residual_cost","translation_RII","AUGCAI","tRNA_adaptation_index","mean_nTE","est.mRNA_abundance","log2.mRNA_abundance",
#                  "delta_initiation","delta_elongation","ratio_elongation","ratio_total.time", "faster.translation", "global.protein_synthesis.rate",
#                  "length_CDS","overal_charge"
#   )
#   # table for early responses
#   master.table.020 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
#   
#   # table for late responses
#   master.table.120 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
#   
#   # rbind the two tables
#   master.table.reduced <- unique(rbind(master.table.020, master.table.120)[, c("ORF","experiment","time",additions) ])#
#   # add stress labels (experiment+time)
#   master.table.reduced$stress <- paste(master.table.reduced$experiment, master.table.reduced$time,sep="_")
#   # add quantile levels of log2-change in mRNA abundance plus an indicator of whether or not a transcript is in the top-20% up-regulated quantile
#   master.table.reduced <- ddply(master.table.reduced, .(stress), function(x) {x$quant.mRNAab = quantile.it(x$log2.mRNA_abundance); x$up_reg=ifelse(x$quant.mRNAab == "0-20", T,F); return(x) } )
#   
#   # Add infos on TRPs
#   require(StingRay); data(GEN)
#   GEN$tandem.repeats$TRPs <- ifelse(GEN$tandem.repeats$class == "homo_TRP", T,F) 
#   TRPs <- subset(GEN$tandem.repeats, TRPs == T)$gene
#   master.table.reduced$TRPs <- ifelse(master.table.reduced$ORF %in% TRPs, T, F)
#   
#   # format as data table
#   master.table.reduced <- data.table(master.table.reduced, key=c("ORF","experiment","time"))
#   # remove gene entries for which no translation event was detected at all
#   master.table.reduced <- subset(master.table.reduced, !is.na(av.initiation_time) & av.initiation_time !=0)
# 
# save(master.table.reduced, file = "results/SMoPT/master.table.reduced.Rda")
# write.table(master.table.reduced, file = "results/SMoPT/master.table.reduced.txt", quote=F, sep="\t",row.names=F)


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- ANTICODON DEMAND AND SUPPLY DATA PROCESSING ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# ----- Prepare data ----

# [1] Correspondance between codons and anticodons
  #   rosetta.stone <- unique(TAB$genetic.code[,-3])
  #   anticodon.types <- melt( list( ANN = unique(rosetta.stone$anticodon[ grep(pattern = "^A.*", x = as.character(rosetta.stone$anticodon), perl = T ) ]), 
  #         CNN = unique(rosetta.stone$anticodon[ grep(pattern = "^C.*", x = as.character(rosetta.stone$anticodon), perl = T ) ]),
  #         GNN = unique(rosetta.stone$anticodon[ grep(pattern = "^G.*", x = as.character(rosetta.stone$anticodon), perl = T ) ]), 
  #         UNN = unique(rosetta.stone$anticodon[ grep(pattern = "^U.*", x = as.character(rosetta.stone$anticodon), perl = T ) ])) )
  #   colnames(anticodon.types) <- c("anticodon","first.pos")
  # 
  #   rosetta.stone <- merge(  rosetta.stone, anticodon.types, by="anticodon")
  #   write.table(rosetta.stone, "data/info/rosetta.stone.txt", sep="\t", quote=F, row.names=F)


##------- Anticodon Supply -----
# [2] Build a table giving the average number of freely available tRNA molecule of each species

# tRNA_molecule.copy.number <- compile.tRNAabundance.master.table()
# write.table(tRNA_molecule.copy.number, file="data/tRNA abundance/tRNA_molecule.copy.number.txt", sep="\t", quote=F, row.names=F)
# 


# #-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# # ---- CODON MASTER TABLE ----
# #-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# 
# 
# #------- Anticodon demand -----
# # DEMAND: Build a table giving the mRNA molecule copy number (estimated from rounded 2^(log2) fold changes from Gasch2000) for each gene in the 4 stresses conditions
# #  log2.mRNA_abundance$genome_normal_conditions <- NA  
# 
# 
# 
# # Compute the number of codon of each sort in every gene of the SMoPT analysis
#   genes <- as.character(unique(master.table.reduced$ORF))
#   codon.count <-  ldply( setNames(genes,genes), 
#                          function(x){ get.uco(SEQ$ORF_CDS[[x]], method="eff") } )
#     colnames(codon.count)[1] <- "ORF"
# 
# # integrating codon occurence, mRNA abundance and inititation rates
# # Compute the codon demand per gene per condition  
# 
#   codon.count.per.gene <- data.frame( merge(master.table.reduced, codon.count, by = c("ORF") ) )
#   #codon.count.per.gene <- codon.count[master.table.reduced,]
# 
#   # (A)
#   codon.count_mRNAab <- codon.count.per.gene # compute the demand based on the abundance of transcripts
#   codon.count_mRNAab[,colnames(codon.count)[-1] ]  <- codon.count_mRNAab$est.mRNA_abundance * as.matrix(codon.count_mRNAab[,colnames(codon.count)[-1] ]  ) 
#   
#   # (B) # Number of events per minute
#   codon.count_events <- codon.count.per.gene # compute the demand based on the number of translation events during the 500s of simulation
#   codon.count_events[,colnames(codon.count)[-1]]  <- codon.count_events$Num_of_events / (500/60) * as.matrix(codon.count_events[,colnames(codon.count)[-1]] )
# 
#   rm(codon.count.per.gene)
# 
#     # Compute the codon demand per condition
#     codon.count_mRNAab.long <- melt(codon.count_mRNAab[,c("ORF","experiment","time", colnames(codon.count)[-1])], id.vars= c("ORF","experiment","time") )
#     colnames(codon.count_mRNAab.long)[(ncol(codon.count_mRNAab.long)-1):ncol(codon.count_mRNAab.long)] <- c("codon","occurrence.mRNAab")
#     codon.count_mRNAab.long <- data.table(codon.count_mRNAab.long, key=c("ORF","experiment","time","codon") )
#   
#     codon.count_events.long <- melt(codon.count_events[,c("ORF","experiment","time", colnames(codon.count)[-1])], id.vars= c("ORF","experiment","time") )
#     colnames(codon.count_events.long)[(ncol(codon.count_events.long)-1):ncol(codon.count_events.long)] <- c("codon","occurrence.events")
#     codon.count_events.long <- data.table(codon.count_events.long, key=c("ORF","experiment","time","codon"))
# 
#     codon.count.long <- merge( codon.count_mRNAab.long, codon.count_events.long, by=c("ORF","experiment","time","codon"), allow.cartesian = T)
#     #codon.count_mRNAab.long[codon.count_events.long, ]
#     # Specific subsets 
#     codon.count.up_reg.long <- subset( merge( codon.count.long, master.table.reduced[,c("ORF","experiment","time","up_reg"), with=F], by =c("ORF","experiment","time"), allow.cartesian = T ), up_reg == T  )
# 
#     rm(codon.count_mRNAab.long, codon.count_events.long)
# 
#       # Sum up the number of expressed codons per condition and codon type to obtain the anticodon demand
#       expressed.codons <- ddply(codon.count.long, .(codon, experiment, time), function(x){ c(demand.mRNAab = sum(x$occurrence.mRNAab), demand.events= sum(x$occurrence.events) ) } ) 
#       expressed.codons.up_reg <- ddply(codon.count.up_reg.long, .(codon, experiment, time), function(x){ c(demand.mRNAab.up = sum(x$occurrence.mRNAab), demand.events.up = sum(x$occurrence.events) ) } ) 
#       
#       # here, when merging, "normal" conditions won't be in because by definitions genes cannot be "up-regulated" in normal conditions. so use all.x=T
#       expressed.codons <- merge(expressed.codons, expressed.codons.up_reg, by=c("codon","experiment","time"), all.x = T)
#       
#       stress.codon.adaptiveness <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1)
#       expressed.codons$stress <- paste(expressed.codons$experiment, expressed.codons$time, sep="_")
# 
#       expressed.codons_vs_adaptiveness <- merge(expressed.codons, stress.codon.adaptiveness, by=c("codon","stress"))
#       expressed.codons_vs_adaptiveness <- expressed.codons_vs_adaptiveness[,-2]
#       expressed.codons_vs_adaptiveness <- ddply( expressed.codons_vs_adaptiveness, .(experiment, time),mutate, 
#                                                  rank.demand.mRNAab = rank(demand.mRNAab), 
#                                                  rank.demand.mRNAab.up = rank(demand.mRNAab.up), 
#                                                  rank.demand.events = rank(demand.events),
#                                                  rank.demand.events.up = rank(demand.events.up),
#                                                  relative.demand.mRNAab = demand.mRNAab / sum(demand.mRNAab) ,
#                                                  relative.demand.mRNAab.up = demand.mRNAab.up / sum(demand.mRNAab.up)
#                                                  )
# 
#       codon.elongation.kinetics <- read.table("results/SMoPT/codon.elongation.kinetics.txt",sep="\t", header=1)
#       expressed.codons_vs_adaptiveness <- merge( expressed.codons_vs_adaptiveness, codon.elongation.kinetics, by=c("codon","experiment","time"))
#       expressed.codons_vs_adaptiveness <- merge( expressed.codons_vs_adaptiveness, rosetta.codons, by="codon")
#   write.table(expressed.codons_vs_adaptiveness, file = "results/stress tAI/stress.codon.adaptiveness_vs_expression.txt",sep="\t", quote=F, row.names=F)
# 


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# the same anticodon can recognise several codons: need to sum up these values
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#   anticodon.demand.codons         <- subset( merge(expressed.codons,        subset(rosetta.codons, select=c(codon, anticodon,aa)), by="codon", all.x=TRUE), !codon %in% c("taa","tag","tga") )
#   anticodon.demand.codons.up_reg  <- subset( merge(expressed.codons.up_reg, subset(rosetta.codons, select=c(codon, anticodon,aa)), by="codon", all.x=TRUE), !codon %in% c("taa","tag","tga") )
#   
#   
#   anticodon.demand.all    <- ddply( anticodon.demand.codons,  .(anticodon, experiment, time), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), 
#                                                                                                                       demand.mRNAab = sum(x$demand.mRNAab), 
#                                                                                                                       demand.events = sum(x$demand.events) ) } ) 
# 
#   anticodon.demand.up_reg <- ddply( anticodon.demand.codons.up_reg,  .(anticodon, experiment, time), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), 
#                                                                                                                               demand.mRNAab.up = sum(x$demand.mRNAab), 
#                                                                                                                               demand.events.up = sum(x$demand.events) ) } ) 
# 
#   # all.x = T because anticodon.demand.all also contains "normal_0", which is by definition absent from the up_reg subset
#   anticodon.demand <- merge(anticodon.demand.all, anticodon.demand.up_reg, by=c("anticodon","experiment","time", "recognised"), all.x=T)
# 
# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# # merge data on anticodon demand and supply
# #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 
#   anticodon.economy <- merge( anticodon.demand, tRNA_molecule.copy.number, by=c("anticodon","experiment","time") )
#   anticodon.economy <- anticodon.economy[, !colnames(anticodon.economy) %in% c("codon","label","copy.number") ]
#   
# write.table(anticodon.economy, file="data/anticodon-codon balance/anticodon.economy.txt",sep="\t", quote=F, row.names=F)

