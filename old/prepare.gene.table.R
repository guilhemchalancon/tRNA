#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- GENE MASTER TABLE ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
redo_everything = T

setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(plyr)
require(reshape)
source("scripts/commands.R")
source("scripts/SMoT.R")





 #### STOCHASTIC SIMULATION OF TRANSLATION ####
# A <- subset(master.table, experiment =="normal", select=c(ORF,global.protein_synthesis.rate, est.mRNA_abundance, log2.mRNA_abundance, av.initiation_time, av.elongation_time, ratio_total.time))
# B <- subset(master.table.120, experiment =="normal" , select=c(ORF,global.protein_synthesis.rate, est.mRNA_abundance, log2.mRNA_abundance, av.initiation_time, av.elongation_time,ratio_total.time))
# 
# A$ORF == B$ORF
# head(A)
# head(B)


# Prepare a master table that summarises the translation dynamics observed in the stochastic simulation of translation for diverse stress conditions

if (redo_everything == T){
    # reduced set of features to include in the table
    additions <- c("valid",
                   "faster.translation", "faster.translation.2",
                   "n.events","n.events.normal","ratio_events",
                   "delta_initiation","delta_elongation",
                   "av.initiation_time","av.initiation_time.normal","ratio_initiation",
                   "av.elongation_time","av.elongation_time.normal","ratio_elongation",
                   "elongation.speed","elongation.speed.normal",
                   "expected.translation_rate","expected.translation_rate.normal", # freq initiation x elongation speed
                   "apparent.translation.rate","apparent.translation.rate.normal", # n.events/n.mRNAs
                   "total.translation_time","total.translation_time.normal", # sum initiation and elongation 
                   "ratio_total.time", 
                   "ratio_expected.translation.rate",
                   "ratio_apparent.translation.rate",
                   #"change.translation.rate",
                   "global.protein_synthesis.rate",
                   "IniProb",
                   "tRNA_adaptation_index","mean_nTE","mRNA_abundance","mRNA_abundance.normal","log2.mRNA_abundance",
                   "residual_cost","translation_RII","AUGCAI",
                   "length_CDS","n.codons","overal_charge"
    )
#     # table for early responses
#     master.table.020 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
#     
#     # table for late responses
#     master.table.120 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
    
    # new type of data processing (for the simulations I perform myself)
    master.table <- compile.master.table.new( path = "data/SMoPT/batch/output/")
    
    # rbind the two tables
    #SMoPT.master.table <- unique(rbind(master.table.020, master.table.120)[, c("ORF","experiment","time",additions) ])#
    SMoPT.master.table <- master.table[, c("ORF","stress","experiment","time",additions) ]
    #colnames(SMoPT.master.table)[grep("tRNA_adaptation_index",colnames(SMoPT.master.table))] <- "tAI"
    # add stress labels (experiment+time)
    #SMoPT.master.table$stress <- paste(SMoPT.master.table$experiment, SMoPT.master.table$time,sep="_")
    # add quantile levels of log2-change in mRNA abundance plus an indicator of whether or not a transcript is in the top-20% up-regulated quantile
    SMoPT.master.table <- ddply(SMoPT.master.table, .(stress), function(x) {x$quant.mRNAab = quantile.it(x$log2.mRNA_abundance); x$up_reg=ifelse(x$quant.mRNAab == "0-20", T,F); return(x) } )
    
    # Add infos on TRPs
    require(StingRay); data(GEN)
    rm(T)
    GEN$tandem.repeats$TRPs <- ifelse(GEN$tandem.repeats$class == "homo_TRP", T,F) 
    TRPs <- subset(GEN$tandem.repeats, TRPs == T)$gene
    SMoPT.master.table$TRPs <- ifelse(SMoPT.master.table$ORF %in% TRPs, T, F)
    
    # format as data table
    #SMoPT.master.table <- data.table(SMoPT.master.table, key=c("ORF","experiment","time"))
    # remove gene entries for which no translation event was detected at all
  #  SMoPT.master.table <- subset(SMoPT.master.table, !is.na(av.initiation_time) & av.initiation_time !=0)
    #SMoPT.master.table$n.codon <- SMoPT.master.table$length_CDS/3
    save(SMoPT.master.table, file = "results/SMoPT/SMoPT.master.table.Rda")
    write.table(SMoPT.master.table, file = "results/SMoPT/SMoPT.master.table.txt", quote=F, sep="\t",row.names=F)
    # clean-up
    rm(TRPs, master.table, 
       #master.table.020, master.table.120, 
       additions, GEN)

} else {
  SMoPT.master.table <- read.table("results/SMoPT/SMoPT.master.table.txt", header =1 )
}




#### STRESS-ADJUSTED tAI ####
if (redo_everything == T){

  RAMP.LENGTH  <-  10
# computed in prepare.codon.table.R
stress.codon.adaptiveness  <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1)

# check correlation of w.2 scores
cor.w <- cor(matrixify(reshape2::dcast(stress.codon.adaptiveness, codon ~ experiment + time, value.var = "w" )))
fields::image.plot(cor.w, zlim=c(0,1))

cor.w2 <- cor(matrixify(reshape2::dcast(stress.codon.adaptiveness, codon ~ experiment + time, value.var = "w.2" )))
fields::image.plot(cor.w2, zlim=c(0,1))

cor.w_FC <- cor(matrixify(reshape2::dcast(stress.codon.adaptiveness, codon ~ experiment + time, value.var = "w_FC" )))
fields::image.plot(cor.w_FC, zlim=c(0,1))

# mRNA abundance quantified from Gasch2000 data for experimental conditions matching those Marc used.  
      #mRNA_molecule.copy.number <- read.table("data/mRNA_molecule.copy.number.txt", sep="\t", header=1)
#GENES = as.character(unique(mRNA_molecule.copy.number$ORF));
require(StingRay); data(NAMES); data(SEQ)
#GENES = intersect(NAMES$systematic.name, names(SEQ$ORF_CDS[ width(SEQ$ORF_CDS) > 120 ]) )
GENES =   intersect( names(SEQ$ORF_CDS[ width(SEQ$ORF_CDS) > 100 ]), subset(NAMES, feature.type == "ORF")$systematic.name   )
GENES <- GENES[grep("^Y", GENES, perl=T)]
# Check the frequency of usage of synonymous anticodon
cat("computing anticodon frequencies...\n")
rosetta.codons <- read.table("data/info/rosetta.codons.txt",sep="\t",header=1)


  
  # anticodon frequencies
  anticodon.freq <- compute.anticodon.enrichment(gene.names = GENES, start =2, codon.info = rosetta.codons)
  write.table(anticodon.freq, file="results/stress tAI/anticodon.freq.txt",sep="\t",quote=F, row.names=F)
  
  
  # First 20 codons (first 10 for initiability, next 10 for elongability)
  anticodon.freq.20 <- compute.anticodon.enrichment(gene.names = GENES, start =2, l = RAMP.LENGTH)
  write.table(anticodon.freq.20, file="results/stress tAI/anticodon.freq.txt",sep="\t",quote=F, row.names=F)
  
  # Construct summary table
  anticodon.usage <- list()
  
  anticodon.usage$preference_in.CDS <- reshape2::dcast( data = subset(anticodon.freq, n.syn.anticodon > 1, select = c(ORF, anticodon, freq.anticodon)), 
                                                        formula = ORF ~ anticodon, value.var = "freq.anticodon" )
  
  anticodon.usage$preference_ramp <- reshape2::dcast( data = subset(anticodon.freq.20, n.syn.anticodon > 1, select = c(ORF, anticodon, freq.anticodon)), 
                                                                formula = ORF ~ anticodon, value.var = "freq.anticodon" )
  
  anticodon.usage$n.in_CDS <- reshape2::dcast( data = subset(anticodon.freq, n.syn.anticodon > 0, select = c(ORF, anticodon, N, n.anticodon)), 
                                               formula = ORF + N ~ anticodon, value.var = "n.anticodon" )
  
  anticodon.usage$n.in_ramp <- reshape2::dcast( data = subset(anticodon.freq.20, n.syn.anticodon > 0, select = c(ORF, anticodon, N, n.anticodon)), 
                                                          formula = ORF + N ~ anticodon, value.var = "n.anticodon" )
  
  anticodon.usage[1:2] <- lapply( names(anticodon.usage)[1:2], function(x){
    d <- melt(anticodon.usage[[x]], id.vars = c("ORF") )
    colnames(d) <- c("ORF","anticodon",x)
    return(d)
  })
  
  anticodon.usage[3:4] <- lapply( names(anticodon.usage)[3:4], function(x){
    d <- melt(anticodon.usage[[x]], id.vars = c("ORF","N") )
    tag  <- unlist(strsplit(x,split = "_"))[2]
    colnames(d) <- c("ORF",paste0("L.", tag ),"anticodon",x)
    return(d)
  })
  
  anticodon.usage.table <- merge_recurse( anticodon.usage, by = c("ORF","anticodon") )  
  anticodon.usage.table$freq.in_CDS       <-  anticodon.usage.table$n.in_CDS / anticodon.usage.table$L.CDS
  anticodon.usage.table$freq.in_ramp       <-  anticodon.usage.table$n.in_ramp / anticodon.usage.table$L.ramp
  anticodon.usage.table$freq.outside_ramp  <-  (anticodon.usage.table$n.in_CDS - anticodon.usage.table$n.in_ramp)/(anticodon.usage.table$L.CDS - anticodon.usage.table$L.ramp ) 
  
  anticodon.frequencies.transcriptome <- ddply(anticodon.usage.table, .(anticodon), summarise, 
                                               agg.freq.in_transciptome = sum(n.in_CDS, na.rm=T)/sum(anticodon.usage.table$n.in_CDS, na.rm=T),
                                               agg.freq.in_av.ramp =  sum(n.in_ramp, na.rm=T)/sum(anticodon.usage.table$n.in_ramp, na.rm=T),
                                               agg.freq.outside_ramp = (sum(n.in_CDS, na.rm=T) - sum(n.in_ramp, na.rm=T)) / (sum(anticodon.usage.table$n.in_CDS, na.rm=T) - sum(anticodon.usage.table$n.in_ramp, na.rm=T) ),
                                               n.in_all.ramps = sum(n.in_ramp, na.rm=T),
                                               n.not_all.ramps = sum(L.ramp, na.rm=T) - sum(n.in_ramp, na.rm=T),
                                               n.not_all.CDS = sum(L.CDS, na.rm=T) - sum(n.in_CDS, na.rm=T),
                                               n.in_all.CDSs  = sum(n.in_CDS, na.rm=T)
                                               )
  anticodon.usage.table <- merge(anticodon.usage.table, anticodon.frequencies.transcriptome, by = "anticodon" )
  
  
  # Measures the odds to find the number of observed anticodons in the ramp compared to the transcriptomic background  
  anticodon.usage.table$OR.ramp  <- ( anticodon.usage.table$n.in_ramp / (anticodon.usage.table$L.ramp - anticodon.usage.table$n.in_ramp) ) / ( anticodon.usage.table$n.in_all.ramps / anticodon.usage.table$n.not_all.ramps )
  anticodon.usage.table$SE.ramp   <- sqrt( 1/anticodon.usage.table$n.in_ramp + 1/(anticodon.usage.table$L.ramp - anticodon.usage.table$n.in_ramp) + 1/anticodon.usage.table$n.in_all.ramps + 1/anticodon.usage.table$n.not_all.ramps ) 
  anticodon.usage.table$OR.ramp_CI.low <- exp( log(anticodon.usage.table$OR.ramp) - 1.96* anticodon.usage.table$SE.ramp )
  anticodon.usage.table$OR.ramp_CI.high <- exp( log(anticodon.usage.table$OR.ramp) + 1.96* anticodon.usage.table$SE.ramp )
  anticodon.usage.table$OR.ramp_valid.interval  <- sign( log ( anticodon.usage.table$OR.ramp_CI.high ) ) * sign( log ( anticodon.usage.table$OR.ramp_CI.low ) )
  
  # Measures the odds to find the number of observed anticodons in the CDS compared to the transcriptomic background
  anticodon.usage.table$OR.CDS    <- ( anticodon.usage.table$n.in_CDS / (anticodon.usage.table$L.CDS - anticodon.usage.table$n.in_CDS) ) / ( anticodon.usage.table$n.in_all.CDSs / anticodon.usage.table$n.not_all.CDS )
  anticodon.usage.table$SE.CDS    <- sqrt( 1/anticodon.usage.table$n.in_CDS + 1/(anticodon.usage.table$L.CDS - anticodon.usage.table$n.in_CDS) + 1/anticodon.usage.table$n.in_all.CDSs + 1/anticodon.usage.table$n.not_all.CDS ) 
  anticodon.usage.table$OR.CDS_CI.low  <- exp( log(anticodon.usage.table$OR.CDS) - 1.96* anticodon.usage.table$SE.CDS )
  anticodon.usage.table$OR.CDS_CI.high <- exp( log(anticodon.usage.table$OR.CDS) + 1.96* anticodon.usage.table$SE.CDS )
  anticodon.usage.table$OR.CDS_valid.interval  <- sign( log ( anticodon.usage.table$OR.CDS_CI.high ) ) * sign( log ( anticodon.usage.table$OR.CDS_CI.low ) )
  
  
  # Measures the odds to find the number of observed anticodons in the ramp  compared to the odds to find it in the CDS
  anticodon.usage.table$OR.ramp.vs.CDS <- ( anticodon.usage.table$n.in_ramp / (anticodon.usage.table$L.ramp - anticodon.usage.table$n.in_ramp) ) / ( anticodon.usage.table$n.in_CDS / (anticodon.usage.table$L.CDS - anticodon.usage.table$n.in_CDS) )
  anticodon.usage.table$SE.ramp.vs.CDS    <- sqrt( 1/anticodon.usage.table$n.in_ramp + 1/(anticodon.usage.table$L.ramp - anticodon.usage.table$n.in_ramp) + 1/anticodon.usage.table$n.in_CDS + 1/(anticodon.usage.table$L.CDS - anticodon.usage.table$n.in_CDS) ) 
  anticodon.usage.table$OR.ramp.vs.CDS_CI.low  <- exp( log(anticodon.usage.table$OR.CDS) - 1.96* anticodon.usage.table$SE.CDS )
  anticodon.usage.table$OR.ramp.vs.CDS_CI.high <- exp( log(anticodon.usage.table$OR.CDS) + 1.96* anticodon.usage.table$SE.CDS )
  anticodon.usage.table$OR.ramp.vs.CDS_valid.interval  <- sign( log ( anticodon.usage.table$OR.ramp.vs.CDS_CI.high ) ) * sign( log ( anticodon.usage.table$OR.ramp.vs.CDS_CI.low ) )
  
  anticodon.usage.table$name <- fNAME(anticodon.usage.table$ORF)
  

#   anticodon.usage.table$rate.ratio_ramp          <-   anticodon.usage.table$freq.in_ramp /  anticodon.usage.table$agg.freq.in_av.ramp
#   anticodon.usage.table$rate.ratio_ramp.overall  <-   anticodon.usage.table$freq.in_ramp /  anticodon.usage.table$agg.freq.in_transciptome
#   anticodon.usage.table$rate.ratio_in.CDS        <-   anticodon.usage.table$freq.in_CDS /  anticodon.usage.table$agg.freq.in_transciptome
  
  
  write.table( anticodon.usage.table, file = "results/master tables/anticodon.usage.table_10.txt",quote=F, row.names=F, sep="\t")
  
  # odd ratio: (p11 x p00) / (p10 x p01) (useful it turns out. way too few events to compute realistically adequate Conf. Intervals )
  #p11: p of having anticodon i in the ramp
  #p00: p of not having anticodon i outside the ramp
  #p10: p of having anticodon i outside the ramp
  #p01: p of not having anticodon i in the ramp
  
  
#   sum( ddply(anticodon.usage.table, .(anticodon), summarise, sum(n.in_CDS) )$..1 ) == 
#     sum(anticodon.usage.table$n.in_CDS)
#     sum(unique( subset(anticodon.usage.table, select=c(ORF, L.CDS)) )$L.CDS)
#   
  
  #subset(SMoPT.master.table, select=c(experiment, time, ORF, est.mRNA_abundance, length_CDS))
  
#             write.table( anticodon.usage$preference.in.CDS,  file = "results/master tables/anticodon.preference.table.txt",  sep="\t", quote=F, row.names=F )
#             write.table( anticodon.usage$preference.20first.codons,  file = "results/master tables/anticodon.preference.20.table.txt",  sep="\t", quote=F, row.names=F )
#             write.table( anticodon.usage$n.in.CDS, file = "results/master tables/anticodon.enumeration.table.txt", sep="\t", quote=F, row.names=F )
#             write.table( anticodon.usage$n.in.20first.codons, file = "results/master tables/anticodon.enumeration.20.table.txt", sep="\t", quote=F, row.names=F )

  
  
  
  # codon frequencies / gene
  cat("computing codon frequencies of each gene (ignoring first codon)...\n")
  # Check the codon frequency for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
  #codon.freq.matrix <- read.table("results/stress tAI/codon.freq.matrix.mat") # Check the codon frequency in entire CDS 
  codon.freq.matrix <- compute.codon.frequency.genes(gene.names = GENES, l = "all", start = 2) 
  write.table(codon.freq.matrix, file="results/stress tAI/codon.freq.mat",sep="\t",quote=F)
  

  cat("computing codon frequencies of each gene outside of the coodn rampe...\n")
  # Check the codon frequency for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
  #codon.freq.matrix <- read.table("results/stress tAI/codon.freq.matrix.mat") # Check the codon frequency in entire CDS 
  codon.freq.matrix.outside.ramp <- compute.codon.frequency.genes(gene.names = GENES, l = "all", start = RAMP.LENGTH+1) 
  write.table(codon.freq.matrix, file="results/stress tAI/codon.freq_outside.ramp.mat",sep="\t",quote=F)


  # Check the codon frequency for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
  cat("computing codon frequencies in the first 20 codons of each gene...\n")
  codon.freq.matrix.20codons <- compute.codon.frequency.genes(gene.names = GENES, l = RAMP.LENGTH-1, start = 2) 
  write.table(codon.freq.matrix.20codons, file="results/stress tAI/codon.freq.matrix_first10codons.mat",sep="\t",quote=F)
  



# codon counts / gene
cat("computing codon frequencies of each gene (ignoring first codon)...\n")
# Check the codon frequency for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
#codon.freq.matrix <- read.table("results/stress tAI/codon.freq.matrix.mat") # Check the codon frequency in entire CDS 
codon.eff.matrix <- compute.codon.frequency.genes(gene.names = GENES, l = "all", start = 2, method = "eff") 
write.table(codon.eff.matrix, file="results/stress tAI/codon.eff.mat",sep="\t",quote=F)


cat("computing codon frequencies of each gene outside of the coodn rampe...\n")
# Check the codon count for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
#codon.freq.matrix <- read.table("results/stress tAI/codon.freq.matrix.mat") # Check the codon frequency in entire CDS 
codon.eff.matrix.outside.ramp <- compute.codon.frequency.genes(gene.names = GENES, l = "all", start = RAMP.LENGTH+1,method = "eff") 
write.table(codon.eff.matrix.outside.ramp, file="results/stress tAI/codon.eff_outside.ramp.mat",sep="\t",quote=F)


# Check the codon count for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
cat("computing codon frequencies in the first 20 codons of each gene...\n")
codon.eff.matrix.20codons <- compute.codon.frequency.genes(gene.names = GENES, l = RAMP.LENGTH-1, start = 2, method = "eff") 
write.table(codon.eff.matrix.20codons, file="results/stress tAI/codon.eff.matrix_first10codons.mat",sep="\t",quote=F)





# grrrr
  
  # Compute s-tAI for whole CDS (minus stop codon)
  cat("computing stress-tAI for whole CDSs...\n")
  stress.tAI.genes <- unique(
    # exclude stop codons and methionine from calculation (As in dosReis2004)
    arrange( ddply(subset(stress.codon.adaptiveness, !codon %in% c("taa","tga","tag") ),  .(experiment, time), 
                   function(d){ 
                     if( unique(as.character(d$experiment)) == "normal" ){
                       data.frame( ORF = rownames(codon.freq.matrix), 
                                   stAI   =   get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w.2" ),
                                   stAIFC =   get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" ),
                                   stAI.old = get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" )
                       ) 
                     } else {
                       data.frame( ORF = rownames(codon.freq.matrix), 
                                   stAI   = get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w.2" ),
                                   stAIFC = get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w_FC"  ),
                                   stAI.old = get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" )
                       )                        
                     }
                   }), 
    ORF, experiment, time)
  )
  stress.tAI.genes <- stress.tAI.genes[,c("ORF","experiment","time","stAI","stAIFC","stAI.old")]
  
  cat("computing stress-tAI for the first 10 codons...\n")
  # s-tAI only for the first N codons  
  stress.tAI.first20codons <- unique(
    arrange( ddply(subset(stress.codon.adaptiveness, !codon %in% c("taa","tga","tag") ),  .(experiment, time), 
                   function(d){ 
                     if( unique(as.character(d$experiment)) == "normal" ){
                     data.frame( ORF = rownames(codon.freq.matrix.20codons), 
                                 stAI.20codons   = get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w.data = d, score = "w.2" ),
                                 stAIFC.20codons = get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w.data = d, score = "w" ),
                                 stAI.20codons.old = get.tai(x=codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" )
                                 ) 
                     } else {
                      data.frame( ORF = rownames(codon.freq.matrix.20codons), 
                                   stAI.20codons   = get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w.data = d, score = "w" ),
                                   stAIFC.20codons = get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w.data = d, score = "w_FC" ),
                                   stAI.20codons.old = get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w.data = d, score = "w.2" )
                                 )   
                     } 
                   } 
    ), 
    ORF, experiment, time)
  )
  stress.tAI.first20codons <- stress.tAI.first20codons[,c("ORF","experiment","time","stAI.20codons","stAIFC.20codons","stAI.20codons.old")]
  
  
  # s-tAI only for the full CDS
stress.tAI.data <- merge( stress.tAI.genes, stress.tAI.first20codons, by=c("ORF","experiment","time"))
#   stress.tAI.normal <- subset(stress.tAI.data, time == 0)[,-c(2,3)]
#   colnames(stress.tAI.normal)[-1] <- paste(colnames(stress.tAI.normal)[-1], ".normal", sep="")
# stress.tAI.data <- merge( stress.tAI.data, stress.tAI.normal , by="ORF",all.x=T)
# stress.tAI.data <- with(stress.tAI.data, data.frame(stress.tAI.data, gain.tAI.20codons = stAI.20codons/stAI.20codons.normal, 
#                                                      gain.tAIFC.20codons =  stAIFC.20codons/stAIFC.20codons.normal ))
tAI.data <- subset(stress.tAI.data, time==0, select=c(ORF,stAI))
colnames(tAI.data)[2] <- "tAI"
stress.tAI.data <- merge(stress.tAI.data, tAI.data, by="ORF")
stress.tAI.data <- stress.tAI.data[,setdiff(colnames(stress.tAI.data), c("stAI.normal","stAIFC.normal","stAI.20codons.normal","stAIFC.20codons.normal"))]
#rm(stress.tAI.normal)       

              
  write.table(stress.tAI.data, file="results/stress tAI/stress.tAI_per.gene.long.txt", sep="\t", row.names=F)
  
  
  stress.tAI.genes.wide <- reshape2::dcast( data = stress.tAI.genes, formula = ORF ~ experiment + time, value.var = 'stAI' )
  colnames(stress.tAI.genes.wide)[-1] <- paste("tAI",colnames(stress.tAI.genes.wide)[-1],sep=".")
  write.table( format(stress.tAI.genes.wide, digits=3), file="results/stress tAI/stress.tAI_per.gene.txt",quote=F, sep="\t", row.names=F)

  stress.tAIFC.genes.wide <- reshape2::dcast( data = stress.tAI.genes, formula = ORF ~ experiment + time, value.var = 'stAIFC' )
  colnames(stress.tAIFC.genes.wide)[-1] <- paste("tAIFC",colnames(stress.tAI.genes.wide)[-1],sep=".")
  write.table( format(stress.tAIFC.genes.wide, digits=3), file="results/stress tAI/stress.tAIFC_per.gene.txt",quote=F, sep="\t", row.names=F)  

  rm(stress.tAI.genes.wide, stress.tAI.genes, stress.tAI.first20codons)#, mRNA_molecule.copy.number)
  
} else {
  stress.tAI.data <- read.table("results/stress tAI/stress.tAI_per.gene.long.txt", sep="\t", header=1)
  # stress.tAI.genes.wide  <- read.table("results/stress tAI/stress.tAI_per.gene.txt", sep="\t", header=1)
}

require(StingRay)
gene.master.table <- merge( SMoPT.master.table, stress.tAI.data, by = c("ORF", "experiment", "time"), all.x=T)
gene.master.table$name  <- fNAME(gene.master.table$ORF)

data.explore.ranking <- subset(gene.master.table, select= c(name, experiment, time, mRNA_abundance, global.protein_synthesis.rate,
                                                            tAI, stAI, stAIFC, stAI.20codons, stAIFC.20codons, log2.mRNA_abundance) )

# [Ok solved]
data.explore.ranking <- arrange( ddply(data.explore.ranking, .(experiment, time), mutate, 
                                       rank.tAI  = rank(tAI),
                                       rank.stAI = rank(stAI), 
                                       rank.stAIFC = rank(stAIFC),
                                       rank.stAI.20codons = rank(stAI.20codons),
                                       rank.stAIFC.20codons = rank(stAIFC.20codons)), experiment, time, rank.stAI )

 
require(reshape2)
require(GGally)
matrices <- list(
  tAI = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "tAI"),
  rank.tAI = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.tAI"),
  # tAI with w computed based on adjusted tGCN
  stAI = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "stAI"),
  rank.stAI = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.stAI"),
  # tAIFC: w computed based on fold chance FC 
  stAIFC = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "stAIFC"),
  rank.stAIFC = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.stAIFC"),
  ### for codon ramp
  # based on adjusted tGCN
  stAI.20codons = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "stAI.20codons"),
  rank.stAI.20codons = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.stAI.20codons"),
  # based on fold change (FC)  
  stAIFC.20codons = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "stAIFC.20codons"),
  rank.stAIFC.20codons = reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.stAIFC.20codons")
  
)

# relative changes in tAI normalised by the maximum (to range up to 1)
matrices$stAI.scaled <- matrices$stAI; 
        matrices$stAI.scaled[,-1] <- sapply( matrices$stAI.scaled[,-1], function(x){x/max(x,na.rm=T)})
matrices$FS.stAI <- matrices$stAI; 
        matrices$FS.stAI[,-1] <- sapply( matrices$FS.stAI[,-1], function(x){ (x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T)) })
matrices$FS.tAI <- matrices$tAI; 
         matrices$FS.tAI[,-1] <- sapply( matrices$FS.tAI[,-1], function(x){ (x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T)) })
matrices$stAI.20codons.scaled <- matrices$stAI.20codons; 
        matrices$stAI.20codons.scaled[,-1] <- sapply( matrices$stAI.20codons.scaled[,-1], function(x){x/max(x,na.rm=T)})
matrices$stAIFC.scaled <- matrices$stAIFC; 
        matrices$stAIFC.scaled[,-1] <- sapply( matrices$stAIFC.scaled[,-1], function(x){x/max(x,na.rm=T)})
matrices$stAI.20codons.scaled <- matrices$stAI.20codons; 
        matrices$stAI.20codons.scaled[,-1] <- sapply( matrices$stAI.20codons.scaled[,-1], function(x){x/max(x,na.rm=T)})
matrices$stAIFC.20codons.scaled <- matrices$stAIFC.20codons; 
        matrices$stAIFC.20codons.scaled[,-1] <- sapply( matrices$stAIFC.20codons.scaled[,-1], function(x){x/max(x,na.rm=T)})
  

against.normal <- function( matrix, tag="" ){
  matrix[,-1] <- sapply(matrix[,-1], function(x){ x/matrix[,"normal_0",drop=F] })
  matrix <- melt(matrix, id.vars=c("name"))
  colnames(matrix)[2] <- "stress"
  colnames(matrix)[3] <- paste("gain",tag,sep=".") 
  return(matrix)  
}

melt.scores <- function(matrix, tag="."){
  matrix <- melt(matrix, id.vars=c("name"))
  colnames(matrix)[2] <- "stress"
  colnames(matrix)[3] <- tag
  return(matrix)
}

score.data <- merge_recurse( lapply(setNames(names(matrices),names(matrices)), function(x) melt.scores(matrices[[x]],tag = x) ), by = c("name","stress"))
gain.data  <- merge_recurse( lapply(setNames(names(matrices)[-c(1:2)],names(matrices)[-c(1:2)]), function(x) against.normal(matrices[[x]],tag = x) ), by=c("name","stress"))

to.add.to.gene.master <- merge(score.data, gain.data, by=c("name","stress"))



# # stress-adjusted tAI master table ------ 
# stAI.data <- melt(stAI.matrix, id.vars = c("name","normal_0"))
# additions <- cbind(stAI.matrix[,c("name","normal_0")],"normal_0",stAI.matrix[,"normal_0"])
# colnames(additions)[3:4] <- c("variable","value")
# stAI.data <- rbind(stAI.data, additions)
# stAI.data$gain.tAI <- log2( stAI.data$value / stAI.data$normal)
# colnames(stAI.data) <- c("name","tAI","experiment","stAI","gain.tAI")


gene.master.table <- merge(gene.master.table, to.add.to.gene.master[, c("name","stress", setdiff(colnames(to.add.to.gene.master), colnames(gene.master.table)))], by = c("name","stress"))
gene.master.table$tAI.scaled <- gene.master.table$tAI / max(gene.master.table$tAI, na.rm=T)

write.table(gene.master.table, file = "results/master tables/gene.master.table.txt", quote=F, sep="\t",row.names=F)

rm(SEQ)

