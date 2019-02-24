# ------------------------------------------------------------------------------------ #
# STRESS-ADJUSTED tRNA ADAPTATION INDEX ANALYSIS
# ------------------------------------------------------------------------------------ #
    # The goal of this script is to compute stress-adjusted codon relative adaptiveness (w)i and 
    # stress-adjusted tRNA Adaptation Indexes for yeast gene in the different stress conditions used in 
    # our study.

    # References
  
    # [1] Evolutionary conservation of codon optimality reveals hidden signatures of cotranslational folding
    # Pechmann 2013
    
    # [2] Solving the riddle of codon usage preferences: a test for translational selection
    # Reis2004

# The selective constraint on codon-anticodon interactions sij is 0 for cognate tRNAs and small for wobble interactions

# requirements
  require(data.table)
  require(Biostrings)
  require(seqinr)
  require(StingRay)
  require(reshape2)

# In-house commands
  setwd("~/Documents/MRC16")
  source("scripts/commands.R")
  source("scripts/SMoT.R")



#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- LOAD DATA ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
master.table.reduced <- read.table("results/SMoPT/master.table.reduced.txt",header=1)


# ---- Codon-based data ----
data(TAB)
TAB$nTE # normalised TE
TAB$CAI # Codon Adaptation Index
TAB$tAI # tRNA Adaptation Index

# Correspondance between codons and anticodons
#   rosetta.stone.stAI <- unique(TAB$genetic.code[,-3])
#   rosetta.stone.stAI$pos.1_anticodon <- sapply( as.character(rosetta.stone.stAI$anticodon), function(x) {s2c(x)[1]} )
#   rosetta.stone.stAI$pos.2.3_anticodon <- sapply( as.character(rosetta.stone.stAI$anticodon), function(x) {c2s(s2c(x)[2:3])} )
#   rosetta.stone.stAI$pos.3_codon <- sapply( as.character(rosetta.stone.stAI$codon), function(x) {s2c(x)[3]} )
#   write.table(rosetta.stone.stAI,"data/info/rosetta.stone2.txt",sep="\t",row.name=F,quote=F)


 ####
# 1 # Wobbling pairs
####

      
      
      # wobble.scores.table <- data.frame( 
      #     pairing = c("I:U","G:C","U:A","C:G","G:U","I:C","I:A","U:G"), 
      #     pos.1_anticodon = c("A", "G", "U", "C", "G","A","A","U"),
      #     pos.3_codon = c("t","c","a","g","t","c","a","g"),
      #     s = 1 - c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68) 
      #     )
      
      # table to compute the relative adaptiveness
      # w <- merge( wobble.scores.table, rosetta.stone, by=c("pos.1_anticodon","pos.3_codon"), all.y=T )
      # write.table(w, file="~/Documents/StingRay/Source/data_sources/tables/wobble_scores.txt",quote=F, sep="\t", row.names=F)


 ####
# 2 # Gene Copy Numbers adjusted to stress
####

#     tRNA_molecule.copy.number <- compile.tRNAabundance.master.table()
#     write.table(tRNA_molecule.copy.number, file="data/tRNA abundance/tRNA_molecule.copy.number.txt", sep="\t", quote=F, row.names=F)
#   
  # check impact of using adjusted.tGCN in the calculation of total.tRNAab.marc
          ggplot(tRNA_molecule.copy.number, aes(x=total.tRNAab, y=total.tRNAab.marc)) +
            geom_point(size=3.5, alpha=0.5, aes(color=CAI.O) ) + geom_abline( slope=1, intercept=0 ) + 
            geom_text(data= subset(tRNA_molecule.copy.number, total.tRNAab > 1e6), aes(label=anticodon), size=4) + 
            scale_y_log10(limits = c(100, 4.5e6)) + scale_x_log10(limits = c(100, 5e6)) +
            scale_color_manual( values= c("skyblue3","gray40","orange") ) +
            facet_wrap( ~ stress) + theme(legend.position="bottom")
          
          ggplot( subset(tRNA_molecule.copy.number, free.tRNAab > total.tRNAab.marc), aes( x=total.tRNAab.marc, y = free.tRNAab) ) + geom_text(aes(color=CAI.O, label=anticodon)) + geom_abline(intercept=0, slope=1) + facet_wrap( experiment ~ time)
    # not convinced this is any useful now (don't I have "data/tRNA abundance/tRNAab_quantified.foldchange.txt" already?)
    #tRNAab.foldchange.simulation <- reshape2::dcast( data = tRNA_molecule.copy.number, formula = anticodon ~ experiment + time, sum, value.var = 'foldchange' )
    #write.table(tRNAab.foldchange.simulation, file="data/tRNAab.fold.change.txt", row.names=F, quote=F, sep="\t")

  
#   # wide format table showing tGCN values for the stress conditions
#   adjusted.tGCN <- ddply( subset(tRNA_molecule.copy.number,time!=0), .(time), function(x){ reshape2::dcast( data =  subset(x, select=c(anticodon, adjusted.tGCN, experiment)), 
#                                                                                                   formula = anticodon ~ experiment, value.var = "adjusted.tGCN")  })
#   adjusted.tGCN <- dcast( data =  subset(tRNA_molecule.copy.number,time!=0, select=c(anticodon, adjusted.tGCN, experiment, time)), formula = anticodon ~ experiment + time, value.var = "adjusted.tGCN")
#     # add normal conditions too
#   adjusted.tGCN$normal_0 <- arrange(unique(tRNA_molecule.copy.number[,c("anticodon","tGCN")]),anticodon)$tGCN
#   #adjusted.tGCN <- cast( data=adjusted.tGCN, formula = time ~ )
# 
# # w_GCN <- merge(w, adjusted.tGCN, by="anticodon")
# # 
# # w_GCN.long <- melt( w_GCN, id.vars = c("anticodon","pos.1_anticodon","pos.3_codon","pairing","s","codon","aa","time") )
# # colnames(w_GCN.long)[(ncol(w_GCN.long)-1):ncol(w_GCN.long)] <- c("stress","tGCN")
# 
#   # ---- Codons ordered as in dosReis2004 (now automatically accounted for in compute.codon.frequency.genes()  )
#   codon.table <- read.table("~/Documents/StingRay/Source/data_sources/tables/wobble_crick.txt",header=1)
#   
#   
#   # same but per-codon and ordered correctly for get.tai()
#   all.tGCN <- arrange(merge(adjusted.tGCN, codon.table, by=c("anticodon"), all.x=T,all.y=T), number)
#   all.tGCN[ is.na(all.tGCN) ]  <- 0 # keeps non-existing codons to conform to the format of get.ws()
#   all.tGCN.long <- melt( all.tGCN, id.vars = c("anticodon","number","codon"))
#   colnames(all.tGCN.long)[(ncol(all.tGCN.long)-1):ncol(all.tGCN.long)] <- c("stress","tGCN") # careful here, tGCN is computed from copy.number (which is the one used in SMoPT)
#   # add values on codon relative adaptiveness
#   all.tGCN.long <- arrange( ddply(all.tGCN.long, .(stress), function(x){ x$w <- get.ws(x$tGCN); return(x) } ), number, plyr::desc(tGCN))
#   all.tGCN.long <- all.tGCN.long[,c("codon","number","anticodon","stress","tGCN","w")]
# 
#   # change in w compaired to normal conditions
#   all.tGCN.long <- ddply( all.tGCN.long, .(stress), function(x) { x$delta_w = round(x$w - subset(all.tGCN.long, stress=="normal_0")$w, 3); return(x)  } )
# 
#   all.tGCN.long <- subset(all.tGCN.long, !(tGCN == 0 | is.na(w)) )
#   
  
 ####
# 3 # Compute codon Frequency
####

# ---- Sequence data ----
data(SEQ)
SEQ$ORF_CDS


# remove 1st codon (Starts from 4th nt position) and stop codon as well
# codon.freq <-  ldply( setNames(unique(mRNA_molecule.copy.number$ORF),unique(mRNA_molecule.copy.number$ORF)), function(x){ get.uco(SEQ$ORF_CDS[[x]][1: (length(SEQ$ORF_CDS[[x]]) - 3)], method="freq") } )
# colnames(codon.freq)[1] <- "ORF"
# codon.freq.matrix <- as.matrix(codon.freq[,-1])
# rownames(codon.freq.matrix) <- codon.freq$ORF
# codon.freq.matrix <- codon.freq.matrix[,as.character(codon.table$codon)]
# write.table(codon.freq.matrix,file="results/stress tAI/codon.freq.matrix.mat",sep="\t")


 #### -------------------------------------- #
# 4 # Compute stress-adjusted tAI scores    #
#### ------------------------------------- #
  
  
#   
# # Check the codon frequency of 
# codon.freq.matrix <- read.table("results/stress tAI/codon.freq.matrix.mat")
# # Compute s-tAI for whole CDS (minus stop codon)
# stress.tAI.genes <- unique(
#   # exclude stop codons and methionine from calculation (As in dosReis2004)
#   arrange( ddply(subset(all.tGCN.long, !codon %in% c("taa","tga","tag") ),  .(stress), 
#             function(d){ 
#               data.frame( ORF = rownames(codon.freq.matrix), 
#                           adj.tAI=get.tai(x=codon.freq.matrix[,-c(11,12,15)], w = d$w )
#                           ) 
#             } 
#            ), 
#            ORF, stress)
#   )
# stress.tAI.genes <- stress.tAI.genes[,c("ORF","stress","adj.tAI")]
# 
#   
# # mRNA abundance quantified from Gasch2000 data for experimental conditions matching those Marc used.  
# mRNA_molecule.copy.number <- read.table("data/mRNA_molecule.copy.number.txt", sep="\t", header=1)
#   
# # Check the codon frequency for the first 20 codons of all genes (used to compute truncacted s-tAI scores)
# codon.freq.matrix.20codons <- compute.codon.frequency.genes(gene.names = as.character(unique(mRNA_molecule.copy.number$ORF)), l = 20)
#   
# # s-tAI only for the first 10 codons  
# stress.tAI.first20codons <- unique(
#                               arrange( ddply(subset(all.tGCN.long, !codon %in% c("taa","tga","tag") ),  .(stress), 
#                                              function(d){ 
#                                                data.frame( ORF = rownames(codon.freq.matrix.20codons), 
#                                                            adj.tAI.20codons=get.tai(x=codon.freq.matrix.20codons[,-c(11,12,15)], w = d$w )
#                                                ) 
#                                              } 
#                               ), 
#                               ORF, stress)
#                             )
# stress.tAI.first20codons <- stress.tAI.first20codons[,c("ORF","stress","adj.tAI.20codons")]
#   
#   
#     
# #stress.tAI.genes$stress <- paste("tAI", stress.tAI.genes$stress, sep=".")
# 
# stress.tAI.data <- merge( stress.tAI.genes, stress.tAI.first20codons, by=c("ORF","stress"))
# 
# stress.tAI.data <- merge( stress.tAI.data, 
#                       subset(master.table.reduced, !is.na(faster.translation) # & as.numeric(time) <120,
#                              , select=c(ORF,stress, experiment,faster.translation,global.protein_synthesis.rate,delta_initiation)),
#                       by=c("ORF","stress")
#                     )
# stress.tAI.data$experiment <- factor(stress.tAI.data$experiment, levels=c("diauxic","ox","osm","temp","normal"))
# write.table(stress.tAI.data, file="results/stress tAI/stress.tAI_per.gene.long.txt", sep="\t", row.names=F)
# 
# stress.tAI.genes.wide <- reshape2::dcast( data = stress.tAI.genes, formula = ORF ~ stress, value.var = 'adj.tAI' )
# colnames(stress.tAI.genes.wide)[-1] <- paste("tAI",colnames(stress.tAI.genes.wide)[-1],sep=".")
# write.table(stress.tAI.genes.wide, file="results/stress tAI/stress.tAI_per.gene.txt",quote=F, sep="\t", row.names=F)
#   
  
  # hist( log2(master.table.reduced$ratio_total.time) )
  
  
head(stress.tAI.data)
  
  
    # Check how adj.tAI varies in genes translated faster or not
    g <- ggplot(stress.tAI.data, aes(x= adj.tAI.20codons, y = adj.tAI, color= faster.translation)) + 
      geom_point() + 
      geom_abline(intercept=0, slope=1, lty=2) +
      labs(title= "early response to stress (20min)", y = "\nadjusted tAI", fill="faster translation") +
      scale_color_manual(values=c("gray60","red")) +
      facet_grid( faster.translation ~ stress ) + 
      them_gul + theme(legend.position="bottom")
    
    ggsave(plot=g, filename = "results/stress tAI/stress.tAI--x.stAI20codons_y.stAI_z.stress_t.faster_translation.png",dpi=300)
    
  
    # Check how adj.tAI varies in genes translated faster or not
    g <- ggplot(stress.tAI.data, aes(x= stress, y = adj.tAI, fill= faster.translation)) + 
      geom_boxplot(notch=T) + 
      coord_flip() + labs(title= "early response to stress (20min)", y = "\nadjusted tAI", fill="faster translation") +
      scale_fill_manual(values=c("gray60","red")) +
      them_gul + theme(legend.position="bottom")
    
    ggsave(plot=g, filename = "results/stress tAI/stress.tAI--x.stAI_y.stress.z_faster_translation.png",dpi=300)

  
  # Check how adj.tAI varies in genes translated faster or not
  g <- ggplot(stress.tAI.data, aes(x= stress, y = adj.tAI.20codons, fill= faster.translation)) + 
    geom_boxplot(notch=T) + 
    coord_flip() + labs(title= "early response to stress (20min)", y = "\nadjusted tAI", fill="faster translation") +
    scale_fill_manual(values=c("gray60","red")) +
    them_gul + theme(legend.position="bottom")
  
  ggsave(plot=g, filename = "results/stress tAI/stress.tAI_20codons_early.response.pdf",dpi=400)
  
  
  ggplot(data=stress.tAI.data, aes(x=adj.tAI, y = log2( global.protein_synthesis.rate ) )) + 
    geom_point() +
    stat_smooth(method="lm") +
    facet_wrap( ~ stress ) + theme_gul
  
  ggplot(data=stress.tAI.data, aes(x= delta_initiation , y = log2( global.protein_synthesis.rate ) )) + 
    geom_point() +
    stat_smooth(method="lm") +
    facet_wrap( ~ stress ) + theme_gul
  
  
  ggplot(data=stress.tAI.data, aes(x= faster.translation, y = delta_initiation  )) + 
    geom_point() +
    geom_hline( yintercept= 0 ) +
    facet_wrap( ~ stress ) + them_gul
  
  
####################################################################################################################
# quick verif 1: do I reproduce the original tAI values I collected from previous publications?

# test <- merge( stress.tAI.genes.wide, unique(subset( master.table.020, select=c(ORF, tRNA_adaptation_index))), by = "ORF" )
# cor.test(test$tAI.normal_0, test$tRNA_adaptation_index) # yes. 0.998 correlation
# 
# ggplot(test, aes(x=tRNA_adaptation_index, y = tAI.normal_0)) + geom_abline(slope=1, intercept=0) + geom_point() + stat_smooth(method="lm")
# 



# ------------------------------------ BIG master table ------------------------------------------------------ #
# build table reporting for each gene and each of the first 100 codons and each stress condition and time point:
# 1 - the specific stress-adjusted value
# 2 - the average elongation time
# 3 - the ration of elongation time compared to normal conditions
# 
# # build codon tables: check the elongation time of individual codons in different conditions and time points
# codon.table.020 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
# codon.table.120 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
# codon.table.controls <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/controls/")
# 
#   
# # a) use codon elongation kinetics
# codon.elongation.kinetics <- subset(rbind( codon.table.020, codon.table.120 ), select=c(codon, experiment, av.elongation_time, num_events, ratio_elongation, ratio_events))
# codon.elongation.kinetics$time       <- sapply( strsplit( as.character(codon.elongation.kinetics$experiment), split = "_"), function(x) x[[3]] )
# codon.elongation.kinetics$experiment <- sapply( strsplit( as.character(codon.elongation.kinetics$experiment), split = "_"), function(x) x[[2]] )
# codon.elongation.kinetics$stress <- paste(codon.elongation.kinetics$experiment, codon.elongation.kinetics$time, sep="_")
# 
#   codon.elongation.kinetics <- subset(codon.elongation.kinetics, select=c(codon,av.elongation_time,num_events,ratio_elongation,ratio_events,stress, experiment,time))
#   #codon.elongation.kinetics <- subset(codon.elongation.kinetics, select=c(codon,experiment, time, av.elongation_time,num_events,ratio_elongation,ratio_events))
#   rm(codon.table.020, codon.table.120)
# 
#   write.table(codon.elongation.kinetics, file = "results/SMoPT/codon.elongation.kinetics.txt",sep="\t",quote=F,row.names=F)
# 
#   
#   
#   
# # b) add stress-adjusted tAI values (per codon)
#   stress.codon.adaptiveness <- all.tGCN.long[,c("codon","stress","w","delta_w")]
#   stress.codon.adaptiveness <- subset(stress.codon.adaptiveness, !is.na(w))
#     write.table(stress.codon.adaptiveness, file="results/stress tAI/stress.codon.adaptiveness.txt",quote=F, sep="\t", row.names=F) # expressed_
# 
#   
# # c) merge a and b together with another master table for codons (processed in tRNA_response_to_stress.prepare_data.R)
#   codon.potential.demand <- read.table("data/tRNA abundance/codon.potential.demand.txt", header=1)
#   
#   codon.master.table <- merge( codon.elongation.kinetics, stress.codon.adaptiveness, by = c("codon","stress"))
#   codon.master.table <- merge( codon.master.table, codon.potential.demand, by = c("codon","stress","experiment","time"))
#   
#   write.table(codon.master.table, file="data/tRNA abundance/codon.master.table.txt",sep="\t", quote=F, row.names=F)
#   
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--

  

  
  
  
  
  
metrics.local.efficiency <- merge( all.tGCN.long.s, codon.elongation.kinetics, by=c("stress","codon") )

merge(TAB$tAI[,1:2], subset(all.tGCN.long.s, stress=="normal_0"), by = "codon")
####################################################################################################################
# quick verif 2
test2 <- merge(TAB$tAI[,1:2], subset(all.tGCN.long.s, stress=="normal_0"), by = "codon")
test2$test  <- round(test2$Scer,3) == round(test2$w,3)
test2


plot(stress.tAI.genes.wide[,-1])


n.positions <- 100
d.ORF.codon.position <-  do.call(rbind, llply( setNames(unique(mRNA_molecule.copy.number$ORF),unique(mRNA_molecule.copy.number$ORF)), 
                     function(x){ 
                       codons <- tolower( splitseq( s2c( as.character( SEQ$ORF_CDS[[x]][1: min(n.positions*3, length(SEQ$ORF_CDS[[x]]) - 3 )] )), word = 3 ) )
                       return( as.data.table(data.frame( ORF=x, codon = codons, position = 1:length(codons), row.names = NULL )) )
                                 } 
                   ) )


d.efficiency.per.codon <- data.table( metrics.local.efficiency, key = c("codon") )
d.ORF.codon.position   <- data.table(d.ORF.codon.position, key=c("codon"))
d1 <- data.table( d.efficiency.per.codon[d.ORF.codon.position, list(position, ORF, stress, w, av.elongation_time, num_events, ratio_elongation, ratio_events), mult="all"], key="ORF")
d.transl.kinetics      <- subset( data.table(subset(master.table.reduced, !is.na(ratio_initiation), select=c(ORF,ratio_initiation,ratio_total.time,faster.translation,up_reg,global.protein_synthesis.rate, stress)), key = c("ORF") ), !is.na(faster.translation) )

setkeyv(d1,c("ORF","stress"))
setkeyv(d.transl.kinetics,c("ORF","stress"))

d2 <- merge(d1, d.transl.kinetics, allow.cartesian=T )
d3 <- d2[, list(mean.w=mean(w),  
                median.w=median(w),  
                mean.ratio_elongation   = mean(ratio_elongation), 
                median.ratio_elongation = median(ratio_elongation),  
                mean.av.elongation.time = mean(av.elongation_time),
                median.av.elongation.time = median(av.elongation_time)
                ),  by=list(ratio_initiation<0.5,stress,position) ]


ggplot(data = d3, aes( x= position, y= median.w, color=ratio_initiation)  ) + geom_line() + facet_grid( stress ~ .)
ggplot(data = d3, aes( x= position, y= mean.ratio_elongation, color=ratio_initiation)  ) + geom_line() + facet_grid( stress ~ .)
ggplot(data = d3, aes( x= position, y= median.ratio_elongation, color=ratio_initiation)  ) + geom_line() + facet_grid( stress ~ .)



ggplot(d.efficiency.per.codon, aes( x=w, y=av.elongation_time) ) + 
  geom_point() + 
#  stat_smooth(method="lm") +
  facet_wrap( ~ stress ) + 
#  scale_x_log10() + scale_y_log10() +
  them_gul

  
  
  
  
  
  
  
  
  


############################################ -_____- ############################################
d4 <- data.table(codon.potential.demand.foldchange.long, key = c("codon","stress") )
d5 <- subset(data.table(d2, key=c("codon","stress"))[d4,], !is.na(ORF))

# missing codons? no
# d5 represents the major table I need for the analysis: each codon of each gene in each conditions (301,469 for genes that go faster, 1,729,222 codons for other genes)

test <- d5[ position < 60 , mean(w), by=list( stress, faster.translation, position ) ]

g <- ggplot(data= test, aes(x=position,y=V1,group=faster.translation)) + 
  geom_line(aes(color=faster.translation)) + 
  them_gul + scale_color_manual(values=c("gray40","red")) +
  facet_grid( stress ~ . ) + labs(y="mean relative adaptiveness (w)\n", x = "\nposition (codon)") +
  theme(legend.position="none")
ggsave(g, filename = "results/stress tAI/codon.adaptiveness_faster.translated.genes.pdf",dpi=400)

test.2 <- d5[ position < 60 , mean(w), by=list( stress, ratio_initiation < 1, position ) ]
g <- ggplot(data= test.2, aes(x=position,y=V1,group=ratio_initiation)) + 
  geom_line(aes(color=ratio_initiation)) + 
  them_gul + scale_color_manual(values=c("gray40","red")) +
  facet_grid( stress ~ . ) + labs(y="mean relative adaptiveness (w)\n", x = "\nposition (codon)") + 
  theme(legend.position="none")
ggsave(g, filename = "results/stress tAI/codon.adaptiveness_faster.translated.genes.pdf",dpi=400)


# cumulated time to read the fist 60 codons, per gene
test.3 <- d5[ position < 20 , list( cumulated.elongation_time=sum(av.elongation_time), initiation = unique(ratio_initiation)), by=list( ORF, stress ) ]
test.3[, mean(cumulated.elongation_time), by = list( initiation < 1, stress ) ]
test.3b <- test.3[, list(cumulated.elongation_time = median(cumulated.elongation_time)), by = list( initiation < 1, stress) ]

 ggplot(data= subset(test.3b, stress %in% c("ox_20","osm_20","temp_20","diauxic_20") ), aes(x=initiation,y=cumulated.elongation_time,group=stress)) + 
  geom_bar(aes(fill=initiation), stat="identity") + 
  them_gul + scale_fill_manual(values=c("gray40","red")) +
  facet_grid( . ~ stress ) + labs(y="time to elongate the first 20 codons\n", x = "\ninitiation faster")


test.3 <- d5[ position < 10 , list( cumulated.elongation_time=sum(av.elongation_time), initiation = unique(ratio_initiation)), by=list( ORF, stress ) ]
test.3[, mean(cumulated.elongation_time), by = list( initiation < 1, stress ) ]
test.3b <- test.3[, list(cumulated.elongation_time = median(cumulated.elongation_time)), by = list( initiation < 1, stress) ]

ggplot(data= subset(test.3b, stress %in% c("ox_20","osm_20","temp_20","diauxic_20") ), aes(x=initiation,y=cumulated.elongation_time,group=stress)) + 
  geom_bar(aes(fill=initiation), stat="identity") + 
  them_gul + scale_fill_manual(values=c("gray40","red")) +
  theme(legend.position="none") +
  facet_grid( . ~ stress ) + labs(y="time to elongate the first 10 codons\n", x = "\ninitiation faster")


# anticodon supply and demand
ggplot( subset(anticodon.economy, anticodon!="CAA") , aes( x = demand.mRNAab, y = total.tRNAab )) + geom_point() + facet_wrap( ~ experiment )

