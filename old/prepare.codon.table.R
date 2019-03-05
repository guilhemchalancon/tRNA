#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- CODON MASTER TABLE ----
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

tRNA_molecule.copy.number <- read.table("data/tRNA abundance/tRNA_molecule.copy.number.txt",sep="\t",header=1)
gene.master.table <- read.table("results/master tables/gene.master.table.txt",header=1, sep="\t")



#------- Anticodon demand -----
# DEMAND: Build a table giving the mRNA molecule copy number (estimated from rounded 2^(log2) fold changes from Gasch2000) for each gene in the 4 stresses conditions
#  log2.mRNA_abundance$genome_normal_conditions <- NA  

####

# Compute the number of codon of each sort in every gene of the SMoPT analysis
require(StingRay)
data(SEQ) # contains sequence data for ORFs
genes <- as.character(unique(gene.master.table$ORF)) # picks up the genes we need (ie those for which the simulation was possible)

if ( redo_everything == T ){
  require(seqinr)
      cat("computing codon.count...\n")
      codon.count <-  ldply( setNames(genes,genes), 
                             function(x){ get.uco(SEQ$ORF_CDS[[x]], method="eff") } )
      colnames(codon.count)[1] <- "ORF"
      write.table(codon.count, file = "results/codon.count.txt", sep="\t",quote=F, row.names=F)  
      
      cat("computing codon.count...\n")
      codon.count.20first <-  ldply( setNames(genes,genes), 
                             function(x){ get.uco(SEQ$ORF_CDS[[x]][1:60], method="eff") } )
      colnames(codon.count.20first)[1] <- "ORF"
      write.table(codon.count.20first, file = "results/codon.count.20first.txt", sep="\t",quote=F, row.names=F)  
      
      rm(codon.count)
} else {
      codon.count <- read.table("results/codon.count.txt",sep="\t",header=1)
}


# integrating codon occurence, mRNA abundance and inititation rates
# Compute the codon demand per gene per condition  
if ( redo_everything == T ){
  cat("computing expressed.codons...\n")
  codon.count <- read.table("results/codon.count.txt",sep="\t",header=1)
  codon.count.per.gene <- data.frame( merge(gene.master.table, codon.count, by = c("ORF") ) )
  #codon.count.per.gene <- codon.count[gene.master.table,]
  
  # (A)
  codon.count_mRNAab <- codon.count.per.gene # compute the demand based on the abundance of transcripts
  codon.count_mRNAab[,colnames(codon.count)[-1] ]  <- codon.count_mRNAab$mRNA_abundance * as.matrix(codon.count_mRNAab[,colnames(codon.count)[-1] ]  ) 
  
  # (B) # Number of events per minute
  codon.count_events <- codon.count.per.gene # compute the demand based on the number of translation events during the 500s of simulation
  codon.count_events[,colnames(codon.count)[-1]]  <- codon.count_events$n.events / (500/60) * as.matrix(codon.count_events[,colnames(codon.count)[-1]] )
  
  
  # Compute the codon demand per condition
  codon.count_mRNAab.long <- melt(codon.count_mRNAab[,c("ORF","experiment","time", colnames(codon.count)[-1])], id.vars= c("ORF","experiment","time") )
  colnames(codon.count_mRNAab.long)[(ncol(codon.count_mRNAab.long)-1):ncol(codon.count_mRNAab.long)] <- c("codon","occurrence.mRNAab")
  codon.count_mRNAab.long <- data.table(codon.count_mRNAab.long, key=c("ORF","experiment","time","codon") )
  
  codon.count_events.long <- melt(codon.count_events[,c("ORF","experiment","time", colnames(codon.count)[-1])], id.vars= c("ORF","experiment","time") )
  colnames(codon.count_events.long)[(ncol(codon.count_events.long)-1):ncol(codon.count_events.long)] <- c("codon","occurrence.events")
  codon.count_events.long <- data.table(codon.count_events.long, key=c("ORF","experiment","time","codon"))
  
  
  # count codons in all genes
  codon.count.long <- merge( codon.count_mRNAab.long, codon.count_events.long, by=c("ORF","experiment","time","codon"), allow.cartesian = T)
  #codon.count_mRNAab.long[codon.count_events.long, ]
  
  
  # Specific subsets for genes defined as up-regulated
  codon.count.up_reg.long <- subset( merge( codon.count.long, subset(gene.master.table[,c("ORF","experiment","time","up_reg")], (up_reg == T )  ), 
                                            by =c("ORF","experiment","time")), allow.cartesian = T )
  
  
  # count codon normal codon occurence (in normal conditions for all transcripts)
  codon.count.normal <- subset(codon.count.long, experiment=="normal", select=c(ORF, codon, occurrence.mRNAab, occurrence.events))
  setnames( codon.count.normal, old = c("occurrence.mRNAab","occurrence.events"), new = c("occurrence.mRNAab.normal", "occurrence.events.normal"))
  
  
  
  # count codon normal codon occurence (in normal conditions for only transcripts that go up in that particular condition)
#   codon.count.up.normal   <- subset( merge( subset(codon.count.long,experiment=="normal", select=c(ORF, codon, occurrence.mRNAab, occurrence.events )), 
#                                             gene.master.table[,c("ORF","experiment","time","up_reg")], by =c("ORF"), allow.cartesian = T ), up_reg == T | (is.na(up_reg) & time == 0)  ) 
  
  codon.count.up.normal <- do.call( rbind, dlply( subset(gene.master.table[c("ORF","experiment","time","up_reg")], up_reg == T  ), .(experiment, time), 
         function(x){ d <- data.frame( experiment = unique(x$experiment), time = unique(x$time) ,subset(codon.count.normal, ORF %in% as.character(x$ORF)) )  
                        } ) )
  rownames(codon.count.up.normal) <- NULL
  setnames( codon.count.up.normal, old = c("occurrence.mRNAab.normal","occurrence.events.normal"), new = c("occurrence.mRNAab.up.normal", "occurrence.events.up.normal"))
  
  
  
  # all observed partitions of up-regulated genes based on all possible experiment x time point combinations
  # eq. to a 8 dimensional Venn Diagramm (196 categories from B_8 ~= 4140 possible partitions [sum of k=1 to 8 of combn(8,k) ] )
  require(reshape2)
  V8 <- split.venn(vennit(dlply(subset(codon.count.up.normal,codon=="aaa"), .(experiment, time), function(x) (unique(x$ORF)) ) ))
  
  
  codon.count.long <- merge(codon.count.long, codon.count.normal, by=c("ORF","codon") )
  codon.count.up_reg.long <- merge(codon.count.up_reg.long, 
                                 codon.count.up.normal, 
                                 by=c("ORF","codon", "experiment","time"), 
                                   all.x=T )
  
#   codon.count.up.normal.sum <- ddply(codon.count.up.normal, .(experiment, time, codon), function(x){  
#     c(demand.mRNAab.up.normal= sum(x$occurrence.mRNAab.up.normal, na.rm=T), 
#     demand.events.up.normal= sum(x$occurrence.events.up.normal, na.rm=T))
#   })
  
  rm(codon.count.per.gene) 
  rm(codon.count_mRNAab.long, codon.count_events.long, codon.count)
  
  
  #### EXPRESSED CODONS ####
 
  # Sum up the number of expressed codons per condition and codon type to obtain the anticodon demand
  expressed.codons <- ddply(codon.count.long, .(codon, experiment, time), function(x){ c(demand.mRNAab = sum(x$occurrence.mRNAab, na.rm=T), 
                                                                                         demand.events= sum(x$occurrence.events, na.rm=T), 
                                                                                         demand.mRNAab.normal= sum(x$occurrence.mRNAab.normal, na.rm=T), 
                                                                                         demand.events.normal= sum(x$occurrence.events.normal, na.rm=T) ) } ) 


  expressed.codons.up_reg <- ddply(codon.count.up_reg.long, .(codon, experiment, time), function(x){ c(demand.mRNAab.up = sum(x$occurrence.mRNAab, na.rm=T), 
                                                                                                       demand.events.up = sum(x$occurrence.events, na.rm=T),  
                                                                                                       demand.mRNAab.up.normal= sum(x$occurrence.mRNAab.up.normal, na.rm=T), 
                                                                                                       demand.events.up.normal= sum(x$occurrence.events.up.normal, na.rm=T)  
                                                                                                       ) } ) 
  
  
  # here, when merging, "normal" conditions won't be in because by definitions genes cannot be "up-regulated" in normal conditions. so use all.x=T
  expressed.codons <- merge(expressed.codons, expressed.codons.up_reg, by=c("codon","experiment","time"), all.x = T)
  expressed.codons$stress <- paste(expressed.codons$experiment, expressed.codons$time, sep="_")
  
  write.table( expressed.codons, file="results/SMoPT/expressed.codons.txt", quote=F, sep="\t", row.names=FALSE)
  write.table( expressed.codons.up_reg, file="results/SMoPT/expressed.codons.up_reg.txt", quote=F, sep="\t", row.names=FALSE)
} else {
  expressed.codons <- read.table("results/SMoPT/expressed.codons.txt", sep="\t",header=1)
}

#### STRESS CODON ADAPTIVENESS ####

# compute codon relative adaptiveness
if (redo_everything == T) {
  cat("computing stress.codon.adaptiveness...\n")
      # ---- Codons ordered as in dosReis2004 (now automatically accounted for in compute.codon.frequency.genes()  )
      adjusted.tGCN  <- read.table("results/stress tAI/adjusted.tGCN.txt", header=1)
      # fold.changes   <- read.table("data/tRNA abundance/tRNAab.fold.change.txt", header=1)
      fold.changes <- read.table("data/tRNA abundance/tRNAab_quantified.foldchange.txt",sep="\t",header=1) # AUGUST 2015: previous table didn't have t=60.
      fold.changes <- subset(fold.changes, anticodon != "CAU2")
      fold.changes$normal_0 <- 1
      # Values are updated here.
      # shouldn't affect analyses though, because this updated table was used in other analyses as well. 
      codon.table    <- read.table("~/Documents/StingRay/Source/data_sources/tables/wobble_crick_get.ws.txt",header=1)
          
      if( "GAU" %in% adjusted.tGCN$anticodon == F ){ 
        # GAU is missing in the quantification, but present in the yeast genome in 1 copy. 
        # I cannot ignore it, even though (as of Sep. 2014) this tRNA is absent from the simulation (thus I wonder how atc codons elongate...)
        # the protocol we used in the calculation of stress-adjusted tGCNs is to assume that existing anticodons must be present in at least 1 copy
        # thus I will keep GAU in 1 copy and 1 copy only in all stress conditions
        adjusted.tGCN$anticodon <- as.character(adjusted.tGCN$anticodon)
         fold.changes$anticodon <- as.character(fold.changes$anticodon)
        adjusted.tGCN <- rbind( adjusted.tGCN, c("GAU", rep(1, ncol(adjusted.tGCN)-1)) )
         fold.changes <- rbind( fold.changes, c("GAU", rep(1, ncol(fold.changes)-1)) )
        #adjusted.tGCN$anticodon <- factor(adjusted.tGCN$anticodon)
      }
        
      # DEGUG #
  ## [Guilhem: dos Reis generates the W scores 4 by 4, concatenating 4 new entries to W, 
  ## each of them being the sum of 1 conventional anticodon-codon pair and one non-conventional one
  ## since only 2 are possible per codon. If the anticodon doesn't exist at all, its contribution to the sum
  ## will be cancelled by having a tGCN of 0. It is this normal to have weird inexisting anticodons in 
  ## Source/data_sources/tables/wobble_crick.txt ]
#       wobble_scores  <- read.table("~/Documents/StingRay/Source/data_sources/tables/wobble_scores.txt",header=1)
#       control <- merge( codon.table, wobble_scores[,c("codon","anticodon")], by = "codon", all=T )
#                           control$anticodon.x  <- as.character(control$anticodon.x)
#                           control$anticodon.y  <- as.character(control$anticodon.y)
#       control <- subset( arrange( control, number ), anticodon.x != anticodon.y)
  
     # if (nrow(control)>0){cat("[!!!] problem with wobble_crick.txt\n")}
  
      # same but per-codon and ordered correctly for get.tai()
      all.tGCN <- arrange(merge(adjusted.tGCN, codon.table, by=c("anticodon"), all.x=T,all.y=T), number)
      all.tGCN[ is.na(all.tGCN) ]  <- 0 # keeps non-existing codon-anticodon relationships... to conform to the format of get.ws()?
      # problem; as of now it drops a lot of codons!
      all.tGCN.long <- melt( all.tGCN, id.vars = c("codon","anticodon","number"))
      colnames(all.tGCN.long)[(ncol(all.tGCN.long)-1):ncol(all.tGCN.long)] <- c("stress","tGCN") # careful here, tGCN is computed from copy.number (which is the one used in SMoPT)

      # why??
      all.tGCN.long$codon <- factor( all.tGCN.long$codon, levels= as.character(codon.table$codon) )    
      all.tGCN.long$tGCN <-  as.numeric( as.character( all.tGCN.long$tGCN )) 



      all.foldchange <- arrange(merge(fold.changes, codon.table, by=c("anticodon"), all.x=T, all.y=T), number)
      all.foldchange[ is.na(all.foldchange) ]  <- 0 # keeps non-existing codon-anticodon relationships... to conform to the format of get.ws()?
      # problem; as of now it drops a lot of codons!
      all.foldchange.long <- melt( all.foldchange, id.vars = c("codon","anticodon","number"))
      colnames(all.foldchange.long)[(ncol(all.foldchange.long)-1):ncol(all.foldchange.long)] <- c("stress","FC") #
      
      # why??
      all.foldchange.long$codon <- factor( all.foldchange.long$codon, levels= as.character(codon.table$codon) )    
      all.foldchange.long$FC <-  as.numeric( as.character( all.foldchange.long$FC )) 


      # add values on codon relative adaptiveness
      all.tGCN.long <- arrange( ddply(all.tGCN.long, .(stress), function(x){ x$w <- get.ws(x$tGCN,
                                                                                          # s =  c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68)
                                                                                           s =  c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
                                                                                             );
                                                                             return(x) } ), number, plyr::desc(tGCN))

      # 23 March 2015: the scaling factor used to compute relative adaptiveness, 1/Wmax, is not effective to compare the relative adaptiveness of a codon in two conditions (only good for 2 codons in one condition).
      test <- get.ws(subset(all.tGCN.long, stress == "normal_0")$tGCN,s=  c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5) )
      Wmax.t_0 <- attributes(test)$Wmax

      require(plyr)
      all.tGCN.long <- arrange( ddply(all.tGCN.long, .(stress), function(x){ x$w.2 <- get.ws(x$tGCN,
                                                                                           s =  c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5))
                                                                              Wcurrent <- attributes(x$w.2)$Wmax
                                                                              x$w.2 <- x$w.2 * Wcurrent / Wmax.t_0 # rescale to compare to normal values
                                                                              x$w.2 <- ifelse( x$w.2 > 2, 2, x$w.2 )
      return(x) } ), number, plyr::desc(tGCN)) 


      # 30 Oct 2014: Madan's suggestion for computing codon relative adaptiveness based on fold changes rather than tGCN
      all.foldchange.long <- arrange( ddply(all.foldchange.long, .(stress), function(x){ x$w_FC <- get.ws(x$FC,
                                                                                           # s =  c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68)
                                                                                           s =  c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
      );
      return(x) } ), number, plyr::desc(FC))
      


      all.metrics.long <- merge( all.tGCN.long, all.foldchange.long, by = c("codon","number","anticodon","stress")
                                 )[,c("codon","number","anticodon","stress","tGCN","w","FC","w_FC","w.2")]

# test <- arrange( ddply(all.tGCN.long, .(stress), function(x){ x$w  <- round(get.ws(x$tGCN, s =  c(0, 0, 0, 0, 0.41, 0.28, 0.9999, 0.68)),2)
#                                                               x$w2 <- round(get.ws(x$tGCN, s =  c(0, 0, 0, 0, 0.5, 0.28, 0.9999, 0.68) ),2)
#                                                               x$w3 <- round(get.ws(x$tGCN, s =  c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5) ),2)               
# return(x) } ), number, plyr::desc(tGCN))
# 
#   #test <- subset(all.tGCN.long, stress=="normal_0")
#   test <- arrange( merge( subset(test, stress=="normal_0"), TAB$tAI[,1:2], by="codon"), number )
#   test$Scer <- round( test$Scer, 2)
#   arrange( subset(test, Scer != w3), anticodon)
#   lapply ( seq(1, 61, by=4), function(x) { test[x:(x+2),] } )

#all.metrics.long <- subset(all.metrics.long, anticodon %in% unique(rosetta.anticodons$anticodon) )

      # change in w compaired to normal conditions
all.metrics.long <- ddply( all.metrics.long, .(stress), function(x) { x$delta_w = round(x$w - subset(all.metrics.long, stress=="normal_0")$w, 3);
                                                                      x$delta_w_FC = round(x$w_FC - subset(all.metrics.long, stress=="normal_0")$w_FC, 3);
                                                                      x$delta_w2 = round(x$w.2 - subset(all.metrics.long, stress=="normal_0")$w.2, 3);
                                                                      return(x)  } )
all.metrics.long$time <- sapply( strsplit( as.character(all.metrics.long$stress), split = "_"), function(x) x[[2]] ) 
all.metrics.long$experiment <- sapply( strsplit( as.character(all.metrics.long$stress), split = "_"), function(x) x[[1]] ) 
all.metrics.long <- subset(all.metrics.long, !is.na(w) )#[,setdiff(colnames(all.tGCN.long),"tGCN")] 
      # mistake in the initial version !(is.na(w) | w==0) ! don't remove the tGCN==0, these values
      # only refer to inexistent anticodons that are accounted for in the calculation. it's simply the tGCN of one of the 
      # 2 possible anticodon sequences that can pair with a given codon. The column tGCN is thus irrelevant here!

      
      # b) add stress-adjusted tAI values (per codon)
      stress.codon.adaptiveness <- all.metrics.long[,c("codon","experiment","time","w","w.2","delta_w","w_FC","delta_w_FC")]
      stress.codon.adaptiveness <- subset(stress.codon.adaptiveness, !is.na(w))
      # save the table
      write.table(stress.codon.adaptiveness, file="results/stress tAI/stress.codon.adaptiveness.txt",quote=F, sep="\t", row.names=F) # expressed_
      
      rm( all.tGCN, all.tGCN.long, adjusted.tGCN, codon.table)
} else {
      # codon   stress   w   delta_w
      stress.codon.adaptiveness <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1) # processed in prepare.anticodon.table.R
}



#### prepare to add elongation kinetics data (from stochastic simulation) ####
if ( redo_everything == T ){
        cat("computing codon.elongation.kinetics ...\n")
        codon.table.020 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
        codon.table.120 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
        codon.table.120 <- subset( codon.table.120, experiment !="genome_normal_0") # avoid duplicates
        # codon.table.controls <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/controls/")
        
        # a) use codon elongation kinetics
        codon.elongation.kinetics <- subset(rbind( codon.table.020, codon.table.120 ), select=c(codon, experiment, av.elongation_time, num_events, delta_elongation, ratio_elongation, ratio_events))
        codon.elongation.kinetics$time       <- sapply( strsplit( as.character(codon.elongation.kinetics$experiment), split = "_"), function(x) x[[3]] )
        codon.elongation.kinetics$experiment <- sapply( strsplit( as.character(codon.elongation.kinetics$experiment), split = "_"), function(x) x[[2]] )
        codon.elongation.kinetics$stress <- paste(codon.elongation.kinetics$experiment, codon.elongation.kinetics$time, sep="_")
        
        codon.elongation.kinetics <- subset(codon.elongation.kinetics, select=c(codon,av.elongation_time,num_events,delta_elongation,ratio_elongation,ratio_events,stress, experiment,time))
        #codon.elongation.kinetics <- subset(codon.elongation.kinetics, select=c(codon,experiment, time, av.elongation_time,num_events,ratio_elongation,ratio_events))
        rm(codon.table.020, codon.table.120)
        write.table(codon.elongation.kinetics, file = "results/SMoPT/codon.elongation.kinetics.txt",sep="\t",quote=F,row.names=F)  
  } else {
        
    codon.elongation.kinetics <- read.table("results/SMoPT/codon.elongation.kinetics.txt",sep="\t", header=1)
        
}



#### prepare to add data on demand potential (from stochastic simulation) ####
if( redo_everything == T ){
            cat("preparing codon.demand.potential...\n")
            tRNA_molecule.copy.number <- read.table("data/tRNA abundance/tRNA_molecule.copy.number.txt",sep="\t",header=1)
            rosetta.codons            <- read.table("data/info/rosetta.codons.txt",header=1)
            codon.demand.potential    <- merge( tRNA_molecule.copy.number[, c("anticodon","experiment","time","stress",
                                                                           "tGCN","adjusted.tGCN",
                                                                           "foldchange","log2.foldchange","free.tRNAab","total.tRNAab","total.tRNAab.marc",
                                                                           "free.tRNAab.normal","total.tRNAab.normal","delta_free.tRNAab",
                                                                           "ratio_free.tRNAab")  
                                                                       ], 
                                             rosetta.codons, by = c("anticodon"), all.y=T  )
            
            codon.demand.potential <- codon.demand.potential[,c("codon",setdiff( colnames(codon.demand.potential), "codon" ) )]
            write.table(codon.demand.potential, file="data/tRNA abundance/codon.demand.potential.txt",sep="\t", row.names=F, quote=F)
            rm(rosetta.codons)
      
} else {  
          codon.demand.potential <- read.table("data/tRNA abundance/codon.demand.potential.txt", header=1)
}


####  prepare to add adaptiveness data (following the algorithm by dosReis2004, processed in prepare.anticodon.table.R) ####
cat("initialising codon.master.table...\n")
codon.master.table <- merge(expressed.codons, stress.codon.adaptiveness, by=c("codon","experiment","time"), all.x=T) # formally expressed.codons_vs_adaptiveness
codon.master.table <- codon.master.table[,setdiff(colnames(codon.master.table),"stress")]
codon.master.table <- ddply( codon.master.table, .(experiment, time),mutate, 
                             rank.demand.mRNAab = rank(demand.mRNAab), 
                             rank.demand.mRNAab.up = rank(demand.mRNAab.up), 
                             rank.demand.events = rank(demand.events),
                             rank.demand.events.up = rank(demand.events.up),
                             relative.demand.mRNAab = demand.mRNAab / sum(demand.mRNAab) ,
                             relative.demand.mRNAab.up = demand.mRNAab.up / sum(demand.mRNAab.up)
)


#### Finalise Merges
codon.master.table <- merge( codon.master.table, codon.elongation.kinetics, by=c("codon","experiment","time"), all.y=T) # add elongation kinetics data (from stochastic simulation)
codon.master.table <- merge( codon.master.table, codon.demand.potential, by=c("codon","experiment","time","stress"), all.x=T) # add potential demand data (from stochastic simulation)
codon.master.table$GC.cod <- gc.content(codon.master.table$codon)
cat("finalising codon.master.table...\n")
#codon.master.table <- merge( codon.master.table, rosetta.codons, by="codon")

# write.table(expressed.codons_vs_adaptiveness, file = "results/stress tAI/stress.codon.adaptiveness_vs_expression.txt",sep="\t", quote=F, row.names=F)

write.table(codon.master.table, file="results/master tables/codon.master.table.txt",sep="\t", quote=F, row.names=F)


et <- function(table){ return( arrange(unique(table[,c("codon","experiment", "time")]), codon, experiment, time) )  } 


rm(stress.codon.adaptiveness, expressed.codons)
rm(codon.elongation.kinetics, codon.demand.potential, 
   codon.count_events, codon.count_mRNAab, codon.count.long, codon.count.up_reg.long)
