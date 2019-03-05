# Analysis of the stochastic simulation by Marc using Josua Plotkin's source code (Cell theory paper of whole-cell modeling of translation)
setwd("~/Documents/Research/MRC16") # change in May 2018 to reflect change of folder location in my drive

      zscore <- function(x){  
        mu <- mean(x[!is.infinite(x)],na.rm=T)
        sigma  <- sd(x[!is.infinite(x)], na.rm=T)
        return( ( x - mu ) / sigma )
      }

code.it <- function( binary.vector ){
  categories <- LETTERS[ 1:length(binary.vector) ]
  paste( categories[ which(binary.vector!=0) ], collapse = "&")
}

get.uco  <- function(biostring.seq,method="freq"){ require(seqinr); require(Biostrings); uco(s2c(Biostrings::toString(biostring.seq)),index=method) }

# 
# compile.master.table  <- function( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/"  ){
#       
#       experiments <- list.files(path = path, full.names = F )
#       files <- lapply( setNames( paste(path, experiments, sep="/"), experiments), list.files)
#       genes <- read.table("~/Documents/MRC16/data/SMoPT/example/input/S.cer.mRNA.abndc.ini.tsv", sep="\t", header=1)
#       genes$Gene <- 0:(nrow(genes)-1)
#       
#       mRNA_abundance <- read.table("~/Documents/MRC16/data/SMoPT/stochastic.simulation/mRNA_abundance.txt", sep="\t", header=1)
#       mRNA_abundance$genome_normal_0  <- genes$rand_mRNA
#       mRNA_abundance_changes <- read.table("~/Documents/MRC16/data/SMoPT/stochastic.simulation/mRNA_abundance_change.txt", sep="\t", header=1)
#       mRNA_abundance_changes$genome_normal_0 <- 0
#       
#       
#       require(StingRay)
#       # collect the unformatted data
#       data <- lapply( setdiff( names(files), "genome_normal_conditions") , 
#                       function(x){ list( initiation = read.table( paste(path, x ,"/", files[[x]][4], sep="" ), sep="\t", header=1),  #,
#                                          elongation = read.table( paste(path, x ,"/", files[[x]][5], sep="" ), sep="\t", header=1)
#                                         ) 
#                                   }
#               )
#       names(data) <- setdiff( names(files), "genome_normal_conditions")
#       data$genome_normal_0 <- list( initiation = read.table( "data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_gene_initimes.out.txt", sep="\t", header=1),  #,
#                                     elongation = read.table( "data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_gene_totetimes.out.txt", sep="\t", header=1) )
#       
#       # format the data into a master table
#       data2 <- ldply( data, function(x){ merge( x$initiation, x$elongation, by=c("Gene","Num_of_events")) } )
#       
#       
#       # reference data set for normal conditions
#       normal.conditions_initiation <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_gene_initimes.out.txt", header=1, sep="\t")
#       colnames(normal.conditions_initiation)[2:3] <- c("Num_of_events.initiation.normal", "av.initiation_time.normal")
#       
#       # reference data set for normal conditions
#       normal.conditions_elongation <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_gene_totetimes.out.txt", header=1, sep="\t")
#       colnames(normal.conditions_elongation)[2:3] <- c("Num_of_events.elongation.normal", "av.elongation_time.normal")
#       
#       data3 <- merge( data2, normal.conditions_initiation, by = "Gene", all.x=T )
#       data3 <- merge( data3, normal.conditions_elongation, by = "Gene", all.x=T )
#       
#       # convert to minutes (will subsequently rename columns)
#       time.variables <- c("Avg_initiation_time.sec.","Avg_total_elong_time.sec.","av.initiation_time.normal","av.elongation_time.normal")
#       data3[,time.variables] <- data3[,time.variables]/60
#       data3$delta_initiation <- data3$Avg_initiation_time.sec. - data3$av.initiation_time.normal
#       data3$delta_elongation <- data3$Avg_total_elong_time.sec. - data3$av.elongation_time.normal
#       data3$ratio_initiation <- data3$Avg_initiation_time.sec. / data3$av.initiation_time.normal
#       data3$ratio_elongation <- data3$Avg_total_elong_time.sec. / data3$av.elongation_time.normal
#       
#       # numeric variables
#       num.var <- names( which( sapply( data3, is.numeric) == T) ) 
#       data3[, num.var] <- do.call(data.frame, lapply(data3[, num.var], function(x) replace(x, is.infinite(x),NA))) 
#       data3[, num.var] <- do.call(data.frame, lapply(data3[, num.var], function(x) replace(x, is.nan(x),NA))) 
#       #function(x){ x[which(is.infinite(x))]  <- NA; x}
#       
#       # check that every gene (3,795) is represented in each experiment
#       # ddply(data3, .(.id), nrow)
#       
#       # master table (including gene names, length and average PARS score - though there are missing atm)
#       master.table <- data3
#       colnames(master.table)[2] <- "experiment"
#       colnames(master.table)[4] <- "av.initiation_time"
#       colnames(master.table)[5] <- "av.elongation_time"
#       master.table <- merge(master.table, genes, by="Gene")
#       master.table$Name <- fNAME(master.table$ORF)
#       master.table <- master.table[,c("ORF","Gene", "Name", setdiff( colnames(master.table), c("ORF","Gene","Name") )) ]
#  
#       
#       # --------------------------------------------------------------------------------------------------------- #
#       
#       # load the master table
# #      master.table <- read.table("data/SMoPT/stochastic.simulation/master.table.txt",header=1,sep="\t")
#         
#       # Add relevant control variables to the master table
#       load("~/Documents/MRC09/data/PLS-PM/df.PLSPM.non_transformed_non_scaled.Rda")
#       
#       ## colnames( df.PLSPM.nTnS )
#       # proto_genes # overexpression_toxicity # average_cost  # total_metabolic_cost
#       # residual_cost # overal_charge # protein_length  # poor_STOP_poor  # AAA_before_STOP
#       # weak_tetra  # NAA_before_STOP # STOP  # AUGCAI  # translation_RII # mean_nTE
#       # tRNA_adaptation_index # ribosome_density_YPD  # ribosome_density_SD
#       # protein_abundance_ypd # protein_abundance_sd  # protein_abundance   # protein_abundance_H2O2
#       
#         features.to.add <- c("proto_genes","overexpression_toxicity","average_cost","total_metabolic_cost","residual_cost","overal_charge","protein_length","length_CDS",
#                              "poor_STOP_poor","AAA_before_STOP","weak_tetra","NAA_before_STOP","STOP","AUGCAI","translation_RII","mean_nTE","tRNA_adaptation_index",
#                              "ribosome_density_YPD","ribosome_density_SD","protein_abundance_ypd","protein_abundance_sd","protein_abundance","protein_abundance_H2O2")
#         
#         Additions <- df.PLSPM.nTnS[ c("gene",features.to.add)]
#         
#         master.table.2 <- merge(master.table, Additions, by.x="Name",by.y="gene", all.x=T)
#       
#         
#         # simply adding the av. intitation and elongation time to get an estime of the expected total translation time
#         master.table.2$expected.translation_time.normal <- master.table.2$av.initiation_time.normal + master.table.2$av.elongation_time.normal
#         master.table.2$expected.translation_time        <- master.table.2$av.initiation_time + master.table.2$av.elongation_time
#         master.table.2$ratio_total.time <- master.table.2$expected.translation_time / master.table.2$expected.translation_time.normal
#         
#         # apply linear regressions of the expected translation time in function of the gene's length - and take residuals to estimate 
#         # the portion of translation time that is indepenent of length (i.e. find how fast genes are translated while accounting for their length)
#         #L <- dlply( master.table.2, .(experiment), function(x) { lm( expected.translation_time ~ length_CDS, data = x ) })
#         conditions <- unique(master.table.2$experiment)
#         L1 <- lapply( conditions, function(x) {  data = subset(master.table.2, experiment==x); rownames(data) <- data$ORF; lm( expected.translation_time ~ length_CDS, data = data )   } )
#         names(L1) <- conditions
#         resid.L1 <- ldply(L1, function(x){ data.frame( ORF = names(x$residuals), residual.translation.time = x$residuals)  } )
#         colnames(resid.L1)[1] <- "experiment"
#         
#         master.table.2 <- merge( master.table.2, resid.L1, by=c("ORF","experiment"),all.x=T)
#       
#         L2 <- lapply( conditions, function(x) {  data = subset(master.table.2, experiment==x); rownames(data) <- data$ORF; lm( delta_elongation ~ length_CDS, data = data )   } )
#         names(L2) <- conditions
#         resid.L2 <- ldply(L2, function(x){ data.frame( ORF = names(x$residuals), residual.delta_elongation = x$residuals)  } )
#         colnames(resid.L2)[1] <- "experiment"
#       
#         master.table.2 <- merge( master.table.2, resid.L2, by=c("ORF","experiment"),all.x=T)
#       
#       
#       # ---- additional processing of the data
#       master.table.2$zscore.delta_initiation          <- zscore( master.table.2$delta_initiation )
#       master.table.2$zscore.delta_elongation          <- zscore( master.table.2$delta_elongation )
#       master.table.2$zscore.residual.delta_elongation <- zscore( master.table.2$residual.delta_elongation )
#       master.table.2$zscore.residual.translation.time <- zscore( master.table.2$residual.translation.time )
# 
#       mRNA_abundance.long <- melt(mRNA_abundance, id.vars = "Gene")
#       colnames(mRNA_abundance.long)[2:3] <- c("experiment","est.mRNA_abundance")
# 
#       mRNA_abundance_changes.long <- melt(mRNA_abundance_changes, id.vars = "ORF")
#       colnames(mRNA_abundance_changes.long)[2:3] <- c("experiment","log2.mRNA_abundance")
#       
#     master.table.2 <-   merge(master.table.2, mRNA_abundance.long, by=c("Gene","experiment"), all.x=T)
#     master.table.2 <-   merge(master.table.2, mRNA_abundance_changes.long, by=c("ORF","experiment"), all.x=T)
# 
#     master.table.2$time       <- sapply( strsplit( as.character(master.table.2$experiment), split = "_"), function(x) x[[3]] )
#     master.table.2$experiment <- sapply( strsplit( as.character(master.table.2$experiment), split = "_"), function(x) x[[2]] )
# 
#     master.table.2$faster.translation <- ifelse( master.table.2$ratio_total.time < 1, T, F)
#     master.table.2$global.protein_synthesis.rate <- 2 ^ (master.table.2$log2.mRNA_abundance ) * (master.table.2$expected.translation_time.normal) / (master.table.2$expected.translation_time) # relative change in mRNA abundance multiplied by the change in translation rate (1/tau)
# 
#     return(master.table.2)
#       
# }

#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

compile.master.table.new <- function( in.path = "data/SMoPT/batch/input/", out.path = "data/SMoPT/batch/output/", ANALYSIS = "GFP"   ){
  
#   require(StingRay)
#   data(SEQ)
#   gene.lengths <- data.frame ( ORF = names(SEQ$ORF_translated), n.codons = width(SEQ$ORF_translated) - 1 )
# 
  
  experiments <- list.files(path = out.path, full.names = F )
  files <- melt(lapply( setNames( paste(out.path, experiments, sep="/"), experiments), list.files))[,c(2,1)]
  colnames(files) <- c("experiment", "time")
 # files$time <- as.numeric(as.character(files$time))
  
#  genes <- read.table( paste0(in.path,"S.cer.mRNA.abndc.ini.tsv"), sep="\t", header=1)
#  genes$Gene <- 0:(nrow(genes)-1)
  
  if(ANALYSIS == "time") {
      require(StingRay)
      load("~/Documents/MRC09/data/PLS-PM/df.PLSPM.non_transformed_non_scaled.Rda")
      df.PLSPM.nTnS$ORF <- fORF(df.PLSPM.nTnS$gene)
      colnames(df.PLSPM.nTnS)[1] <- "name"
      genes <- read.table( paste0(in.path,"S.cer.mRNA.abndc.ini.tsv"), sep="\t", header=1)
      genes$Gene <- 0:(nrow(genes)-1)
      genes <- merge(genes, df.PLSPM.nTnS[,c("ORF","name","length_CDS")], by="ORF", all.y=T)
    
  } else {
    genes <- read.table( paste0(in.path,"S.cer.mRNA.abndc.ini.tsv"), sep="\t", header=1)
    genes$Gene <- 0:(nrow(genes)-1)
  }
  
  
  
  mRNA_abundance <- read.table( paste0(in.path,"mRNA_abundance.txt"), sep="\t", header=1)
  #mRNA_abundance$normal_0  <- genes$rand_mRNA
  mRNA_abundance_changes <- read.table( paste0(in.path,"mRNA_abundance_change.txt"), sep="\t", header=1)
  mRNA_abundance_changes$normal_0 <- 0
  
  if( ANALYSIS == "time"){
  # collect the unformatted data
  require(StingRay)
    data <- apply( files, 1, function(x){ 
      list( 
        initiation = read.table( paste0(out.path, x[["experiment"]],"/",sprintf("%03.0f", as.numeric(x["time"])),"/Scer_gene_initimes.out",  sep="" ), sep="\t", header=1),
        elongation = read.table( paste0(out.path, x[["experiment"]],"/",sprintf("%03.0f", as.numeric(x["time"])),"/Scer_gene_totetimes.out", sep="" ), sep="\t", header=1)
      ) 
    })
    names(data) <- paste( files$experiment, files$time, sep="_" )
  }
  if( ANALYSIS == "GFP"){
    # bug to fix: NORMAL FAILS TO RUN
    FILES <- subset(files, experiment != "normal")
    data <- apply( FILES , 1, function(x){ 
      list( 
        initiation = read.table( paste0(out.path, x[["experiment"]],"/",x["time"],"/Scer_gene_initimes.out",  sep="" ), sep="\t", header=1),
        elongation = read.table( paste0(out.path, x[["experiment"]],"/",x["time"],"/Scer_gene_totetimes.out", sep="" ), sep="\t", header=1)
      ) 
    })
    names(data) <- paste( FILES$experiment, FILES$time, sep="_" )
  }
  
  # format the data into a master table
  data2 <- ldply( data, function(x){ merge( x$initiation, x$elongation, by=c("Gene","Num_of_events")) } )
  data2$Avg_initiation_time.sec.[which( is.nan(data2$Avg_initiation_time.sec.) | data2$Avg_initiation_time.sec.==0)] <- NA
  
  # reference data set for normal conditions
  normal.conditions_initiation <- read.table("data/SMoPT/batch/output/normal/000/Scer_gene_initimes.out", header=1, sep="\t")
  colnames(normal.conditions_initiation)[2:3] <- c("Num_of_events.initiation.normal", "av.initiation_time.normal")
  normal.conditions_initiation$av.initiation_time.normal[ is.nan(normal.conditions_initiation$av.initiation_time.normal) | normal.conditions_initiation$av.initiation_time.normal == 0 ] <- NA
  
  # reference data set for normal conditions
  normal.conditions_elongation <- read.table("data/SMoPT/batch/output/normal/000/Scer_gene_totetimes.out", header=1, sep="\t")
  colnames(normal.conditions_elongation)[2:3] <- c("Num_of_events.elongation.normal", "av.elongation_time.normal")
  
  data3 <- merge( data2, normal.conditions_initiation, by = "Gene", all.x=T )
  data3 <- merge( data3, normal.conditions_elongation, by = "Gene", all.x=T )
  
  # convert to minutes (will subsequently rename columns) FEB2015: OR NOT! codons/sec is easier I tend to think now
  time.variables <- c("Avg_initiation_time.sec.","Avg_total_elong_time.sec.","av.initiation_time.normal","av.elongation_time.normal")
  data3[,time.variables] <- data3[,time.variables] #/60
  data3$delta_initiation <- data3$Avg_initiation_time.sec. - data3$av.initiation_time.normal
  data3$delta_elongation <- data3$Avg_total_elong_time.sec. - data3$av.elongation_time.normal
  data3$ratio_initiation <- data3$Avg_initiation_time.sec. / data3$av.initiation_time.normal
  data3$ratio_elongation <- data3$Avg_total_elong_time.sec. / data3$av.elongation_time.normal
  data3$ratio_events     <- data3$Num_of_events / data3$Num_of_events.initiation.normal
  data3$ratio_events[is.infinite(data3$ratio_events)] <- NA
  data3$ratio_events[is.nan(data3$ratio_events)] <- NA
  # numeric variables
  num.var <- names( which( sapply( data3, is.numeric) == T) ) 
  data3[, num.var] <- do.call(data.frame, lapply(data3[, num.var], function(x) replace(x, is.infinite(x),NA))) 
  data3[, num.var] <- do.call(data.frame, lapply(data3[, num.var], function(x) replace(x, is.nan(x),NA))) 
  # function(x){ x[which(is.infinite(x))]  <- NA; x}
  
  # check that every gene (3,795) is represented in each experiment
  # ddply(data3, .(.id), nrow)
  
  # master table (including gene names, length and average PARS score - though there are missing atm)
  master.table <- data3
  
  colnames(master.table)[2] <- "stress"
  colnames(master.table)[3] <- "n.events"
  colnames(master.table)[4] <- "av.initiation_time"
  colnames(master.table)[5] <- "av.elongation_time"
  colnames(master.table)[6] <- "n.events.normal"
  
  require(StingRay)
  master.table <- merge(master.table, genes, by="Gene", all.y=T)
  if(ANALYSIS == "time"){ 
    master.table$name <- fNAME(master.table$ORF) 
  } else {
    other.genes <- setdiff(as.character(master.table$ORF), as.character(NAMES$systematic.name))  
    master.table$name <- as.character(master.table$ORF)
    master.table$name [which(!master.table$name %in% other.genes)] <- fNAME( as.character(master.table$name[which(!master.table$name %in% other.genes)]) )
    
  }
  master.table <- master.table[,c("ORF","Gene", "name", setdiff( colnames(master.table), c("ORF","Gene","name") )) ]
  colnames(master.table)[grep("rand_mRNA",colnames(master.table))] <- "mRNA_abundance.normal"
  head(master.table)
  
  # variable to subsequently reject entries with wrongly estimated initiation/elongation values 
  
    # mRNA abundance data
    mRNA_abundance.long <- melt(mRNA_abundance, id.vars = "Gene")
    colnames(mRNA_abundance.long)[2:3] <- c("stress","mRNA_abundance")
    
    mRNA_abundance_changes.long <- melt(mRNA_abundance_changes, id.vars = "ORF")
    colnames(mRNA_abundance_changes.long)[2:3] <- c("stress","log2.mRNA_abundance")
  
  master.table <-   merge(master.table, mRNA_abundance.long, by=c("Gene","stress"), all.x=T)
  master.table <-   merge(master.table, mRNA_abundance_changes.long, by=c("ORF","stress"), all.x=T)
  
  # avoid errors for mRNAs with abundance lower than 1. Replace 0 by it's expected value (given the log2). 
  # I think Marc replaced them by 1 by default (he must have done something since there is recorded translation for them, so they couldn't be 0 in abundance)
  recomputed.abundance <- master.table$mRNA_abundance.normal * 2^(master.table$log2.mRNA_abundance)  
  master.table$mRNA_abundance[ !is.na(master.table$mRNA_abundance) & master.table$mRNA_abundance == 0 ]  <- recomputed.abundance[ !is.na(master.table$mRNA_abundance) & master.table$mRNA_abundance == 0 ]
  master.table$valid  <- complete.cases(master.table[,c("av.initiation_time","av.elongation_time","av.initiation_time.normal","av.elongation_time.normal")])
  
  # --------------------------------------------------------------------------------------------------------- #

  # EDIT FEB2015: more correct estimation of translation rates

  # apply linear regressions of the expected translation time in function of the gene's length - and take residuals to estimate 
  # the portion of translation time that is indepenent of length (i.e. find how fast genes are translated while accounting for their length)
  #L <- dlply( master.table.2, .(experiment), function(x) { lm( expected.translation_time ~ length_CDS, data = x ) })
      #   conditions <- unique(master.table$stress)
      #   L1 <- lapply( conditions, function(x) {  data = subset(master.table.2, stress==x); rownames(data) <- data$ORF; lm( total.translation_time ~ length_CDS, data = data )   } )
      #   names(L1) <- conditions
      #   resid.L1 <- ldply(L1, function(x){ data.frame( ORF = names(x$residuals), residual.translation.time = x$residuals)  } )
      #   colnames(resid.L1)[1] <- "stress"
      #   
      #   
      #   L2 <- lapply( conditions, function(x) {  data = subset(master.table.2, stress==x); rownames(data) <- data$ORF; lm( delta_elongation ~ length_CDS, data = data )   } )
      #   names(L2) <- conditions
      #   resid.L2 <- ldply(L2, function(x){ data.frame( ORF = names(x$residuals), residual.delta_elongation = x$residuals)  } )
      #   colnames(resid.L2)[1] <- "stress"
      #   
      #   master.table <- merge( master.table, resid.L1, by=c("ORF","stress"),all.x=T)
      #   master.table <- merge( master.table, resid.L2, by=c("ORF","stress"),all.x=T)
      #   
  
  if (ANALYSIS == "time"){
        master.table$n.codons  <- master.table$length_CDS/3
        master.table$elongation.speed.normal <- (master.table$n.codons)/master.table$av.elongation_time.normal
        master.table$elongation.speed        <- (master.table$n.codons)/master.table$av.elongation_time
        
        # "APPARENT" = n.events of translation for gene g / n. of g-mRNAs 
        master.table$apparent.translation.rate        <- master.table$n.events / master.table$mRNA_abundance
        master.table$apparent.translation.rate.normal <- master.table$n.events.normal / master.table$mRNA_abundance.normal
        master.table$ratio_apparent.translation.rate  <- master.table$apparent.translation.rate / master.table$apparent.translation.rate.normal
        master.table$ratio_apparent.translation.rate[is.infinite(master.table$ratio_apparent.translation.rate)] <- NA
        
        # "EXPECTED" frequency of initiation x the speed it takes to elongate
        master.table$expected.translation_rate.normal <- 1/master.table$av.initiation_time.normal * master.table$elongation.speed.normal
        master.table$expected.translation_rate        <- 1/master.table$av.initiation_time * master.table$elongation.speed
        master.table$ratio_expected.translation.rate <- (master.table$expected.translation_rate ) / (master.table$expected.translation_rate.normal)
        
        # "TOTAL" sum av. waiting time for initiation + av. elongation time
        # simply adding the av. intitation and elongation time to get an estime of the expected total translation time
        master.table$total.translation_time.normal <- master.table$av.initiation_time.normal + master.table$av.elongation_time.normal
        master.table$total.translation_time        <- master.table$av.initiation_time + master.table$av.elongation_time
        master.table$ratio_total.time <- (master.table$total.translation_time ) / (master.table$total.translation_time.normal)
        
    features.to.add <- c("proto_genes","overexpression_toxicity","average_cost","total_metabolic_cost","residual_cost","overal_charge","protein_length","poor_STOP_poor","AAA_before_STOP","weak_tetra","NAA_before_STOP","STOP","AUGCAI","translation_RII","mean_nTE","tRNA_adaptation_index","ribosome_density_YPD","ribosome_density_SD","protein_abundance_ypd","protein_abundance_sd","protein_abundance","protein_abundance_H2O2")
    Additions <- df.PLSPM.nTnS[ c("name",features.to.add)]
    # Merge data into master table
    master.table.2 <- merge(master.table, Additions, by="name", all.x=T)
    master.table.2$experiment <- sapply( strsplit( as.character(master.table.2$stress), split = "_"), function(x) x[[1]] )
    master.table.2$time       <- NA
    master.table.2$time[ !is.na(master.table.2$stress) ] <- sapply( strsplit( as.character(master.table.2$stress[!is.na(master.table.2$stress)] ), split = "_"), function(x) x[[2]] )
    # ---- additional processing of the data
    # Add relevant control variables to the master table
    master.table.2$zscore.delta_initiation          <- zscore( master.table.2$delta_initiation )
    master.table.2$zscore.delta_elongation          <- zscore( master.table.2$delta_elongation )
  #   master.table.2$zscore.residual.delta_elongation <- zscore( master.table.2$residual.delta_elongation )
  #   master.table.2$zscore.residual.translation.time <- zscore( master.table.2$residual.translation.time )
  #   
    
    master.table.2$global.protein_synthesis.rate <- 2 ^ (master.table.2$log2.mRNA_abundance ) * (master.table.2$total.translation_time.normal) / (master.table.2$total.translation_time) # relative change in mRNA abundance multiplied by the change in translation rate (1/tau)
  #  master.table.2$change.translation.rate  <- master.table.2$ratio_events * 2 ^ ( - master.table.2$log2.mRNA_abundance) # n_events(s)/mRNA(s) / n_events(n)/mRNA(n)
      
  # ggplot(master.table.2, aes(x=change.translation.rate, y= ratio_apparent.translation.rate)) + geom_point()
  # with( data = master.table.2, expr = cor.test(change.translation.rate, ratio_apparent.translation.rate))
  
    master.table.2$faster.translation   <- ifelse( master.table.2$ratio_total.time < 1, T, F)
    master.table.2$faster.translation.2 <- ifelse( master.table.2$ratio_apparent.translation.rate > 1, T, F)
    master.table.2$increased.production <- ifelse( master.table.2$ratio_events > 1, T, F)
    
    
  } else {
    master.table.2 <- master.table
  }
  
  return(master.table.2)
}







get.stats <- function(master.table){
  stats  <- list()

  master.table <- subset(master.table, !is.na(faster.translation) )
    
  stats$n_faster <- length(unique(subset(master.table, faster.translation == T)$Name))
  
  EXP1 <- unique( master.table$experiment )[1]
  
  y <- ddply(subset(master.table,!is.na(faster.translation) & experiment == EXP1 ), .(faster.translation), function(x){median(x$length_CDS/3, na.rm=T)}); 
  stats$effect_length = round( (y[2,2]-y[1,2])/median(subset(master.table,!is.na(faster.translation))$length_CDS/3, na.rm=T),2)*100
  
  y <- ddply(subset(master.table,!is.na(faster.translation) & experiment == EXP1 ), .(faster.translation), function(x){median(x$tRNA_adaptation_index, na.rm=T)})           
  stats$effect_tAI = round( (y[2,2]-y[1,2])/median(subset(master.table,!is.na(faster.translation))$tRNA_adaptation_index, na.rm=T),2)*100
  
  y.av_ini <- ddply(subset(master.table,!is.na(faster.translation) & experiment == EXP1 ), .(faster.translation), function(x){median(x$av.initiation_time, na.rm=T)})         
  stats$median_diff_initiation_temp <- round( (y.av_ini[2,2]-y.av_ini[1,2]),2)
  
  y.ratio_ini <- ddply(subset(master.table,!is.na(faster.translation) & experiment == EXP1 ), .(faster.translation), function(x){median(x$ratio_initiation, na.rm=T)})         
  stats$effect_ratio_initiation_temp <- round( (y.ratio_ini[2,2]-y.ratio_ini[1,2])/median(subset(master.table,!is.na(faster.translation))$ratio_initiation, na.rm=T),2)*100
  
  y.median_diff_ini_normal <- y.av_ini[,2]/y.ratio_ini[,2] 
  stats$median_diff_initiation_normal <- round( y.median_diff_ini_normal[2] - y.median_diff_ini_normal[1] ,2)

  y <- ddply( master.table, .(experiment), function(x){ sum(x$faster.translation==T) / nrow(x) } ) ### REMEMBER HERE  
  stats$prop_faster <- y[,2]
  names(stats$prop_faster) <- y[,1]
  
  return(stats)
}

compile.codon.elongation.matrix <- function(  path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/" ){
  
      codons <- read.table("~/Documents/MRC16/data/SMoPT/stochastic.simulation/template.txt", header = 1)[,c("codon","c.number")]
      #codons <- read.table("~/Documents/MRC16/data/SMoPT/stochastic.simulation/codons.txt", header=1) # no, wrong!
      colnames(codons)[2] <- "Codon"
      experiments <- list.files(path = path, full.names = F )
      files <- lapply( setNames( paste(path, experiments, sep="/"), experiments), list.files)
      
      require(StingRay)
      # collect the unformatted data
      data <- lapply( names(files)  , function(x){ codon_elongation = read.table( paste(path, x ,"/", files[[x]][2], sep="" ), sep="\t", header=1) })
      names(data) <- names(files) #, "genome_normal_conditions")
      data$genome_normal_0 <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_etimes.out.txt", header=1)
      colnames(data$genome_normal_0)[3] <- "Avg_elong_time.sec."
      
      # format the data into a master table
      data2 <- ldply( data )
      colnames(data2)[1] <- "experiment"
  
       # reference data set for normal conditions
      normal.conditions_codons <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_etimes.out.txt", header=1, sep="\t")
      
      data3 <- merge( data2, normal.conditions_codons, by = "Codon", all.x=T )
      colnames(data3) <- c("Codon","experiment", "num_events","av.elongation_time","num_events.normal","av.elongation_time.normal")
      
      # convert to minutes (will subsequently rename columns)
      time.variables <- c("av.elongation_time","av.elongation_time.normal")
      data3[,time.variables] <- data3[,time.variables]
      data3$delta_elongation <- data3$av.elongation_time - data3$av.elongation_time.normal
      data3$ratio_elongation <- data3$av.elongation_time / data3$av.elongation_time.normal
      data3$ratio_events <- data3$num_events / data3$num_events.normal
     
      master.table <- merge( codons, data3, by = "Codon", all.y=T)
      
      require(StingRay); data(TAB); tAI <- subset( TAB$tAI[,c("codon","Scer")], !is.nan(Scer) )
      colnames(tAI)[2] <- "tAI"
      
      master.table <- merge( master.table, tAI, by = "codon")
      
      return(master.table)
}


#---------------------------------------------------------------------------------------------------------#


compile.codon.elongation.matrix.new <- function(  path = "data/SMoPT/batch/output/" ){
  
  codons <- read.table("~/Documents/MRC16/data/SMoPT/stochastic.simulation/template.txt", header = 1)[,c("codon","c.number")]
  colnames(codons)[2] <- "Codon"
  
  experiments <- list.files(path = path, full.names = F )
  files <- melt(sapply( setNames( paste(path, experiments, sep="/"), experiments), list.files))[,c(2,1)]
  colnames(files) <- c("experiment", "time")
  files$time <- as.numeric(as.character(files$time))
  
  require(StingRay)
  # collect the unformatted data
#  data <- lapply( names(files)  , function(x){ codon_elongation = read.table( paste(path, x ,"/", files[[x]][2], sep="" ), sep="\t", header=1) })
  
  data <- apply( files, 1, function(x){
    read.table( paste(path, x[["experiment"]],"/",sprintf("%03.0f", as.numeric(x["time"])), "/Scer_etimes.out", sep=""), header=1 )
  })
  
  names(data) <- paste( files$experiment, files$time, sep="_" )
  
  # format the data into a master table
  data2 <- ldply( data )
  colnames(data2)[1] <- "stress"
  
  # reference data set for normal conditions
  normal.conditions_codons <- subset(data2, stress=="normal_0")[,-1]
  
  data3 <- merge( data2, normal.conditions_codons, by = "Codon", all.x=T )
  colnames(data3) <- c("Codon","stress", "num_events","av.elongation_time","num_events.normal","av.elongation_time.normal")
  
  # convert to minutes (will subsequently rename columns)
  time.variables <- c("av.elongation_time","av.elongation_time.normal")
  data3[,time.variables] <- data3[,time.variables]
  data3$delta_elongation <- data3$av.elongation_time - data3$av.elongation_time.normal
  data3$ratio_elongation <- data3$av.elongation_time / data3$av.elongation_time.normal
  data3$ratio_events <- data3$num_events / data3$num_events.normal
  
  master.table <- merge( codons, data3, by = "Codon", all.y=T)
  
  require(StingRay); data(TAB); tAI <- subset( TAB$tAI[,c("codon","Scer")], !is.nan(Scer) )
  colnames(tAI)[2] <- "tAI"
  
  master.table <- merge( master.table, tAI, by = "codon")
  
  return(master.table)
}

# genes that go down in total translation time upon stress
extract.venn.data <- function( dataset= master.table.020 ){
  require(venneuler)
  
  list.genes <- dlply( subset( dataset, faster.translation == T), .(experiment), function(x) x$Name )
  
  venn.genes <- melt( list.genes )
  venn.genes$L2 <- 1
  
  test <- dcast(venn.genes, value ~ L1, value.var = "L2", fill=0 ) # adjacency df
  
  output <- list()
  output$coded <- data.frame( gene = test[,1], code = apply( test[,-1], 1, code.it ) )
  combinations <- as.matrix(test[,-1])
  rownames(combinations) <- test[,1]
  
  # Venn diagram
  output$venn  <- venneuler( combinations = combinations )
  
  output$venn$labels <- paste( output$venn$labels, laply(list.genes, length), sep="\n" ) # show the (total) number of genes per sets
  
  output$in_all  <- Reduce( intersect, list.genes ) # get the intersection of all sets
  output$jaccard <- length(output$in_all) / length(unique(output$coded$gene))
  return( output )   
}

# Please ignore my last email, I was writing from my iPad and thinking out loud, but I should have waited a bit: 
#   
#   The basic ideas of our calculations are (obviously) similar, but differed about precision (especially when and how rounding are used), and here is why: 
#   
#   I first computed the theoretical abundance of each tRNA in normal conditions:
#   N_tRNAn(i) = 3.3x10^6 x tGCNn(i) / sum tGCNn
# 
# Then the N_tRNA under stress:
#   N_tRNAs(i) = N_tRNAn(i) x FC
# 
# Which I then used to compute adjusted tGCN
# tGCNs(i) = round( tGCNn(i) x FC )

compile.tRNAabundance.master.table <- function(total_number_tRNAs.normal = 3.3*10^6, # parameter used by Shah20013 and us
                                               compute.from.tGCN = F,
                                               add.GAU = F
                                               ){
  # Info table
  rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",sep="\t", header=1)
  # Marc's relative quantifications of tRNA abundance
  tRNAab_quantif <-  read.table("data/tRNA abundance/tRNAab_quantified.foldchange.txt",sep="\t", header=1)
  tRNAab_quantif$normal_0 <- 1
  tRNA.abundance <- melt(tRNAab_quantif, id.vars = "anticodon")
  colnames(tRNA.abundance)[2:3] <- c("stress", "foldchange") 
  tRNA.abundance$time       <- sapply( strsplit( as.character(tRNA.abundance$stress), split = "_"), function(x) x[[2]] )
  tRNA.abundance$experiment <- sapply( strsplit( as.character(tRNA.abundance$stress), split = "_"), function(x) x[[1]] )
  tRNA.abundance$log2.foldchange <- log2(tRNA.abundance$foldchange)
  
  # Computing (1) a table to distinguish elongating and intiating Methionine tRNAs and (2) account for their relative contribution to a global 
  # tRNA abundance of CAU (necessary for the SMoPT analysis where no distinction can be made between CAU and CAU2)
  
  CAU.table <- merge( subset(tRNA.abundance, anticodon == "CAU")[,-c(1,6)], subset(tRNA.abundance, anticodon == "CAU2")[,-c(1,6)], by = c("stress","experiment","time") )
  colnames(CAU.table)[c(4:5)] <- c("foldchange.eMet_tRNA","foldchange.iMet_tRNA")
  
  CAU.table$adjusted.tGCN <- round( 5 * CAU.table$foldchange.eMet_tRNA + 5 * CAU.table$foldchange.iMet_tRNA )
  # because the fold changes are defined compare to the same condition where both CAU.1 and CAU.2 have equal tGCN (5)
  CAU.table$foldchange    <- 0.5*( CAU.table$foldchange.eMet_tRNA + CAU.table$foldchange.iMet_tRNA )
  write.table(CAU.table, file="results/master tables/CAU.table.txt",sep="\t",quote=F, row.names=F)
  
  
  # Integrate data on abundance of free tRNAs
  #free.tRNA_t020   <- compile.free.tRNA.copy.numbers(path="data/SMoPT/stochastic.simulation/unformated/20 minutes/")
  #free.tRNA_t120   <- compile.free.tRNA.copy.numbers(path="data/SMoPT/stochastic.simulation/unformated/120 minutes/")
  
  #free.tRNA <- unique( rbind(free.tRNA_t020, free.tRNA_t120) )
  free.tRNA <- compile.free.tRNA.copy.numbers.new()
  
  
  # dcast( free.tRNA, experiment + time ~ anticodon, value.var =  "tGCN" )
  
  # warning: copy.number (data from SMoPT) and tGCN (which I collected from Marc's textbook table) differ! # I will use copy.number to closely
  # match with the conditions of the stochastic simulation
  #rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",sep="\t", header=1)
  # arrange( merge( unique( subset(free.tRNA, select= c(anticodon, copy.number) ) ),
  #         subset(rosetta.anticodons, select=c(anticodon, tGCN)), by="anticodon"), anticodon )
  
  # AUGUST 2015: quick patch. The addition of free.tRNA to the master data frame is a bit obsolete. It is helpful to provide tGCN, but it's additional
  # information relates to simulations of translation (which ideally should be added much later, since this current script is meant to help build a master table with basic information).
  # the consequences of the repurposing of scripts and procedures during the project weren't fully address in the interest of time.
  # the fix consists of pre-merging a subset of free.tRNA without the simulation data, and add these subsequently, so that one can add info on time points
  # for which simulation data may not be available (typically t=60min).
  
  theoretical.abundance <- unique( subset(free.tRNA, select= c(anticodon, tGCN, label) ))
  theoretical.abundance$anticodon <- as.character(theoretical.abundance$anticodon)
  if( add.GAU == T ){
    theoretical.abundance <- rbind(theoretical.abundance, c(anticodon = "GAU", tGCN = 1)) # add the missing tRNA tI(GAU)Q
    theoretical.abundance$tGCN <- as.numeric(theoretical.abundance$tGCN)  
  }
    theoretical.abundance$total.tRNAab.normal <- theoretical.abundance$tGCN / sum(theoretical.abundance$tGCN) * total_number_tRNAs.normal
  
  
  tRNA.abundance.master.table <- merge( tRNA.abundance, theoretical.abundance, by="anticodon")
  tRNA.abundance.master.table <- merge( tRNA.abundance.master.table, free.tRNA, by = c("anticodon","experiment","time", "stress","tGCN","label"), all.x=T )
  # patch APRIL 2015 to discriminate w^s of rare tRNAs
  tRNA.abundance.master.table$adjusted.tGCN  <- with(tRNA.abundance.master.table, round(tGCN * foldchange, 1)  )
  #tRNA.abundance.master.table[tRNA.abundance.master.table$adjusted.tGCN == 0, ]$adjusted.tGCN <- 1 # readjust to 1 tRNA at least 
  
#   if(compute.from.tGCN == T){
#       #tRNA.abundance.master.table$total.tRNAab <- tRNA.abundance.master.table$total.tRNAab.normal * tRNA.abundance.master.table$adjusted.tGCN / tRNA.abundance.master.table$tGCN
#   } else {
#       #tRNA.abundance.master.table$total.tRNAab <- tRNA.abundance.master.table$total.tRNAab.normal * tRNA.abundance.master.table$foldchange  
#   }
#   
  tRNA.abundance.master.table$total.tRNAab      <- tRNA.abundance.master.table$total.tRNAab.normal * tRNA.abundance.master.table$foldchange  
  tRNA.abundance.master.table$total.tRNAab.marc <- tRNA.abundance.master.table$total.tRNAab.normal * tRNA.abundance.master.table$adjusted.tGCN / tRNA.abundance.master.table$tGCN
  
  # merge with the rosetta table
  tRNA.abundance.master.table <- merge( tRNA.abundance.master.table, rosetta.anticodons, by=c("anticodon","tGCN") ) 
  
  # relative availability: how many genes a particular tRNA has in relation to all tRNAs coding for the same amino acid?
  tRNA.abundance.master.table <- ddply( tRNA.abundance.master.table, .(aa, stress), mutate, relative.availability = round(total.tRNAab / sum(total.tRNAab),2)  )
  
  # reorder columns
  tRNA.abundance.master.table <- subset(tRNA.abundance.master.table, 
                                        select=c(anticodon, aa, codon, experiment, time, stress, 
                                                 tGCN, adjusted.tGCN,
                                                 foldchange, log2.foldchange, 
                                                 free.tRNAab, total.tRNAab, total.tRNAab.marc,
                                                 free.tRNAab.normal, total.tRNAab.normal,
                                                 delta_free.tRNAab, ratio_free.tRNAab, 
                                                 relative.availability,relative.availability.normal,
                                                 CAI.O, nTE.O, pos.1_anticodon, pos.2.3_anticodon, pos.3_codon, 
                                                 aa.with.1.anticodon, label,
                                                 topology
                                                 ))

  all.tRNAs <- tRNA.abundance.master.table
  #merge(all.tRNAs, CAU.table[,c("stress","foldchange","adjusted.tGCN")], by = c("stress"))

#   CAU <- subset(CAU.table, time != 60 )
#   all.tRNAs[ which(  as.character(CAU$stress) %in% all.tRNAs$stress & all.tRNAs$anticodon =="CAU" ), 
#              c("stress","adjusted.tGCN","foldchange") ] <- arrange( CAU[,c("stress","adjusted.tGCN","foldchange")], stress )
#   
  all.tRNAs$log2.foldchange <- log2(all.tRNAs$foldchange)

  return(all.tRNAs) 
}




# no no no, currently I look into FREE tRNA abundance, NOT total tRNA abundance!
compile.free.tRNA.copy.numbers  <- function( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/"
                                           ){  
      # RETRIEVE DATA FROM THE SIMULATION OUTPUTS ON THE QUANTITY OF FREE tRNAs
      experiments <- list.files(path = path, full.names = F )
      files <- lapply( setNames( paste(path, experiments, sep="/"), experiments), list.files)
      # anticodons <- read.table("data/SMoPT/stochastic.simulation/anticodons.txt",header=1) # No. Wrong.
      # colnames(anticodons)[3] <- "tGCN"
      
      anticodons <- subset( read.table("data/SMoPT/utilities/S.cer.tRNA_copy_num",header=1), gcn != 0 )
      anticodons$label  <- paste0("Free_tRNA",0:(nrow(anticodons)-1))
      colnames(anticodons)[1:2] <- c("anticodon","tGCN")
      anticodons$anticodon <- gsub( pattern = "T", replacement = "U", x= as.character(anticodons$anticodon))
      
        write.table(anticodons, "data/SMoPT/stochastic.simulation/anticodons.txt", sep="\t", quote = F , row.names=F)
      
      # DEBUG: VERIFY CONSISTENCY OF tGCN SCORES
      anticodons.2 <- read.table("data/SMoPT/utilities/S.cer.tRNA_copy_num",header=1)
      
      anticodons.2$tRNA <- gsub( "T", "U", x = anticodons.2$tRNA)
      control <- subset( merge(anticodons, anticodons.2,by.x="anticodon",by.y="tRNA"), tGCN !=gcn)
      if (nrow(control) > 0 ){ control }
      
      data <- lapply( names(files) , function(x){ read.table( paste(path, x ,"/", files[[x]][1], sep="" ), sep="\t", header=F)
                                                                           }
            )
      names(data) <- names(files)
      data$genome_normal_0 <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_avg_ribo_tRNA.out.txt", header=F)


      # format the data into a master table
      data2 <- ldply( data )
      colnames(data2) <- c("experiment","label","free.tRNAab")

      normal.conditions_tRNA <- read.table("data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/output_avg_ribo_tRNA.out.txt", header=F, sep="\t")
      colnames(normal.conditions_tRNA) <- c("label", "free.tRNAab.normal")
      
      data3 <- merge( data2, normal.conditions_tRNA, by = "label", all.x=T )
      
      data3$delta_free.tRNAab <- data3$free.tRNAab - data3$free.tRNAab.normal
      data3$ratio_free.tRNAab <- data3$free.tRNAab / data3$free.tRNAab.normal

      data3$time       <- sapply( strsplit( as.character(data3$experiment), split = "_"), function(x) x[[3]] )
      data3$experiment <- sapply( strsplit( as.character(data3$experiment), split = "_"), function(x) x[[2]] )

      # merge with anticodon ID
      return( merge(data3, anticodons, by="label") )
}

#--------------------------------------------------------------------------------------------------------------------------------#

# no no no, currently I look into FREE tRNA abundance, NOT total tRNA abundance!
compile.free.tRNA.copy.numbers.new  <- function( path = "data/SMoPT/batch/output/", recreate.anticodons.table = F ){  
  # RETRIEVE DATA FROM THE SIMULATION OUTPUTS ON THE QUANTITY OF FREE tRNAs
  experiments <- list.files(path = path, full.names = F )
  files <- melt(sapply( setNames( paste(path, experiments, sep="/"), experiments), list.files))[,c(2,1)]
  colnames(files) <- c("experiment", "time")
  files$time <- as.numeric(as.character(files$time))
  
  
  if ( recreate.anticodons.table == T){
      anticodons <- subset( read.table("data/SMoPT/utilities/S.cer.tRNA_copy_num",header=1), gcn != 0 )
      anticodons$label  <- paste0("Free_tRNA",0:(nrow(anticodons)-1))
      colnames(anticodons)[1:2] <- c("anticodon","tGCN")
      anticodons$anticodon <- gsub( pattern = "T", replacement = "U", x= as.character(anticodons$anticodon))
      
      write.table(anticodons, "data/SMoPT/stochastic.simulation/anticodons.txt", sep="\t", quote = F , row.names=F)    

      # VERIFY CONSISTENCY OF tGCN SCORES
      anticodons.2 <- read.table("data/SMoPT/utilities/S.cer.tRNA_copy_num",header=1)
      anticodons.2$tRNA <- gsub( "T", "U", x = anticodons.2$tRNA)
      control <- subset( merge(anticodons, anticodons.2,by.x="anticodon",by.y="tRNA"), tGCN !=gcn)
      
      if (nrow(control) > 0 ){ control }
      
  } else {
    anticodons  <- read.table("data/SMoPT/stochastic.simulation/anticodons.txt", header=1)
  }

  data <- apply( files, 1, function(x){
    read.table( paste(path, x[["experiment"]],"/",sprintf("%03.0f", as.numeric(x["time"])), "/Scer_avg_ribo_tRNA.out", sep=""), header=F )
  })
  
  names(data) <- paste( files$experiment, files$time, sep="_" )
  
  # format the data into a master table
  data2 <- ldply( data )
  colnames(data2) <- c("stress","label","free.tRNAab")
    
  normal.conditions_tRNA <- subset(data2, stress == "normal_0")[,-1]; colnames(normal.conditions_tRNA)[2] <- "free.tRNAab.normal"
  
  data3 <- merge( data2, normal.conditions_tRNA, by = "label", all.x=T )
  
  
  data3$delta_free.tRNAab <- data3$free.tRNAab - data3$free.tRNAab.normal
  data3$ratio_free.tRNAab <- data3$free.tRNAab / data3$free.tRNAab.normal
  
  data3$time       <- sapply( strsplit( as.character(data3$stress), split = "_"), function(x) x[[2]] )
  data3$experiment <- sapply( strsplit( as.character(data3$stress), split = "_"), function(x) x[[1]] )
  
  # merge with anticodon ID
  return( merge(data3, anticodons, by="label") )
}




map.dynamics <- function(data=ILoveThatGirl, what = "translation.dynamics", GO.slim=NULL){

  var <- data.frame( x = c("log2.mRNA_abundance.x","ratio_total.time.x"), 
                     y = c("log2.mRNA_abundance.y","ratio_total.time.y"), 
                     desc = c("change in mRNA abundance","relative translation time"),
                     what = c("transcription.dynamics","translation.dynamics"),
                     stringsAsFactors  = F
                     ) 
  
  # convert the chosen variables to log2 ratios
  if(what!="transcription.dynamics"){
    VARs <- unlist(var[ which( var$what == what), ])
    data[, VARs[1:2] ] <- as.data.frame( sapply(  VARs[1:2], function(x){ log2(data[,x]) } ) )    
  }
  
  if(!is.null(GO.slim)){
      GO.slim <- subset(GO.slim, ORF %in% data$ORF )
      data2   <- merge(GO.slim, subset(data, ORF %in% GO.slim$ORF), by="ORF" )
  }
  
  
  # n.cat <- ddply(data, .(experiment, translation.dynamics), nrow ) without argument 'what'
  # n.cat <-  do.call("ddply", list( .data = data, .variables = c( "experiment", what ), nrow ) ) # slower than the melt->table approach
  n.cat <- melt( table( data[,c(what, "experiment")] ) ) 
  n.tot <- ddply(data, .(experiment), nrow )
  counts <- merge(n.cat,n.tot,by="experiment")
  colnames(counts)[3:4] <- c("n.cat","n.tot")
  counts$percentage <- paste0(round( counts$n.cat / counts$n.tot * 100, 1), "%")
  corners <- data.frame(x     = c(-Inf, -Inf, Inf, Inf),   y = c(-Inf, Inf, -Inf, Inf), 
                        hjust = c(-0.1,-0.1,1.05,1.05),vjust = c(-2.65,2.95,-2.65,2.95), 
                        dynamics= levels(counts[,what]) # WARNING, levels are alphabetically ordered, and it just works because up <-> slow and down <-> fast #c("fast->fast","fast->slow","slow->fast","slow->slow")
                        )
  colnames(corners)[5] <- what

  labels  <-  data.frame( x     = c(-Inf, -Inf, Inf, Inf),   y = c(-Inf, Inf, -Inf, Inf), 
                        hjust = c(-0.1,-0.1,1.05,1.05),vjust = c(-1,1.5,-1,1.5), 
                        labels= levels(counts[,what])  )

  n.graphs <- merge(corners, counts, by=what, all.y=T)
  
  g <- ggplot(data, aes_string(x= VARs[1], y= VARs[2] ) ) + 
          geom_rect( xmin=0, xmax=Inf, ymin=0, ymax=Inf ,fill="gray90",color="gray99") + 
          geom_vline(xintercept=0,lty=3, color="gray30",size=0.7)+
          geom_hline(yintercept=0,lty=3, color="gray30",size=0.7)+
          # labels
          geom_text( data = labels, aes(x=x,y=y,label=labels,hjust=hjust,vjust=vjust), color="gray30",size=3.5 ) +
          # n 
          geom_text( data = n.graphs, aes(x=x,y=y,label=percentage,hjust=hjust,vjust=vjust), color="gray30",size=3.5 )
  
  if(is.null(GO.slim)){ 
        g <- g %+% geom_point(aes(color = tRNA_adaptation_index)) +
             labs( x= "t= 20min", y ="t= 120min", color="tAI", title = what )
    } else {
        g <- g %+% geom_point( data = subset(data, !ORF %in% GO.slim$ORF ), color = "gray60",pch=21, size=1.5, fill="white" ) +
                   geom_point( data = data2, aes(fill =  GO.description ), pch=21  ) +
                  # geom_text(  data = subset(data, abs(log2(ratio_total.time.x)) > 1.6 & abs(log2(ratio_total.time.y)) > 1 & translation.dynamics != "slow->slow" ) , aes(label=Name), size=2, color="black" ) +
                   labs( x= "t= 20min", y ="t= 120min", fill="" ) + guides(fill = guide_legend(ncol = 3))
    }      
  g <- g %+% facet_wrap( ~ experiment ) + labs(title= VARs[3]) +
          theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position="bottom")
  return(g)
}

quantile.it <- function(x, N=5){
  require(gtools)
  LABELS <- paste( round((1 -seq(1:N)/N)*100 ), round(100/N) + round((1 -seq(1:N)/N)*100 ),sep="-" )  
  if( length(table(x)) >1 ) { 
    quantcut( x, q = seq(0,1,by=1/N), na.rm=T, labels=LABELS, ordered = T )
  } else {
    rep(NA,times = length(x))
  }
}



automate.chi.squared <- function( data = master.table.dynamics, # data set containing the response variable
                                  expression = "global.protein_synthesis.rate.x > 1",
                                  keys = "cell budding",  # category (or categories) of choice
                                  classification.table = go.slim.ORFs # table classifying genes into categories
                                ){
  
                tests <- do.call( rbind, lapply( setNames(keys, keys), function(x){
                # genes in the category
#                selection <- as.character( classification.table[ which(classification.table[,2] == x ),"ORF" ]) 
                 selection <- as.character( subset(classification.table, GO.description %in% x , select=c(ORF)))  
                if(length(selection) > 5 ){                  
                    x_test <- ifelse( eval( parse(text=expression), envir =  data ), T, F ) 
                    y_target <- ifelse( data[["ORF"]] %in% selection, T, F) 
                    chisq <- chisq.test( x = x_test, y = y_target )
                                     
                    table  = data.frame( x = expression, y = x, 
                                         n = sum(rowSums(chisq$observed)),
                                         n_xy.TT = chisq$observed["TRUE","TRUE"],
                                         xy.TT = round( ( chisq$observed["TRUE","TRUE"] / chisq$expected["TRUE","TRUE"] ), 2),
                                         xy.TF = round( ( chisq$observed["TRUE","FALSE"] / chisq$expected["TRUE","FALSE"] ), 2),
                                         xy.FT = round( ( chisq$observed["FALSE","TRUE"] / chisq$expected["FALSE","TRUE"] ), 2),
                                         xy.FF = round( ( chisq$observed["FALSE","FALSE"] / chisq$expected["FALSE","FALSE"] ), 2),
                                         chi2 = round(chisq$statistic,2), 
                                         p.value = chisq$p.value, row.names=NULL )
                  } else { NULL }
                  return(table)
                }
           ) )
  return(tests)
}


#####
automate.fisherexact <- function( data = gene.master.table, # data set containing the response variable
                                  expression = "global.protein_synthesis.rate > 1",
                                  id = "name", keys,
                                  #keys = "cell budding",  # category (or categories) of choice
                                  classification.table = go.slim.ORFs # table classifying genes into categories
){
  
  require(scales)
    
  classification.table <- as.data.frame(classification.table)
  tests <- do.call( rbind, apply(keys, 1, function(k){
    # genes in the category
    #print(x)
    #selection <- as.character( classification.table[ which(classification.table[,"GO.description"] == x ), id ]) 
    print(k[["GO.description"]])
    
    selection <- as.character( subset(classification.table, GO.description %in% k[["GO.description"]] ,
                                      select=c(eval(parse(text=id) )),drop = T))
    
   k_test <- ifelse( eval( parse(text=expression), envir =  data ), T, F ) 
   y_target <- ifelse( data[[id]] %in% selection, T, F)
    
   t <- table(k_test, y_target)
   min.cell <- min(t)
  #  print(table(k_test, y_target)) 
  #  cat( paste( dim(table(k_test)), dim(table(y_target)), min(table(k_test, y_target)), sep="\t"), "\n")
  
  # of interest = black ball 
  p.hyper = phyper( q = t["TRUE", "TRUE"]-1, # number of black balls drawn from the urn
                    k = (t["TRUE","TRUE"]+t["TRUE","FALSE"]), # number of balls drawn from the urn
                    m = (t["TRUE","TRUE"]+t["FALSE","TRUE"]), # total number of black balls in the urn
                    n = (t["TRUE","FALSE"]+t["FALSE","FALSE"]), # total number of white balls in the urn
                    lower.tail = F # P(X>q) 
  )
  
  
      if( 
          dim(table(k_test)) == 2 & dim(table(y_target)) == 2 & min.cell > 5 
        ){
        #cat("ouais, TRUE!!!!\n")
        chisq <- chisq.test( x = k_test, y = y_target, simulate.p.value = T )
        fisher <- fisher.test( x = k_test, y = y_target)
        
     
        table  = data.frame( x = expression, y = as.character(k[["GO.description"]]),
                             type = as.character(k[["GO.type"]]), 
                             n = sum(rowSums(chisq$observed)),
                             min.cell = min.cell,
                             # TP/ (TP+FN)
                             sensitivity = gul_percent(chisq$observed["TRUE","TRUE"] / (chisq$observed["TRUE","TRUE"] + chisq$observed["FALSE","TRUE"])),
                             # TN / (TN + FP)
                             specificity = gul_percent(chisq$observed["FALSE","FALSE"] / (chisq$observed["FALSE","FALSE"] + chisq$observed["TRUE","FALSE"])),
                             # TP / (TP + FP)
                             precision   = gul_percent(chisq$observed["TRUE","TRUE"] / (chisq$observed["TRUE","TRUE"] + chisq$observed["TRUE","FALSE"])),
                             # (TP + TN) / (P+N)
                             accuracy    = gul_percent( (chisq$observed["TRUE","TRUE"] + chisq$observed["FALSE","FALSE"])/sum(rowSums(chisq$observed)) ),
                             n.TT = chisq$observed["TRUE","TRUE"],
                             n.TF = chisq$observed["TRUE","FALSE"],
                             n.FT = chisq$observed["FALSE","TRUE"],
                             n.FF = chisq$observed["FALSE","FALSE"],
                             r.TT = round( ( chisq$observed["TRUE","TRUE"] / chisq$expected["TRUE","TRUE"] ), 2),
                             r.TF = round( ( chisq$observed["TRUE","FALSE"] / chisq$expected["TRUE","FALSE"] ), 2),
                             r.FT = round( ( chisq$observed["FALSE","TRUE"] / chisq$expected["FALSE","TRUE"] ), 2),
                             r.FF = round( ( chisq$observed["FALSE","FALSE"] / chisq$expected["FALSE","FALSE"] ), 2),
                             chi2 = round(chisq$statistic,2), 
                             risk = round( (chisq$observed["TRUE","TRUE"] * sum(chisq$observed["FALSE",]))/(chisq$observed["FALSE","FALSE"] * sum(chisq$observed["TRUE",])) ,2),
                             odds = round(fisher$estimate,2),
                             CI.low = round(fisher$conf.int[1],2),
                             CI.high = round(fisher$conf.int[2],2),
                             p.value = fisher$p.value,
                             p.hyper = p.hyper
        )
        rownames(table) <- NULL
        return(table)
      } else{ 
        #cat("Sample too small.\n");  
              table  = data.frame( x = expression, y = as.character(k[["GO.description"]]),
                                   type = as.character(k[["GO.type"]]), 
                                   n = sum(t),
                                   min.cell = min.cell,
                                   # TP/ (TP+FN)
                                   sensitivity = gul_percent(t["TRUE","TRUE"] / (t["TRUE","TRUE"] + t["FALSE","TRUE"])),
                                   # TN / (TN + FP)
                                   specificity = gul_percent(t["FALSE","FALSE"] / (t["FALSE","FALSE"] + t["TRUE","FALSE"])),
                                   # TP / (TP + FP)
                                   precision   = gul_percent(t["TRUE","TRUE"] / (t["TRUE","TRUE"] + t["TRUE","FALSE"])),
                                   # (TP + TN) / (P+N)
                                   accuracy    = gul_percent( (t["TRUE","TRUE"] + t["FALSE","FALSE"])/sum(t) ),
                                   n.TT = t["TRUE","TRUE"],
                                   n.TF = t["TRUE","FALSE"],
                                   n.FT = t["FALSE","TRUE"],
                                   n.FF = t["FALSE","FALSE"],
                                   r.TT = NA,
                                   r.TF = NA,
                                   r.FT = NA,
                                   r.FF = NA,
                                   chi2 = NA, 
                                   risk = NA,
                                   odds = NA,
                                   CI.low = NA,
                                   CI.high = NA,
                                   p.value = NA,
                                   p.hyper = p.hyper
              )
              rownames(table) <- NULL
              
              return(table)
              
      }
  
  }
  ) )
  
  if( nrow(tests) > 0 ){ 
    tests$adjusted.p <- p.adjust(tests$p.value, method="BH")
    tests$consistent.CI <- ifelse( sign( log(tests$CI.low) * log(tests$CI.high)) == -1, "no", "yes")
    tests$diagnostic <- ifelse( tests$odds > 1, "enriched", "depleted" )
    tests$high.confidence <- ifelse( tests$consistent.CI == "yes" & tests$adjusted.p < 0.05, "yes", "no" )
#     tests$p.value <- noquote(format.pval(tests$p.value, digits = 3)) 
#     tests$p.hyper <- noquote(format.pval(tests$p.hyper, digits = 3)) 
#     tests$adjusted.p <- noquote(format.pval(tests$adjusted.p, digits = 3)) 
    return(arrange(tests, plyr::desc(diagnostic), type, adjusted.p))
  } else {return(NULL)}
}

