setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)

source("scripts/commands.R")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- COMMANDS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#


get.uco  <- function(biostring.seq,method="freq"){ uco(s2c(toString(biostring.seq)),index=method) }


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
# ---- ANTICODON DEMAND AND SUPPLY ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# ----- Prepare data ----
# Build a table giving the average number of freely available tRNA molecule of each species
  tRNA_molecule.copy.number_t020   <- compile.tRNA.copy.numbers(path="data/SMoPT/stochastic.simulation/unformated/20 minutes/")
  tRNA_molecule.copy.number_t120   <- compile.tRNA.copy.numbers(path="data/SMoPT/stochastic.simulation/unformated/120 minutes/")
  tRNA_molecule.copy.number_normal <- compile.tRNA.copy.numbers(path="data/SMoPT/stochastic.simulation/unformated/genome_normal_conditions/")
  anticodon.supply <- arrange( Reduce( rbind, list( tRNA_molecule.copy.number_t020, tRNA_molecule.copy.number_t120, tRNA_molecule.copy.number_normal)), anticodon, label, experiment)
  colnames(anticodon.supply)[3] <- "supply"

# Correspondance between codons and anticodons
  rosetta.stone <- unique(TAB$genetic.code[,-3])
  
# Build a table giving the mRNA molecule copy number (estimated from rounded 2^(log2) fold changes from Gasch2000) for each gene in the 4 stresses conditions
  mRNA_abundance <- read.table("data/SMoPT/stochastic.simulation/mRNA_abundance.txt", header=1, stringsAsFactors=F)
  log2.mRNA_abundance <- read.table("data/SMoPT/stochastic.simulation/mRNA_abundance_change.txt", header=1, stringsAsFactors=F)
  log2.mRNA_abundance$genome_normal_conditions <- NA  
  genes <- read.table("~/Documents/MRC16/data/SMoPT/example/input/S.cer.mRNA.abndc.ini.tsv", sep="\t", header=1,stringsAsFactors=F)
  genes$Gene <- 0:(nrow(genes)-1)

  mRNA_molecule.copy.number <- melt( merge( genes[,c("Gene","ORF","IniProb")], mRNA_abundance, by="Gene" ), id.vars = c("Gene","ORF","IniProb") )[,-1]
  colnames(mRNA_molecule.copy.number)[3:4] <- c("experiment","mRNA_abundance")
  mRNA_molecule.copy.number <- rbind( mRNA_molecule.copy.number, data.frame( ORF = genes$ORF, experiment = "genome_normal_conditions", IniProb = genes$IniProb, mRNA_abundance = genes$rand_mRNA) )


  mRNA_log2.change <- melt( log2.mRNA_abundance, id.vars = "ORF" )
  colnames(mRNA_log2.change)[2:3] <- c("experiment","log2.mRNA_abundance")

  mRNA_molecule.copy.number <-  merge( mRNA_molecule.copy.number, mRNA_log2.change, by=c("ORF","experiment") )

  write.table(mRNA_molecule.copy.number, file="data/mRNA_molecule.copy.number.txt",sep="\t", quote=F, row.names=F)

# Compute the number of codon of each sort in every gene of the SMoPT analysis
  codon.count <-  ldply( setNames(unique(mRNA_molecule.copy.number$ORF),unique(mRNA_molecule.copy.number$ORF)), function(x){ get.uco(SEQ$ORF_CDS[[x]], method="eff") } )
  colnames(codon.count)[1] <- "ORF"

# Compute the codon demand per gene per condition  
  codon.count.mRNAab <- merge( mRNA_molecule.copy.number, codon.count, by="ORF" )
  codon.count.mRNAab[,colnames(codon.count)[-1]]  <- codon.count.mRNAab$mRNA_abundance * as.matrix(codon.count.mRNAab[,colnames(codon.count)[-1]] )

# Compute the codon demand per condition
  codon.count.mRNAab.long <- melt(codon.count.mRNAab, id.vars= c("ORF","IniProb","experiment","mRNA_abundance","log2.mRNA_abundance") )
  colnames(codon.count.mRNAab.long)[6:7] <- c("codon","occurrence")

  codon.count.mRNAab.long2 <- ddply( subset(codon.count.mRNAab.long,!is.na(log2.mRNA_abundance)), .(experiment), function(x){ x$quant.log2mRNA <- quantile.it(x$log2.mRNA_abundance); 
                                                                                                                             x$up_reg <- ifelse(x$quant.log2mRNA == "0-20", T,F); 
                                                                                                                             return(x) })
  
  codon.count.mRNAab.long2 <- rbind(  codon.count.mRNAab.long2, 
                                      data.frame( subset(codon.count.mRNAab.long,is.na(log2.mRNA_abundance)), quant.log2mRNA = NA, up_reg = NA)
                                  )

# Sum up the number of expressed codons per condition and codon type to obtain the anticodon demand
  expressed.codons <- ddply(codon.count.mRNAab.long2, .(codon, experiment), function(x){ sum(x$occurrence) } )
  colnames(expressed.codons)[3] <- "total.number"

  expressed.codons.up_reg <- ddply( subset(codon.count.mRNAab.long2,up_reg==T), .(codon, experiment), function(x){ sum(x$occurrence) } )
  colnames(expressed.codons.up_reg)[3] <- "total.number"

# the same anticodon can recognise several codons: need to sum up these values
  anticodon.demand.codons  <- merge(expressed.codons, rosetta.stone, by="codon", all.x=TRUE)
  anticodon.demand.codons.up_reg  <- merge(expressed.codons.up_reg, rosetta.stone, by="codon", all.x=TRUE)
  
  anticodon.demand        <- ddply( anticodon.demand.codons, .(anticodon, experiment), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), demand = sum(x$total.number) ) } ) 
  anticodon.demand.up_reg <- ddply( anticodon.demand.codons.up_reg, .(anticodon, experiment), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), demand = sum(x$total.number) ) } ) 

# merge data on anticodon demand and supply
  anticodon.economy <- merge( anticodon.demand, anticodon.supply[,-1], by=c("anticodon","experiment") )
  anticodon.economy <- merge( anticodon.economy, unique(rosetta.stone[,c("anticodon","aa")]), by="anticodon", all.x=T )

  anticodon.economy.up_reg <- merge( anticodon.demand.up_reg, anticodon.supply[,-1], by=c("anticodon","experiment") )
  anticodon.economy.up_reg <- merge( anticodon.economy.up_reg, unique(rosetta.stone[,c("anticodon","aa")]), by="anticodon", all.x=T )


data <- subset(anticodon.economy, experiment %in% c("genome_ox_20","genome_ox_120","genome_osm_20rep") )
data <- anticodon.economy.up_reg
data <- anticodon.economy


# arrange( subset(rosetta.stone, aa == "Val"), anticodon, codon )

ggplot( data = data, aes( x = demand, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=4,alpha=0.4, aes(color = copy.number) ) + 
  geom_text(aes(label=recognised), size=3) + 
  stat_smooth( method="lm" ) + 
  scale_x_sqrt() + scale_y_sqrt() +
  theme_bw() + coord_fixed() +
  facet_wrap( ~ experiment, scales = "free" )









#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- FRACTION OF TRANSCRIPTOME ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# percentage.transcriptome <- ddply( mRNA_molecule.copy.number, .(experiment), function(x){ x$percentage <-  100*x$mRNA_abundance / sum(x$mRNA_abundance); return(x) })
 
percentage.transcriptome.020 <- ddply( master.table.020, .(experiment), function(x){ x$percentage <-  100*x$est.mRNA_abundance / sum(x$est.mRNA_abundance); return(x) })
percentage.transcriptome.120 <- ddply( master.table.120, .(experiment), function(x){ x$percentage <-  100*x$est.mRNA_abundance / sum(x$est.mRNA_abundance); return(x) })

percentage.transcriptome.020 <- ddply( percentage.transcriptome.020, .(experiment), function(x){x$quant.log2mRNA <- quantile.it(x$log2.mRNA_abundance); x$up_reg <- ifelse(x$quant.log2mRNA == "0-20", T,F); return(x)} )
percentage.transcriptome.120 <- ddply( percentage.transcriptome.120, .(experiment), function(x){x$quant.log2mRNA <- quantile.it(x$log2.mRNA_abundance); x$up_reg <- ifelse(x$quant.log2mRNA == "0-20", T,F); return(x)} )


ggplot( data = percentage.transcriptome.020, aes(x=percentage,y=av.initiation_time) ) + 
  geom_point() +
  #scale_x_log10() + scale_y_log10() +
  theme_bw() + coord_fixed() +
  facet_wrap( ~ experiment, scales = "free" )


ggplot( data = percentage.transcriptome.020, aes(x=up_reg,y=tRNA_adaptation_index) ) + 
  geom_boxplot(notch=T) +
  #scale_x_log10() + scale_y_log10() +
  theme_bw() + coord_fixed() +
  facet_wrap( ~ experiment, scales = "free" )

###
  percentage.transcriptome.020.up_reg  <- subset(percentage.transcriptome.020, up_reg==T)

  percentage.transcriptome.020.up_reg$experiment

# Sum up the number of expressed codons per condition and codon type to obtain the anticodon demand
  expressed.codons <- ddply(codon.count.mRNAab.long, .(experiment, codon), function(x){ sum(x$occurrence) } )
  colnames(expressed.codons)[3] <- "total.number"

  # the same anticodon can recognise several codons: need to sum up these values
  anticodon.demand.codons  <- merge(expressed.codons, rosetta.stone, by="codon", all.x=TRUE) 

  anticodon.demand <- ddply( anticodon.demand.codons, .(experiment, anticodon), function(x){ data.frame( recognised = paste(x$codon,collapse="+"), demand = sum(x$total.number) ) } ) 

# merge data on anticodon demand and supply
  anticodon.economy <- merge( anticodon.demand, anticodon.supply[,-1], by=c("anticodon","experiment") )
  anticodon.economy <- merge( anticodon.economy, rosetta.stone[,c("anticodon","aa")], by="anticodon", all.x=T )
 
