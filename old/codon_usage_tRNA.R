setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)

source("scripts/commands.R")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- LOAD DATA ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# ---- Codon-based data ----
data(TAB)
TAB$nTE # normalised TE
TAB$CAI # Codon Adaptation Index
TAB$tAI # tRNA Adaptation Index

# ---- Transcriptomic data ----
data(TRS)
TRS$Gasch2000_response.stress

# ---- Sequence data ----
data(SEQ)
SEQ$ORF_CDS

# ---- tRNA abundance data ----
# Data generated and collected by Marc and processed here:scripts/tRNA_abundance_Autumn2013.R
load("data/tRNA abundance/tM.Rda")
tM # log2 scaled fold change of tRNA abundance between stress condition and wt normal condition

tRNAab$conditions$label <- factor(tRNAab$conditions$treatment, levels=unique(tRNAab$conditions$treatment), labels=c("peroxide","heat_shock","ethanol","sorbitol"))
tRNAab$conditions$experiment <- factor(tRNAab$conditions$treatment, levels=unique(tRNAab$conditions$treatment), labels=c("ox","temp","diauxic","osm"))
tRNAab$conditions$stress <- paste(tRNAab$conditions$experiment, tRNAab$conditions$time, sep="_")
colnames(tRNAab$abundance)[1] <- "anticodon"
save(tRNAab, file="data/tRNA abundance/tRNAab.Rda")

load("data/tRNA abundance/tRNAab.Rda")

tRNAab.foldchange <-  tRNAab$abundance[ , c("anticodon", subset(tRNAab$conditions, type == "average")$reference ) ]
colnames(tRNAab.foldchange)[-1] <- as.character(factor( colnames(tRNAab.foldchange)[-1],  levels = tRNAab$conditions$reference, labels = tRNAab$conditions$stress))
write.table(tRNAab.foldchange, "data/tRNA abundance/tRNAab_quantified.foldchange.txt", sep="\t", quote=F, row.names=F)

# tRNAab and tRNAab.foldchange are important datasets to keep!

# ---- Genetic code -----

TAB$genetic.code


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- APPROACH ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# Single example
uco(s2c(toString(SEQ$ORF_CDS[[1]])),index="freq")


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- COMMANDS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

get.uco  <- function(biostring.seq,method="freq"){ uco(s2c(toString(biostring.seq)),index=method) }


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- DATA PROCESSING ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# 1 - Compute codon usage per gene per codon (3 methods)
uco.freq <- sapply(SEQ$ORF_CDS, function(x) get.uco(x, method="freq") ) 
uco.freq <- as.data.table(data.frame( gene = names(SEQ$ORF_CDS), t(uco.freq)) )
save(uco.freq,file="data/Lister/orf_all.codon.usage_freq.Rda")

uco.eff <- sapply(SEQ$ORF_CDS, function(x) get.uco(x, method="eff") ) 
uco.eff <- as.data.table(data.frame( t(uco.eff)) )
save(uco.eff,file="data/Lister/orf_all.codon.usage_eff.Rda")

uco.rscu <- sapply(SEQ$ORF_CDS, function(x) get.uco(x, method="rscu") ) 
uco.rscu <- as.data.table(data.frame(  gene = names(SEQ$ORF_CDS), t(uco.rscu)) );
save(uco.rscu,file="data/Lister/orf_all.codon.usage_rscu.Rda")

      # RSCU is a simple measure of non-uniform usage of synonymous codons in a coding sequence (Sharp et al. 1986). 
      # RSCU values are the number of times a particular codon is observed, relative to the number of times that the 
      # codon would be observed for a uniform synonymous codon usage (i.e. all the codons for a given amino-acid have
      # the same probability). In the absence of any codon usage bias, the RSCU values would be 1.00 (this is the case
      # for sequence cds in the exemple thereafter). A codon that is used less frequently than expected will have an 
      # RSCU value of less than 1.00 and vice versa for a codon that is used more frequently than expected.
      # Do not use correspondence analysis on RSCU tables as this is a source of artifacts (Perriere and Thioulouse 2002, Suzuki et al. 2008). 
      # Within-aminoacid correspondence analysis is a simple way to study synonymous codon usage (Charif et al. 2005). 
      # For an introduction to correspondence analysis and within-aminoacid correspondence analysis see the chapter 
      # titled Multivariate analyses in the seqinR manual that ships with the seqinR package in the doc folder. 
      # You can also use internal correspondence analysis if you want to analyze simultaneously a row-block structure
      # such as the within and between species variability (Lobry and Chessel 2003).
      
# 2 - Define up-regulated transcripts
Gash_quantilized <- data.table(data.frame( gene = rownames(TRS$Gasch2000_response.stress), quantilize(data=TRS$Gasch2000_response.stress,
                                                                                                      bin.number=20, name="gash") )) 
Gash_quantilized <- data.table(data.frame( gene = rownames(TRS$Gasch2000_response.stress), quantilize(data=TRS$Gasch2000_response.stress,
                                                                                                      bin.number=5, name="gash") )) 
setnames(Gash_quantilized, colnames(Gash_quantilized)[2:ncol(Gash_quantilized)], colnames(TRS$Gasch2000_response.stress))
sum(is.na(TRS$Gasch2000_response.stress)) == sum(is.na(Gash_quantilized))  # verify that the number of NAs is exactly identical prior and post quantile-transformation

Gene_categories <- data.table(melt(Gash_quantilized,id.vars="gene"))
setnames(Gene_categories, old=c("variable","value"), c("condition","bin"))
Gene_categories$expression_change <- apply(Gene_categories, 1, function(x){ TRS$Gasch2000_response.stress[ x[[1]], x[[2]] ] })
bin.number <- 20; LABELS <- paste( round((1 -seq(1:bin.number)/bin.number)*100 ), round(100/bin.number) + round((1 -seq(1:bin.number)/bin.number)*100 ),sep="-" )

Gene_categories$bin <- factor(Gene_categories$bin, levels=LABELS)
Gene_categories <- subset(Gene_categories, !is.na(bin))

# 3 - Define selections
read.table("data/gash_conditions.txt")
selection.bins <- c("0-5","5-10","5-15","15-20","80-85","85-90","90-95","95-100")
selection.expression  <- read.table("data/Lister/selected_conditions.txt", header=T)
colnames(tM)  
  
# 4 - Subset genes
ggplot(data=subset(Gene_categories,condition %in% selection.expression$condition[2]),aes(x=bin, y=expression_change)) + geom_boxplot(notch=T) + facet_grid( condition ~ .)
ggsave(filename="results/lister/binning.pdf",dpi=300)
Gene_selected <- subset(Gene_categories, bin %in% selection.bins & condition %in% selection.expression$condition )
Gene_selected$change <- ifelse(Gene_selected$bin %in% c("0-5","5-10","5-15","15-20"), "up", "down")

# 5 - RNA Seq data
data(XPR)
XPR$mRNA.abundance
Gene_mRNAab <- subset(merge(Gene_categories, XPR$mRNA.abundance, by="gene"), !is.na(mRNA.ab_YPD.rpkm) & condition %in% selection.expression$condition & mRNA.ab_YPD.rpkm < 1000 ) # proxy for removing far outliers
ggplot(data=subset(Gene_mRNAab,condition %in% selection.expression$condition[2]),aes(x=bin, y=mRNA.ab_YPD.rpkm)) + geom_boxplot(notch=T) + facet_grid( condition ~ .)

Gene_mRNAab_exp <-  Gene_mRNAab[ ,putative_expression := 2^(expression_change)*mRNA.ab_YPD.rpkm ]

# 6 - Integrate at codon usage
Gene_eff <- data.table( melt( subset( uco.eff, gene %in% Gene_mRNAab_exp$gene ), id.vars="gene"), key="gene" )
setnames(Gene_eff, c("variable","value"),c("codon","eff"))
setkey(Gene_mRNAab_exp, key="gene")

# ! memory problem here ! 
Gene_eff_mRNA_codons <- merge( subset(Gene_eff, !codon %in% c("atg","tga", "taa", "tag","tgg")), Gene_mRNAab_exp, by="gene", allow.cartesian=TRUE )

save(Gene_eff_mRNA_codons, file="data/Lister/Gene_eff_mRNA_codons.Rda")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- PRELIMINARY ANALYSIS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

# tRNAab_data <- subset(tRNAab$abundance, select= c("tRNA_anticodon", subset(tRNAab$conditions, type=="average" & time == 20)$reference) )
# colnames(tRNAab_data)<- c("anticodon",subset(tRNAab$conditions, type=="average" & time == 20)$label)
# 
 upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)
 rosetta.stone <- merge(upstress_codons, unique(TAB$genetic.code[,c("codon","anticodon")]), by="codon")
# 
# preliminary.data <- melt( merge(tRNAab_data, rosetta.stone, by="anticodon"), id.vars=c("codon","up_in","anticodon"))
# preliminary.data$uniq_pair  <- paste(preliminary.data$codon, preliminary.data$anticodon, sep="_")
# 
# tAI <- TAB$tAI
# #tAI$codon <- sapply(tAI$codon, function(x){ tolower(toString(transcribe(DNAString(x)))) })
# tAI$codon <-  tolower(tAI$codon)
# 
# preliminary.data <- merge(preliminary.data, tAI[,1:2], by="codon")
# colnames(preliminary.data)[c(5,7)] <- c("tRNAab","tAI")
# 
# preliminary.data <- arrange(preliminary.data, tAI)
# preliminary.data$uniq_pair <- factor(preliminary.data$uniq_pair, levels= preliminary.data$uniq_pair )
# 
# ggplot(data=preliminary.data, aes(x=uniq_pair, y=log2(tRNAab))  ) + geom_bar(aes(fill=up_in), stat="identity", position="dodge") + coord_flip() + facet_grid(  . ~ variable ) + labs(fill="codon up in", y="log2 change in tRNA ab for corresponding anticodon", x="codon-anticodon pairs")


#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- GLOBAL ANALYSIS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

anticodon.demand  <- data.table( merge(subset(Gene_eff, gene %in% RNAseq$gene), rosetta.stone, by="codon"), key = c("gene", "anticodon"))
anticodon.demand  <- anticodon.demand[ , sum(eff), by= c("gene","anticodon") ]
setnames(anticodon.demand, "V1", "n.anticodon")

tRNAab_data <- subset(tRNAab$abundance, select= c("tRNA_anticodon", subset(tRNAab$conditions, type=="average" )$reference) )
colnames(tRNAab_data)<- c("anticodon",subset(tRNAab$conditions, type=="average" )$label)

# Merge the number of anticodons / gene + tRNA abundance per condition (not necessary yet, see below)
#merge( anticodon.demand, tRNAab_data, by="anticodon") 
anticodon.demand
# Select up-regulated genes for the 4 experiments marc did, based on consistency of belonging to 0-20 across time (maybe just make 5 bins and proceed further)
Gene_categories_sel <- merge(subset(Gene_categories, condition %in% selection.expression$condition & bin %in% c("0-5","5-10","5-15","15-20","80-85","85-90","90-95","95-100") ), selection.expression, by="condition")
Gene_categories_sel$diff_exp <- ifelse(Gene_categories_sel$bin %in% c("80-85","85-90","90-95","95-100"), "down","up" )

Gene_categories_sel_codons <- merge(Gene_categories_sel, anticodon.demand, by="gene", allow.cartesian=T)

# Now have a list of up-regulated genes For each experiment as well as normal conditions (hey that is top 20% highly expressed genes? )... 
# melt it in a table with genes in columns, and the experiments in which they are up=regulated, relabel the experiments with simple terms (peroxide, heat, osmotic, diauxid_shift)
# Relabel experiments in Marc's data with simple temrs (peroxide, heat, osmotic, diauxid_shift) that match the expression data ones.

# 1 - tRNA supply
tRNAab_global <- subset(tRNAab$abundance, select= c("tRNA_anticodon", subset(tRNAab$conditions, type=="average")$reference) )
global_tRNA_supply <- melt(tRNAab_global, id.vars=c("tRNA_anticodon") )
global_tRNA_supply$label <- factor(global_tRNA_supply$variable, levels=subset(tRNAab$conditions, type=="average")$reference, labels=subset(tRNAab$conditions, type=="average")$label )
global_tRNA_supply$time <- factor(global_tRNA_supply$variable, levels=subset(tRNAab$conditions, type=="average")$reference, labels=subset(tRNAab$conditions, type=="average")$time )
colnames(global_tRNA_supply)[c(1,3)] <- c("anticodon","change_tRNAabundance")

# 2 - anticodon demand
global_anticodon_demand <- ddply(Gene_categories_sel_codons, .(anticodon, diff_exp, label, condition, time), summarise, sum(n.anticodon),median(n.anticodon, na.rm=T),mean(n.anticodon, na.rm=T))
colnames(global_anticodon_demand)[6:8] <- c("N.anticodon", "median.anticodon","mean.anticodon")

data(NAMES)
RNAseq <- subset(XPR$mRNA.abundance, gene %in% subset(NAMES, feature.type == "ORF" )$systematic.name )
up.regulated.genes_ypd <- as.character(subset( data.frame( gene = RNAseq$gene, quantilize(data=RNAseq[,-c(1,3),drop=F],name="rna",bin.number=10) ), mRNA.ab_YPD.rpkm %in% c("0-10","10-20") )$gene)

anticodon_demand_ypd <- ddply( subset(anticodon.demand, gene %in% up.regulated.genes_ypd), .(anticodon), summarise, sum(n.anticodon), median(n.anticodon,na.rm=T), mean(n.anticodon,na.rm=T))
colnames(anticodon_demand_ypd)[2:4] <- c("N.anticodon_ypd", "median.anticodon_ypd","mean.anticodon_ypd")

global_anticodon_demand_d <-  merge(global_anticodon_demand, anticodon_demand_ypd, by="anticodon")

global_anticodon_demand_d$rel.N.anticodon <- scale(global_anticodon_demand_d$N.anticodon/global_anticodon_demand_d$N.anticodon_ypd)
#global_anticodon_demand_d$rel.median.anticodon  <- scale(global_anticodon_demand_d$median.anticodon/global_anticodon_demand_d$median.anticodon_ypd)
global_anticodon_demand_d$rel.mean.anticodon  <- scale(global_anticodon_demand_d$mean.anticodon/global_anticodon_demand_d$mean.anticodon_ypd)

# 3 - merge the datasets
data <- merge( global_tRNA_supply[,-2], global_anticodon_demand_d, by=c("anticodon","label","time") )
data <- merge( data, unique(TAB$genetic.code[,c(2,4)]), by="anticodon" )

data$log2_change_tRNAabundance <- log2(data$change_tRNAabundance)
data$exclude  <- F
data[ which(data$time == 20 & data$label == "sorbitol" & data$anticodon %in% c("CAA","UAU","GCU","UCC") ), ]$exclude <- T

ddply(subset(Gene_categories_sel, diff_exp=="up"), .(condition), nrow) # how many genes are selected?
length(up.regulated.genes_ypd) # 1,061 genes

# 4 - Thiel-Sen regression
require(zyp)
test <- ddply( subset(data, diff_exp=="up" & time == 20 & exclude ==F & rel.mean.anticodon < 3 & log2(change_tRNAabundance) < 3 ), .(label), 
               function(x) zyp.sen(x,formula= rel.mean.anticodon ~ log2_change_tRNAabundance)$coefficients)
test$time  <- 20

# 5 - Plot results
ggplot(subset(data, diff_exp=="up" & time == 20 & anticodon !="CCG" ), aes(x=log2(change_tRNAabundance), y=rel.mean.anticodon )) + 
  geom_point(pch=16,alpha=0.3) + 
  #stat_smooth(aes(group=diff_exp), method="lm") + 
  geom_abline(data=test, colour="red", aes(intercept=Intercept, slope=log2_change_tRNAabundance, data=test, group=label)) +
  geom_text(aes(label=anticodon),size=3,alpha=0.75) + 
  facet_grid( time ~ label, scale="free") + theme_bw() 

ggsave("results/lister/adjustment.pdf",dpi=300)



# tRNAab_melt <- melt(tRNAab_data, id.vars="anticodon")
# write.table(tRNAab_melt[,"variable",drop=F],file="~/Desktop/tRNA_conditions")

# merge the tables and subset so that you end up with n.anticodon / up-regulated gene / experiment
# sum the anticodons overall per experiment (across the genes)

# make the ratios n.anticodon_upregulated_condition / n.anticodon_high_expressed_normal --> anticodon demand, row = anticodons
# merge that with tRNA abundance data





#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- GLOBAL ANALYSIS ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
test <- merge(Gene_rscu, subset(Gene_selected, condition %in% selection.expression[2]), by="gene" )
RSCU_per_bin <- ddply(test, .(codon, change, bin), summarize, median(rscu,na.rm=T) )
colnames(RSCU_per_bin)[4] <- "median.RSCU"

Delta_RSCU_per_bin  <- cast(RSCU_per_bin, formula= codon ~ bin )
heatmap(as.matrix(Delta_RSCU_per_bin), scale="row")
ggplot(test, aes(x=codon, y= rscu ) ) + geom_point(aes(colour=change), stat="identity") 


