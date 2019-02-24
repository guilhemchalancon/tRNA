# 
setwd("~/Documents/MRC16/")

require(plyr)
require(Biostrings)
require(ggplot2)

##### In-House Commands ----
translation.efficiency.scoring.per.codon <- function(codon_seq_list = all_CDS, aa_seq_list= all_AA, CAI= CAI.table, tAI = tAI.table, nTE = nTE.table, optimality = optimality.table, species= "Scer"){
  genes <- names(codon_seq_list)
  require(seqinr)
  scores <- lapply( genes, function(x){ data.frame( codon = codon_seq_list[[x]], aa = s2c(as.character(aa_seq_list[[x]])), CAI = round(CAI[ codon_seq_list[[x]], species],3), tAI = round(tAI[ codon_seq_list[[x]], species],3), nTE = nTE[ codon_seq_list[[x]], species], optimality = optimality[ codon_seq_list[[x]], species] ) } )
  names(scores) <- genes
  return(scores)
} 

nTE.scoring.per.gene <- function(x) {
  require(entropy)
  data.frame( mean.nTE = round(mean(x$nTE, na.rm=T),3), 
              fraction.optimal.codons = round(length( which(x$optimality == "O") )/(nrow(x)-1), 3) ,
              median.nTE = round(median(x$nTE, na.rm=T), 3), 
              ramp.nTE = round(mean(x$nTE[1:(min(nrow(x), 45))], na.rm=T),3), 
              dip.nTE= round(mean(x$nTE[1:(min(nrow(x), 15))], na.rm=T),3),
              plateau.nTE= round(mean( x$nTE[ (min(nrow(x), 46):nrow(x))], na.rm=T),3),
              too.short = ifelse( nrow(x) < 46, T, F),
              entropy.codons = round(entropy(table(x$codon)),3),
              entropy.rle.optimality = round(entropy(table(rle(as.character(x$optimality))[[1]])),3),
              seq_length=nrow(x)
  )}

check.seq <- function(genes){ lapply( genes, function(x){  y <- all_nTE[[x]]$optimality
                                                           y <- factor(y, levels=c("N","O"),labels=c("_","O"))
                                                           return(c2s(y)) }  ) }


##### Construct data ----

# Table for calculating tRNA Adaptation Index
CAI.table <- read.table("~/Documents/StingRay/Source/data_sources/expression/CAI_table.txt",sep="\t",header=1, row.names=1)
# Table for calculating tRNA Adaptation Index
tAI.table <- read.table("~/Documents/StingRay/Source/data_sources/expression/tAI_table.txt",sep="\t",header=1, row.names=1)
# Table giving the normalised translation efficiency scores per codon and per species
nTE.table <- read.table("~/Documents/StingRay/Source/data_sources/expression/frydman_translation_efficiency/normalizedTE.txt", header=1,row.names=1)
optimality.table <- read.table("~/Documents/StingRay/Source/data_sources/expression/frydman_translation_efficiency/normalizedcodonoptimality.txt", header=1,row.names=1)


ura.constructs <- readDNAStringSet("data/GFP/Ura_constructs.fasta")
ura.AA <- readAAStringSet("data/GFP/Ura_original_translated.fasta")  

# Verify that all constructs lead to the same protein
t <- translate(ura.constructs)

# Compute normalised translation efficiency per codon, per construct
ura_codons <- lapply( ura.constructs, function(x) as.character(codons(x)) )
ura_AA <- lapply(ura.constructs, function(x) ura.AA)
ura_TE <- translation.efficiency.scoring.per.codon(codon_seq_list= ura_codons, aa_seq_list=ura_AA )
ura_TE <- lapply(ura_TE, function(x){ data.frame( position=1:nrow(x), x) })

save(ura_TE, file="results/gfp sequence optimization/ura_TE.Rda")

ura_TE.df <- ldply(ura_TE)
colnames(ura_TE.df)[1] <- "construct"

##### Analysis ----

ggplot(data=ura_TE.df, aes(x=position, y=CAI)) + geom_line(aes(colour=construct)) + scale_color_brewer(palette="Paired") + facet_grid( construct ~ .) + theme_bw()
ggplot(data=ura_TE.df, aes(x=position, y=tAI)) + geom_line(aes(colour=construct)) + scale_color_brewer(palette="Paired") + facet_grid( construct ~ .) + theme_bw()
ggplot(data=ura_TE.df, aes(x=position, y=nTE)) + geom_line(aes(colour=construct)) + scale_color_brewer(palette="Paired") + facet_grid( construct ~ .) + theme_bw()

ggplot(data=ura_TE.df, aes(x=construct, y=CAI)) + geom_violin(aes(fill=construct)) + geom_boxplot(alpha=0.5,width=0.5,notch=T) + scale_fill_brewer(palette="Paired") + theme_bw() + stat_summary(aes(group=construct), fun.y=mean, geom="point", colour="white", size=5, alpha=0.85)
ggplot(data=ura_TE.df, aes(x=construct, y=nTE)) + geom_violin(aes(fill=construct)) + geom_boxplot(alpha=0.5,width=0.5,notch=T) + scale_fill_brewer(palette="Paired") + theme_bw() + stat_summary(aes(group=construct), fun.y=mean, geom="point", colour="white", size=5, alpha=0.85)
ggplot(data=ura_TE.df, aes(x=construct, y=tAI)) + geom_violin(aes(fill=construct)) + geom_boxplot(alpha=0.5,width=0.5,notch=T) + scale_fill_brewer(palette="Paired") + theme_bw() + stat_summary(aes(group=construct), fun.y=mean, geom="point", colour="white", size=5, alpha=0.85)

sapply(ura_TE, function(x){nTE.scoring.per.gene(x)});

