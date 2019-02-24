#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#     OBJECTIVES    #--------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# Here I try to identify the ESR from Gasch et al. 2000
PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"

require(data.table)
require(reshape)
require(reshape2)
require(ggplot2)
require(plyr)
require(seqinr)
require(Biostrings)
require(XLConnect)
require(pheatmap)
require(StingRay)
require(grDevices) # for the colorRampPalette
require(scales)
data(SEQ);
setwd("~/Documents/MRC16")
source("scripts/commands.R")
source("scripts/SMoT.R")

gene.master.table       <- read.table("results/master tables/gene.master.table.txt", header=1) # especially s-tAI scores and translation dynamics data from stochastic simulation

require(StingRay)
data(TRS)
Gasch.data <- as.matrix(TRS$Gasch2000_response.stress)
conditions <- read.table("~/Documents/MRC16/data/Gash_and_Chen/gash_conditions.txt")

Gasch.fig1 <- read.csv("~/Documents/MRC16/data/Gash_and_Chen/gasch.fig1_all.txt",header=1, sep="\t", comment.char = "#")
Gasch.fig3 <- read.csv("~/Documents/MRC16/data/Gash_and_Chen/gasch.fig3_ESR.txt",header=1, sep="\t", comment.char = "#")

# Remove NAs
Gasch.fig1[,2:ncol(Gasch.fig1)] <- sapply(Gasch.fig1[,-1], function(x){ x[which(is.na(x))] <- 0; x }) 
Gasch.fig3[,2:ncol(Gasch.fig3)] <- sapply(Gasch.fig3[,-1], function(x){ x[which(is.na(x))] <- 0; x }) 


# Check the data file matches - Yes they do!
# merge(Gasch.data[ , "heat_shock_33to37C_20min", drop=F], Gasch.fig1[, c("heat.shock.33.to.37..20.minutes", "name")], by.x=0, by.y="name")
# merge(Gasch.data[ , "1M_sorbitol_60_min", drop=F], Gasch.fig1[, c("X1M.sorbitol...60.min", "name")], by.x=0, by.y="name")


#--------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------#
# Retrieve the Environmental Response Cluster # remove conditions like YPD and steady state which aren't environmental stresses per se
pdf(file =  paste0(PATH,"part1/ESR.clusters.pdf"), width=10, height=10 )
h.ESR <- gulmap(matrixify(Gasch.fig3[, - c(grep(pattern = "steady", x = colnames(Gasch.fig3)), grep(pattern = "YP", x = colnames(Gasch.fig3))) ]), 
                cluster_cols = F, treeheight_row = 15, fontsize_row = 0, cellwidth=6, cellheight=0.5) # hclust(dist(matrixify(Gasch.fig3)))
dev.off()
table(clusters)

clusters <- cutree( h.ESR$tree_row, k =2 )
ESR.clusters   <- data.frame( ORF = names(clusters), cluster = factor(clusters, levels=c(1,2), labels=c("ESR.down","ESR.up")))
# gene.master.table or SMoPT
ESR.table <-      merge(subset(SMoPT.data, time > 0 ), subset(ESR.clusters,ORF %in% gene.master.table$ORF), by="ORF",  
                        all.x=T)
ESR.table$time <- factor(ESR.table$time, levels = c(20,120), labels = c("20 min", "120 min"))

# stAI.additions <- subset(SMoPT.data, time != "0") 
# stAI.additions$time <- factor(stAI.additions$time, levels = c(20,120), labels = c("20 min", "120 min"))
# ESR.table <- merge(ESR.table, stAI.additions , by = c("name","experiment","time"), all.y=T)

ESR.table$cluster <- as.character(ESR.table$cluster)
ESR.table$cluster[is.na(ESR.table$cluster)] <- "other"

ESR.table$label <- paste(ESR.table$experiment,ESR.table$time,sep="\n")
ESR.table$label <- factor(ESR.table$label, 
                           levels= c(paste( rep(c("diauxic","ox","osm","temp"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min"),
                           labels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min")
)

save(ESR.table, file = "results/master tables/ESR.table.Rd")

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------#

load("results/master tables/ESR.table.Rd")
dim(ESR.table)
head(ESR.table)


# Check variations in mRNA abundance in each condition of our study (box plots are consistent with what we see on the pheatmap just above)
g <- ggplot(ESR.table, aes(x=cluster , y=log2.mRNA_abundance )) + geom_boxplot(aes(fill=cluster), notch=T) + 
     geom_hline(yintercept=0,lty=3) +
     scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
     labs(y="mRNA abundance fold change (log2)", x="") +
     facet_wrap( ~ label, nrow=2) + theme(legend.position="none")

ggsave(plot=g, filename = paste0(PATH,"part1/ESR.log2.mRNAab.pdf"), width = 4.55, height=11.4, useDingbats=F )


# Check variations in tAI in each condition of our study (box plots are consistent with what we see on the pheatmap just above)
g <- ggplot(ESR.table, aes(x=cluster , y=stAI )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="s-tAI (unscaled)", x="") +
  facet_wrap( ~ label, nrow=2) + theme(legend.position="none")

ggsave(plot=g, filename = paste0(PATH,"part1/ESR.stAI.pdf"), width = 4.55, height=11.4, useDingbats=F )

ggplot(ESR.table, aes(x=cluster , y=FS.tAI )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="s-tAI (scaled)", x="") +
  facet_wrap( ~ label, nrow=2) + theme(legend.position="none")


ggplot(ESR.table, aes(x=cluster , y=gain.stAI )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="s-tAI (scaled)", x="") + 
  #  ylim(0,2) +
  facet_wrap( ~ label, nrow=2, scales = "free_y") + theme(legend.position="none")


ggplot(ESR.table, aes(x=cluster , y=tAI )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="s-tAI (scaled)", x="") + 
#  ylim(0,2) +
  facet_wrap( ~ label, nrow=2) + theme(legend.position="none")


ggplot(ESR.table, aes(x=tAI , y=FS.stAI )) + geom_point(pch=21,aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="s-tAI (scaled)", x="") +
  facet_wrap( ~ label, nrow=2) + theme(legend.position="none")

# gule
wc.ESR.tAI <- ddply(ESR.table, .(experiment,time), function(y) wilcox.test.batch(x = y, grouping.vars = "cluster", value = "stAI"))
wc.ESR.gain.tAI <- ddply(ESR.table, .(experiment,time), function(y) wilcox.test.batch(x = y, grouping.vars = "cluster", value = "gain.stAI"))

subset(wc.ESR.tAI, X1 == "other" & X2 == "ESR.up")
subset(wc.ESR.tAI, X1 == "ESR.down" & X2 == "ESR.up")

subset(wc.ESR.gain.tAI, X1 == "other" & X2 == "ESR.up")
subset(wc.ESR.gain.tAI, X1 == "ESR.down" & X2 == "ESR.up")


# Check variations in translation rate in each condition of our study (box plots are consistent with what we see on the pheatmap just above)
g <-  ggplot(ESR.table, aes(x=cluster , y= log2(ratio_events) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="translation rate (log2)", x="") +
  facet_wrap( ~ label, nrow=2) + theme(legend.position="none")

  ggplot(ESR.table, aes(x=cluster , y= n.events )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_y_log10() +
  scale_fill_manual(values=c("#D6604D","#4393C3")) + 
  #labs(y="translation events", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

  ggplot(ESR.table, aes(x=cluster , y= ratio_events * 2 ^( - log2.mRNA_abundance ) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_y_log10() +
  scale_fill_manual(values=c("#D6604D","#4393C3")) + 
  #labs(y="translation events", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

  ggplot(ESR.table, aes(x=cluster , y= n.events/est.mRNA_abundance )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  scale_y_log10() +
  scale_fill_manual(values=c("#D6604D","#4393C3")) + 
  #labs(y="translation events", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

cbind(ESR.table$est.mRNA_abundance * 2^(-ESR.table$log2.mRNA_abundance),
ESR.table$est.mRNA_abundance / 2^(ESR.table$log2.mRNA_abundance))

summary(ESR.table$est.mRNA_abundance / 2^(ESR.table$log2.mRNA_abundance))
summary(ESR.table$est.mRNA_abundance)


  ggplot(subset(ESR.table, !is.na(log2.mRNA_abundance)), 
         aes(x= n.events.normal/(est.mRNA_abundance * 2^(-log2.mRNA_abundance) ) , y= n.events / est.mRNA_abundance )) + 
  geom_point(aes(fill=cluster), notch=T, pch=21) + 
  geom_abline(slope=1, intercept=0, lty=3) +
  #scale_y_log10() +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="translation/mRNA (stress)", x="translation/mRNA (non-stress)") +
  facet_grid( experiment ~ time) + theme(legend.position="none")


  ggplot(subset(ESR.table, !is.na(log2.mRNA_abundance)), 
       aes(x = AUGCAI , y= av.initiation_time )) + 
  geom_point(data = subset(ESR.table, !is.na(log2.mRNA_abundance) & cluster == "other"), alpha=0.2,fill="gray", pch=21) +
  geom_point(data = subset(ESR.table, !is.na(log2.mRNA_abundance) & cluster != "other"), aes(color=cluster), size =2.6) + 
  geom_abline(slope=1, intercept=0, lty=3) +
  #scale_y_log10() +
  scale_color_manual(values=c("#D6604D","#4393C3")) + 
  labs(y="translation/mRNA (stress)", x="translation/mRNA (non-stress)") +
  facet_grid( experiment ~ time) + theme(legend.position="top")




g <- ggplot(ESR.table, aes(x=cluster , y= log2(ratio_total.time) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="change in total time of translation (log2)", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

ggplot(ESR.table, aes(x=cluster , y= (delta_elongation)/n.codon )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3")) + 
  labs(y=expression(Delta[e]), x="") + ylim(c(-0.2,0.2)) +
  facet_grid( experiment ~ time) + theme(legend.position="none")

ggplot(ESR.table, aes(x=cluster , y= log2(1/ratio_initiation) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y=expression(Delta[i]), x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")



ggplot(ESR.table, aes(x=cluster , y= (elongation.speed) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  #geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="elongation speed", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

ggplot(ESR.table, aes(x=cluster , y= (n.codon) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  #geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="length", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

ggplot(ESR.table, aes(x=cluster , y= (n.codon/av.elongation_time) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  #geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="elongation speed", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

ggplot(ESR.table, aes(x=cluster , y= (av.initiation_time) )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  #geom_hline(yintercept=0,lty=3) +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y="waiting time initiation", x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")


# Check variations in global protein production rates
g <- ggplot(ESR.table, aes(x=cluster , y=global.protein_synthesis.rate )) + geom_boxplot(aes(fill=cluster), notch=T) + 
  geom_hline(yintercept=1,lty=3) +
  scale_y_sqrt() +
  scale_fill_manual(values=c("#D6604D","#4393C3","gray")) + 
  labs(y=expression(Pi), x="") +
  facet_grid( experiment ~ time) + theme(legend.position="none")

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------#
# Codon usage
upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)

load("data/Lister/orf_all.codon.usage_eff.Rda")
load("data/Lister/orf_all.codon.usage_rscu.Rda")
load("data/Lister/orf_all.codon.usage_freq.Rda")

require(reshape2)
RSCU.transcriptome <- ddply( arrange(melt(uco.rscu, id.var="gene"), gene), .(variable), summarise, mean.RSCU = mean(value,na.rm=T), SD.RSCU = sd(value,na.rm=T) )
colnames(RSCU.transcriptome)[1] <- c("codon")
FREQ.transcriptome <- ddply( arrange(melt(uco.freq, id.var="gene"), gene), .(variable), summarise, mean.FREQ = mean(value,na.rm=T), SD.FREQ = sd(value,na.rm=T) )
colnames(FREQ.transcriptome)[1] <- c("codon")
FREQ.transcriptome <- data.table(FREQ.transcriptome, key = "codon")

tmp <- function(x) { 
  mm <- colMeans(x, na.rm=T)
  ss = sapply(x, sd, na.rm=T)
  list(names=names(x), mean=mm, sd=ss)
}


codon.usage.ESR <- ddply(ESR.clusters,
                            .(cluster), function(x){ 
                              genes <- as.character(x$ORF)
                              # Relative Synonymous Codon Usage
                              rscu <- uco.rscu[ gene %in% genes,] ;
                              rscu.stats <- rscu[,-1, with=FALSE][, tmp(.SD)]
                              # Effective counts
                              eff <- uco.eff[ gene %in% genes,] ;
                              eff.stats <- eff[,-1, with=FALSE][, tmp(.SD)]
                              # Codon Frequency
                              freq <- uco.freq[ gene %in% genes,] ;
                              freq.stats <- freq[,-1, with=FALSE][, tmp(.SD)]
                              setkey(rscu.stats,"names")
                              setkey(eff.stats,"names")
                              setkey(freq.stats,"names")
                              # Merge stat tables
                              output <- rscu.stats[freq.stats[eff.stats]]
                              # Add total counts per (as opposed to average counts)
                              output[, eff.sum := colSums(eff[,-1,with=FALSE])]
                              return(output) })
setnames(codon.usage.ESR, old = c("mean","sd","i.mean","i.sd","i.mean.1","i.sd.1","names"), new=
           c(paste(rep(c("rscu","freq","eff"),each=2),rep(c("mean","sd"),times = 2), sep="."),"codon"))
         # Add info on clusters determined with the secondary analysis of the transcriptomics data from Gasch et al. 2000 
         codon.usage.ESR <- data.table(merge(codon.usage.ESR, upstress_codons, by="codon"),key="codon")
         

codon.usage.ESR <- merge(codon.usage.ESR, FREQ.transcriptome, by="codon")
codon.usage.ESR[ , rel.freq := freq.mean/mean.FREQ]


#--------------------------------------------------------------------------------------------------------------------------------------------------#








         # Compute Z-scores per codon [options] per group (down, stable or up)
          codon.usage.ESR[, scaled.rscu := zscore(rscu.mean) , by=list(codon)]
          codon.usage.ESR[, scaled.freq := zscore(freq.mean) , by=list(codon)]
         
         write.table(codon.usage.ESR, file="results/codon.usage.ESR.txt", sep="\t", quote=F,row.names=F)
         
         # --------------------------------------------------------------------------------------------------- #
          codon.usage.ESR <- fread("results/codon.usage.ESR.txt")
         
         
         # Relative RSCU (scaled by codon across conditions)
         RSCU.heatmaps <- list(
           up = matrixify(reshape2::dcast(subset(codon.usage.ESR, !codon %in% c("tag","taa","tga") & cluster %in% "ESR.up" ), 
                                          cluster ~ codon, value.var="rel.freq"), 
                          ref.n = 1),
           down =  matrixify(reshape2::dcast(subset(codon.usage.ESR, !codon %in% c("tag","taa","tga") & cluster %in% "ESR.down" ), 
                                             cluster ~ codon, value.var="rel.freq"), 
                             ref.n = 1),
           all = matrixify(reshape2::dcast(subset(codon.usage.ESR, !codon %in% c("tag","taa","tga") ), 
                                           cluster ~ codon, value.var="rel.freq"), 
                           ref.n = 1)
         )
         
         clusters.RSCU.stable <- heatmap(RSCU.heatmaps$all)
         codons.up <- colnames(RSCU.heatmaps$up)[clusters.RSCU.stable$colInd]
         

pdf(paste0(PATH,"part1/ESR.codon.usage_RSCU.scaled.pdf"), width = 15, height = 7)
lapply( names(RSCU.heatmaps), function(x){
  gulmap(RSCU.heatmaps[[x]][,codons.up], cellwidth=12, cellheight=12, 
         cluster_cols = F, centeredOn = 1,
         treeheight_row=10, treeheight_col=10, main = x, 
         annotation= data.frame(matrixify(rosetta.codons[,c("codon","CAI.O")])), #"up_in"
         annotation_colors = list(    'CAI.O' = c(O="orange",N="skyblue") ) #,
         #'up_in' = c(mild="green", severe="red")
         #)
  )  
})
dev.off()

# colnames(Gasch.fig1)
# Heat.Shock.015.minutes
# diauxic.shift.timecourse
# 1M.sorbitol...15.min
# constant.0.32.mM.H2O2..20.min..redo
# constant.0.32.mM.H2O2..120.min..redo
# X1M.sorbitol...120.min

ggplot(data=ESR.table, aes( x = tAI.profile, y = cluster ) ) + geom_jitter(binaxis = "x",position = "dodge" , stackdir = "center")

table( ESR.table$tAI.profile, ESR.table$cluster )



chi <- chisq.test( x = ESR.table$tAI.profile, y= ESR.table$cluster)
chi$observed
chi$expected

# renew the test focussing only on similar vs increased adaptation
d <- subset(ESR.table, tAI.profile !=c("decreased\nadaptation"))
d$set <- ifelse(d$tAI.profile == "increased\nadaptation", "increased", "not-increased" )
chi.d <- ddply( subset(d, cluster!="other"), .(experiment,time), function(x){ 
                                                  chi <- chisq.test( x = x$cluster, y= x$set)
#                                                   as.numeric(chi$observed)[4]
#                                                   as.numeric(chi$observed)[2]
#                                                   as.numeric(chi$expected)
                                                  data.frame( chi2 = round(chi$statistic, 2),
                                                              p.value = chi$p.value,
                                                              n = sum(chi$observed),
                                                              phi = round(sqrt(chi$statistic/sum(chi$observed)),2)
                                                              )
                                                  } )


#chi$observed
chi.d


ggplot(data=ESR.table, aes(x=tAI.profile)) + geom_histogram( aes(fill=cluster), position="fill" ) + 
  scale_fill_manual(values=c("red","blue","gray"))  + facet_wrap( time ~ experiment)

ggplot(data=ESR.table, aes(x=tAI.profile)) + geom_histogram( aes(fill=cluster) ) + 
  scale_fill_manual(values=c("red","blue","gray"))

ggplot(data=subset(ESR.table, cluster !="other" ), aes(x=cluster)) + geom_histogram( aes(fill=tAI.profile), position="fill" ) + 
  scale_fill_manual(values=c("red","gray","orange"))

ggplot(data=subset(ESR.table, cluster !="other" ), aes(x=cluster)) + geom_histogram( aes(fill=tAI.profile), position="fill" ) + 
  scale_fill_manual(values=c("red","gray","orange"))


ggplot(data=subset(ESR.table, cluster!="other" & tAI.profile != "decreased\nadaptation" ) , aes(x=tAI.profile)) + 
  geom_histogram( aes(fill=cluster), position="fill" ) + 
  scale_fill_manual(values=c("red","blue","gray"))


ggplot(data=subset(ESR.table, cluster!="other" & tAI.profile != "decreased\nadaptation" ) , aes(x=tAI.profile)) + 
  geom_histogram( aes(fill=cluster), position="fill" ) + 
  scale_fill_manual(values=c("red","blue","gray"))


ggplot(data=subset(ESR.table, cluster!="other" ) , aes(x=cluster)) + 
  geom_histogram( aes(fill=tAI.profile), position="fill" ) + 
  scale_fill_manual(values=c("red","gray","orange"))


# -------------------------------------------------------------------------

ggplot(data=subset(ESR.table, cluster!="blah" ) , aes(x=gain.stAI>1)) + 
  geom_histogram( aes(fill=cluster), position="fill" ) + 
  scale_fill_manual(values=c("red","orange","gray")) + facet_wrap( ~ experiment + time)

ggplot(data=subset(ESR.table, cluster!="blag" ) , aes(x=gain.stAI>1)) + 
  geom_histogram( aes(fill=cluster) ) + 
  scale_fill_manual(values=c("red","orange","gray")) + facet_wrap( ~ experiment + time)


