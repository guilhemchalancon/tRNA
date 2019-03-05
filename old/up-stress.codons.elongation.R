setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(plyr)
require(reshape)
require(RColorBrewer)
require(pheatmap)
source("scripts/commands.R")


# >> correspondence between codon, anticodons and amino acids
rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1,sep="\t")


upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)

df.upstress_codons <- merge(codon.master.table, upstress_codons, by= "codon")


df.upstress_codons <- ddply(df.upstress_codons, .(experiment, time), mutate, 
                            zscore.change.demand = zscore(demand.mRNAab.up/demand.mRNAab.up.normal),
                            zscore.change.supply = zscore(total.tRNAab/total.tRNAab.normal),
                            log2.change.demand = log2(demand.mRNAab.up/demand.mRNAab.up.normal),
                            log2.change.supply = log2(total.tRNAab/total.tRNAab.normal)
                            ) 

head(df.upstress_codons)

dlply(df.upstress_codons, .(experiment, time), summarise, summary(zscore.change.demand) )

ggplot(data= subset(df.upstress_codons, time>0), aes( x = factor(time), y = as.character(codon))) + 
  geom_tile(aes(fill= zscore.change.demand)) +
  facet_grid( up_in ~ experiment, drop = T, scales = "free_y") + scale_fill_gradient2(low = "skyblue", mid="black", high="yellow")

ggplot(data= subset(df.upstress_codons, time>0), aes( x = factor(time), y = as.character(codon))) + 
  geom_tile(aes(fill= zscore.change.supply)) +
  facet_grid( up_in ~ experiment, drop = T, scales = "free_y") + scale_fill_gradient2(low = "skyblue", mid="black", high="yellow")


m <- matrixify(reshape2::dcast( data = subset(df.upstress_codons, time > 0 ), value.var="zscore.change.demand", formula = codon ~ experiment + time))
pheatmap(m, cellwidth=12, cellheight=12, clustering_method = "centroid", cluster_cols = F, treeheight_row=20)



# Feb 2015: check again the top-regulated genes in each of the conditions
# CODON.USAGE.STRESS: table that gives RSCU, Frequency, mean and total occurrence of each codon, crossed by STRESS x TRANSCRIPTION BEHAVIOUR (based on log2 fold change in mRNA abundance) 
# Distinguish between quant.mRNAab (which is about expression levels) and quant.log2mRNAab (which is about changes in expression levels)

# Visualise quantile distributions
pdf(paste0(PATH,"part3/codon.usage_quantiles.pdf"), width = 15, height = 7)
ggplot(subset(ddply(subset(gene.master.table, time>0), .(experiment, time), function(x){x$quant.log2mRNAab <- quantile.it(x$log2.mRNA_abundance, N = 10); return(x) }),!is.na(quant.log2mRNAab)), aes(x=quant.log2mRNAab, y=log2.mRNA_abundance)) + 
  coord_flip() +
  geom_hline(yintercept=0, lty=3) +
  geom_jitter(alpha=0.25, size=1) +
  geom_violin() +
  #geom_boxplot(notch=F) + 
  labs(y="fold change in mRNA abundance (log2)", x="quantiles") +
  facet_wrap(  ~ experiment + time, nrow=2)
dev.off()


pdf(paste0(PATH,"part3/codon.usage_quantiles.ox20.pdf"), width = 5, height = 5)
ggplot(subset(ddply(subset(gene.master.table, time == 20 & experiment == "ox"), .(experiment, time), function(x){x$quant.log2mRNAab <- quantile.it(x$log2.mRNA_abundance, N = 10); return(x) }),!is.na(quant.log2mRNAab)), aes(x=quant.log2mRNAab, y=log2.mRNA_abundance)) + 
  coord_flip() +
  geom_hline(yintercept=0, lty=3) +
  geom_jitter(alpha=0.25, size=1) +
  geom_violin() +
  #geom_boxplot(notch=F) + 
  labs(y="fold change in mRNA abundance (log2)", x="quantiles") +
  facet_wrap(  ~ experiment + time, nrow=2)
dev.off()
 

load("data/Lister/orf_all.codon.usage_eff.Rda")
load("data/Lister/orf_all.codon.usage_rscu.Rda")
load("data/Lister/orf_all.codon.usage_freq.Rda")

RSCU.transcriptome <- ddply( arrange(melt(uco.rscu), gene), .(variable), summarise, mean.RSCU = mean(value,na.rm=T), SD.RSCU = sd(value,na.rm=T) )
colnames(RSCU.transcriptome)[1] <- c("codon")

# Table that merges the gene.master.table with quantile information about the log2 fold change (could be integrated directly into gene.master.table, if time permits)
transcription.response.table <- ddply(subset(gene.master.table, time>0), .(experiment, time), function(x){x$quant.log2mRNAab <- quantile.it(x$log2.mRNA_abundance, N = 10); return(x) })
transcription.response.table$transcription.behaviour  <- ifelse( transcription.response.table$quant.log2mRNAab %in% c("0-10","10-20"), "up", 
                                                                 ifelse( transcription.response.table$quant.log2mRNAab %in% c("90-100","80-90"), "down", "stable" )  )


# How similar are the sets of up and down-regulated genes?
gene.sets <- dlply( transcription.response.table, .(experiment, time, transcription.behaviour), function(x) as.character(x$ORF) )

gene.sets.similarity <- data.frame(t(combn( names(gene.sets), 2)))
  colnames(gene.sets.similarity)[1:2] <- c("id1","id2")
  gene.sets.similarity$group1  <- gsub(gene.sets.similarity$id1, pattern="^(.*)\\.(.*)\\.(.*)$", replacement="\\3")
  gene.sets.similarity$group2  <- gsub(gene.sets.similarity$id2, pattern="^(.*)\\.(.*)\\.(.*)$", replacement="\\3")
  gene.sets.similarity$Jaccard <- apply(gene.sets.similarity, 1, function(x){ 
  require(sets)
  set_similarity(x= as.set(gene.sets[[ x[[1]] ]]), y= as.set(gene.sets[[ x[[2]] ]]) 
)
})

# Similarity of up and down regulated genes
# transform into a matrix for visualisation with Pheatmap
m <- get.sym.matrix(df=gene.sets.similarity, value.var="Jaccard")
diag(m) <- 0

# Within group similarity
within.group.Jaccard <- subset(gene.sets.similarity, as.character(id1) != as.character(id2) & as.character(group1) == as.character(group2))
ddply(within.group.Jaccard, .(group1), summarise, median(Jaccard) )


require(pheatmap)
require(RColorBrewer)
require(reshape2)
pdf(paste0(PATH,"part3/similarity_up_down.genes.pdf"), width = 11, height = 11)
pheatmap( m, cellwidth = 12, cellheight= 12, 
          color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100))[1:80],
          cluster_rows= T, cluster_cols=T, treeheight_row=15, treeheight_col=15,
          main = "Similarity of transcriptional responses: Jaccard index analysis")
dev.off()


tmp <- function(x) { 
  mm <- colMeans(x, na.rm=T)
  ss = sapply(x, sd, na.rm=T)
  list(names=names(x), mean=mm, sd=ss)
}


# Compute codon usage metrics for (up, down or stable) genes across stress conditions and time points
codon.usage.stress <- ddply(subset(transcription.response.table, !is.na(log2.mRNA_abundance) ), 
      .(experiment, time, transcription.behaviour), function(x){ 
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
        return(output ) }))
setnames(codon.usage.stress, old = c("mean","sd","i.mean","i.sd","i.mean.1","i.sd.1","names"), new=
           c(paste(rep(c("rscu","freq","eff"),each=2),rep(c("mean","sd"),times = 2), sep="."),"codon")
# Add info on clusters determined with the secondary analysis of the transcriptomics data from Gasch et al. 2000 
codon.usage.stress <- data.table(merge(codon.usage.stress, upstress_codons, by="codon"))


# Compute Z-scores per codon [options] per group (down, stable or up)
codon.usage.stress[, scaled.rscu := zscore(rscu.mean) , by=list(codon)]
codon.usage.stress[, scaled.freq := zscore(freq.mean) , by=list(codon)]

stables <- subset(codon.usage.stress, transcription.behaviour == "stable", select=c(codon, experiment, time,freq.mean))
setnames( stables, "freq.mean", "freq.mean.stable")
setkey(stables, codon, experiment, time)
setkey(codon.usage.stress, codon, experiment, time, transcription.behaviour)
codon.usage.stress <- merge( codon.usage.stress, stables, by= c("experiment","time","codon"))
codon.usage.stress[, rel.freq := zscore(freq.mean/freq.mean.stable) ]


write.table(codon.usage.stress, file="results/codon.usage.stress.txt", sep="\t", quote=F,row.names=F)

# --------------------------------------------------------------------------------------------------- #
codon.usage.stress <- fread("results/codon.usage.stress.txt")


# Relative RSCU (scaled by codon across conditions)
RSCU.heatmaps <- list(
  up = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "up" ), 
                                 experiment + time + transcription.behaviour ~ codon, value.var="scaled.rscu"), 
                 ref.n = 1:3),
  down =  matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "down" ), 
                                    experiment + time + transcription.behaviour ~ codon, value.var="scaled.rscu"), 
                    ref.n = 1:3),
  stable = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "stable" ), 
                                    experiment + time + transcription.behaviour ~ codon, value.var="scaled.rscu"), 
                    ref.n = 1:3) ,
  all = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") ), 
                                  experiment + time + transcription.behaviour ~ codon, value.var="scaled.rscu"), 
                  ref.n = 1:3),
  updown = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga")  & transcription.behaviour %in% c("up", "down") ), 
                                  experiment + time + transcription.behaviour ~ codon, value.var="scaled.rscu"), 
                  ref.n = 1:3)
  )

clusters.RSCU.stable <- heatmap(RSCU.heatmaps$stable)
codons.stable <- colnames(RSCU.heatmaps$stable)[clusters.RSCU.stable$colInd]

pdf(paste0(PATH,"part3/codon.usage_RSCU.scaled.pdf"), width = 15, height = 7)
lapply( names(RSCU.heatmaps), function(x){
  gulmap(RSCU.heatmaps[[x]][,codons.stable], cellwidth=12, cellheight=12, 
         cluster_cols = F,
         treeheight_row=10, treeheight_col=10, main = x, 
         annotation= data.frame(matrixify(rosetta.codons[,c("codon","CAI.O")])), #"up_in"
         annotation_colors = list(    'CAI.O' = c(O="orange",N="skyblue") ) #,
                                      #'up_in' = c(mild="green", severe="red")
         #)
  )  
})
dev.off()

FREQ.heatmaps <- list(
  up = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "up" ), 
                                 experiment + time + transcription.behaviour ~ codon, value.var="scaled.freq"), 
                 ref.n = 1:3),
  down =  matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "down" ), 
                                    experiment + time + transcription.behaviour ~ codon, value.var="scaled.freq"), 
                    ref.n = 1:3),
  stable = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour %in% "stable" ), 
                                     experiment + time + transcription.behaviour ~ codon, value.var="scaled.freq"), 
                     ref.n = 1:3) ,
  all = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") ), 
                                  experiment + time + transcription.behaviour ~ codon, value.var="scaled.freq"), 
                  ref.n = 1:3),
  updown = matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga")  & transcription.behaviour %in% c("up", "down") ), 
                                     experiment + time + transcription.behaviour ~ codon, value.var="scaled.freq"), 
                     ref.n = 1:3)
)

pdf(paste0(PATH,"part3/codon.usage_FREQ.scaled.pdf"), width = 15, height = 7)
lapply( FREQ.heatmaps, function(x){
  gulmap(x, cellwidth=12, cellheight=12, treeheight_row=10, treeheight_col=10, 
         clustering_method = "ward",
         annotation= data.frame(matrixify(rosetta.codons[,c("codon","up_in","CAI.O")])),
         annotation_colors = list(    'up_in' =     c(mild="skyblue",severe="orange")  
         )
  )  
})
dev.off()

codon.usage.stress <- read.table("results/codon.usage.stress.txt",header=1)
head(codon.usage.stress)
# don't forget that plot
ggplot(data= subset(codon.usage.stress), aes( x = factor(transcription.behaviour), y = as.character(codon))) + 
  geom_tile(aes(fill= scaled.freq ), width=0.8) + coord_fixed() +
  facet_grid( up_in ~ experiment + time, drop = T, scales = "free_y") + 
  scale_fill_gradient2(low = "skyblue", mid="black", high="orange")

# Side note on STOP codon usage
rel.use.STOP <- ddply( subset(codon.usage.stress, codon %in% c("tag","taa","tga")), .(experiment,time,transcription.behaviour), 
       summarise, codon = codon, rel.use = (freq.mean)/sum(freq.mean), freq = freq.mean )

rel.use.STOP <- ddply(rel.use.STOP, .(codon), transform, scaled.use = scale(rel.use), scaled.use.2 = scale(freq) )

m.STOP <- matrixify(reshape2::dcast( rel.use.STOP, experiment + time + transcription.behaviour ~ codon, value.var="freq"), ref.n = 1:3)

require(pheatmap)
require(RColorBrewer)
pdf(paste0(PATH,"part3/STOP_codons.pdf"), width = 11, height = 11)
pheatmap(m.STOP, cellwidth=12, cellheight=12, treeheight_row=10, treeheight_col=10, 
         clustering_method = "mcquitty",
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100))[1:80]
       )
dev.off()


# Total demand (not accounting for abundance)
ggplot(data=codon.usage.stress, aes(x= time, y=freq.mean, fill = transcription.behaviour)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid( experiment ~ codon)  

ggplot(data=codon.usage.stress, aes(x= codon, y=freq.mean, fill = up_in )) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid( transcription.behaviour ~ time + experiment)  



m.freq <- matrixify(reshape2::dcast(subset(codon.usage.stress, !codon %in% c("tag","taa","tga") & transcription.behaviour !="stable" ), 
                          experiment + time + transcription.behaviour ~ codon, value.var="rel.freq"), 
          ref.n = 1:3)

pdf(paste0(PATH,"part3/codon.usage_scaled.freq.pdf"), width = 15, height = 7)
gulmap(m.freq, centeredOn = 1, cellwidth=12, cellheight=12, 
       treeheight_row=10, treeheight_col=10, 
       annotation= data.frame(matrixify(rosetta.codons[,c("codon","up_in")])),
       annotation_colors = list(    'up_in' =     c(mild="skyblue",severe="orange") )
)
dev.off()