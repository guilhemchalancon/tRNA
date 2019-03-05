source("scripts/commands.R"); source("scripts/SMoT.R")
PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
setwd("~/Documents/Research//MRC16/") # May 2018 change in folder location
require(reshape2)
require(GGally)
require(ggplot2)
require(plyr)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
#                                                                               DATA
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #

anticodon.master.table  <- read.table("results/master tables/anticodon.master.table.txt", header=1) # especially adjusted tGCN, anticodon supply, relative availability
codon.master.table      <- read.table("results/master tables/codon.master.table.txt", header=1) # especially relative adaptiveness and anticodon deman

codon.master.table$label <- with(codon.master.table, paste(experiment, paste(time,"min"),sep="\n")) 
codon.master.table$label  <- factor(codon.master.table$label, 
         levels= c(paste( rep(c("diauxic","ox","osm","temp"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min"),
         labels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min")
  )

# gene.master.table       <- read.table("results/master tables/gene.master.table.txt", header=1) # especially s-tAI scores and translation dynamics data from stochastic simulation
# gene.master.table$experiment <- factor(gene.master.table$experiment, levels = c("normal","diauxic","ox","osm","temp"))
# gene.master.table$stress <- factor(gene.master.table$stress, 
#                                                  levels = c("normal_0",paste(rep(c("diauxic","ox","osm","temp"), each=2), rep(c(20,120),times = 4), sep="_")), 
#                                                  labels = c("Non-stress",
#                                                             "Diauxic\n20min", "Diauxic\n120min",
#                                                             "Oxidative\n20min","Oxidative\n120min",
#                                                             "Osmotic\n20min","Osmotic\n120min",
#                                                             "Temperature\n20min","Temperature\n120min")  )


load("results/master tables/SMoPT.table.Rd")

# ddply(subset(gene.master.table, is.na(quant.mRNAab)), .(experiment, time), nrow)

# data.explore.ranking <- subset(gene.master.table, select= c(name, experiment, time, mRNA_abundance, global.protein_synthesis.rate,
#                                                             tRNA_adaptation_index, adj.tAI, adj.tAIFC, adj.tAI.20codons, adj.tAIFC.20codons, log2.mRNA_abundance) )
# 
# # [Ok solved]
# data.explore.ranking <- arrange( ddply(data.explore.ranking, .(experiment, time), mutate, 
#                                        rank = rank(adj.tAI), 
#                                        rank.FC = rank(adj.tAIFC),
#                                        rank.20codons = rank(adj.tAI.20codons),
#                                        rank.FC20codons = rank(adj.tAIFC.20codons)
#                                        ), 
#                                  experiment, time, rank )

      # cor( subset(data.explore.ranking, experiment == "normal")[, c("tRNA_adaptation_index","adj.tAI")], use="pairwise.complete.obs" )
      # 
      # require(reshape2)
      # require(GGally)
      # stAI.matrix <- reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "adj.tAIFC")
      # rank.matrix <- reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.FC")
      # 
      # stAI.20codons.matrix <- reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "adj.tAIFC.20codons")
      # rank.20codons.matrix <- reshape2::dcast( subset(data.explore.ranking, time %in% c(0,20, 120) ), name ~ experiment + time, value.var = "rank.20codons")
      # 
      # 
      # # stAI.matrix <- dcast( subset(data.explore.ranking, time %in% c(0,20) ), name ~ experiment, value.var = "adj.tAI")
      # # rank.matrix <- dcast( subset(data.explore.ranking, time %in% c(0,20) ), name ~ experiment, value.var = "rank")
      # 
      # correlations.rank.stAIFC <- cor(matrixify(rank.matrix))
      # 
      # # stress-adjusted tAI master table ------ 
      # stAI.data <- melt(stAI.matrix, id.vars = c("name","normal_0"))
      # additions <- cbind(stAI.matrix[,c("name","normal_0")],"normal_0",stAI.matrix[,"normal_0"])
      # colnames(additions)[3:4] <- c("variable","value")
      # stAI.data <- rbind(stAI.data, additions)
      # stAI.data$gain.tAI <- log2( stAI.data$value / stAI.data$normal)
      # colnames(stAI.data) <- c("name","tAI","experiment","stAI","gain.tAI")

#rank.20codons
# Compute scrores for rankings and Feature-Scaled stress-adjusted tAI
# stAI.data <- ddply(stAI.data, .(experiment), mutate, 
#                    rank.stAI = rank(stAI), 
#                    rank.tAI = rank(tAI), 
#                    rank.stAI.20codons = rank(stAI), 
#                    rank.tAI.20codons = rank(tAI), 
#                   scaled.stAI = scale(stAI), 
#                   scaled.tAI  = scale(tAI) ) # EDIT MARCH 2015!
# stAI.data <- ddply(stAI.data, .(experiment), mutate, FS.stAI = (scaled.stAI - min(scaled.stAI,na.rm=T)) / (max(scaled.stAI,na.rm=T) - min(scaled.stAI,na.rm=T) ),
#                                                      FS.tAI  = ( scaled.tAI - min(scaled.tAI,na.rm=T) ) / (max( scaled.tAI,na.rm=T) - min( scaled.tAI,na.rm=T) )
#                       ) # for feature-scaled
# stAI.data <- ddply(stAI.data, .(experiment), mutate, FS.stAI = stAI/max(stAI,na.rm=T),
#                                                       FS.tAI = tAI/max(tAI,na.rm=T)
#                       ) # for feature-scaled
# stAI.data$FS.gain.tAI = stAI.data$FS.stAI / stAI.data$FS.tAI

# 
# stAI.data<- ddply(stAI.data, .(experiment), mutate, gain.rank = rank.stAI / max(rank.stAI) / (rank.tAI / max(rank.tAI)), 
#                    quant.gain = ifelse(experiment!="normal_0", quantile.it(gain.rank), NA),
#                    quant.tAI  = ifelse(experiment!="normal_0", quantile.it(stAI), NA)
#                    )
# stAI.data <- stAI.data[,c("name","experiment","tAI","rank.tAI","stAI","rank.stAI", "gain.tAI", "gain.rank","quant.gain","quant.tAI","FS.stAI","FS.tAI","FS.gain.tAI")]
# 
# require(dplyr)
# require(tidyr) # split columns neatly
# stAI.data <- separate(stAI.data, col = "experiment", into = c("experiment","time"), sep="_" )
# 
# stAI.data <- merge(stAI.data, subset(data.explore.ranking, select= c("name","experiment", "time", "mRNA_abundance", "log2.mRNA_abundance", "global.protein_synthesis.rate"),
#                                      by=c("name","experiment", "time")))
# 
# # Rank-based
# stAI.data$tAI.profile  <- factor( ifelse(stAI.data$gain.rank < 0.5, "a", ifelse(stAI.data$gain.rank <2, "b", "c")), 
#                                   levels = letters[1:3], labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
# # Gain-based (change of +/- 30% in tAI) -- least discriminatory (in stress tAI decreases for most genes, but it doesn't reflect relative adapt.)
# stAI.data$tAI.profile2  <- factor( ifelse(stAI.data$gain.tAI < log2(0.7), "a", ifelse(stAI.data$gain.tAI < log2(1.3), "b", "c")), 
#                                    levels = letters[1:3], labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
# # Feature-Scaled-based
# stAI.data$tAI.profile3  <- factor( ifelse(stAI.data$FS.gain.tAI < 0.7, "a", ifelse(stAI.data$FS.gain.tAI < 1.3, "b", "c")), 
#                                    levels = letters[1:3], 
#                                    labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
# stAI.data$experiment <- factor(stAI.data$experiment, levels=c("normal", "diauxic","osm","ox","temp"))
# stAI.data$time <- factor(stAI.data$time, levels=c(0,20,120))
# stAI.data$label <- paste(stAI.data$experiment, paste(stAI.data$time,"min"), sep="\n")
# stAI.data$label <- factor(stAI.data$label, levels= 
#                             c("normal\n0 min",paste( rep(c("diauxic","osm","ox","temp"), times=2), rep(paste(c(20,120), "min"),each=4), sep="\n"))
#                           )
# 
# save(stAI.data, file="results/master tables/stAI.data.R")
#   load("results/master tables/stAI.data.R")






# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #


  # stAI.plot <- ggpairs( stAI.matrix, columns = 2:ncol(stAI.matrix), title = "Stress-adjusted tAI" ) 
  # rep_cells<-c(6,11,12,16:18,21:24)
  # add_params<-"alpha=0.4, size=0.9, color='skyblue3'"
  # 
  # ggally_data<-stAI.plot$data                         # makes sure that the internal parameter picks up your data (it always calls it's data 'ggally_data'
  # calls<-add_p(stAI.plot,rep_cells,params=add_params)  #call the function
  # for(i in 1:length(calls)){stAI.plot<-putPlot(stAI.plot,calls[[i]][1],as.numeric(calls[[i]][2]),as.numeric(calls[[i]][3]))}
  # 
  # 
  # rank.plot <- ggpairs( rank.matrix, columns = 2:ncol(rank.matrix), title = "Rank in stress-adjusted tAI" ) 
  # 
  # pdf(file = paste0(PATH,"part3/adjusted.tAIFC--correlation.stAI.pdf"), width=8.6, height=8.31, useDingbats=FALSE )
  # stAI.plot
  # dev.off()

# ----- [] Principal Component Analysis of ranks in s-tAI
    # require(graphics)
    # require(MASS)
    # 
    # d <- matrixify(rank.matrix)
    # PCA <- princomp( x = d )
    # loadings(PCA)
    # plot(PCA)
    # #biplot(PCA,)
    # 
    # princomp( d, covmat = MASS::cov.rob(d))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #

UP <- "#185665"
DOWN <- "#A36A23"

# [] tAI in general ----- ok
ggplot(subset(gene.master.table, time == 0 ), aes(x=tRNA_adaptation_index, y = stAI) ) + geom_point()

g <- ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = FS.tAI, y = FS.stAI ) ) + 
  geom_point(aes(color=tAI.profile)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  stat_smooth( method = "lm", color = "white", alpha=0.5) + coord_fixed() + scale_color_manual(values=c(DOWN,UP,"gray60")) +
  scale_y_continuous(labels=prettyNum) + scale_x_continuous(labels=prettyNum) +
  facet_wrap( ~ label, nrow=2) + ggtitle("Variations in tRNA Adaptation Index during stress")

ggsave(plot = g, filename = paste0(PATH,"part3/adjusted.tAIFC--variations.pdf"), 
       dpi=250, width=14.1, height=8.1, useDingbats=FALSE )

ggsave(plot = g, filename = paste0(PATH,"part3/adjusted.tAIFC--variations.png"), 
       dpi=250, width=14.1, height=8.1 )



g <- ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = FS.tAI, y = FS.stAI ) ) + 
  geom_point(aes(color=tAI.profile)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  stat_smooth( method = "lm", color = "white", alpha=0.5) + coord_fixed() + scale_color_manual(values=c(DOWN,UP,"gray60")) +
  scale_x_continuous(labels=prettyNum) + scale_y_continuous(labels=prettyNum) +
  facet_wrap(  ~ label, nrow = 2) + ggtitle("Variations in tRNA Adaptation Index during stress")

# ok. Johann.
ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = tAI, y = stAIFC ) ) + 
  geom_point(aes(color=tAI.profile1)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  stat_smooth( method = "lm", color = "white", alpha=0.5)  + scale_color_manual(values=c("red","gray40","orange")) +
  facet_wrap( time ~ experiment, nrow = 2) + ggtitle("Normalisation of tAI")

ggsave(plot = g, filename = paste0(PATH,"part3/adjusted.FS.tAIFC--variations.pdf"), 
       dpi=250, width=14.1, height=8.1, useDingbats=FALSE )


ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = tAI, y = stAIFC ) ) + 
  geom_point(aes(color=tAI.profile3)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  stat_smooth( method = "lm", color = "white", alpha=0.5)  + scale_color_manual(values=c("red","gray40","orange")) +
  facet_wrap( time ~ experiment, nrow = 2) + ggtitle("Normalisation of tAI")



# Rank of tAI
g <- ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = rank.tAI, y = rank.stAI ) ) + 
  geom_point(aes(color=tAI.profile)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  geom_abline( intercept=0, slope =1, lty=3, size = 0.6, color="black", alpha=0.95) +
  geom_abline( intercept=0, slope =2, lty=1, size = 0.6, color="black", alpha=0.95) +
  geom_abline( intercept=0, slope =0.5, lty=1, size = 0.6, color="black", alpha=0.95) +
  coord_fixed() + 
  scale_color_manual(values=c(DOWN,UP,"gray60")) +
  facet_wrap(  ~ label, nrow =2 ) + ggtitle("How the tAIFC-rank of genes changes during stress")

ggsave(plot = g, filename = paste0(PATH,"part3/adjusted.tAIFC--rank.variations.pdf"), 
       dpi=250, width=14.1, height=8.1, useDingbats=FALSE )

ggsave(plot = g, filename = paste0(PATH,"part3/adjusted.tAIFC--rank.variations.png"), 
       dpi=250, width=14.1, height=8.1 )




# ------------------------------------------------ mRNA abundance of genes with biggest changes in adaptability 

g <- ggplot( data = subset(SMoPT.data,experiment!="normal"), aes(x= tAI.profile, y=log2.mRNA_abundance, fill = tAI.profile) ) +
  geom_boxplot() + scale_fill_manual(values=c(DOWN,UP,"gray60")) +
  labs(y="change in mRNA abundance (log2)", x="", fill="tAI profile")+
  facet_wrap(  ~ label, nrow=2) + ggtitle("Changes in mRNA abundance in function of s-tAI")

ggplot( data = subset(SMoPT.data,experiment!="normal"), aes(x= tAI.profile, y=mRNA_abundance, fill = tAI.profile) ) +
  geom_boxplot() + scale_fill_manual(values=c(DOWN,UP,"gray60")) + scale_y_sqrt() +
  facet_wrap(  ~ experiment) + ggtitle("Changes in mRNA abundance in function of s-tAI")

ggplot( data = subset(gene.master.table, time != 0 & !is.na(up_reg)), aes(x= stress, y=stAIFC, fill = up_reg) ) +
  geom_boxplot() + scale_fill_manual(values=c(DOWN,UP,"gray60")) +
   ggtitle("Changes in mRNA abundance in function of s-tAI")

ggplot( data = subset(SMoPT.data,experiment!="normal"), aes(x=log2.mRNA_abundance, y= gain.stAI, color = tAI.profile) ) +
  geom_point() + scale_color_manual(values=c(DOWN,UP,"gray60")) +
  geom_hline(yintercept=0, color = "gray30") + 
  geom_vline(xintercept=0, color = "gray30") + 
  labs(y="gain in tAI (log2)", x="change in mRNA abundance (log2)", fill="tAI profile")+
  facet_wrap(  ~ label, nrow=2) + ggtitle("Changes in translation efficiency in function of changes in mRNA abundance")



# ------------------------------------------------ Gain and loss in adaptiveness ------------------------ #

require(reshape2)
require(venneuler)


gene.sets.table <- reshape2::dcast( subset(SMoPT.data, experiment != "normal"), name ~ stress, value.var = "tAI.profile1" )

# How similar are the sets of genes with increase/decrease in tAI?
gene.sets <- dlply( subset(SMoPT.data, experiment != "normal"), .(experiment, time, tAI.profile1), function(x) as.character(x$name) )

gene.sets.similarity <- data.frame(t(combn( names(gene.sets), 2)))
colnames(gene.sets.similarity)[1:2] <- c("id1","id2")
gene.sets.similarity$group1  <- gsub(gene.sets.similarity$id1, pattern="^(.*)\\.(.*)\n.*$", replacement="\\2")
gene.sets.similarity$group2  <- gsub(gene.sets.similarity$id2, pattern="^(.*)\\.(.*)\n.*$", replacement="\\2")
gene.sets.similarity$Jaccard <- apply(gene.sets.similarity, 1, function(x){ 
  require(sets)
  set_similarity(x= as.set(gene.sets[[ x[[1]] ]]), y= as.set(gene.sets[[ x[[2]] ]]) 
  )
})


# within-group similarity
within.group.Jaccard <- subset(gene.sets.similarity, as.character(id1) != as.character(id2) & as.character(group1) == as.character(group2))
mean(within.group.Jaccard$Jaccard)
median(within.group.Jaccard$Jaccard)


outside.group.Jaccard <- subset(gene.sets.similarity, as.character(id1) != as.character(id2) & as.character(group1) != as.character(group2))
mean(outside.group.Jaccard$Jaccard)

require(venneuler)
# get genes whose rank has at increased most
d.rankup <- dlply( subset(SMoPT.data, experiment != "normal" & gain.rank.stAIFC > 2 & time == 20), .(experiment), function(x) as.character(x$name))
venn.up <- split.venn(vennit(d.rankup))
vd.up <- venneuler(combinations = melt(d.rankup) )

# get genes whose rank has most decreased
d.rankdown <- dlply( subset(SMoPT.data, experiment != "normal" &  gain.rank.stAIFC < 0.5 & time == 20), .(experiment), function(x) as.character(x$name))
venn.down <- split.venn(vennit(d.rankdown))
vd.down <- venneuler(combinations = melt(d.rankdown) )



#         # get genes whose rank has at increased most
#         d.rankhigh <- dlply( subset(SMoPT.data, experiment != "normal" &  quant.tAI == "0-20" & time == 20 ), .(experiment), function(x) as.character(x$name))
#         venn.high <- split.venn(vennit(d.rankhigh))
#         vd.high <- venneuler(combinations = melt(d.rankhigh) )
#         
#         # get genes whose rank has most decreased
#         d.ranklow <- dlply( subset(SMoPT.data, experiment != "normal" &  quant.tAI == "80-100" & time == 20), .(experiment), function(x) as.character(x$name))
#         venn.low <- split.venn(vennit(d.ranklow))
#         vd.low <- venneuler(combinations = melt(d.ranklow) )
#         


# Similarity of up and down regulated genes
# transform into a matrix for visualisation with Pheatmap
m <- get.sym.matrix(df=gene.sets.similarity, value.var="Jaccard")
diag(m) <- 0

require(pheatmap)
require(RColorBrewer)
require(reshape2)
pdf(paste0(PATH,"part3/similarity_stAI_groups.genes.pdf"), width = 11, height = 11)
pheatmap( m, cellwidth = 20, cellheight= 20, 
          color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100))[1:80],
          cluster_rows= T, cluster_cols=T, treeheight_row=15, treeheight_col=15,
          main = "Patterns of tAI variation during stress are similar"
)
dev.off()


pdf(file = paste0(PATH,"part3/adjusted.tAI--sets_with_high_or_low_stAIFC.pdf"))
par(mfrow=c(2,2))
    plot(vd.high, main = "high s-tAI")
    plot(vd.low,  main = "low s-tAI")
    plot(vd.up, main = "s-tAI going increasing")
    plot(vd.down, main = "s-tAI decreasing")
dev.off()


# 
# subset(go.slim.ORFs, Name %in% venn.up$'diauxic&osm&ox&temp')
# subset(go.slim.ORFs, Name %in% venn.down$'diauxic&osm&ox&temp')

go.slim.ORFs <- read.table("data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
  colnames(go.slim.ORFs)[2] <- "name"

  condensed.GO.P <- ddply(subset(go.slim.ORFs, GO.type == "P" & !GO.description %in% c("biological_process","other") ), .(name), function(x){ paste(x$GO.description, collapse = " | ") } )
  colnames(condensed.GO.P)[2] <- "GOs"
  
O <- unique(subset(go.slim.ORFs, GO.type == "P", select=c(GO.type,GO.ID,GO.description)))


# treat conditions separately
enrichment.tests.up <- dlply( subset(SMoPT.data, time == 20 & experiment!="normal"), .(experiment), function(x) {
        subset( automate.fisherexact( data = x, expression = 'gain.rank.stAIFC > 2', keys = O), adjusted.p < 0.01 )
})



enrichment.tests.down <- dlply( subset(SMoPT.data, time == 20 & experiment!="normal"), .(experiment), function(x) {
  #O <- unique(as.character(subset(go.slim.ORFs, GO.type == "P")$GO.description))
  subset( automate.fisherexact( data = x, expression = 'gain.rank.stAIFC < 0.5', keys = O), adjusted.p < 0.01 )
})


enrichment.tests.translation.rate.up <- dlply( subset(SMoPT.data, time == 20 & experiment!="normal"), .(experiment, time), function(x) {
  #O <- unique(as.character(subset(go.slim.ORFs, GO.type == "P")$GO.description))
  subset( automate.fisherexact( data = x, expression = 'ratio_events > 1.5', keys = O), adjusted.p < 0.01 )
})

enrichment.tests.translation.rate.down <- dlply( subset(SMoPT.data, time == 20 & experiment!="normal"), .(experiment, time), function(x) {
  #O <- unique(as.character(subset(go.slim.ORFs, GO.type == "P")$GO.description))
  subset( automate.fisherexact( data = x, expression = 'ratio_events < 0.5', keys = O), adjusted.p < 0.01 )
})


save( enrichment.tests.down, enrichment.tests.up, enrichment.tests.translation.rate.down, enrichment.tests.translation.rate.up, file = "results/stress tAI/GO_enrichments.Rda")

# take all "red" genes and all "orange" genes together (if they are orange or red in at least 2 conditions)
venn.up$`diauxic&osm&ox&temp` 
venn.down$`diauxic&osm&ox&temp` 


# JUNE 2015
# # Check GO Slim enrichment for genes that are consistently red (or consistently orange) in at least 5 conditions (Chapter 5)
# d.rankup.early   <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 & gain.rank.stAIFC > 2), .(experiment), function(x) as.character(x$name))
# d.rankdown.early <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 &  gain.rank.stAIFC < 0.5), .(experiment), function(x) as.character(x$name))
# d.rankup.late   <- dlply( subset(SMoPT.data, experiment != "normal" & time == 120 & gain.rank.stAIFC > 2), .(experiment), function(x) as.character(x$name))
# d.rankdown.late <- dlply( subset(SMoPT.data, experiment != "normal" & time == 120 &  gain.rank.stAIFC < 0.5), .(experiment), function(x) as.character(x$name))

d.rankup.early   <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 & gain.rank.stAIFC > 2), .(experiment), function(x) as.character(x$name))
d.rankdown.early <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 & gain.rank.stAIFC < 0.5), .(experiment), function(x) as.character(x$name))
d.rankup.late   <- dlply( subset(SMoPT.data, experiment != "normal" & time == 120 & gain.rank.stAIFC > 2), .(experiment), function(x) as.character(x$name))
d.rankdown.late <- dlply( subset(SMoPT.data, experiment != "normal" & time == 120 & gain.rank.stAIFC < 0.5), .(experiment), function(x) as.character(x$name))

counts.venn.up.early   <-  dcast(melt(d.rankup.early),  value ~  ., value.var = "L1")  
  colnames(counts.venn.up.early)   <- c("name","n_stAIFC.up")
counts.venn.down.early <-  dcast(melt(d.rankdown.early),  value ~  ., value.var = "L1")  
  colnames(counts.venn.down.early)   <- c("name","n_stAIFC.down")
counts.venn.up.late   <-  dcast(melt(d.rankup.late),  value ~  ., value.var = "L1")  
  colnames(counts.venn.up.late)   <- c("name","n_stAIFC.up")
counts.venn.down.late <-  dcast(melt(d.rankdown.late),  value ~  ., value.var = "L1")  
  colnames(counts.venn.down.late)   <- c("name","n_stAIFC.down")

  
  
d.ranksky.early  <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 & gain.rank.stAIFC > 10), .(experiment), function(x) as.character(x$name))
d.rankdrop.early <- dlply( subset(SMoPT.data, experiment != "normal" & time == 20 & gain.rank.stAIFC < 0.1), .(experiment), function(x) as.character(x$name))

counts.venn.sky.early   <-  dcast(melt(d.ranksky.early),  value ~  ., value.var = "L1")  
  colnames(counts.venn.sky.early)   <- c("name","n_stAIFC.up")
counts.venn.drop.early <-  dcast(melt(d.rankdrop.early),  value ~  ., value.var = "L1")  
  colnames(counts.venn.drop.early)   <- c("name","n_stAIFC.down")





d.initiation.early  <- dlply( subset(study.translation.kinetics, experiment != "normal" & time == 20 & variation_initiation_frequency > 0), .(experiment), function(x) as.character(x$name))

counts.venn.initiation.early <-  dcast(melt(d.initiation.early),  value ~  ., value.var = "L1")  
colnames(counts.venn.initiation.early)   <- c("name","n_initiation.up")


genes.systematically.up.initiation <- as.character(subset(counts.venn.initiation.early, n_initiation.up > 3)$name)
  

table.initiation <- subset(study.translation.kinetics, name %in% genes.systematically.up.initiation & experiment == "temp" & time == 20, select=c(name, av.initiation_time.normal, av.initiation_time))
View( merge(table.initiation, condensed.GO.P, by = "name", all.x =T ) )

  
# intersect(counts.venn.up$name, counts.venn.down$name)
# intersect(subset(counts.venn.up, n_stAIFC.up > 2)$name, subset(counts.venn.down, n_stAIFC.down > 2)$name)
# # no intersection anymore when selecting genes whose rank in stAI goes up at least 3/8 times or goes down at least 3/8 times
# intersect(subset(counts.venn.up, n_stAIFC.up > 3)$name, subset(counts.venn.down, n_stAIFC.down > 3)$name)
# 
# # 788 "orange" genes 
# nrow(subset(counts.venn.up, n_stAIFC.up >= 5))
# # 676 "red" genes
# nrow(subset(counts.venn.down, n_stAIFC.down >= 5))

SMoPT.data$stAI.orange.early <- SMoPT.data$name %in% as.character(subset(counts.venn.up.early, n_stAIFC.up >= 3)$name)
SMoPT.data$stAI.red.early    <- SMoPT.data$name %in% as.character(subset(counts.venn.down.early, n_stAIFC.down >= 3)$name)
SMoPT.data$stAI.orange.late <- SMoPT.data$name %in% as.character(subset(counts.venn.up.late, n_stAIFC.up >= 3)$name)
SMoPT.data$stAI.red.late    <- SMoPT.data$name %in% as.character(subset(counts.venn.down.late, n_stAIFC.down >= 3)$name)
SMoPT.data$strong.gain.rank <- SMoPT.data$name %in% as.character(subset(counts.venn.sky.early, n_stAIFC.up >= 3)$name)
SMoPT.data$strong.loss.rank <- SMoPT.data$name %in% as.character(subset(counts.venn.drop.early, n_stAIFC.down >= 3)$name)

#         complexes <- unique(as.character(go.slim.ORFs$GO.description[ grep("complex", go.slim.ORFs$GO.description) ]))
#         subunit <- unique(as.character(go.slim.ORFs$GO.description[ grep("subunit", go.slim.ORFs$GO.description) ]))
#         ribosom <- unique(as.character(go.slim.ORFs$GO.description[ grep("ribosom", go.slim.ORFs$GO.description) ]))
#         
#         go.slim.ORFs$GO.type <- as.character(go.slim.ORFs$GO.type)
#         go.slim.ORFs[ go.slim.ORFs$GO.description %in% unique(c(complexes,subunit,ribosom)) & go.slim.ORFs$GO.type == "C", "GO.type"] <- "CX" 
#         write.table(go.slim.ORFs, file="data/annotations/go_slim_ORFs.txt", sep="t",quote = F, row.names=F)

# go.slim.ORFs$GO.type <- as.character(go.slim.ORFs$GO.type)
# go.slim.ORFs[ go.slim.ORFs$GO.type == "C", "GO.type"] <- "CX"
# C <- c("membrane","cellular_component","cytoplasm","mitochondrial envelope","mitochondrion","nucleus","chromosome","endomembrane system","cell wall","plasma membrane","vacuole","cytoplasmic membrane-bounded vesicle","endoplasmic reticulum","cellular bud","cytoskeleton","site of polarized growth","spindle","chromosome, centromeric region","nucleolus","Golgi apparatus","condensed nuclear chromosome kinetochore","condensed nuclear chromosome, centromeric region","peroxisome","extracellular region","cell cortex")         
# go.slim.ORFs[ go.slim.ORFs$GO.description %in% C,"GO.type"] <- "C"


        
O <- unique(subset(go.slim.ORFs, GO.type == "C", select=c(GO.type,GO.ID,GO.description)))
#write.table(O[,3,drop=F], file = "~/Desktop/truc.text",quote=F, row.names=F)

enrichment.stAI.orange.early <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'stAI.orange.early == TRUE', keys = O), adjusted.p < 0.01 )
enrichment.stAI.red.early    <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'stAI.red.early    == TRUE', keys = O), adjusted.p < 0.01 )
enrichment.stAI.orange.late <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'stAI.orange.late == TRUE', keys = O), adjusted.p < 0.01 )
enrichment.stAI.red.late    <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'stAI.red.late    == TRUE', keys = O), adjusted.p < 0.01 )

enrichment.stAI.sky.early   <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'strong.gain.rank == TRUE', keys = O), adjusted.p < 0.01 )
enrichment.stAI.drop.early  <- subset( automate.fisherexact( data = subset(SMoPT.data, experiment=="normal"), expression = 'strong.loss.rank == TRUE', keys = O), adjusted.p < 0.01 )


subset(SMoPT.data, name == "SKP2" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "MMS4" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "ILV6" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "ILV6" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, gain.rank.FC, stAI, stAIFC))
subset(SMoPT.data, name == "RPL4A" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "RPL13A" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "RPS16B" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "CYS4" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
subset(SMoPT.data, name == "HIS4" & experiment == "ox", select = c(name, time, gain.rank.stAI, gain.rank.stAIFC, stAI, stAIFC))
unique(SMoPT.data$name)[ grepl(pattern = "^R", unique(gene.master.table$name)) ]

#unique(stAI.data.ox$name)[ grepl(pattern = "^R", unique(stAI.data.ox$name)) ]



# Genes whose rank in stAI increases at least 10 times
sky.stAIFC.rank.genes <- as.character(subset(counts.venn.sky.early, n_stAIFC.up >= 3)$name)
drop.stAIFC.rank.genes <- as.character(subset(counts.venn.drop.early, n_stAIFC.down >= 3)$name)

table.sky <- subset(SMoPT.data, name %in% sky.stAIFC.rank.genes & experiment == "temp" & time == 20, select=c(name, tAI, stAIFC, stAI, gain.rank.stAIFC, rank.tAI, rank.stAIFC ))
table.drop <- subset(SMoPT.data, name %in% drop.stAIFC.rank.genes & experiment == "temp" & time == 20, select=c(name, tAI, stAIFC, stAI, gain.rank.stAIFC, rank.tAI, rank.stAIFC ))
View( merge(table.sky, condensed.GO.P, by = "name", all.x =T ) )
View( merge(table.drop, condensed.GO.P, by = "name", all.x =T ) )






# ------------------------------------------------ High and low in adaptiveness ------------------------ #




# table for most adapted
table.most.adapted <- do.call( cbind, dlply(SMoPT.data, .(experiment), function(x){ y <- head( subset(arrange( x, plyr::desc(stAI)), select = c(name, experiment, stAI, gain.rank,tAI) ), 100 );
                                             y[,3:5] <- round(y[,3:5],2)
                                             return(y)
                                             } ) )


table.least.adapted <- do.call( cbind, dlply(stAI.data, .(experiment), function(x){ y <- head( subset(arrange( x, stAI), select = c(name, experiment, stAI, gain.rank,tAI) ), 100 );
                                                                                   y[,3:5] <- round(y[,3:5],2)
                                                                                   return(y)
} ) )

# table for highest increase in adaptation
table.most.up <- do.call( cbind, dlply(stAI.data, .(experiment), function(x){ y <- head( subset(arrange( x, plyr::desc(gain.rank)), 
                                                                                                     select = c(name,tAI, stAI, gain.rank, quant.tAI, quant.gain ) ), 100 );
                                                                                   y[,2:4] <- round(y[,2:4],2)
                                                                                   return(y)
} ) )


# table for highest increase in adaptation
table.most.down <- do.call( cbind, dlply(stAI.data, .(experiment), function(x){ y <- head( subset(arrange( x, gain.rank), 
                                                                                                select = c(name,tAI, stAI, gain.rank, quant.tAI, quant.gain ) ), 100 );
                                                                              y[,2:4] <- round(y[,2:4],2)
                                                                              return(y)
} ) )


copy2clipboard(table.least.adapted)
copy2clipboard(table.most.adapted)
copy2clipboard(table.most.up)
copy2clipboard(table.most.down)

ddply(stAI.data, .(experiment), nrow)

# Genes whose stAI increases most during stress

copy2clipboard( as.character(subset(counts.venn.up.early, n_stAIFC.up >= 3)$name))
copy2clipboard( as.character(subset(counts.venn.down.early, n_stAIFC.down >= 3)$name) )

copy2clipboard( venn.up$'diauxic&osm&ox&temp' )
copy2clipboard( venn.up$'diauxic' )
copy2clipboard( venn.up$'osm' ) 
copy2clipboard( venn.up$'ox' )
copy2clipboard( venn.up$'temp' )

# Genes whose stAI decreases most during stress
copy2clipboard( venn.down$'diauxic&osm&ox&temp' )
copy2clipboard( venn.down$'diauxic' )
copy2clipboard( venn.down$'osm' ) 
copy2clipboard( venn.down$'ox' )
copy2clipboard( venn.down$'temp' )

# Genes with highest stAI during stress
copy2clipboard( venn.high$'diauxic&osm&ox&temp' )
copy2clipboard( venn.high$'diauxic' )
copy2clipboard( venn.high$'osm' ) 
copy2clipboard( venn.high$'ox' )
copy2clipboard( venn.high$'temp' )

# Genes with lowest stAI
copy2clipboard( venn.low$'diauxic&osm&ox&temp' )
copy2clipboard( venn.low$'diauxic' )
copy2clipboard( venn.low$'osm' ) 
copy2clipboard( venn.low$'ox' )
copy2clipboard( venn.low$'temp' )


# FOCUS
copy2clipboard( venn.high$'diauxic&osm&ox&temp' )
copy2clipboard( venn.low$'diauxic&osm&ox&temp' )
copy2clipboard( venn.up$'diauxic&osm&ox&temp' )
copy2clipboard( venn.down$'diauxic&osm&ox&temp' )


# FOCUS on ORANGE
copy2clipboard( venn.up$diauxic )
copy2clipboard( venn.up$osm )
copy2clipboard( venn.up$ox )
copy2clipboard( venn.up$temp )

# FOCUS on RED
copy2clipboard( venn.down$diauxic )
copy2clipboard( venn.down$osm )
copy2clipboard( venn.down$ox )
copy2clipboard( venn.down$temp )

ggplot( data = subset(study.translation.kinetics, experiment == "ox" ), 
        aes( 
          x = variation_speed,
          colour = name %in% as.character(subset(go.slim.ORFs, GO.ID ==  "GO:0006979" )$name),
          y = variation_initiation_frequency
        )) + geom_point() + facet_wrap( ~ time)



venn.up.stAI.early <- read.table(file = "results/GO_yeastmine/venn.up.stAI.early.tsv",header=F, sep="\t")
colnames(venn.up.stAI.early) <- c("GO","p.value","genes")
venn.up.stAI.early$matches <- sapply(as.list(venn.up.stAI.early$genes), function(x) length(unlist(strsplit(as.character(x), fixed = T,split = ","))) )
venn.up.stAI.early$example <- sapply(as.list(venn.up.stAI.early$genes), function(x) paste(fNAME(unlist(strsplit(as.character(x), fixed = T,split = ",")))[1:5], collapse="|") )
require(stargazer)

stargazer( arrange( venn.up.stAI.early[,-3], p.value ), summary = F,
           out = "~/Dropbox/PhD/thesis/TeX/items/tables/appendix/table_stAIup.GO.draft.tex"
)







#----- [] Check Metabolic enzymes ------
require(StingRay)
data(GEN)
Pleio.Met.Enz <- as.character(unique(GEN$metabolic.pleiotropic.enzymes$orf))
PTM.Enz <- as.character( unique(GEN$ptm.enzymes$orf) )
Riboproteins <- unique(GEN$riboproteins$ORFA, GEN$riboproteins$ORFB)
translation <- unique(do.call(c, lapply(GEN$translation, function(x) as.character(x$gene))))

# --- PTM enzymes

#### Global protein synthesis rate
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = global.protein_synthesis.rate, 
                                                                          fill = ifelse(ORF %in% PTM.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + scale_y_log10() + labs( fill = "PTM enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Transcriptional response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2.mRNA_abundance, 
                                                                          fill = ifelse(ORF %in% PTM.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "PTM enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Adjusted tAI (estimated elongation response)
ggplot( data = subset(gene.master.table, !is.na(stress) ), aes( x = stress, y = adj.tAI, fill = ifelse(ORF %in% PTM.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "PTM enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Elongation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_elongation), 
                                                                          fill = ifelse(ORF %in% PTM.Enz, T,F) )
) + geom_hline( yintercept = 0, color = "gray") +
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "PTM enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Initiation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_total.time), 
                                                                          fill = ifelse(ORF %in% PTM.Enz, T,F) )
) + geom_hline( yintercept = 0, color = "gray") + 
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "PTM enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))



# --- Pleiotropic metabolic enzymes

#### Global protein synthesis rate
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = global.protein_synthesis.rate, 
                                                                          fill = ifelse(ORF %in% Pleio.Met.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + scale_y_log10() + labs( fill = "Metabolic enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Transcriptional response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2.mRNA_abundance, 
                                                                          fill = ifelse(ORF %in% Pleio.Met.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "Metabolic enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Adjusted tAI (estimated elongation response)
ggplot( data = subset(gene.master.table, !is.na(stress) ), aes( x = stress, y = adj.tAI, 
                                                                fill = ifelse(ORF %in% Pleio.Met.Enz, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "Metabolic enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Elongation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_elongation), 
                                                                          fill = ifelse(ORF %in% Pleio.Met.Enz, T,F) )
) + geom_hline( yintercept = 0, color = "gray") +
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "Metabolic enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Initiation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_total.time), 
                                                                          fill = ifelse(ORF %in% Pleio.Met.Enz, T,F) )
) + geom_hline( yintercept = 0, color = "gray") + 
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "Metabolic enzymes") + scale_fill_manual( values = c("gray40","#A4C663"))


# --- Riboproteins 


#### Global protein synthesis rate
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = global.protein_synthesis.rate, 
                                                                          fill = ifelse(ORF %in% Riboproteins, T,F) )
) + geom_boxplot(notch = F, width= 0.8) + scale_y_log10() + labs( fill = "Riboproteins") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Transcriptional response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2.mRNA_abundance, 
                                                                          fill = ifelse(ORF %in% Riboproteins, T,F) )
) + geom_boxplot(notch = F, width= 0.8) + labs( fill = "Riboproteins") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Adjusted tAI (estimated elongation response)
ggplot( data = subset(gene.master.table, !is.na(stress) ), aes( x = stress, y = adj.tAI, 
                                                                fill = ifelse(ORF %in% Riboproteins, T,F) )
) + geom_boxplot(notch = F, width= 0.8) + labs( fill = "Riboproteins") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Elongation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_elongation), 
                                                                          fill = ifelse(ORF %in% Riboproteins, T,F) )
) + geom_hline( yintercept = 0, color = "gray") +
  geom_boxplot(notch = F, width= 0.8) + labs( fill = "Riboproteins") + scale_fill_manual( values = c("gray40","#A4C663"))

#### Initiation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_initiation), 
                                                                          fill = ifelse(ORF %in% Riboproteins, T,F) )
) + geom_hline( yintercept = 0, color = "gray") + 
  geom_boxplot(notch = F, width= 0.8) + labs( fill = "Riboproteins") + scale_fill_manual( values = c("gray40","#A4C663"))


# --- Translation genes

#### Global protein synthesis rate
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = global.protein_synthesis.rate, 
                                                                          fill = ifelse(ORF %in% translation, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + scale_y_log10() + labs( fill = "translation") + scale_fill_manual( values = c("gray40","orange"))

#### Transcriptional response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2.mRNA_abundance, 
                                                                          fill = ifelse(ORF %in% translation, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "translation") + scale_fill_manual( values = c("gray40","orange"))

#### Adjusted tAI (estimated elongation response)
ggplot( data = subset(gene.master.table, !is.na(stress) ), aes( x = stress, y = adj.tAI, 
                                                                fill = ifelse(ORF %in% translation, T,F) )
) + geom_boxplot(notch = T, width= 0.8) + labs( fill = "translation") + scale_fill_manual( values = c("gray40","orange"))

#### Elongation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_elongation), 
                                                                          fill = ifelse(ORF %in% translation, T,F) )
) + geom_hline( yintercept = 0, color = "gray") +
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "translation") + scale_fill_manual( values = c("gray40","orange"))

#### Initiation response
ggplot( data = subset(gene.master.table, !is.na(stress) & time >0 ), aes( x = stress, y = log2(ratio_total.time), 
                                                                          fill = ifelse(ORF %in% translation, T,F) )
) + geom_hline( yintercept = 0, color = "gray") + 
  geom_boxplot(notch = T, width= 0.8) + labs( fill = "translation") + scale_fill_manual( values = c("gray40","orange"))

