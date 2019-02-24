#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#     OBJECTIVES    #--------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# Here I try to identify biological factors that explain variations in tRNA abundance, be it by stress condition or at a multidimensional level
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


# plot_labeller <- function(variable,value){
#   experiment_names = list( 'diauxic' = "diauxic shift", 'ox' = "oxidative stress", 'osm' = "osmotic shock", 'temp' = "temperature shock" )
#   time_names = list( '0' = "t=0", '20' = "t=20",'60' = "t=60",'120' = "t=120")  
#   if (variable=='experiment') {
#     return(experiment_names[value])
#   } else {
#     return(time_names[value])
#   }
# }

#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#|||    Data        #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# >> correspondence between codon, anticodons and amino acids
rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1,sep="\t")
rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",header=1,sep="\t")
    
    nonunique.anticodons <- as.character(subset(rosetta.anticodons, aa.with.1.anticodon == F )$anticodon)

# >> orgininal quantification (raw data)
      # fold change data with standard errors for the triplicates (aim is to evaluate the quality of the quantification)
marc.data <- read.table("data/tRNA abundance/tRNAab_quantified.mean_sd.txt",header=1)
marc.data <- merge(marc.data, rosetta.anticodons, by="anticodon", all.x=T)
marc.data$anticodon <- factor( marc.data$anticodon, levels = as.character(unique(arrange(marc.data, aa)$anticodon)) )
marc.data$GC.anti  <- gc.content(marc.data$anticodon)

cross.condition.SD <- ddply(marc.data, .(anticodon), summarise, cross.cond.sd = sd(average), cross.cond.disp = cross.cond.sd/mean(average) )
marc.data <- merge(marc.data, cross.condition.SD, by="anticodon")
marc.data$sample.disp <- marc.data$std^2/marc.data$average
marc.data$label <- paste(marc.data$anticodon, marc.data$aa, sep="-")
marc.data$CE <- apply(marc.data, 1, function(x){ conf.int(mu = as.numeric(x[["average"]]), sigma = as.numeric(x[["std"]]), conf = 0.1) })
marc.data$sem <- with( marc.data, sqrt( std^2/3  ))

##write.table(marc.data, "results/tRNA.abundance/marc.data.txt",sep="\t", row.names=F)
marc.data <- read.csv("results/tRNA.abundance/marc.data.txt",sep="\t")

# >> case of Met-tRNAs -----

methionine.table        <- read.table("results/master tables/CAU.table.txt",header=1)  

# >> master tables (processed data) -----
            # complex summary tables that integrate many calculations, summary statistics and theoretical information
anticodon.master.table  <- read.table("results/master tables/anticodon.master.table.txt", header=1) # especially adjusted tGCN, anticodon supply, relative availability
codon.master.table      <- read.table("results/master tables/codon.master.table.txt", header=1) # especially relative adaptiveness and anticodon deman
gene.master.table       <- read.table("results/master tables/gene.master.table.txt", header=1) # especially s-tAI scores and translation dynamics data from stochastic simulation

# >> master table for the simulation of protein translation
SMoPT.data <- subset(gene.master.table, !is.na(experiment) & !is.na(faster.translation) ) # version of the gene.master.table restricted to genes (in a given condition and time point) that have successfully ran through the stochastic simulation (ie used in the simulation and had a measurement associated)
SMoPT.data <- ddply(SMoPT.data, .(experiment, time, stress), function(x){  x$ratio_total.time.q <- quantile.it( x$ratio_total.time, N = 7 ); return(x) })
SMoPT.data <- ddply(SMoPT.data, .(experiment, time, stress), function(x){  x$ratio_expected.translation.rate.q <- quantile.it( x$ratio_expected.translation.rate, N = 7 ); return(x) })
SMoPT.data <- ddply(SMoPT.data, .(experiment, time, stress), function(x){  x$total.time.q <- quantile.it( x$av.initiation_time + x$av.elongation_time, N = 10 ); return(x) })
SMoPT.data <- ddply(SMoPT.data, .(experiment, time, stress), function(x){  x$ratio_initiation.q <- quantile.it( x$ratio_initiation , N = 10 ); return(x) })

SMoPT.data$tAI.scaled # ex FS.tAI
SMoPT.data$stAI.scaled # ex FS.stAI
SMoPT.data$gain.stAI.scaled # ex. FS.gain.tAI

SMoPT.data<- ddply(SMoPT.data, .(experiment, time), mutate, 
                  gain.rank = rank.stAI / max(rank.stAI) / (rank.tAI / max(rank.tAI)), 
                  gain.rank.FC = rank(stAIFC) / max(rank(stAIFC)) / (rank.tAI / max(rank.tAI)),
                  quant.gain = ifelse(experiment!="normal", quantile.it(gain.rank), NA),
                  quant.tAI  = ifelse(experiment!="normal", as.character(quantile.it(stAI)), NA)
)

# Rank-based
SMoPT.data$tAI.profile  <- factor( ifelse(SMoPT.data$gain.rank < 0.5, "a", ifelse(SMoPT.data$gain.rank <2, "b", "c")), 
                                  levels = letters[1:3], labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))

SMoPT.data$tAI.profile1  <- factor( ifelse(SMoPT.data$gain.rank.FC < 0.5, "a", ifelse(SMoPT.data$gain.rank.FC <2, "b", "c")), 
                                   levels = letters[1:3], labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))


# Gain-based (change of +/- 30% in tAI) -- least discriminatory (in stress tAI decreases for most genes, but it doesn't reflect relative adapt.)
SMoPT.data$tAI.profile2  <- factor( ifelse(SMoPT.data$gain.stAI < log2(0.7), "a", ifelse(SMoPT.data$gain.stAI < log2(1.3), "b", "c")), 
                                   levels = letters[1:3], labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
# Feature-Scaled-based
SMoPT.data$tAI.profile3  <- factor( ifelse(SMoPT.data$gain.stAI.scaled < 0.7, "a", ifelse(SMoPT.data$gain.stAI.scaled < 1.3, "b", "c")), 
                                   levels = letters[1:3], 
                                   labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))

# Feature-Scaled-based
SMoPT.data$tAI.profile3  <- factor( ifelse(SMoPT.data$gain.stAI.scaled < 0.7, "a", ifelse(SMoPT.data$gain.stAI.scaled < 1.3, "b", "c")), 
                                    levels = letters[1:3], 
                                    labels = c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))




SMoPT.data$label <- paste(SMoPT.data$experiment,SMoPT.data$time,sep="\n")
SMoPT.data$label <- factor(SMoPT.data$label, 
                           levels= c(paste( rep(c("diauxic","ox","osm","temp"), times = 2 ), rep(c(20,120), each = 4 ), sep="\n"),"normal\n0 min"),
                           labels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min")
)




      # load("results/master tables/stAI.data.R")
      # SMoPT.data <- merge(SMoPT.data, stAI.data[,setdiff(colnames(stAI.data), c("est.mRNA_abundance","log2.mRNA_abundance", "global.protein_synthesis.rate") )], by=c("name","experiment","time"))
write.table(SMoPT.data, file = "results/master tables/SMoPT.table.txt",sep="\t", quote=T, row.names=F)
save(SMoPT.data, file = "results/master tables/SMoPT.table.Rd")

# May 2018  ---- Science Signalling 
load( file = "results/master tables/SMoPT.table.Rd")
head(SMoPT.data)

# >> Cool == verify that the binning sort of works :) ----- 
ddply(SMoPT.data, .(experiment, time, ratio_total.time.q), summarise, max = max(ratio_total.time), n= length(ratio_total.time) )


# >> tRNA response to stress over time
            # variant of the anticodon master table data that contain t=60min (not used in the stochastic simulation), but no simulation data
tRNAab_quantified.foldchange_long <- read.table("data/tRNA abundance/tRNAab_quantified.foldchange_long_plusCAU2.txt",header=1,sep="\t")
  tRNAab.heatmap  <- subset( tRNAab_quantified.foldchange_long, time != 0)
  tRNAab.heatmap$stress  <- factor( tRNAab.heatmap$stress, 
                                    levels = paste( rep(c("diauxic","ox","osm","temp"), times = 3 ), rep(c(20,60,120), each = 4 ), sep="_")
                                  )
  tRNAab.heatmap$time   <- factor(tRNAab.heatmap$time, levels = c(20,60,120), labels=c("20min","60min","120min") )  
  tRNAab.heatmap$label  <-  paste(tRNAab.heatmap$anticodon, tRNAab.heatmap$aa, sep="-")

  head(tRNAab.heatmap)
  head(subset(marc.data, anticodon == "CAU2"))

#   # September 2015
#   tRNAab.heatmap.update <- subset(marc.data, anticodon == "CAU2", select = c(anticodon, experiment, time, average ))
#   tRNAab.heatmap.update$log2.FC <- round(log2(tRNAab.heatmap.update$average),2)
#     
#   arrange( subset(tRNAab.heatmap.update, anticodon == "CAU2"), experiment, time)
#   
# ------ transform that latter dataset (tRNAab_quantified.foldchange_long) into a list of matrices to create heat maps of fold change in tRNA abundance
# Visualise the changes in abundance over time by using heatmaps
require(reshape2)
      M.tRNAab <- list()
      M.tRNAab$all.together   <- reshape2::dcast(subset(tRNAab.heatmap), formula = label ~ stress, value.var = "log2.foldchange" )
      M.tRNAab$all.but.temp   <- reshape2::dcast(subset(tRNAab.heatmap, experiment != "temp"),    formula = anticodon ~ stress,          value.var = "log2.foldchange" )
      M.tRNAab$diauxic        <- reshape2::dcast(subset(tRNAab.heatmap, experiment == "diauxic"), formula = label ~ time,value.var = "log2.foldchange" )
      M.tRNAab$osm            <- reshape2::dcast(subset(tRNAab.heatmap, experiment == "osm"),     formula = label ~ time,value.var = "log2.foldchange" )
      M.tRNAab$ox             <- reshape2::dcast(subset(tRNAab.heatmap, experiment == "ox"),      formula = label ~ time,value.var = "log2.foldchange" )
      M.tRNAab$temp           <- reshape2::dcast(subset(tRNAab.heatmap, experiment == "temp"),    formula = label ~ time,value.var = "log2.foldchange" )
      M.tRNAab$all.020        <- reshape2::dcast(subset(tRNAab.heatmap, time == "20min"), formula = label ~ stress, value.var = "log2.foldchange" )
      M.tRNAab$all.060        <- reshape2::dcast(subset(tRNAab.heatmap, time == "60min"), formula = label ~ stress, value.var = "log2.foldchange" )
      M.tRNAab$all.120        <- reshape2::dcast(subset(tRNAab.heatmap, time == "120min"), formula = label ~ stress, value.var = "log2.foldchange" )
      M.tRNAab <- lapply(M.tRNAab, matrixify )
      M.tRNAab$all.together.2 <- M.tRNAab$all.together[, paste( rep(c("diauxic","ox","osm","temp"), each = 3 ), rep(c(20,60,120), time = 4 ), sep="_") ]

  # define colors for heatmaps
  # coolColors <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100)





# Correlations between stress conditions over time 
COR <- list( t020 = cor(M.tRNAab$all.020), 
             t060 = cor(M.tRNAab$all.060), 
             t120 = cor(M.tRNAab$all.120))

COR <- lapply( COR, function(x){ rownames(x) <- gsub(rownames(x), pattern="^(.*)_(.*)$", replacement="\\1")
                            colnames(x) <- gsub(colnames(x), pattern="^(.*)_(.*)$", replacement="\\1");
                            x
                            } )



pdf(paste0(PATH,"part1/tRNAab--correlations.per.stress.pdf"), width = 4, height = 7)
gulmap( COR$t020, clustering_method = "complete",scale = "none", 
        treeheight_row=20, treeheight_col=20,
        cluster_rows = T, main = "20 min", display_numbers = T, legend = F,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        cluster_cols = T, cellwidth=25, cellheight=25 )

gulmap( COR$t060, clustering_method = "complete",scale = "none",
        treeheight_row=20, treeheight_col=20,
        cluster_rows = T, main = "60 min", display_numbers = T, legend = F,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        cluster_cols = T, cellwidth=25, cellheight=25 )

gulmap( COR$t120, clustering_method = "complete",scale = "none", 
        treeheight_row=20, treeheight_col=20,
        cluster_rows = T, main = "120 min", display_numbers = T, 
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10", notecol="white",
        cluster_cols = T, cellwidth=25, cellheight=25 )
dev.off()

# Correlation between tRNAs across stress and over time
annotations <- anticodon.master.table[, c("anticodon", "aa","stress", "adjusted.tGCN", "relative.availability", 
                                          "CAI.O", "nTE.O","demand.mRNAab.up", "demand.mRNAab") ]

annotations <- reshape(annotations, idvar=c('anticodon',"aa", "CAI.O","nTE.O"), timevar='stress', direction='wide')

# AAAAGH!!!!
annotations.2 <- matrixify(annotations)
rownames(annotations) <- paste(annotations$anticodon, annotations$aa, sep="-")
annotations <- annotations[,-c(1:2)]

pdf(paste0(PATH,"part1/tRNAab--correlations.per.tRNAs.pdf"), width = 11, height = 11)
gulmap( cor(t(M.tRNAab$all.020)), clustering_method = "complete",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations at 20 min", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

gul <- gulmap( cor(t(M.tRNAab$all.120)), clustering_method = "complete",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations at 120 min", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )


gulmap( cor(t(M.tRNAab$all.together)), clustering_method = "complete",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_120"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Cross-dataset correlations", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

gulmap( cor(t(M.tRNAab$diauxic)), clustering_method = "mcquitty",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations diauxic shift", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

gulmap( cor(t(M.tRNAab$ox)), clustering_method = "mcquitty",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations oxidative stress", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

gulmap( cor(t(M.tRNAab$osm)), clustering_method = "mcquitty",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations osmotic stress", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

gulmap( cor(t(M.tRNAab$temp)), clustering_method = "mcquitty",scale = "none",
        annotation = annotations[,c("CAI.O", "adjusted.tGCN.normal_0","demand.mRNAab.up.ox_20"),drop=F], 
        annotation_colors= list(CAI.O = c(N="skyblue",N.O="gray",O="orange"),
                                adjusted.tGCN.normal_0  = c("white","purple")), 
        cluster_rows = T, main = "Correlations temperature stress", display_numbers = F, legend = T,
        Min = -1, centeredOn=0, Max = 1, border_color = "gray10",
        treeheight_row =25, treeheight_col=25,
        cluster_cols = T, cellwidth=12, cellheight=12 )

dev.off()

cluster.labels <- ddply(melt(cutree(gul$tree_row, k=4)), .(value), nrow)
colnames(cluster.labels) <- c("cluster","n")
clusters <- cutree( gul$tree_row, k =4 )
clusters <- data.frame( anticodon = gsub(names(clusters), pattern="^(.*)-(.*)$", replacement="\\1"), cluster = clusters, nametag= names(clusters))
clusters <- merge(subset(marc.data), clusters, by = "anticodon")
clusters <- merge(cluster.labels, clusters, by="cluster", all.y=T)
clusters$n <- factor(clusters$n, levels=c(8,14,3,16))


cluster.profiles <- ddply( clusters, .(cluster, time), summarise, av.fold.change = log2(mean(average)) )
cluster.profiles$cluster <- factor(cluster.profiles$cluster, levels=cluster.labels$cluster, labels= cluster.labels$n )

g <- ggplot( data = cluster.profiles, aes( x = factor(time), y = av.fold.change, shape = factor(cluster)) ) + 
  geom_line(aes(group=cluster)) + theme(legend.position="bottom") +
  geom_point()   
ggsave(plot=g, filename = paste0(PATH,"part1/tRNAab--correlations.clusters.pdf"), width = 3, height=2.5, useDingbats=F) 


g <- ggplot( data = subset(clusters, time==120 & anticodon != "CAA"), aes( x = experiment, y = (average)) ) +
  geom_line( aes(group=anticodon), alpha=0.5 ) + 
  geom_hline(yintercept=1, color = "gray50") +
  geom_line( data = ddply(subset(clusters, time==120 & anticodon !="CAA"), .(n, experiment), summarise,  average = mean(average)), lwd=1.3, aes(group="ah") ) +
  geom_point( data = ddply(subset(clusters, time==120 & anticodon !="CAA"), .(n, experiment), summarise,  average = mean(average)), aes(shape=factor(n)), size = 3 ) +
  geom_text( data = subset(clusters, time==120 & experiment=="temp" & anticodon !="CAA" ), aes(label=nametag), size = 2, hjust=-0.1, angle=0 ) +
  labs(x="", y="log2") +
  facet_wrap( ~ n) + theme(legend.position="none") 
ggsave(plot=g, filename = paste0(PATH,"part1/tRNAab--correlations.clusters_stress.pdf"), width = 4, height=3.5, useDingbats=F) 

# Figure "quality control" 
ggplot( data = marc.data, aes( x = (average), y = std/average, group = anticodon) ) + 
  geom_point() +  ggtitle("quality of tRNA abundance measurements") +
  facet_wrap( ~ experiment + time, nrow=2)



# >> focus on specific time points -----------
    # t= 20min ----------------- #  
    data.20        <- subset( anticodon.master.table, time %in% c(20,0))
    data.20$stress <- factor(data.20$stress, levels=c("normal_0","ox_20","diauxic_20","osm_20","temp_20"),
                              labels=c("Non-stress", "20' Oxidative", "20' Diauxic", "20' Osmotic", "20' Temperature"))
    data.20$jitter.x <- jitter(as.numeric(data.20$stress), factor = 0.5 )

    # t=120min ----------------- #
    data.120        <- subset( anticodon.master.table, time %in% c(120,0))
    data.120$stress <- factor(data.120$stress, levels=c("normal_0","ox_120","diauxic_120","osm_120","temp_120"),
                          labels=c("Non-stress", "120' Oxidative", "120' Diauxic", "120' Osmotic", "120' Temperature"))
    data.120$jitter.x <- jitter(as.numeric(data.120$stress), factor = 0.5 )



# >> codon frequency table (10 first codons) ----------------- 
    codon.freq.matrix.20codons <- read.table("results/stress tAI/codon.freq.matrix_first20codons.mat",sep="\t",header=1)


# >> amino acid demand ----------------- 

    aa.master.table <- ddply( anticodon.master.table, .(aa, experiment, time), 
                              summarise, 
                              aa.demand = sum(demand.mRNAab),   aa.demand.up = sum(demand.mRNAab.up), 
                              aa.supply = sum(total.tRNAab),  aa.supply.free = sum(free.tRNAab)
                            )

    aa.master.table <- ddply( aa.master.table, .(experiment, time), mutate,  
                         rel.aa.demand    = round( aa.demand / sum(aa.demand), 3),
                         rel.aa.demand.up = round( aa.demand.up / sum(aa.demand.up), 3),
                         rel.aa.supply    = round( aa.supply / sum(aa.supply), 3),
                       )

    aa.master.table <- merge(aa.master.table, TAB$aminoacid.properties, by = "aa" )

# >> anticodon enrichment in genes -----

    enrichment.table  <- read.table("results/master tables/anticodon.enrichment.table.txt",  sep="\t", header=1)
    enrichment.table  <- merge(enrichment.table, 
                               subset( study.translation.kinetics, select=c(ORF, experiment, time, stress, ratio_initiation, ratio_elongation, faster.translation, global.protein_synthesis.rate  )),
                               by = "ORF")
#    enumeration.table <- read.table("results/master tables/anticodon.enumeration.table.txt", sep="\t", header=1)
#     enumeration.table <- merge(enumeration.table, 
#                                subset( SMoPT.data, select=c(ORF, experiment, time, stress, ratio_initiation, ratio_elongation, faster.translation, global.protein_synthesis.rate  )),
#                                by = "ORF")    
#     # should be optional
    anticodon.enrichment <- melt(enrichment.table,
                                 id.vars = c("ORF", "experiment", "time", "stress", "ratio_initiation", "ratio_elongation", "faster.translation", "global.protein_synthesis.rate") )

#     CCG.enrichment.table <- subset(unique(anticodon.enrichment), variable == "CCG" & !is.na(ratio_initiation) )   
#     CCG.enrichment.table <- ddply(CCG.enrichment.table, .(experiment, time, stress), function(x){  
#         x$ratio_initiation.q <- as.numeric(quantile.it( x$ratio_initiation, N = 7 )); 
#         return(x) })

# >> post-loading processing -------

#     # get box plot statistics for everyone -- t = 20min
#     anticodon.enrichment.stats <- do.call(rbind, dlply( subset(anticodon.enrichment, time == 20), .( stress, variable, faster.translation ), function(x){
#       data.frame( t(c( stress = unique(as.character(x$stress)), 
#                        anticodon = unique(as.character(x$variable)), 
#                        class = unique(as.character(x$faster.translation))
#       )),
#       t( as.data.frame( do.call(c, boxplot.stats(x$value)[1:3] ) ) )
#       )
#     }) 
#     )
#     rownames(anticodon.enrichment.stats) <- NULL
#     
#     anticodon.enrichment.FTS <- cast( data = melt( subset(anticodon.enrichment.stats, select = c(stress, anticodon, class, stats3, conf1, conf2)), id.vars = c("stress","anticodon","class")) , 
#                                       formula = stress + anticodon ~ class + variable
#     )
#     anticodon.enrichment.FTS$median.diff  <- anticodon.enrichment.FTS$TRUE_stats3 - anticodon.enrichment.FTS$FALSE_stats3
#     anticodon.enrichment.FTS <- arrange( anticodon.enrichment.FTS, median.diff)
# 
# 
#     # Anticodon relative frequency in the first 20 codon positions
#     enumeration.20.table <- read.table("results/master tables/anticodon.enumeration.20.table.txt", sep="\t", header=1)
#     enumeration.20.table <- merge(enumeration.20.table, 
#                                   subset( SMoPT.data, select=c(ORF, experiment, time, stress, faster.translation )),
#                                   by = "ORF")
#     
#     anticodons <- unique(as.character(anticodon.enrichment$variable)) 
#     anticodon.counts.20 <- ddply( enumeration.20.table, .(experiment, time, stress, faster.translation), function(x){ colSums(x[,anticodons]) } )
#     
#     # percentage of presence of anticodons in the first 20 codon positions of the CDS
#     anticodon.percent.first.20 <- ddply( anticodon.counts.20, .(experiment, time, stress, faster.translation), 
#                                         function(x) { round( x[,anticodons] / rowSums(x[,anticodons]), 4)*100 } )
# 
#     test <- subset( cast( data = melt(anticodon.percent.first.20, id.vars = c("experiment","time","stress", "faster.translation")), 
#       formula = experiment + time + stress + variable ~ faster.translation ), time != 0 )
#     test$diff <- test$`TRUE` - test$`FALSE`
#     test$enrichment  <- test$`TRUE` / test$`FALSE`
# 
#     anticodon.percent.first.20 <-  merge( anticodon.master.table, test, 
#            by.x = c("anticodon","stress","experiment","time"), 
#            by.y = c("variable", "stress","experiment","time") )


# >> order factors ----------------- 

anticodon.enrichment$experiment   <- factor(anticodon.enrichment$experiment, levels=c("normal","diauxic","ox","osm","temp")  )
anticodon.master.table$experiment <- factor(anticodon.master.table$experiment, levels=c("normal","diauxic","ox","osm","temp")  )
codon.master.table$experiment     <- factor(codon.master.table$experiment, levels=c("normal","diauxic","ox","osm","temp")  )
gene.master.table$experiment      <- factor(gene.master.table$experiment, levels=c("normal","diauxic","ox","osm","temp")  )
SMoPT.data$experiment             <- factor(SMoPT.data$experiment, levels=c("normal","diauxic","ox","osm","temp")  )

 
# first estimate change demand.mRNAab in  conditions for up-regulated genes 

up.reg.genes <- dlply( gene.master.table, .(experiment,time), function(x){ subset(x, up_reg == T, select=c(experiment, time, name)) })

#subset(anticodon.master.table, ORF %in%  )



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------#||| Analysis Part I  #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# I want to present clearly the variations of tRNA abundance in response to stress, in particular identify which tRNAs (anticodon) go up/down when, 
# and how this affect the balance of supply and demand for codons and amino acids

# In this section, the main research objective is to characterise the variations in tRNA abundance measured in different stress conditions.
  # Are the changes in abudnance signficant? 
  # Can we expect them to have consequences on global protein synthesis rates? 
  # Can we expect these changes to differentially impact genes in their translation efficiency?
  # Which tRNA go up/down? 
  # Can we find determinants / explanatory factors for why some tRNA go up and some go down?



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------# Global tRNA response #-------------------------------------------------------------------------------------------------#
#-----------------------------# Multivariate Analysis #-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# How to characterise the changes in tRNA abundance - multidimensional analysis
# Here the goal is to find a way to distinguish patterns of changes in abundance that 
# distinguish sets of tRNAs


# For this I will try affinity propagation-based clustering
          #require("apcluster")
          
          #tRNAab.affinity.prop <- build.ap.clusters(data = M.tRNAab$all.together )
          #tRNAab.affinity.prop <- build.ap.clusters(data = M.tRNAab$all.but.temp[setdiff( nonunique.anticodons, "GAU"),] )
          #tRNAab.affinity.prop <- build.ap.clusters(data = M.tRNAab$diauxic )
            # 
            # plot( tRNAab.affinity.prop$negDistMat, M.tRNAab$all.but.temp )
            # plot( tRNAab.affinity.prop$expSimMat,  M.tRNAab$all.but.temp )
            # plot( tRNAab.affinity.prop$linKernel,  M.tRNAab$all.but.temp )
            # 
            # plot( tRNAab.affinity.prop$negDistMat, M.tRNAab$diauxic )
            # plot( tRNAab.affinity.prop$expSimMat,  M.tRNAab$diauxic )
            # plot( tRNAab.affinity.prop$linKernel,  M.tRNAab$diauxic )
          
            # heatmap( negDistMat(M.tRNAab$all.but.temp) )
            # heatmap( negDistMat(M.tRNAab$diauxic) )
            # heatmap( negDistMat(M.tRNAab$osm) )
            # heatmap( negDistMat(M.tRNAab$ox) )
          
          # dlply( identify.clusters( tRNAab.affinity.prop$negDistMat ), .(cluster.id) )
          # dlply( identify.clusters( tRNAab.affinity.prop$expSimMat  ), .(cluster.id) )
          # dlply( identify.clusters( tRNAab.affinity.prop$linKernel  ), .(cluster.id) )


# ---- t-distributed Stochastic Neighbour Embedding ----
# Another way to reduce dimensions
require("tsne")
data(TAB)
rosetta.anticodons$heatmap.label <- paste(rosetta.anticodons$anticodon, rosetta.anticodons$aa, sep="-")
aa.with.multiple.tRNAs <- names( which( table(rosetta.anticodons$aa) > 1 ) )
table <- subset(rosetta.anticodons, aa %in% aa.with.multiple.tRNAs)
isoacceptor.anticodons <- intersect( unique(table$heatmap.label), rownames(M.tRNAab$all.together))

M <- M.tRNAab$all.together[isoacceptor.anticodons,]
tRNAab.tsne <- data.frame(tsne( dist(M), perplexity = 15, k = 2 ))
colnames(tRNAab.tsne)[1:2] <- c("Dim1","Dim2")
tRNAab.tsne$heatmap.label <- rownames(M)
tRNAab.tsne <- merge(tRNAab.tsne, rosetta.anticodons,by="heatmap.label")

assignment <- ldply( dlply(tRNAab.tsne, .(aa), function(x){
  m <- x[,c("Dim1", "Dim2")]
  rownames(m) <- x$heatmap.label
  d <- dist(m)
  h <- hclust(d)
  c <- melt(cutree(h, k = ifelse(nrow(m) == 2, 1, 2 )  ))
  colnames(c)[1] <- "cluster"
  c$heatmap.label  <- rownames(c)
  c
}) )

tRNAab.tsne <- merge(tRNAab.tsne, assignment, by = c("heatmap.label","aa"))
tRNAab.tsne <- merge(tRNAab.tsne, TAB$aminoacid.properties, by = "aa")
write.table(tRNAab.tsne, file="results/tRNA.abundance/tRNAab.tSNE_clusters.txt", sep="\t", quote=F, row.names=F)

tRNAab.tsne <- read.table("results/tRNA.abundance/tRNAab.tSNE_clusters.txt", sep="\t", header=1)


# [] tsne, colored by aa -----
# plot the 2 dimensions and show how anticodon are arranged with regard to the amino acid they code for
g <- ggplot( data = tRNAab.tsne, aes(x=Dim1, y=Dim2)) + 
  #geom_abline( intercept= 0, slope=1, color="gray") +
  #geom_line( aes(group= paste(cluster, AA), color = AA ), lty=2 ) +
  geom_point( data = subset(tRNAab.tsne, nTE.O != "O"), color="black", size=6, pch = 21, fill="white" ) +
  geom_point( aes(color=factor(cluster)), size=4 ) +
  geom_text(aes(label=AA),size=3, vjust=2) + coord_equal() +
  labs(title="multidimensional reduction of\nsynonymous anticodons' response to stress") +
  theme_gul + 
  guides( color = guide_legend(ncol = 2) )

ggsave(g,filename = paste0(PATH,"part1/tRNAab.tSNE_clusters.pdf"),
       dpi = 250, width=8.25, height=6.72)

# [] tsne, split by topology -----
# ---- plot the 2 dimensions and show how anticodon are arranged with regard to the amino acid they code for
  g <- ggplot( data = tRNAab.tsne, aes(x=Dim1, y=Dim2)) + 
    geom_line( aes(group= paste(cluster, AA) ), lty=2 ) +
    geom_point( data = subset(tRNAab.tsne, nTE.O != "O"), color="black", size=4.8 ) +
    geom_point( aes(color=topology), size=4 ) +
    geom_text(aes(label=AA),size=3, vjust=1.6) + 
    labs(title="multidimensional reduction of\nsynonymous anticodons' response to stress") +
    facet_wrap( ~ topology) +
    theme_gul + guides( color = guide_legend(ncol = 2), guide.position="bottom" ) + theme(legend.position="bottom")


# [] tsne, split by metabolic cost and topology -----
# ---- plot the 2 dimensions and show how anticodon are arranged with regard to the amino acid they code for
  # ggplot( data = tRNAab.tsne, aes(x=Dim1, y=Dim2)) + 
  #   geom_line( aes(group= paste(cluster, aa) ), lty=2 ) +
  #   geom_point( data = subset(tRNAab.tsne, nTE.O != "O"), color="black", size=4.8 ) +
  #   geom_point( aes(color=metab.cost.ferm), size=4 ) +
  #   geom_text(aes(label=aa),size=3, vjust=1.6) + 
  #   labs(title="multidimensional reduction of\nsynonymous anticodons' response to stress") +
  #   facet_wrap( ~ topology) +
  #   theme_gul + guides( color = guide_legend(ncol = 2), guide.position="bottom" ) + theme(legend.position="bottom")
  
# Redundancy analysis OR Multiple Correspondance Analysis



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------# Global tRNA response #-------------------------------------------------------------------------------------------------#
#-------------------------------#     Heat maps     #----------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# [] pheatmaps ------ 
pdf(paste0(PATH,"part1/tRNAab--x.time_y.foldchange_z.anticodon.per.stress.pdf"), width = 4, height = 7)
require(pheatmap)
require(RColorBrewer)
require(classInt)
  gulmap(M.tRNAab$all.together, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols = F, 
       annotation = rosetta.anticodons[,-c(1:3,5,7)], main="stress responses",cellwidth = 9, cellheight = 9 )

  gulmap(M.tRNAab$all.together.2, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols = F, 
       annotation = rosetta.anticodons[,-c(1:3,5,7)], main="stress responses",cellwidth = 9, cellheight = 9 )

  gulmap(M.tRNAab$diauxic, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols=F, 
           annotation = rosetta.anticodons[,-c(1:3,5,7)], main="diauxic shift",cellwidth = 9, cellheight = 9 )
  
  gulmap(M.tRNAab$osm, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols=F, 
           annotation = rosetta.anticodons[,-c(1:3,5,7)], main="osmotic shock",cellwidth = 9, cellheight = 9  )
  
  gulmap(M.tRNAab$ox, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols=F, 
           annotation = rosetta.anticodons[,-c(1:3,5,7)], main="oxidative stress",cellwidth = 9, cellheight = 9  )
  
  gulmap(M.tRNAab$temp, clustering_method = "complete",scale = "none",cluster_rows = T, cluster_cols=F, 
           annotation = rosetta.anticodons[,-c(1:3,5,7)], main="temperature stress",cellwidth = 9, cellheight = 9  )
dev.off()



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------# Global tRNA response #-------------------------------------------------------------------------------------------------#
#-------------------------------#     Fold change    #---------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# [] Quality assessment ----
# Warning: slow to run, plus needs a big display
        dodge = position_dodge(width=0.9)
        g <- ggplot( data = subset(marc.data), aes(x=factor(time), y=average, group = anticodon, fill=aa) ) + 
          geom_bar(stat = "identity", position= position_dodge(), width=0.4 ) + 
          geom_hline(yintercept=1, color="gray40") +
          geom_errorbar( aes(ymin = average - std, ymax= average + std), position=dodge, width=0.1, size=0.3, color="gray40" ) +
          scale_y_continuous( labels = prettyNum ) +
          facet_grid( anticodon ~ experiment , scales = "free") + 
          labs( x="time", y="fold change in tRNA abundance") +
          theme_gul + 
          theme(legend.position="none", strip.text.y=element_text(angle = 0), axis.text.x=element_text( size=rel(.65), angle=0 ),  axis.text.y=element_text( size=rel(.65), angle=0 ) )
        
        ggsave(plot = g, filename = paste0(PATH,"part1/tRNAab--x.time_y.foldchange_z.experiment_t.anticodon.png"), dpi = 250, 
               width = 4.45*1.15, height = 10*1.15 )

ggplot( data = marc.data, aes(x=time, y=average) ) +
  geom_line( aes( group = anticodon ) ) +
  geom_errorbar( aes(ymin = average - std, ymax= average + std), position=dodge, width=0.1, size=0.3, color="gray40" ) +
  geom_point() + 
  facet_wrap( ~ experiment )



#### [] Outlier detection ####
library(fitdistrplus)
library("extremevalues")
require(gridExtra)
#test <- dlply(marc.data, .(experiment,time), summarise, fit =  fitdist(average, "lnorm") )

boxcox_skewness_minimisation(data = marc.data, v = "average")
marc.data$boxcox.average <- boxcoxtransform(x = marc.data$average, lambda = 0.2)
hist(marc.data$boxcox.average)

# Fitting mean change in tRNA abundance
# "norm", "lnorm", "pois", "exp", "gamma", "nbinom", "geom", "beta", "unif" and "logis"
fitln    <- fitdist(subset(marc.data, average < max(average) & average > min(average) )$average, "lnorm", method='mle')
fitexp   <- fitdist(subset(marc.data, average < max(average) & average > min(average) )$average, "exp", method='mle')
fitnorm  <- fitdist(subset(marc.data, boxcox.average < max(average) & average > min(average) )$average, "norm", method='mle')
fitgamma <- fitdist(subset(marc.data, average < max(average) & average > min(average) )$average, "gamma", method='mle')
fitweibull    <- fitdist(subset(marc.data, average < max(average) & average > min(average) )$average, "weibull", method='mle')
gofstat(list(fitln, fitexp,fitgamma, fitweibull), fitnames=c("lognormal", "exponential","gamma","weibull"))
qqcomp(list(fitln, fitexp, fitgamma, fitweibull), legendtext=c("lognormal", "exponential","gamma","weibull"))
ppcomp(list(fitln, fitexp, fitgamma, fitweibull), legendtext=c("lognormal", "exponential","gamma","weibull"))

plot(fitln)
plot(fitexp)
plot(fitgamma)
plot(fitnorm)
plot(fitweibull)

gofstat(fitln)$chisqpvalue
gofstat(fitexp)$chisqpvalue
gofstat(fitgamma)$chisqpvalue
gofstat(fitweibull)$chisqpvalue

# reviens
fit.sample.disp <- fitdist(subset(marc.data, average < max(average) & average > min(average) )$sample.disp, "lnorm")
fit.sample.disp$gof <- gofstat( fit.sample.disp)
fit.sample.disp$gof$chisqpvalue
plot(fit.sample.disp)

L <- getOutliers(marc.data$sample.disp, distribution="lognormal", method = "II")
L <- getOutliers(marc.data$CE, distribution="lognormal", method = "II")
#L <- getOutliers(marc.data$sample.disp, distribution="weibull")
#outlierPlot(marc.data$average,L,mode="qq")
outlierPlot(marc.data$sample.disp,L,mode="qq")
outlierPlot(marc.data$CE,L,mode="qq")
outlierPlot(marc.data$sample.disp,L,mode="residual")

Lav <- getOutliers(marc.data$average, distribution="lognormal", method = "II")
outlierPlot(marc.data$average,Lav,mode="qq")
outlierPlot(marc.data$average,Lav,mode="residual")

Lstd <- getOutliers(marc.data$std, distribution="lognormal", method = "II")
outlierPlot(marc.data$std,Lstd,mode="qq")
outlierPlot(marc.data$std,Lstd,mode="residual")



tRNA.qqplots <- dlply( marc.data, .(experiment, time), function(x){ gulqqFitPlot(df=x, distribution = "lognormal", ref="label")}  )

pdf( paste0(PATH,"part1/tRNAab--outliers.pdf"), useDingbats = F, height=27, width=15)
  grid.arrange( tRNA.qqplots[[1]], 
                tRNA.qqplots[[2]], 
                tRNA.qqplots[[3]], 
                tRNA.qqplots[[4]], 
                tRNA.qqplots[[5]], 
                tRNA.qqplots[[6]], 
                tRNA.qqplots[[7]], 
                tRNA.qqplots[[8]], 
                tRNA.qqplots[[9]],
                ncol=3 ) # 
dev.off()

g.biological.outliers <- gulqqFitPlot(df=marc.data, y = "average", distribution = "lognormal", ref = "label")
colnames(g.biological.outliers$data)[grep("outlier", colnames(g.biological.outliers$data))]  <- "biological.outlier"

# plot biological outliers all together
pdf( paste0(PATH,"part1/tRNAab--outliers_alltogether.pdf"), useDingbats = F, height=7.5, width=7.5)
  g.biological.outliers %+% 
    geom_point( data = subset(g.biological.outliers$data, biological.outlier == T ), 
                aes(fill= paste(experiment, paste("(",time," min)",sep="") )), pch=21,size=3.5 ) + 
    scale_fill_brewer(guide= guide_legend(title="stress condition", nrow=2),palette="Set2") +
    theme(legend.position="bottom")
dev.off()

g.measurement.outliers <- gulqqFitPlot(df=marc.data, y = "sample.disp", distribution = "lognormal", ref = "label")
colnames(g.measurement.outliers$data)[grep("outlier", colnames(g.measurement.outliers$data))] <- "measurement.outlier"

g.measurement.outliers %+% 
  geom_point( data = subset(g.measurement.outliers$data, measurement.outlier == T ), 
              aes(fill= paste(experiment, paste("(",time," min)",sep="") )), pch=21,size=3.5 ) + 
  scale_fill_brewer(guide= guide_legend(title="stress condition", nrow=2),palette="Set2") +
  theme(legend.position="bottom")


# do not categorise measurements with very low noise as "outliers" (these are the most reproducible)
marc.data.outliers <- merge(g.biological.outliers$data, subset(g.measurement.outliers$data, select=c(anticodon,experiment,time,measurement.outlier)), by=c("anticodon","experiment","time"))
marc.data.outliers$measurement.outlier[ marc.data.outliers$measurement.outlier == T & marc.data.outliers$std < mean(marc.data.outliers$std) ] <- F
marc.data.outliers$outlier <- paste(marc.data.outliers$biological.outlier, marc.data.outliers$measurement.outlier)
marc.data.outliers$outlier <- factor(marc.data.outliers$outlier, labels=c("inlier","measurement outlier?","biological outlier?"))

pdf(  paste0(PATH,"part1/tRNAab--outliers_counts.pdf") )
require(RColorBrewer)
pheatmap(cluster_rows = F, color = brewer.pal(n = 7, "Greys"), treeheight_col = 10, annotation = annotations,
  table(subset(marc.data.outliers, outlier !="inlier")$outlier, as.character(subset(marc.data.outliers, outlier !="inlier")$label) )[2:3,], 
         cellwidth=10, cellheight=10, )
dev.off()


require(scales)
require(RColorBrewer)
require(gridExtra)

g1 <- ggplot( data = marc.data.outliers, aes(x=average, y=sample.disp) ) + 
  geom_point( aes(fill=outlier), pch=21, size=4) + 
  scale_fill_manual(values=c("white", "red", "#0571B0","red4")) +
  stat_smooth(method="lm", color="gray50") +
  labs(x="Sample Mean\n(3 replicates)",y="Sample Index of Dispersion\n(3 replicates)") +
  theme(legend.position="none") +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  

g2 <- ggplot( data = marc.data.outliers, aes(x=cross.cond.disp, y=sample.disp) ) + 
  geom_point( aes(fill=outlier), pch=21, size=4) + 
  scale_fill_manual(values=c("white", "red", "#0571B0","red4")) +
  stat_smooth(method="lm", color="gray50") +
  labs(x="Cross-condition Index of Dispersion\n(4x3 conditions)",y="Sample Index of Dispersion\n(3 replicates)") +
  theme(legend.position="none") +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  

pdf(  paste0(PATH,"part1/tRNAab--outliers_sampling.pdf"), height= 5, width=10, useDingbats=F )
grid.arrange(g1,g2,nrow=1 )
dev.off()



### GC content of anticodon
ggplot( data = subset( anticodon.master.table, time > 0 ) ) + geom_boxplot( aes(x= factor(GC.anti), y= log2.foldchange )  ) +
  facet_grid( experiment ~ time)



# [] description change in free tRNA concentration (20min)----------------- #
g.a <- ggplot(data=  data.20, 
            aes(y=free.tRNAab/100000)) +
  geom_boxplot( aes(x= as.numeric(stress), group = stress), 
                notch=T, fill="white", color = "gray80", outlier.size = 0) + 
  geom_point( aes(x = jitter.x, color=CAI.O), size=3.75, alpha=0.8 ) +
  geom_line( aes(x= jitter.x ,group = anticodon, color = CAI.O), size=0.5, alpha=0.15 ) +
  geom_text(data = subset(data.20, free.tRNAab > 2.5e5 | (CAI.O == "N" & free.tRNAab > 1.25e5)), 
            aes(x = jitter.x, label = anticodon), size = 2, vjust = -1.25 ) +
  scale_color_manual(values=c( "skyblue3","gray40","orange" ) ) +
  scale_x_continuous( labels = levels(data.20$stress), breaks = 1:5 ) +
  scale_y_continuous( labels = prettyNum ) +
  labs (x="\nstress", y="Free abundance/tRNA\n(x 100,000)\n", title="Abundance of free tRNAs 20min after stress induction") + 
  theme_gul + theme( legend.position ="none" )


# [] description change in free tRNA concentration (120min)----------------- #
g.b <- ggplot(data=  data.120, 
       aes(y=free.tRNAab/100000)) +
  geom_boxplot( aes(x= as.numeric(stress), group = stress), 
                notch=T, fill="white", color = "gray80", outlier.size = 0) + 
  geom_point( aes(x = jitter.x, color=CAI.O), size=3.75, alpha=0.8 ) +
  geom_line( aes(x= jitter.x ,group = anticodon, color = CAI.O), size=0.5, alpha=0.15 ) +
  geom_text(data = subset(data.120, free.tRNAab > 2.5e5 | (CAI.O == "N" & free.tRNAab > 1.25e5)), 
            aes(x = jitter.x, label = anticodon), size = 2, vjust = -1.25 ) +
  scale_color_manual(values=c( "skyblue3","gray40","orange" ) ) +
  scale_x_continuous( labels = levels(data.120$stress), breaks = 1:5 ) +
  scale_y_continuous( labels = prettyNum ) +
  labs (x="\nstress", y="Free abundance/tRNA\n(x 100,000)\n", title="Abundance of free tRNAs 120min after stress induction") + 
  theme_gul + theme( legend.position ="none" )

png( paste0(PATH,"part1/free.tRNAab--x.stress_y.free.tRNAab.png"), res = 250, width = 2*6.75, height = 6.75, units = "in" )
grid.arrange( g.a, g.b, ncol = 2)
dev.off()



# [] free tRNA abundance in function of CAI optimality and stress  ----------------- 
g <- ggplot(data=  ddply(data.120, .(stress, CAI.O), summarise, free.tRNAab=sum(free.tRNAab) ), 
       aes(x=stress, y=free.tRNAab/100000, fill=CAI.O, group=CAI.O)) +
  geom_bar(stat="identity", position="dodge" ) + 
  scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
  labs (x="\nstress", y="Free tRNAs\n(x 100,000)\n", title="Abundance of free tRNAs\n120min after stress induction") + 
  theme_gul


# [] free tRNA abundance in function of total tRNA abundance and stress  ----------------- 
require(scales)
g <- ggplot(data= anticodon.master.table, aes( x = total.tRNAab.marc, y = free.tRNAab) )  +
  geom_abline( slope = 1, intercept = 0, color = "gray", size = 2, alpha=0.25) +
  stat_smooth( method = "lm", color = "skyblue2", size = 1.5, se = F ) +
  geom_point( aes( fill = tGCN <= 3), pch = 21) +
  geom_text( data = subset(anticodon.master.table, free.tRNAab > total.tRNAab.marc), size = 3, vjust = -0.5, aes( label = anticodon) ) + 
  facet_wrap(time ~ experiment) +
  scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 4),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  + coord_fixed() +
  scale_fill_manual( values = c("black","white") ) +
  labs( x="tRNA abundance", y = "free tRNA abundance") +
  theme_gul


ggsave(plot = g, filename = paste0(PATH,"part1/free.tRNAab--x.total.tRNAab_y.free.tRNAab_z.stress.png"), 
       dpi = 250, width=5.82, height=5.82)






# ggplot(data=  data , 
#        aes(x=stress, y=free.tRNAab/100000)) +
#   geom_boxplot(notch=T) + 
#   #scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
#   labs (x="\nstress", y="Free tRNAs\n(x 100,000)\n", title="Abundance of free tRNAs\n120min after stress induction") + 
#   theme_gul



# [] description change in total tRNA concentration ----------------- 
g <- ggplot(data=  data.120, 
       aes(x=stress, y=total.tRNAab/100000)) +
  geom_boxplot(notch=T) + 
#  geom_hex() +
  labs (x="\nstress", y="Total abundance/tRNA\n(x 100,000)\n", title="Abundance of total tRNAs\n120min after stress induction") + 
  theme_gul

# [] description change in total tRNA abundance ----------------- 
g <- ggplot(data=  ddply(data.120, .(stress, CAI.O), summarise, total.tRNAab=sum(total.tRNAab) ), 
       aes(x=stress, y=total.tRNAab/100000,  fill=CAI.O)) +
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
  labs (x="\nstress", y="Total tRNAs\n(x 100,000)\n", title="Abundance of tRNAs 120min after stress induction\n") + 
  theme_gul

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.stress_y.total.tRNAab_z.CAI.O.pdf"), 
       dpi = 250, width=8.18, height=7.24)


g <- ggplot(data=  data.120, 
       aes(x=stress, y=total.tRNAab/100000)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
  labs (x="\nstress", y="Total tRNAs\n(x 100,000)\n", title="Abundance of tRNAs 120min after stress induction") + 
  theme_gul

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.stress.120min_y.total.tRNAab.pdf"), 
       dpi = 250, width=6.24, height=4.61)

#-----------------------------------------------------------------#
       
# [] description change in bound tRNA abundance ----------------- #
dodge = position_dodge(width=0.9)
g <- ggplot(data=  subset( anticodon.master.table, (time == 120 ) | (experiment=="normal") ), 
            aes(x=experiment, y=(total.tRNAab.marc - free.tRNAab)/100000,  fill=CAI.O)) +
  geom_bar(stat="identity", aes(group=anticodon), position=position_dodge() ) + 
  geom_text( aes(label = anticodon, group=anticodon, ymax = (total.tRNAab.marc - free.tRNAab)/100000), 
             position=dodge, size=2.25, hjust=-0.2, color="gray30", angle = 90 ) +
  scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
  labs (x="\nstress", y="Bound tRNAs\n(x 100,000)\n", 
        title="Abundance of bound tRNAs 120min after stress induction\n",
        fill="anticodon") + 
  facet_wrap( ~ aa) +
  theme_gul + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="bottom") 

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.stress_y.bound.tRNAab_z.CAI.O.t_aa.png"), 
       dpi=250, width=9.1, height=8.88)
#-----------------------------------------------------------------#



# [] description change in bound tRNA abundance ----------------- #
g <- ggplot(data=  data.120, 
            aes(x=stress, y=(1 - free.tRNAab/total.tRNAab),  fill=CAI.O)) +
  geom_bar(stat="identity", position="dodge") + 
  scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
  labs (x="\nstress", y="Fraction of bound tRNAs\n", title="Abundance of tRNAs 120min after stress induction\n") + 
  facet_wrap( ~ aa) + coord_fixed() +
  theme_gul + theme(axis.text.x=element_text(angle=90, hjust=1))


#-----------------------------------------------------------------#


# [] description change in free tRNA concentration ----------------- #
g <- ggplot(data= anticodon.master.table, 
       aes(x=stress, y=total.tRNAab/100000, fill = (anticodon == "CAA")) ) +
  geom_bar(stat="identity") + 
  scale_fill_manual( values = c("gray40", "orange") ) +
  labs (x="\nstress", y="Total tRNAs\n(x 100,000)\n", 
        fill = "CAA-tRNA",
        title="Total abundance of tRNAs") + 
  theme_gul

total.tRNAab <- ddply(anticodon.master.table, .(experiment,time), summarise, total.tRNAab = sum(total.tRNAab) )
total.tRNAab$label <- paste(total.tRNAab$experiment, total.tRNAab$time, sep="\n")
total.tRNAab$experiment <- factor(total.tRNAab$experiment, levels = c("diauxic","ox","osm","temp"), labels=c("diauxic\nshift","oxidative\nstress","osmotic\nstress","temperature\nstress"))


require(gridExtra)

g1 <- ggplot( data=subset(total.tRNAab,time!=0), aes( x = time, y= total.tRNAab )) + 
  geom_hline(yintercept= subset(total.tRNAab, time==0)$total.tRNAab , lty=3) +
  geom_bar(stat="identity") + 
  scale_x_continuous(breaks=c(20,120)) + 
  scale_y_continuous( label=scientific_10) +
  facet_wrap( ~ experiment, scales = "free_x", nrow=1 ) + labs(y="total tRNA abundance")

pie.chart.data <- subset(anticodon.master.table, time !=0)
pie.chart.data <- ddply(pie.chart.data, .(experiment,time), mutate, rank = rank(foldchange), alpha= rescale( sqrt(abs(log2.foldchange)),to = c(0.75,1)) )
pie.chart.data$label <- paste(pie.chart.data$experiment, pie.chart.data$time, sep="\n")
pie.chart.data$label <- factor(pie.chart.data$label, levels=  paste( rep(c("diauxic","ox","osm","temp"), times = 2 ), rep(c(20,120), each = 4 ), sep="\n"),
                               labels = paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), rep(c("20min","120min"), each = 4 ), sep="\n") )
pie.chart.data$behaviour <- ifelse( pie.chart.data$log2.foldchange > log2(1.5), "up", ifelse( pie.chart.data$log2.foldchange < log2(0.5),"down","steady") )

pie.props.up   <- ddply(pie.chart.data, .(experiment,time), summarise, prop = round( sum(log2.foldchange>= log2(1.5)) / length(anticodon) ,2))
pie.props.down <- ddply(pie.chart.data, .(experiment,time), summarise, prop = round( sum(log2.foldchange<= log2(0.5)) / length(anticodon) ,2))

g2 <- ggplot( data = pie.chart.data, aes(x=factor(1), group = rank, color = behaviour, fill = behaviour, alpha=alpha),  ) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") + 
  #geom_text( data = pie.props.up, aes(x=-Inf, y=+Inf, label= prop) ) + #percent_format()(pie.props.up$prop)
  labs(x="",y="") + 
  scale_alpha(range = c(0.5,1)) + 
  scale_fill_manual(values=c("#CC3333","gray","#336699")) + 
  scale_color_manual(values=c("#CC3333","gray","#336699")) + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none") +
  facet_wrap( ~ label , nrow =2) 


grid.arrange(g1,g2, nrow=2, heights=c(1/3,2/3) )
dev.copy2pdf(device = quartz, file = paste0(PATH,"part1/fig_results_total.tRNAab.pdf") )


# Voronoi tree map
require(treemap)
treemap(subset(anticodon.master.table, experiment=="normal" & time == 0), 
        index=c("label"), 
        vSize="total.tRNAab", 
        vColor="CAI.O",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=T,
        range=c(0.2,0.7),
        palette = c("skyblue","gray","orange"),
        type="categorical",
        position.legend = "none"
        )
dev.copy2pdf(device = quartz, file = paste0(PATH,"part1/tRNAab_normal.voronoi.pdf") )



treemap(subset(anticodon.master.table, experiment=="osm" & time == 20), 
        index=c("label"), 
        vSize="total.tRNAab", 
        vColor="CAI.O",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=T,
        range=c(0.2,0.7),
        palette = c("skyblue","gray","orange"),
        type="categorical",
        position.legend = "none"
)


# fold changes

treemap(subset(anticodon.master.table, experiment=="osm" & time == 20), 
        index=c("label"), 
        vSize="foldchange", 
        vColor="CAI.O",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=T,
        range=c(0.2,0.7),
        palette = c("skyblue","gray","orange"),
        type="categorical",
        position.legend = "none"
)


treemap(subset(anticodon.master.table, experiment=="osm" & time == 120), 
        index=c("label"), 
        vSize="foldchange", 
        vColor="CAI.O",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=T,
        range=c(0.2,0.7),
        palette = c("skyblue","gray","orange"),
        type="categorical",
        position.legend = "none"
)

treemap(subset(anticodon.master.table, experiment=="ox" & time == 120), 
        index=c("label"), 
        vSize="foldchange", 
        vColor="CAI.O",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=T,
        range=c(0.2,0.7),
        palette = c("skyblue","gray","orange"),
        type="categorical",
        position.legend = "none"
)

g <- ggplot(data= subset(anticodon.master.table, anticodon == "GCC"), 
       aes(x=stress, y=total.tRNAab/100000, fill = (anticodon == "GCC")) ) +
  geom_bar(stat="identity") + 
  scale_fill_manual( values = c("gray60") ) +
  labs (x="\nstress", y="Total tRNAs\n(x 100,000)\n", 
        fill = "tRNA(Gly-GCC)",
        title="Total abundance of tRNA(GCC)") + 
  theme_gul

ggsave(g, filename = paste0(PATH,"part1/tRNApool--x.stress_y.total.tRNAab(GCC).pdf"), dpi =250, useDingbats=F )

#-----------------------------------------------------------------#
# ---- changes of tRNA abundance over time, grouped per amino acid and experiment
  # g <- ggplot( data = subset(tRNAab_quantified.foldchange_long, experiment != "normal" ), 
  #              aes( x = as.numeric(time), y = total.tRNAab/100000, group = anticodon ) ) + 
  #   geom_line(aes(color=CAI.O), size=1) +
  #   scale_y_log10( labels = prettyNum ) +
  #   facet_grid( experiment ~ aa) +
  #   labs( x= "time\n(0,20,60,120min)", y="tRNA abundance\n(x100,000)", color="recognised codon(s)" ) +
  #   scale_color_manual( values=c("skyblue3","gray","orange") ) +
  #   theme_gul + theme( axis.text.x = element_blank( ), legend.position="bottom" ) 
  # 
  # ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.time_y.total.tRNAab_z.aa.t_stress.png"), 
  #        dpi=250, width=9.1, height=6.31)

#-----------------------------------------------------------------#

# ---- variation in free tRNA abundance (compared to normal conditions) in function of total tRNA abudance
# ---- no apparent scaling (low or high abundance doesn't associate with delta free.tRNAab)

R2.data <- ddply( subset(anticodon.master.table, experiment != "normal" ), .(experiment, time), function(x){ y <- lm( data = x, formula = log2(ratio_free.tRNAab)  ~ total.tRNAab );
                                                                 as.character(as.expression(substitute(italic(R)^2~"="~r2 , list( r2= format(summary(y)$r.squared, digits = 2) ))))
} )


g <- ggplot( data = subset(anticodon.master.table, experiment != "normal" ), 
        aes( x = total.tRNAab, y = log2( ratio_free.tRNAab) ) ) + 
  geom_hline(yintercept=0, lty=1, size=0.5, color="gray70") +
  stat_smooth(method="lm", color = "gray30", size = 0.75) +
  geom_point(aes(color=CAI.O), size = 2.5 ) +
  geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  facet_grid( experiment ~ time ) +
  labs( color="recognised codon(s)", 
        y = "variation in free tRNA abundance (log2)", 
        x = "tRNA abundance" ) +
  scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  scale_color_manual( values=c("skyblue3","gray","orange") ) +
  theme_gul + theme( legend.position="bottom" ) 

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.total.tRNAab_y.free.tRNAab_z.time_t.experiment.pdf"),
       dpi = 250,
       width = 4.33,
       height = 6.47, useDingbats=FALSE 
)
#-----------------------------------------------------------------#








# QUESTION # How do tRNA vary in abundance over time for distinct stresses?
data.mean <- ddply( subset(tRNAab_quantified.foldchange_long, experiment != "normal"), .(experiment, time), function(x){ log2.foldchange = mean(x$log2.foldchange) } )
colnames(data.mean)[3] <- "log2.foldchange"

g <- ggplot( data = subset(tRNAab_quantified.foldchange_long, experiment != "normal"), aes( x = (time), y = log2.foldchange ) ) + 
  geom_hline(yintercept=0, lty=1, color="gray80",size=2.5) +
  geom_line(aes(color=tGCN, group = anticodon)) +
  geom_line(data = data.mean, color = "red", aes(group=NA), size = 2 ) +
  geom_text( data =  subset(tRNAab_quantified.foldchange_long, experiment != "normal" & time == "t=120"), 
             aes(label = anticodon, color = tGCN), hjust = -0.25, size = 2, alpha=0.75 ) +
  facet_wrap( ~ experiment ) +
  labs( x= "time", y="tRNA abundance\n(fold-change log2)", color="tRNA Gene Copy Number" ) +
  theme_gul + theme( legend.position="none" )

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.time_y.log2.foldchange_z.stress.pdf"), 
       dpi=250, width=5.92, height=5.89, useDingbats=FALSE  )

rm(data.mean)

#-----------------------------------------------------------------------------------------------------------------------------------------------#



# What is the variation in tRNAab over time for individual amino acids?

# [] log2 fold change 
g <- ggplot( data = subset(tRNAab_quantified.foldchange_long,experiment!="normal"), 
             aes( x = as.numeric(as.factor(time)), y = log2.foldchange, group = anticodon ) ) + 
  geom_rect( xmin = -Inf, xmax = +Inf, ymin = -1, ymax = 1, fill = "gray80", alpha = 0.7 ) +
  geom_line(aes(color=CAI.O), size=0.75) +
  geom_hline(yintercept=0, lty=1, size=0.5) +
  geom_text( data = subset( tRNAab_quantified.foldchange_long,experiment!="normal" & time == 120), 
             x = 5, hjust=0.15, aes(label = anticodon, color = CAI.O), size = 1.75 ) + 
  facet_grid( experiment ~ aa) +
  labs( x= "time\n(t=0,20,60,120min)", y="tRNA abundance (fold-change log2)\n", color="recognised codon(s)" ) +
  scale_color_manual( values=c("skyblue3","gray40","orange") ) +
  scale_x_continuous( breaks = 1:7, limits = c(1,7) ) + 
  theme_gul + theme( axis.text.x = element_blank( ), 
                     axis.ticks.x = element_blank(),
                     legend.position="bottom" ) 

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.time_y.log2.foldchange.abundance_z.aa_t.stress.pdf"), 
       dpi=250,  width=9.1, height=6.31, useDingbats=FALSE )
#-----------------------------------------------------------------#


# [] changes in tRNA abundance per amino acid ------
# (assuming 3.3M tRNA molecules in normal conditions, using the tGCN to get an estimate/anticodon during normal conditions, and adjusting with fold change for stress conditions)
g <- ggplot( data = subset(tRNAab_quantified.foldchange_long,experiment!="normal" ), 
             aes( x = as.numeric(as.factor(time)), y = total.tRNAab, group = anticodon ) ) + 
  geom_line(aes(color=CAI.O), size=0.75) +
  geom_text( data = subset( tRNAab_quantified.foldchange_long,experiment!="normal" & time == 120), 
             x = 5, hjust=0.15, aes(label = anticodon, color = CAI.O), size = 1.75 ) + 
  #geom_hline(yintercept=0, lty=1, size=0.5) +
  facet_grid( experiment ~ aa) +
  labs( x= "time\n(t=0,20,60,120min)", y="tRNA abundance\n", color="recognised codon(s)" ) +
  scale_color_manual( values=c("skyblue3","gray40","orange") ) +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x, n = 2),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  scale_x_continuous( breaks = 1:7, limits = c(1,7) ) + 
  theme_gul + theme( axis.text.x = element_blank( ), 
                     axis.ticks.x = element_blank(),
                     legend.position="bottom" ) 

ggsave(plot = g, filename = paste0(PATH,"part1/tRNApool--x.time_y.total.tRNAab_z.aa_t.stress.pdf"), 
       dpi=250, width=9.1, height=6.31, useDingbats=FALSE )



#-----------------------------------------------------------------#


# QUESTION # Do the optimality of codons recognised by anticodons associate with specific changes in tRNA abundance?

        # # [] t=20min
        # g3 <- ggplot( data = subset(tRNAab_quantified.foldchange_long,time=="t=20" & anticodon != "CAA"), aes( x = CAI.O, y = log2.foldchange ) ) + 
        #   geom_hline(yintercept=0, lty=3) +
        #   geom_boxplot(width=0.5, notch=F) +
        #   facet_wrap( ~ experiment ) +
        #   labs( y="tRNA abundance\n(fold-change log2)", 
        #         x="\nType of codons recognised by tRNAs",
        #         title="20min after induction\n") +
        #   theme_gul + theme( legend.position="bottom" ) 
        # 
        # # t=60min
        # g4 <- ggplot( data = subset(tRNAab_quantified.foldchange_long,time==60), aes( x = CAI.O, y = log2.foldchange ) ) + 
        #   geom_hline(yintercept=0, lty=3) +
        #   geom_boxplot(width=0.5, notch=F) +
        #   facet_wrap( ~ experiment ) +
        #   labs( y="tRNA abundance\n(fold-change log2)", 
        #         x="\nType of codons recognised by tRNAs",
        #         title="60min after induction\n") +
        #   theme_gul + theme( legend.position="bottom" ) 
        # 
        # # t=120min
        # g5 <- ggplot( data = subset(tRNAab_quantified.foldchange_long,time==120), aes( x = CAI.O, y = log2.foldchange ) ) + 
        #   geom_hline(yintercept=0, lty=3) +
        #   geom_boxplot(width=0.5, notch=F) +
        #   facet_wrap( ~ experiment ) +
        #   labs( y="tRNA abundance\n(fold-change log2)", 
        #         x="\nType of codons recognised by tRNAs",
        #         title="120min after induction\n") +
        #   theme_gul + theme( legend.position="bottom" ) 
        # 
        # 
        # ggsave(plot = g3, filename = paste0(PATH,"tRNA.abundance/foldchange.abundance.vs.stress.vs.optimality.t020.pdf"), dpi = 250 )
        # ggsave(plot = g4, filename = paste0(PATH,"tRNA.abundance/foldchange.abundance.vs.stress.vs.optimality.t060.pdf"), dpi = 250 )
        # ggsave(plot = g5, filename = paste0(PATH,"tRNA.abundance/foldchange.abundance.vs.stress.vs.optimality.t120.pdf"), dpi = 250 )
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------# Global tRNA response #-------------------------------------------------------------------------------------------------#
#-----------------------------# * The case of Met-tRNAs #--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

t0 <- subset(marc.data, time == 120) # the time here is not important, all values will be overriden
t0$time <- 0
t0$average <- 1
t0$std <- 0
t0$CE <- 0
t0$up <- 0
t0$down <- 0
marc.data <- rbind(t0, marc.data)
marc.data$up <- marc.data$average + marc.data$CE
marc.data$down <- marc.data$average - marc.data$CE
marc.data$up[which(marc.data$up > 9)] <- 9
marc.data$down[which(marc.data$down < 0)]  <- 0

Met.data <-  subset(marc.data, anticodon %in% c("CAU","CAU2")) # CAU2 is now absent...
rm(t0)

# Met.data <- methionine.table
        #   reshape2::dcast( data = melt( subset(marc.data, anticodon %in% c("CAU","CAU2")), id.vars=c("anticodon","experiment","time")) , 
        #                    formula =  experiment + time ~ anticodon + variable, 
        #                    value.var = "value")
Met.data$anticodon <- factor(Met.data$anticodon, levels = c("CAU2","CAU"), labels= c("iMet-tRNA","eMet-tRNA" ))
Met.data$experiment <- factor(Met.data$experiment, levels=c("diauxic","ox","osm","temp"), 
                              labels=c("diauxic\nshift","oxidative\nstress","osmotic\nstress","temperature\nstress"))

# ---- > foldchange elongation/initiation Met-tRNA ----
dodge = position_dodge(width=0.9)
g <- ggplot( data = subset(Met.data,time>0), aes( x=factor(time), y= average, group = anticodon) ) +
  geom_hline( yintercept = 1) +
  geom_bar( stat = "identity", aes(fill = anticodon), position = position_dodge(), width = 0.75 ) +
  geom_errorbar( aes(ymin = average - CE, ymax= average + CE), position=dodge, width=0.1, size=0.3, color="gray20" ) +
  scale_fill_manual( values = c("gray70","gray40")) +
  scale_y_continuous( labels = prettyNum ) +
  facet_wrap( ~ experiment) +
  labs( x = "time", y="fold change") +
  theme_gul + theme(legend.position= "bottom")

ggsave( plot=g, filename = paste0(PATH,"part1/Met-tRNA--x.time_y.foldchange_z.stress.pdf"), 
             dpi=250,  width=5.92, height=5.89)


#reviens le gul!

g <- ggplot(data=subset(Met.data, anticodon=="iMet-tRNA"), aes(x=time,y=average)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_line(aes(group=experiment))+
  geom_point(aes(shape=experiment), pch=21, fill="white",size=3) +
  geom_errorbar( aes(ymin = average - CE, ymax= average + CE), 
                 position=dodge, size=0.3, color="gray40", width=6) +
  labs(y=expression("i-tRNA"^"Met"*" fold change")) +
  facet_wrap( ~ experiment, nrow=1 ) + scale_x_continuous(breaks=c(0,20,60,120))

ggsave( plot=g, filename = paste0(PATH,"part1/fig_results_initiation_tRNA.pdf"), 
        dpi=250, useDingbats=F, width=7.28, height=2.43 )


d.Met <- reshape(  data = melt(Met.data, 
                               id.vars=c("experiment","anticodon","time"), 
                               measure.vars = "average")[,-4], 
                   timevar="anticodon",
                   idvar = c("experiment", "time"), 
                   direction = "wide")

d.Met$i_e <- d.Met$'value.iMet-tRNA' / d.Met$'value.eMet-tRNA'
                   
ggplot( data = d.Met, 
        aes(x= experiment, y = i_e, group=time )) + geom_bar(position = "dodge", stat="identity", aes(fill=time))


# Long to run....
g <- ggplot(data=subset(marc.data), aes(x=time,y=average)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_errorbar( aes(ymin = down, ymax= up), position=dodge, width=0.1, size=0.3, color="gray20" ) +
  geom_line(aes(group=experiment))+
  geom_point(aes(shape=experiment), pch=21, fill="white",size=2) +
  facet_wrap( ~ experiment + anticodon ) + ylim(0,9)
ggsave(plot=g,filename="~/Desktop/measurements.pdf", useDingbats = F, width=12, height=12)


# ------ essential tRNAs -----
essential.single.gene.coded.tRNAs <- c("CCU", "CCG", "CUG","GAG", "CGU")

ggplot(data=subset(marc.data, anticodon %in% essential.single.gene.coded.tRNAs), 
       aes(x=time,y=log2(average))) + 
  geom_hline(yintercept=0, lty=3) +
  geom_line(aes(group=anticodon))+
  geom_point(aes(shape=anticodon, fill=aa), pch=21, size=3) +
#    geom_errorbar( aes(ymin = average - CE, ymax= average + CE), 
#                   position=dodge, size=0.3, color="gray40", width=6) +
  labs(y="fold change\n(log2)", title="essential single-gene encoded tRNAs") +
  scale_x_continuous(breaks=c(0,20,60,120)) +
  facet_wrap( ~ experiment, nrow=1 ) 


# induced tRNAs at ox.20
osm20 <- subset(marc.data, experiment == "osm")
# from heatmap:
osm20.induced.tRNAs <- c("CAA","UAG","CCU","CGU","GUC","CCG","UUG","UGC","AGA","UGA","CUG","UAA","GUA","CCA","GAG")
osm20$induced  <- ifelse( as.character(osm20$anticodon) %in% osm20.induced.tRNAs, T,F  )
osm20$essential  <- ifelse( as.character(osm20$anticodon) %in% c("CCU", "CCG", "CUG","GAG", "CGU"), T,F  )

head(osm20)
m.osm20 <- melt(subset(osm20,anticodon!="CAU2"), id.vars = c("anticodon","experiment","time","aa","induced","essential","CAI.O"), 
                measure.vars = c("average","GC.anti","tGCN"))


ggplot( subset(m.osm20, time >0), aes(x=essential, y=value) ) + 
  geom_dotplot(position = "identity", stackdir="center", binaxis="y", aes(fill=CAI.O)) +
  facet_wrap( time ~ variable, scales = "free_y") + scale_fill_manual(values=c("skyblue3","gray","orange"))


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------# Global tRNA response #-------------------------------------------------------------------------------------------------#
#-----------------------------# The case of Leu-tRNAs(CAA) #--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
CAA.table  <- subset(marc.data, anticodon %in% "CAA")


dodge = position_dodge(width=0.9)
g <- ggplot( data = CAA.table, aes( x=factor(time), y= average, group = anticodon) ) +
  geom_hline( yintercept = 1) +
  geom_bar( stat = "identity", fill = "white", color = "gray30", 
            position = position_dodge(), width = 0.75 ) +
  geom_errorbar( aes(ymin = average - std, ymax= average + std), position=dodge, width=0.1, size=0.3, color="gray20" ) +
  scale_y_continuous( labels = prettyNum ) +
  facet_wrap( ~ experiment) +
  labs( x = "time", y="fold change", title = "Relative abundance of tRNA-Leu(CAA) during stress") +
  theme_gul

ggsave( plot=g, filename = paste0(PATH,"part1/Leu-tRNA(CAA)--x.time_y.foldchange_z.stress.png"), 
             dpi=250,  width=5.92, height=5.89)


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#------------------------------#||| Analysis Part II  #--------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

      
      
      # In this section, the main research objective is to investigate whether or not the changes observed in tRNA abundance affect translation. 
      # Addressing this requires to also account for the changes in the demand (i.e. the transcriptional response). This is what justifies the use of a stochastic model.  
      
      # Are anticodon demand and supply covarying?
      # What is the combined effect of changes in transcription and translation during stress on the production of individual proteins?
      
      
      
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------# Anticodon demand #--------------------------------------------------------------------------------------------------#
#-------------------------------# Relative availability  #-----------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# [] balance: early response to stress -----
data.figure  <- subset(anticodon.master.table)
data.figure$experiment <- factor( data.figure$experiment, levels = c("normal","diauxic","ox","osm","temp"), labels = c("No stress","Diauxic shift","Oxidative stress","Osmotic shock","Temperature stress") )
data.figure$time <- factor( data.figure$time, labels = c("0min","20min","120min"))
#R2.data <- ddply( data.figure, .(experiment, time), function(x){ y <- lm( data = x, formula = total.tRNAab ~ demand.mRNAab );
#                                                                 as.character(as.expression(substitute(italic(R)^2~"="~r2 , list( r2= format(summary(y)$r.squared, digits = 2) ))))
#                                                                           } )

lm_eqn_special = function(df){
  m1 = lm( total.tRNAab ~ demand.mRNAab, df);
  m2 = lm( total.tRNAab ~ demand.mRNAab, subset(df, anticodon != "CAA"));
  pears = cor.test( formula = ~ total.tRNAab + demand.mRNAab, data = subset(df, anticodon !="CAA"), method = "pearson",exact = F);
  line1 <- substitute(italic(r)~"="~pears~" | "~alpha~"="~coeff, 
                      list( pears = format(pears$estimate, digits = 2),
                            coeff = format(summary(m2)$coefficients[[2]], digits = 2)
                            )
                      )
  line2 <- substitute(italic(R)^2~"="~r2~" | "~italic(R)[+CAA]^2~"="~r2.caa,
                       list( r2.caa     = format(summary(m1)$r.squared, digits = 2), 
                             r2         = format(summary(m2)$r.squared, digits = 2)
                             )
                      )
  as.character( as.expression( paste( c('atop(', line1, ',', line2,')'), collapse=" ") ) )
}

R2.data <- ddply( data.figure, .(experiment, time), function(x){ lm_eqn_special(x) })

regression.data <- dlply( data.figure, .(experiment,time), function(x){ y <- data.frame(matrixify(x[,c("anticodon","total.tRNAab","demand.mRNAab")])); lm(total.tRNAab ~ demand.mRNAab, data = y) })

lapply( names(regression.data), function(X){ ggplot( data = with( regression.data[[X]], 
                data.frame( anticodon = names(residuals), residuals = residuals, y = fitted.values ) 
                ), aes(x=rank(y), y= abs(residuals/y) ) ) + geom_line() + labs(title = X) })



plots_residuals <- lapply( names(regression.data), function(X){ ggplot( data = with( regression.data[[X]], 
                                                                  data.frame( anticodon = names(residuals), residuals = residuals, fit = fitted.values ) 
), aes(x=rank(fit), y= abs(residuals) ) ) + 
    geom_line(color="gray") +
    geom_point( color = "gray60" ) + scale_y_sqrt(labels=scientific_10) +
    labs(title = X, y="abs residuals", x="tRNA abundance rank") + geom_text(aes(label=anticodon), size = 3, vjust = -1.4) })


# check that residuals increase with tRNA abudnance to make the point that abundant tRNAs are the most off-balanc
res_ab <- do.call(rbind, lapply( names(regression.data), function(X){  
  anticodons <- names(residuals)
  data = with( regression.data[[X]], data.frame( anticodon = names(residuals), residuals = residuals, fit = fitted.values )) 
  d <- with(data, cor.test( abs(residuals), fit, method="spearman"));
  d <- data.frame( r=round(d$estimate,2), p = d$p.value)
  d$condition <- X
  return(d)
            }) )[-1,] 



#R2.data$caption1 <- paste0("R^2==", R2.data$V1)

# [] Does the total tRNA abundance of synonymous anticodons scale with the anticodon demand? ----
g <- ggplot(  data = data.figure, aes(x= demand.mRNAab/100000, y = total.tRNAab/100000 ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_sqrt( labels = prettyNum, limits =c(0,15), breaks = c(0,1,3,6,9,12,16) ) + coord_fixed() +
  scale_y_sqrt( labels = prettyNum, breaks = c(0,1,3,6,9,12) ) +
  labs( x = expression(atop("anticodon demand",atop("x100,000",""))), 
        y = expression(atop("anticodon supply",atop("tRNA abundance, x100,000",""))), 
        color = "recognised codons",
        title="Balance between anticodon demand and supply") +
  facet_wrap( time ~ experiment, scale="free_x", nrow =2 ) + 
  theme_gul + theme(legend.position="bottom")

ggsave( plot = g, filename = paste0(PATH,"part2/fig_results_balance_pre.pdf"), 
        dpi=250, width = 10.9, height=6.7, useDingbats=FALSE)

#------------------------#

### [] Variation in demand and variation in supply #####
dddd <- subset(anticodon.master.table,time!=0)
dddd$var.demand.up <- (dddd$demand.mRNAab.up - dddd$demand.mRNAab.up.normal)
dddd$var.suppl <- (dddd$total.tRNAab - dddd$total.tRNAab.normal)/dddd$total.tRNAab.normal
  
ggplot(data = dddd, aes(x = demand.mRNAab.up, y = demand.mRNAab.up.normal)) + 
  geom_point() + stat_smooth(method="lm") + facet_grid( experiment ~ time)

ggplot(data = subset(anticodon.master.table, time != 0), aes(x = (demand.mRNAab - demand.mRNAab.normal) , y = total.tRNAab - total.tRNAab.normal )) + 
  geom_point() + stat_smooth(method="lm") + facet_grid( experiment ~ time)


ggplot(data = subset(anticodon.master.table, time != 0), aes(x = (demand.mRNAab) , y = total.tRNAab )) + 
  geom_point() + stat_smooth(method="lm") + facet_grid( experiment ~ time)




#------------------------#

# [] Does the total tRNA abundance of synonymous anticodons scale with the anticodon demand in up regulated genes? ----
data.figure  <- subset(anticodon.master.table, time == 20)
data.figure$experiment <- factor( data.figure$experiment, levels = c("normal","diauxic","ox","osm","temp"), labels = c("No stress","Diauxic shift","Oxidative stress","Osmotic shock","Temperature stress") )
data.figure$time <- factor( data.figure$time, labels = c("20min"))


        lm_eqn_special = function(df, y = "total.tRNAab", x ="demand.mRNAab.up"){
          f <- as.formula( paste(y, x, sep = " ~ ") )
          m1 = lm( formula = f, df);
          m2 = lm( formula = f, subset(df, anticodon != "CAA"));
          spear = cor.test( formula = as.formula( paste( "~", y, "+", x)), data = df, method = "spearman",exact = F);
          eq.1 <- substitute(italic(R)^2~"="~r2 ,          list( r2     = format(summary(m1)$r.squared, digits = 2) ))
          eq.2 <- substitute(italic(R)[-CAA]^2~"="~r2.caa, list( r2.caa = format(summary(m2)$r.squared, digits = 2) ))
          eq.3 <- substitute(italic(rho)~"="~spear, list( spear = format(spear$estimate, digits = 2) ))
          as.character( as.expression( paste( c('atop(', eq.2, ',', eq.3, ')'), collapse=" ") ) )
        }



R2.data <- ddply( data.figure, .(experiment, time), function(c){ lm_eqn_special(df = c) })
#R2.data$caption1 <- paste0("R^2==", R2.data$V1)

g <- ggplot(  data = data.figure, aes(x= demand.mRNAab.up/100000, y = total.tRNAab/100000 ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  labs( x = expression(atop("anticodon demand",atop("x100,000",""))), 
        y = expression(atop("anticodon supply",atop("tRNA abundance, x100,000",""))), 
        color = "recognised codons",
        title=expression(atop("Balance between anticodon demand and supply",atop("in early response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")

ggsave( plot = g, filename = paste0(PATH,"part2/upreg--x.anticodon_demand_y.total.tRNAab_z.stress_t.20min.pdf"), 
        dpi=250, width = 10.9, height=4.35, useDingbats=FALSE)

#-------------------------#




#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# [] balance: late response to stress -----
data.figure  <- subset(anticodon.master.table, time %in% c(0,120))
data.figure$experiment <- factor( data.figure$experiment, levels = c("normal","diauxic","ox","osm","temp"), labels = c("No stress","Diauxic shift","Oxidative stress","Osmotic shock","Temperature stress") )
data.figure$time <- factor( data.figure$time, labels = c("0min","120min"))
                                                                 #y1 <- lm( data = x, formula = total.tRNAab ~ demand.mRNAab );
                                                                 #y2 <- lm( data = subset(x,anticodon!="CAA"), formula = total.tRNAab ~ demand.mRNAab )
                                                                 #c( round( summary(y1)$r.squared, 2), round( summary(y2)$r.squared, 2))
##} )

R2.data <- debug(  ddply( data.figure, .(experiment, time), function(c){ lm_eqn_special(df = c) }) )
#R2.data$caption1 <- paste0("R^2==", R2.data$V1)
#R2.data$caption2 <- paste0("R[-CAA]^2==", R2.data$V2)

# [] Does the total tRNA abundance of synonymous anticodons scale with the anticodon demand? ----
g <- ggplot(  data = data.figure, aes(x= demand.mRNAab/100000, y = total.tRNAab/100000 ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  labs( x = expression(atop("anticodon demand",atop("x100,000",""))), 
        y = expression(atop("anticodon supply",atop("tRNA abundance, x100,000",""))), 
        color = "recognised codons",
        title=expression(atop("Balance between anticodon demand and supply",atop("in late response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")

ggsave( plot = g, filename = paste0(PATH,"part2/tRNAab--x.anticodon_demand_y.total.tRNAab_z.stress_t.120min.pdf"), 
        dpi=250, width = 10.9, height=4.35,  useDingbats=FALSE)





#-------------------------#


data.figure  <- subset(anticodon.master.table, time %in% c(0,20))
data.figure$experiment <- factor( data.figure$experiment, levels = c("normal","diauxic","ox","osm","temp"), labels = c("No stress","Diauxic shift","Oxidative stress","Osmotic shock","Temperature stress") )
data.figure$time <- factor( data.figure$time, labels = c("0min","20min"))
data.figure$variation.demand  <- (data.figure$demand.mRNAab - data.figure$demand.mRNAab.normal)/data.figure$demand.mRNAab.normal
data.figure$variation.supply  <- (data.figure$total.tRNAab - data.figure$total.tRNAab.normal)/data.figure$total.tRNAab.normal
data.figure$variation.demand.up  <- (data.figure$demand.mRNAab.up - data.figure$demand.mRNAab.up.normal)/data.figure$demand.mRNAab.up.normal

### BUG!!! R2.data <- ddply( data.figure, .(experiment, time), function(c){ lm_eqn_special(df = c) })

## [] Do variation in demand corresponding to up-regulated transcripts match with anticodon supply? #####
g <- ggplot(  data = subset(data.figure, time !="0min"), aes(x= variation.demand, y = total.tRNAab/100000 ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = percent) +
  scale_y_continuous( labels = prettyNum ) +
  labs( x = expression(atop("variation in anticodon demand",atop("(in %)",""))), 
        y = expression(atop("anticodon supply",atop("tRNA abundance, x100,000",""))), 
        color = "recognised codons",
        title=expression(atop("Balance between anticodon demand and supply",atop("in early response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")


## [] Do variation in demand corresponding to up-regulated transcripts match with anticodon supply? #####
g <- ggplot(  data = subset(data.figure, time !="0min"), aes(x= variation.demand.up, y = total.tRNAab/100000 ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = percent) +
  scale_y_continuous( labels = prettyNum ) +
  labs( x = expression(atop("variation in anticodon demand",atop("of up-regulated mRNAs (top-20%)",""))), 
        y = expression(atop("anticodon supply",atop("tRNA abundance, x100,000",""))), 
        color = "recognised codons",
        title=expression(atop("Balance between anticodon demand and supply",atop("in early response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")


## [] Do variation in demand correlate with variation in supply? #####
g <- ggplot(  data = subset(data.figure, time !="0min"), aes(x= variation.demand, y = variation.supply ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = percent) +
  scale_y_continuous( labels = percent ) +
  labs( x = expression(atop("variation in anticodon demand",atop("(in %)",""))), 
        y = expression(atop("variation in anticodon supply",atop("(in %)",""))), 
        color = "recognised codons",
        title=expression(atop("Balance between anticodon demand and supply",atop("in early response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")


## [] Do variation in demand from up-regulated transcripts correlate with variation in supply? #####
g <- ggplot(  data = subset(data.figure, time !="0min"), aes(x= variation.demand.up, y = variation.supply ) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  #scale_x_continuous( labels = percent) +
  #scale_y_continuous( labels = percent ) +
  labs( x = expression(atop("variation in anticodon demand",atop("(in %)",""))), 
        y = expression(atop("variation in anticodon supply",atop("(in %)",""))), 
        color = "recognised codons",
        title=expression(atop("Variation in anticodon demand and supply",atop("in early response to stress","")))) +
  facet_wrap( time ~ experiment, scale="free", nrow =1 ) + 
  theme_gul + theme(legend.position="bottom")

## 
p1 <- ggplot(  data = subset(data.figure, time !="0min"), aes(x= as.numeric(experiment), y = variation.demand.up ) ) + 
  geom_line( aes(group=anticodon), color = "gray", size=0.5, alpha=0.7) +
  geom_point( size = 3.25, aes(color=CAI.O), alpha=0.9 ) + 
  #geom_text( size = 1.5, aes(label=anticodon) ) +
  #geom_text( data = R2.data, aes_string( label = expression( V1 ) ), parse = T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  scale_x_continuous( labels = levels(data.figure$experiment)[-1]) +
  scale_y_continuous( labels = percent, limits = c(0,1.80)  ) +
  labs( x = "", 
        y = expression(atop("variation in anticodon supply",atop("(in %)",""))), 
        color = "recognised codons",
        title=expression(atop("Variation in anticodon demand",atop("in early response to stress","")))) +
  theme_gul + theme(legend.position="none")

temp.data <- subset(gene.master.table, time == 20 )
temp.data$experiment <- factor( as.character(temp.data$experiment), levels = c("normal","diauxic","ox","osm","temp"), labels = c("No stress","Diauxic shift","Oxidative stress","Osmotic shock","Temperature stress") ) 
p2 <- ggplot(  data = temp.data, aes(x= experiment, y = log2.mRNA_abundance ) ) +  
  #geom_jitter(aes(color=up_reg),stat = "identity", height =0 ) + 
  geom_dotplot(data = subset(temp.data, abs(log2.mRNA_abundance) > 0.01), aes(fill=up_reg, color =up_reg), binaxis="y", 
               stackdir="center", binwidth=0.03, alpha=1 ) +
  geom_hline( yintercept = 0, color = "gray40", lty=2) +
  geom_boxplot(data = temp.data, aes(group=paste(experiment,up_reg)), notch=T, outlier.size = 0, width = 0.5 ) +
  scale_y_continuous( labels = prettyNum) +
  scale_fill_manual( values = c("gray60","skyblue3") ) + 
  scale_color_manual( values = c("gray60","skyblue3") ) + 
  labs( x = "", 
        y = expression(atop("variation in mRNA abundance (log2)",atop("Top-20% genes",""))), 
        color = "recognised codons",
        title=expression(atop("Increase in mRNA abundance from up-regulated genes",atop("in early response to stress","")))) +
  theme_gul + theme(legend.position="none")

rm(temp.data)
require(gridExtra)

pdf( file = paste0(PATH,"part2/upreg--impact_of_mRNAab_on_demand.pdf"), width = 9, height = 9)
grid.arrange(p1, p2, nrow = 2, heights=c(4/9, 5/9) )
dev.off()


#--------------------------#

# [] how similar the anticodon demand is across stress conditions?

anticodon.demand.comparison <- list(
  m.20  = matrixify( cast( subset(anticodon.master.table, time !=120) , formula = anticodon ~ experiment, value = "demand.mRNAab" ) ),
  m.120 = matrixify( cast( subset(anticodon.master.table, time != 20) , formula = anticodon ~ experiment, value = "demand.mRNAab" ) ) 
)

anticodon.demand.correlation     <- ldply( anticodon.demand.comparison, function(x){ melt(cor(x[complete.cases(x),])) })
anticodon.demand.correlation$X1  <- factor(anticodon.demand.correlation$X1, levels = c("normal","diauxic","ox","osm","temp"))
anticodon.demand.correlation$X2  <- factor(anticodon.demand.correlation$X2, levels = c("normal","diauxic","ox","osm","temp"))
anticodon.demand.correlation$.id <- factor(anticodon.demand.correlation$.id, levels = c("m.20","m.120"), labels = c("t=20","t=120"))


g <- ggplot( data = melt(anticodon.demand.correlation), aes(x=X1, y=X2, fill=value)) + 
  geom_tile() + 
  geom_text( aes(label = round(value,2) ), size = 3, alpha=0.5, color = "white" ) +
  coord_fixed() + labs(x="", y = "", title = "Correlation of anticodon demand") +
  theme_gul +
  theme( axis.text.x = element_text(angle=90)) +
  facet_wrap( ~ .id , drop = T)

ggsave(plot=g, filename = paste0(PATH,"part2/mRNAab.correlation.heatmaps.png"),dpi = 250, width=6.31, height=6.38 )

# ------------------------- #
# [] how similar the mRNA abundance is across stress conditions? -------
mRNAab.comparison <- list(
      m.20  = matrixify( cast( subset(gene.master.table, !is.na(est.mRNA_abundance) & time !=120) , formula = ORF ~ experiment, value = "est.mRNA_abundance" ) ),
      m.120 = matrixify( cast( subset(gene.master.table, !is.na(est.mRNA_abundance) & time != 20) , formula = ORF ~ experiment, value = "est.mRNA_abundance" ) ) 
  )


mRNAab.correlation     <- ldply( mRNAab.comparison, function(x){ melt(cor(x[complete.cases(x),])) })
mRNAab.correlation$X1  <- factor(mRNAab.correlation$X1, levels = c("normal","diauxic","ox","osm","temp"))
mRNAab.correlation$X2  <- factor(mRNAab.correlation$X2, levels = c("normal","diauxic","ox","osm","temp"))
mRNAab.correlation$.id <- factor(mRNAab.correlation$.id, levels = c("m.20","m.120"), labels = c("t=20","t=120"))


g <- ggplot( data = melt(mRNAab.correlation), aes(x=X1, y=X2, fill=value)) + 
  geom_tile() + 
  geom_text( aes(label = round(value,2) ), size = 3, alpha=0.5, color = "white" ) +
  coord_fixed() + labs(x="", y = "", title = "Correlation of mRNA abundance") +
  theme_gul +
  theme( axis.text.x = element_text(angle=90)) +
  facet_wrap( ~ .id , drop = T)

ggsave(plot=g, filename = paste0(PATH,"part2/mRNAab.correlation.heatmaps.png"),dpi = 250, width=6.31, height=6.38 )

#-------------------------#



# [] Does the relative availability of synonymous anticodons scale with the anticodon demand? ----
g <- ggplot(  data = anticodon.master.table, aes(x= demand.mRNAab, y = relative.availability ) ) + 
  geom_point( size = 3, aes(color=CAI.O) ) + 
  stat_smooth( method="lm", color="gray30" ) + 
  facet_wrap( time ~ experiment ) + 
  scale_color_manual( values = c("skyblue3","gray","orange")) + 
  labs( x = "\nanticodon demand", 
        y = "relative availability of tRNA\n", 
        color = "recognised codons",
        title="The relative availability of synonymous anticodons scale with anticodon demand\n") +
  theme_gul + theme(legend.position="bottom")


ggsave(plot = g, file=paste0(PATH,"part2/anticodon.economy--x.anticodon.demand_y.relative.availability_z.stress.png"),dpi = 250, width=7.14, height=5.11)




#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------# Elongation kinetics #-------------------------------------------------------------------------------------------------#
#-------------------------------# Relative adaptiveness  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# ---- are codon relative adaptiveness values reflected by average elongation times measured in the stochastic model?
# in fact, not much. but this points out that the readjustments in the mRNA pool are at least equally important to explain variation in elongation kinetics
# - or that the method to compute ws is obsolote
# - or that there are descripencies between my calculations of tGCNs and the numbers Marc used. 

# wait ---- what? the tAI of the of some codons is totally underestimated! in 0, normal it should be as close as possible to the line 
# well, it could also be that the tGCN values used back in the days by dosReis2004 don't match with the ones Marc or I used later
# (since tAI is the mere reflect of dos Reis' calculations )
g <- ggplot(data = subset(codon.master.table,!is.na(tAI)), aes(x= tAI, y= w, color=CAI.O) ) +
  geom_abline(slope=1, intercept=0, color="skyblue2") +
  geom_point() +
  geom_text(data= subset(codon.master.table, CAI.O == "O" & tAI < 0.2), 
            aes(label=codon), size = 4, color = "gray20" ) +
  geom_text(data= subset(codon.master.table, CAI.O == "N" & tAI > 0.75), 
            aes(label=codon), size = 4, color = "gray20" ) +
  labs(y="codon relative adaptiveness", x="tAI") +
  scale_color_manual( values = c("skyblue3","orange")) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  facet_wrap( time ~ experiment ) +
  theme_gul

ggsave(plot=g, filename=paste0(PATH,"part2/adaptiveness--x.tAI_y.cra_z.stress_t.CAI.O.png"), dpi=250,
       width=5.92, height=5.89)

subset(codon.master.table, experiment == "normal", select= c(codon, tAI, w))

# update April 2015, to account for the new scaling of w scores

# [] w^n vs w^s --------
ww <- as.data.frame( cast( subset(codon.master.table,!is.na(tAI), select = c(codon, experiment, time, w.2, CAI.O)), codon + CAI.O ~ experiment + time, value = "w.2"   ) )
w.table <- data.frame( melt( ww, id.vars=c("codon","CAI.O","normal_0")), row.names=NULL)
colnames(w.table) <- c("codon","CAI.O","wn","stress", "ws")
w.table$label <- factor( w.table$stress, 
        levels= c("diauxic_20","ox_20","osm_20","temp_20","diauxic_120","ox_120","osm_120","temp_120"),
        labels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"))
)




ww <- arrange( ww, CAI.O, codon)[,c("codon","CAI.O","normal_0","diauxic_20","diauxic_120","ox_20","ox_120","osm_20","osm_120","temp_20","temp_120")]
write.table(ww, "~/Dropbox/PhD/thesis/TeX/items/tables/chapter5/table_adaptiveness.txt", row.names=F, quote=F, sep="\t")

# rosetta.codons[,c("codon", "anticodon","aa", "tAI", "nTE", "CAI.O")]

lm_eqn_special2 = function(df, exclude = "ttg"){
  m1 = lm( ws ~ wn, df);
  m2 = lm( ws ~ wn, subset(df, !codon %in% exclude));
  eq.1 <- substitute(italic(R)^2~"="~r2 ,          list( r2     = format(summary(m1)$r.squared, digits = 2) ))
  #eq.2 <- substitute(italic(R)[-ttg]^2~"="~r2.exclude, list( r2.exclude = format(summary(m2)$r.squared, digits = 2) ))
  #as.character( as.expression( paste( c('atop(', eq.1, ',', eq.2, ')'), collapse=" ") ) )
  as.character( as.expression( paste( c('atop(', eq.1, ')'), collapse=" ") ) )
}

# grume de grutte
R2.data <- ddply( w.table, .(label), function(x){ lm_eqn_special2(x) }  )

  g <- ggplot(data = w.table, aes(x= wn, y= ws) ) +
  stat_smooth(method = "lm", color = "gray60") +
  geom_point(aes(fill=CAI.O), color="gray30",pch=21, size = 3) +
  geom_hline(yintercept=1, color="gray30", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray30", lty=2) +
  geom_text( data = subset(w.table, wn * ws > 0.16 ), aes(label=codon), size = 2.5, alpha=0.6, vjust = 1.65) +
  geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  labs(y="codon relative adaptiveness during stress", x="tAI") +
  labs( x = expression(atop(w^n, atop("codon relative adaptiveness in normal conditions",""))), 
        y = expression(atop(w^s, atop("codon relative adaptiveness during stress",""))), 
        fill = "optimality"
  ) +
  coord_fixed(ratio = 0.5) +
  scale_fill_manual( values = c("skyblue3","orange")) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  facet_wrap( ~ label, nrow = 2 ) +
  theme_gul + theme(legend.position = "none")

ggsave(plot=g, filename=paste0(PATH,"fig_results_relative_adaptiveness.png"), dpi=250,
       width=7, height=5.5)

ggsave(plot=g, filename=paste0(PATH,"part2/adaptiveness--x.wn_y.ws.stress_t.CAI.O.png"), dpi=250,
       width=10.8, height=6.38)


ggsave(plot=g, filename=paste0(PATH,"part2/adaptiveness--x.wn_y.ws.stress_t.CAI.O.pdf"), dpi=250,
       width=10.8, height=6.38, useDingbats=F)

rm(R2.data)



# ---------------------------------- #

ddply(study.translation.kinetics, .(experiment, time), summarise, av.speed = mean(elongation.speed), av.normal = mean(elongation.speed.normal) )

ddply(codon.master.table, .(experiment, time), summarise, min.speed = min(1/av.elongation_time), max.speed = max(1/av.elongation_time), dnr = max(1/av.elongation_time)/min(1/av.elongation_time) )

# ---------------------------------- #


# [] how do elongation kinetics scale with relative adaptiveness? ----
g <- ggplot(data = subset(codon.master.table, experiment!="normal" & topology != "aa(tRNA(c))" ), 
            aes(x= w, y= av.elongation_time) ) +
  #stat_smooth( method="gam", formula = y ~ s(x, bs = "ts"), se = F ) +
  geom_point(color="gray30", size=3.1) +
  geom_point(aes(color=log2.foldchange), size =2.5) +
  labs( x = expression(atop(w^s, atop("(codon relative adaptiveness)",""))), 
        y = "average codon\nelongation time (s)", 
        color = "log2 foldchange\ntRNA abundance"
  ) +
  stat_smooth( method="gam", formula = (y) ~ log(x)*x, se = F, color = "gray30") + # , aes(outfit = fit <<- ..y..) ) +
  facet_wrap( ~ label, nrow = 2 ) +
  #scale_color_gradient(low = "skyblue",  high="orange") +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_color_gradient2(low = "skyblue", mid="white", high="orange",midpoint = 0) +
  theme_gul + theme(legend.position="bottom" )  

ggsave(plot = g, filename = paste0(PATH,"fig_results_relative_adaptiveness_speed.png"), 
       dpi=250, width=9, height=5.5 )

ggsave(plot = g, filename = paste0(PATH,"part2/tRNApool--x.codon.rel.adaptiveness_y.av.elongation.time_z.stress_t.time.png"), 
       dpi=250, width=7, height=5 )


#--------------------------------------------------------------#

# I will need to add a delta.av.elongation_time to the codon.master.table and a ratio as well, and use it here as "y"
g <- ggplot(data = subset(codon.master.table, experiment!="normal" ), aes(x= delta_w, y= av.elongation_time) ) +
  geom_point(color="gray30", size=3.1) +
  geom_point(aes(color=w), size =2.5) +
  geom_vline(xintercept=0, size=rel(0.25)) +
  labs( x = expression(atop(Delta~w^s, atop("(change in codon relative adaptiveness)",""))), 
        y = "average elongation time (s)", 
        color = expression(w^s)
  ) +
  facet_grid( time ~ experiment ) +
  scale_color_gradient(low = "white",  high="orange") +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  #scale_color_gradient2(low = "skyblue", mid="white", high="orange",midpoint = 0.25) +
  theme_gul + theme(legend.position="bottom" )  

ggsave(plot = g, filename = paste0(PATH,"part2/tRNApool--x.delta.w_y.av.elongation.time_z.stress_t.time.png"), 
       dpi=250, width=7, height=5 )


##

# [] delta_elongation in function of delta_w ----
g <- ggplot(data = subset(codon.master.table, experiment!="normal" ), aes(x= delta_w, y= delta_elongation) ) +
  geom_point(color="gray30", size=3.1) +
  geom_point(aes(color=w), size =2.6) +
  geom_hline(yintercept=0, size=rel(0.25)) +
  geom_vline(xintercept=0, size=rel(0.25)) +
  labs( x = expression(atop(Delta~w^s, atop("(change in codon relative adaptiveness)",""))), 
        y = expression(atop(Delta~e, atop("change in average elongation time",""))), 
        color = expression(w^s)
  ) +
  facet_grid( time ~ experiment ) +
  scale_color_gradient(low = "white",  high="orange") +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  #scale_color_gradient2(low = "skyblue", mid="white", high="orange",midpoint = 0.25) +
  theme_gul + theme(legend.position="bottom" )  

ggsave(plot = g, filename = paste0(PATH,"part2/tRNApool--x.delta.w_y.delta_elongation_z.stress_t.time.png"), 
       dpi=250, width=7, height=5 )

#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#  Anticodon demand  #-------------------------------------------------------------------------------------------------#
#--------------------------------# Up-regulated genes ---------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#



# [] in top-20% up-regulated genes, blue codons "invade" the right region (increase in relative demand) ----
  # g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= rank.demand.mRNAab.up, y= w )) + 
  #   stat_smooth(method="lm",color="gray40") +
  #   geom_point(alpha=0.8, aes(color=CAI.O), size=5) + 
  #   geom_text(aes(label=codon), size=2.25, color="gray20") +
  #   facet_grid(experiment ~ time, scale="free") + 
  #   scale_color_manual(values=c("skyblue3","orange")) +
  #   labs( x = expression(atop("ranked demand",atop("in top-20% up-regulated mRNAs",""))), 
  #         y = expression(atop("w",atop("codon relative adaptiveness",""))), 
  #         color = "recognised codons"
  #       ) +
  #   theme_gul + theme(legend.position="bottom")
  # 
  # ggsave( plot = g, filename = paste0(PATH,"part2/up_reg--x.ranked.demand_y.cra_z.time_t.stress.png"), 
  #         dpi=250, width=6.5, height=6 )




# [] gain in ranks for codon demand (mRNA abundance) ----
    # g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= rank.demand.mRNAab.up - rank.demand.mRNAab, y=w)) + 
    #   geom_vline(xintercept=0, lty=2) +
    #   geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
    #   geom_text(aes(label=codon), size=2.25, color="gray20") +
    #   facet_grid(experiment ~ time, scale="free") + 
    #   scale_color_manual(values=c("skyblue3","orange")) +
    #   labs(x="\ndifferential ranked demand\n(codon occurence x mRNA abundance summed over genes)\nup-regulated mRNAs - all mRNAs",
    #        y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
    #   theme_gul + theme(legend.position="bottom")
    # 
    # ggsave( plot = g, filename = paste0(PATH,"part2/up_reg--x.diff.ranked.demand_y.cra_z.time_t.stress.png"), 
    #         dpi=250, width=6.5, height=6 )





# [] same, based on actual demand (not rank) -----
g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= demand.mRNAab.up - demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x=expression(Delta[demand]), y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_gul + theme(legend.position="bottom")

ggsave( plot = g, filename = paste0(PATH,"part2/up_reg--x.diff.demand_y.cra_z.time_t.stress.png"), 
        dpi=250, width=6.5, height=6 )






# [] relative change in demand of up-reg mRNAs compared to all mRNAs ----
  # g <- ggplot( data = subset(codon.master.table, experiment !="normal"), 
  #              aes(x= (demand.mRNAab.up - demand.mRNAab)/demand.mRNAab, y=w)) + 
  #   #geom_vline(xintercept=0, lty=2) +
  #   geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  #   geom_text(aes(label=codon), size=2.25, color="gray20") +
  #   facet_grid(experiment ~ time, scale="free") + 
  #   scale_color_manual(values=c("skyblue3","orange")) +
  #   labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  #   theme_gul + theme(legend.position="bottom")
  # 
  # ggsave( plot = g, filename = paste0(PATH,"part2/up_reg--x.rel.diff.demand_y.cra_z.time_t.stress.png"), 
  #         dpi=250, width=7, height=5 )





# [] relative change (compared to all mRNAs) in relative demand (compared to other codons) ----
g <- ggplot( data = subset(codon.master.table, experiment !="normal"), 
             aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2, color ="gray") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20",alpha=0.9) +
  facet_grid(experiment ~ time, scale="free") + 
  scale_y_continuous( labels = prettyNum ) +
  scale_color_manual(values=c("skyblue3","orange")) +
  labs( x = expression(atop("relative difference in demand", atop("(in top-20% up-regulated mRNAs vs all mRNAs)",""))), 
        y = expression(atop("w", atop("codon relative adaptiveness",""))), 
        color = "optimality\n(normal cond.)"
  ) +
  theme_gul + theme(legend.position="bottom")

ggsave(plot=g, filename = paste0(PATH,"part2/up_reg--x.rel.diff_y.cra_z.time_t.stress.png"),
       dpi=250, width=6.5, height=6)





# [] relative diff. in demand vs change in codon relative adaptiveness ------
g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y=delta_w)) + 
  geom_hline(yintercept=0, lty=2, color = "gray") +
  geom_vline(xintercept=0, lty=2, color = "gray") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free" ) + 
  scale_color_manual(values=c("skyblue3","orange")) +
  scale_y_continuous( labels = prettyNum ) +
  labs( x = expression(atop("relative difference in demand", atop("(in top-20% up-regulated mRNAs vs all mRNAs)",""))), 
        y = expression(atop(Delta~w, atop("change in codon relative adaptiveness",""))), 
        color = "optimality\n(normal cond.)"
  ) +
  #labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="w\nchange in relative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_gul + theme(legend.position="bottom")

ggsave(plot=g, filename = paste0(PATH,"part2/up_reg--x.rel.diff_y.delta_w_z.time_t.stress.png"),
       dpi=250, width=6.5, height=6)






# ---- relative diff. in mRNA vs elongation time
g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y= av.elongation_time )) + 
  geom_vline(xintercept=0, lty=2, color = "gray") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs( x = expression(atop("relative difference in demand", atop("(in top-20% up-regulated mRNAs vs all mRNAs)",""))), 
        y = "elongation time (s)", 
        color = "optimality\n(normal cond.)"
  ) +
  theme_gul + theme(legend.position="bottom")

ggsave(plot=g, filename = paste0(PATH,"part2/up_reg--x.rel.diff_y.elongation.time_z.time_t.stress.png"),
       dpi=250, width = 6.5, height =6 )


# [] change in number of elongation events in function of codon relative adaptiveness
g <- ggplot( data = subset(codon.master.table, experiment !="normal"),
             aes(x=w, y= demand.events.up - demand.events)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(y=expression(Delta[events]), x="\nw\nrelative adaptiveness",color="optimality\n(normal cond.)") +
  theme_gul + theme(legend.position="bottom")


ggsave(plot=g, filename = paste0(PATH,"part2/up_reg--x.cra_y.delta_events_z.time_t.stress.png"),dpi=250)



g <- ggplot( data = subset(codon.master.table, experiment !="normal"), 
        aes(x= round(w,2) - round(delta_w,2), y=w)) + 
  geom_abline(intercept=0, slope=1, size = 0.25, color="gray30") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  scale_color_manual(values=c("skyblue3","orange")) +
  scale_y_continuous(labels=prettyNum, limits= c(0,1), breaks = c(0.25,0.75)) +
  scale_x_continuous(labels=prettyNum, limits= c(0,1), breaks = c(0.25,0.75)) +
  coord_fixed() +
  facet_grid(experiment + time ~ aa) +
  labs(x=expression(w^n), y=expression(w^s),
       color="optimality\n(normal cond.)") +
  theme_gul +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=90, hjust=1))

ggsave(plot=g, filename = paste0(PATH,"part2/adaptiveness--x.wn_y.ws_z.aa_t.stress.png"),
       dpi=200, width = 9.75, height = 5)


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------# Elongation kinetics  #----------------------------------------------------------------------------------------------#
#---------------------------------#  Amino acid demand  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


lm_eqn_special3 = function(df, exclude = "Leu"){
  m1 = lm( aa.supply ~ aa.demand, df);
  m2 = lm( aa.supply ~ aa.demand, subset(df, !aa %in% exclude));
  eq.1 <- substitute(italic(R)^2~"="~r2 ,              list( r2     = format(summary(m1)$r.squared, digits = 2) ))
  eq.2 <- substitute(italic(R)[-Leu]^2~"="~r2.exclude, list( r2.exclude = format(summary(m2)$r.squared, digits = 2) ))
  as.character( as.expression( paste( c('atop(', eq.1, ',', eq.2, ')'), collapse=" ") ) )
}

R2.data <- ddply( aa.master.table, .(experiment, time), function(x){ lm_eqn_special3(x) }  )


# [] how the amino acid demand scales with total tRNA abundance ----
g <- ggplot( data = aa.master.table, aes(x= aa.demand/100000, y= aa.supply/100000 )) + 
  stat_smooth( method="lm",color="gray40") +
  geom_point(alpha=0.8, size=5, aes(color=metab.cost.resp)) + 
  geom_text(aes(label=aa), size=2.25, color="white") +
  geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  #scale_color_manual(values=c("skyblue3","orange")) +
  labs( x = expression(atop("amino acid demand (x100,000)", atop("(pooled anticodon demand)",""))), 
        y = "total tRNA abundance\n(x100,000)", 
        color = "metabolic\ncost"
  ) +
  facet_wrap(experiment ~ time, scale="free", ncol=3) + 
  theme_gul + theme(legend.position="bottom")

ggsave(plot=g, filename = paste0(PATH,"part1/aa.economy--x.aa.demand_y.total.tRNAab_z.stress.png"),
       dpi=250, width = 5.5, height = 6.5)

rm(R2.data)

# compare correlations between changes in supply ( s - n )/n and metabolic cost. ok?


aa <- as.data.frame( cast( subset(aa.master.table, select = c(aa, experiment, time, aa.supply, aa.demand, metab.cost.ferm)), aa + metab.cost.ferm ~ experiment + time, value = "aa.supply"   ) )
aa.delta.table <- data.frame( melt( aa, id.vars=c("aa","metab.cost.ferm","normal_0")), row.names=NULL)
colnames(aa.delta.table) <- c("aa","metab.cost.ferm","aa.supply.n","stress", "aa.supply")
aa.delta.table$change.supply <- (aa.delta.table$aa.supply - aa.delta.table$aa.supply.n) / aa.delta.table$aa.supply.n


g <- ggplot( data = aa.delta.table, aes(x= metab.cost.ferm, y= change.supply )) + 
  geom_hline( yintercept = 0, color = "gray") +
  stat_smooth( method="lm",color="gray40") +
  geom_point(alpha=0.8, size=5, aes(color=metab.cost.ferm)) + 
  geom_text(aes(label=aa), size=2.25, color="white") +
  #geom_text( data = R2.data, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  labs( 
    x = expression(atop("metabolic cost", atop("(n. active phosphate needed in fermentation)",""))),
    y = expression(atop("change in supply", atop("(pooled anticodon supply)",""))), 
        color = "metabolic\ncost"
  ) +
  facet_wrap( ~ stress, scale="free", ncol=4) + 
  theme_gul + theme(legend.position="bottom")


ggsave(plot = g, filename = paste0(PATH,"part1/aa.economy--x.metab.cost_y.delta_elongation_z.experiment_z.time.png"))


# compare how the overal demand for amino acid and the demand in up=regulated transcripts fit with the actual tRNA supply
    # data.comp <- melt( subset(aa.master.table, experiment !="normal"), id.vars = c("aa","experiment","time", "aa.supply", "aa.supply.free") )
    # 
    # g <- ggplot( data = subset(data.comp, variable %in% c("rel.aa.demand","rel.aa.demand.up")) , aes(x= value * 100, y= aa.supply.free/100000, color = variable, group=variable )) + 
    #   geom_line( aes(group=aa), color="orange", alpha=0.4) +
    #   stat_smooth(method="lm") +
    #   geom_point(alpha=0.8, size=5) + 
    #   geom_text(aes(label=aa), size=2.45, color="white") +
    #   scale_color_manual(values=c("gray20","#D46A6A50")) +
    #   labs(x="\nrelative amino acid demand(%)", y="total tRNA abundance\n(x100,000)\n") +
    #   facet_wrap(experiment ~ time, scale="free", ncol=2) + 
    #   theme_gul + theme(legend.position="bottom")


# codon competition
codon.competition <- read.table("results/master tables/codon.competition.txt",sep="\t",header=1)
codon.competition$label <- factor( codon.competition$stress, 
                                   levels = c("diauxic_20","ox_20","osm_20","temp_20","diauxic_120","ox_120","osm_120","temp_120","normal_0"),
                                   labels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min")
                                   )

codon.competition$percent.var.competition <- with(codon.competition, (competition - competition.normal)/competition.normal  )
codon.competition$percent.var.elongation  <- with(codon.competition, (delta_elongation)/(av.elongation_time/ratio_elongation)  )

regressions.competition   <- ddply( subset( codon.competition, !is.na(competition.change) & time >0 & ratio_elongation < 15 ), .(label), function(x){ batch.regressions(x, x = "competition.change",y = "ratio_elongation") })
regressions.competition.2 <- ddply( subset( codon.competition, !is.na(percent.var.competition) & time >0 & ratio_elongation < 15 ), .(label), function(x){ batch.regressions(x, x = "percent.var.competition",y = "percent.var.elongation") })

# ggplot( data = codon.competition, aes( x = competition, y = av.elongation_time) ) + 
#   stat_smooth(method="lm",color="white") +
#   geom_point( aes(color = CAI.O ), size = 2.5 ) + scale_x_log10() +
#   facet_wrap( ~ experiment + time ) + scale_color_manual(values=c("skyblue","orange"))

require(scales)
g <- ggplot( data = subset( codon.competition, time >0 ) , aes( x = competition.change, y = ratio_elongation) ) + 
  geom_vline( xintercept = 1, lty =3, color="gray40" ) +
  geom_hline( yintercept = 1, lty =3, color="gray40" ) +
  stat_smooth(method="lm",color="#20222D") +
  geom_text( data = subset(codon.competition, competition.change > 5 | (competition.change > 3 & ratio_elongation > 3) ), aes(label=codon), size = 2.5, alpha=0.9, vjust = 1.65) +
  geom_text( data = subset(codon.competition, competition.change < 0.3 ), aes(label=codon), size = 2.5, alpha=0.9, vjust = 1.65) +
  geom_text( data = subset(codon.competition, codon == "cgg" & time > 0), aes(label=codon), size = 2.5, alpha=0.9, vjust = 1.65) +
  geom_text( data = regressions.competition, aes_string( label = expression(V1) ), parse =T, x=Inf, y=-Inf, hjust = 1.05, vjust = -0.5, size = 2.75 ) +
  geom_point( aes(fill = CAI.O ), pch=21, size = 2.5 ) +
  scale_x_log10(labels=prettyNum) +
  scale_y_log10(labels=prettyNum) +
  coord_fixed() +
  labs(x="competition strength: increase in the relative abundance of near-cognate tRNAs\n(log scale)", y= "change in codon elongation time\n(log scale)", fill = "codon type") +
  facet_wrap( ~ label, ncol= 4) + scale_fill_manual(values=c("skyblue3","orange")) + theme(legend.position="none")

ggsave(plot=g, filename = paste0(PATH,"fig_results_codon_elongation_competition_new.png"), dpi=200, width=9, height=5.5 )


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#------------------------------#||| Analysis Part III #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


    # In this section, the main research objective is to investigate how individual proteins are affected by changes in tRNA abundance and in mRNA abundance
    # We will thus exploit the stochastic model of protein translation (SMoPT) published by Shah et al. 2013 to analyse changes in initiation frequency and 
    # in elongation speed. 
    
    # Are up-regulated genes (either at the transcriptional or translational levels, or both) somehow better adapted to the current tRNA pool in a given condition?
    # Are down-regulated genes less adapted 
        # (note that this question and the previous are completely independent, one could be true without the need for the other to be true or false)
        # If such relationship exist, one should observe a correlation between the codon usage of these genes and the concentrations in the tRNA pool, 
        # or indirectly between codon usage and changes in tRNA abundance
        # this might be observed for genes that are clearly undergoing some form of regulation (ie where the cell actively spends energy to up or down regulate their expression)




#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------# Elongation kinetics #-----------------------------------------------------------------------------------------------#
#--------------------------------#  Defining TFS genes  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# [] variations in initiation vs variations in elongation, and hihglight FTS genes ----


g <- ggplot(data= subset(SMoPT.data, experiment !="normal"), 
       aes( x = delta_initiation, y = delta_elongation, 
            #color = tAI.profile
            color = faster.translation.2
            #color = delta_initiation < - delta_elongation
            )  ) + 
  geom_point()  +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  geom_abline(intercept=0, slope=-1) +
  labs(x=expression(paste(Delta["i"]~"\n(change in waiting time for initiation)")), 
       y=expression(paste(Delta["e"]~"\n(change in average elongation time)"))
       ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_color_manual(values = c("gray20","red")) +
  facet_grid( time + tAI.profile ~ experiment) + 
  theme_gul + 
  theme(legend.position="none")

ggsave( plot = g, filename = paste0(PATH,"part3/SMoPT--x.delta_initiation_y.delta_elongation_z.experiment_z.time.pdf"), useDingbats =F,
        dpi = 250, width = 9.75, height = 5)



# [] rel. translation time in function of relative w time for initiation ----
# --- there is a clear linear dependency between changes in waiting time and changes in total translation time 
# which highlights the strong contribution of initiation to variation in translation time during stress
g <- ggplot( data=subset(SMoPT.data, experiment !="normal" & ratio_initiation > 0 ), 
        aes( x = log2(ratio_initiation), y = log2(ratio_total.time), color = faster.translation )  ) + 
  geom_point()  +
  stat_smooth( method = "lm", color ="gray70",size=1, alpha=0.8) +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  labs( x="relative waiting time\nbetween initation events (log2)", y="relative\ntranslation time (log2)" ) +
  scale_color_manual(values = c("gray20","red")) +
  facet_grid( time ~ experiment) + 
  theme_gul + 
  theme(legend.position="none")

ggsave( plot = g, filename = paste0(PATH,"part3/SMoPT--x.relative_initiation_y.relative_translation.time_z.experiment_z.time.pdf"), useDingbats=F,
       dpi = 250, width = 9.75, height = 6)

# [] rel. translation time in function of relative elongation time ----
# [] such relationship doesn't exist for elongation time, except for genes that have a very high tAI (>0.5) (shown in triangles) in normal conditions -----
g <- ggplot( data=subset(SMoPT.data, experiment !="normal" & ratio_initiation > 0 ), 
        aes( x = log2(ratio_elongation), y = log2(ratio_total.time), color = faster.translation, shape = tRNA_adaptation_index >=0.5 )  ) + 
  geom_point(  data = subset(SMoPT.data, experiment !="normal" & ratio_initiation > 0 & tRNA_adaptation_index < 0.5))  +
  stat_smooth( method = "lm", color ="gray70",size=1, alpha=0.8) +
  geom_point(  data = subset(SMoPT.data, experiment !="normal" & ratio_initiation > 0 & tRNA_adaptation_index >=0.5), fill="white")  +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  labs( x="relative\ntime of elongation (log2)", y="relative\ntranslation time (log2)" ) +
  scale_color_manual(values = c("gray20","red")) +
  scale_shape_manual(values = c(16,24) ) +
  facet_wrap( ~ label) + 
  theme_gul + 
  theme(legend.position="none")

ggsave( plot = g, filename = paste0(PATH,"part3/SMoPT--x.relative_elongation_y.relative_translation.time_z.experiment_z.time.pdf"), useDingbats=F, 
        dpi = 250, width = 9.75, height = 6)


# [] quantiles translation time vs elongation vs changes in stAI
g <- ggplot(data=SMoPT.data, aes(x=total.time.q, y= av.elongation_time )) + 
  geom_boxplot(aes(fill=tAI.profile), outlier.colour = "gray60", notch=F) + 
  facet_wrap( ~ label) + 
  labs(x="average translation time (quantiles)", y="average elongation time", fill="stAI profile") +
  scale_fill_manual(values=c("red", "orange", "gray")) +
  theme(axis.text.x = element_text(angle=90, hjust=0))

ggsave( plot = g, filename = paste0(PATH,"part3/SMoPT-quantiles.elongation_adaptation.pdf"), useDingbats=F, 
        dpi = 250, width = 9.75, height = 9.75)

# ------ technical controls/assessments

# [] is the frequency of initiation dependent on initiation probability? ----
g <- ggplot( data = subset(SMoPT.data, experiment !="normal" & av.initiation_time > 0), 
        aes(x= IniProb, y = 1/av.initiation_time ) ) + 
  geom_point(alpha=0.25) +
  scale_x_sqrt( breaks = c(5e-3, 5e-2, 5e-1, 1e-1), labels = prettyNum ) +
  scale_y_sqrt( breaks = c(0.1, 1, 5, 15, 30, 45),  labels = prettyNum ) +
  facet_grid( time ~ experiment ) + 
  #stat_smooth( method="lm" ) +
  stat_smooth( method="gam", formula = y ~ s(x, bs = "ts") ) +
  scale_color_manual(values = c("gray20","red")) +
  labs(x=expression( paste(p["i"]~"(probability of initiation)") ), 
       y = expression( paste(1/tau["i"]~~~"(frequency of initiation, in"~s^{-1},")") ),
       color = "FTS genes") +
  theme_gul + #scale_x_reverse() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.8), legend.position= "bottom")


# [] is the frequency of initiation dependent on mRNA abundance? ----
g <- ggplot( data = subset(SMoPT.data, experiment !="normal" & av.initiation_time > 0), 
             aes(x= log2.mRNA_abundance, y = log2(ratio_initiation) ) ) + 
  geom_vline( xintercept = 0, color ="gray") + geom_hline( yintercept = 0, color ="gray") +
  geom_point(alpha=0.25) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  facet_grid( time ~ experiment ) + 
  #stat_smooth( method="lm" ) +
  stat_smooth( method="gam", formula = y ~ s(x, bs = "ts") ) +
  scale_color_manual(values = c("gray20","red")) +
  labs(x = "Change in mRNA abundance (log2)" , 
       y = "Change in waiting time for intiation (log2)",
       color = "FTS genes") +
  theme_gul + #scale_x_reverse() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.8), legend.position= "bottom")

ggsave(plot=g, filename = paste0(PATH,"part3/SMoPT--x.log2.mRNAab_y.log2.ration.initiation.z_experiment_t.time.png"),
       dpi=250, width = 9.61, height = 5.82 )


# [] length CDS vs n. events -----
g <- ggplot( data = subset(SMoPT.data, experiment !="normal" & av.initiation_time > 0), 
        aes(x= length_CDS, y = n.events ) ) + 
  geom_vline( xintercept = 0, color ="gray") + geom_hline( yintercept = 0, color ="gray") +
  geom_point(alpha=0.25) +
  scale_x_log10( labels = prettyNum ) +
  scale_y_log10( labels = prettyNum ) +
  facet_grid( time ~ experiment ) + 
  stat_smooth( method="lm" ) +
  #stat_smooth( method="gam", formula = y ~ s(x, bs = "ts") ) +
  scale_color_manual(values = c("gray20","red")) +
  labs(x = "CDS length" , 
       y = "Num_of_events") +
  theme_gul + #scale_x_reverse() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.8), legend.position= "bottom")


ggsave(plot=g, filename = paste0(PATH,"part3/SMoPT--x.length_CDS_y.n_events.z_experiment_t.time.png"),
       dpi=250, width = 9.61, height = 5.82 )



# [] length CDS vs initiation time -----
g <- ggplot( data = subset(SMoPT.data, experiment !="normal" & av.initiation_time > 0), 
             aes(x= length_CDS, y = av.initiation_time ) ) + 
  geom_vline( xintercept = 0, color ="gray") + 
  geom_hline( yintercept = 0, color ="gray") +
  geom_point( alpha=0.25 ) +
  scale_x_log10( labels = prettyNum ) +
  scale_y_log10( labels = prettyNum ) +
  facet_wrap(  ~ label, nrow=2 ) + 
  stat_smooth( method="lm" ) +
  scale_color_manual(values = c("gray20","red")) +
  labs(x = "CDS length" , 
       y = "Initiation timing") +
  theme_gul + #scale_x_reverse() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.8), legend.position= "bottom")

ggsave(plot=g, filename = paste0(PATH,"part3/SMoPT--x.length_CDS_y.av.initiation_time.z_experiment_t.time.png"),
       dpi=250, width = 9.61, height = 5.82 )


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------# Adaptation to tRNA pool #---------------------------------------------------------------------------------------------#
#--------------------------------#  s-tAI of FTS genes  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#



# [] Check how adj.tAI varies in genes translated faster or not ----
  # g <- ggplot( data = subset(SMoPT.data, experiment != "normal"), aes(x= adj.tAI.20codons, y = adj.tAI, color= faster.translation)) + 
  #   geom_point() + 
  #   geom_abline(intercept=0, slope=1, lty=2) +
  #   labs(title= "The adaptation index of the first 10 codons\n does not predict the overal adaptation (stAI)", 
  #        x= "stress-adjusted tAI\n(first 10 codons)", 
  #        y = "\nstress-adjusted tAI\n(whole CDS)", fill="FTS genes") +
  #   scale_color_manual(values=c("gray60","red")) +
  #   scale_x_continuous( labels = prettyNum ) +
  #   scale_y_continuous( labels = prettyNum ) +
  #   facet_wrap(  stress ~ faster.translation, ncol = 4 ) + 
  #   theme_gul + theme(legend.position="bottom")
  # 
  # ggsave(plot=g, filename = paste0(PATH,"stress tAI/stress.tAI--x.stAI20codons_y.stAI_z.stress_t.faster_translation.png"),dpi=300)


# [] Check how adj.tAI varies in genes translated faster or not -----
g <- ggplot( subset( SMoPT.data, experiment != "normal"), 
             aes(x= experiment , y = adj.tAI, fill= faster.translation)) + 
  #geom_violin() + 
  geom_boxplot(notch=T) + 
  coord_flip() + labs(title= "stress-adjusted tAI", x ="stress", y = "\nadjusted tAI", fill="FTS genes") +
  scale_y_continuous( labels = prettyNum ) +
  scale_fill_manual(values=c("gray60","red")) +
  facet_wrap( ~ time ) +
  theme_gul + theme(legend.position="bottom")



g <- ggplot( data = ddply( subset(SMoPT.data, experiment!="normal"),
                           .(experiment, time, faster.translation), summarise, adj.tAI = mean(adj.tAI)), 
             aes(x= factor(time), y = adj.tAI, fill= faster.translation)) + 
  geom_bar(stat="identity",position="dodge") + 
  coord_flip() + labs(title= "stress-adjusted tAI", x= "time", y = "\nadjusted tAI", fill="FTS genes") +
  scale_y_continuous( labels = prettyNum ) +
  scale_fill_manual(values=c("gray60","red")) +
  facet_wrap( ~ experiment ) +
  theme_gul + theme(legend.position="bottom")


ggsave(plot=g, filename = paste0(PATH,"part3/adjusted.tAI--x.stAI_y.stress.z_faster_translation.png"),dpi=300)


# [] Check how adj.tAI varies in genes translated faster or not ----
g <- ggplot( data = subset(SMoPT.data, time != 0), 
             aes(x= factor(ratio_total.time.q, levels = 7:1), y = adj.tAI.20codons, fill = factor(ratio_total.time.q) )) + 
  geom_boxplot(notch=T) + 
  labs(title= "Contribution of the 20 first codons to translation speed", 
       x = "relative increase in translation time",
       y = "\nadjusted tAI\n(in 20 first codons)", fill="FTS genes") +
  facet_wrap( time ~ experiment,  ncol = 4 ) +
  scale_x_discrete( labels = c("slowest","","","intermediate","","","fastest")) +
  scale_fill_brewer( "clarity" ) +
  theme_gul + theme(legend.position="none")

ggsave(plot=g, filename = paste0(PATH,"part3/adjusted.tAI--x.stAI20codon_y.stress_early.response.pdf"),
       dpi = 250, width = 11, height = 6)




# [] stress-adjusted tAI in the first 20 codons ----
g <- ggplot( data = ddply(SMoPT.data, .(experiment, time, faster.translation), 
                          summarise, 
                          adj.tAI.20codons = mean(adj.tAI.20codons)), 
             aes(x= paste(experiment,time, sep="\n") , y = adj.tAI.20codons, fill= faster.translation)) + 
  geom_bar(stat="identity",position="dodge") + 
  coord_flip() + labs(title= "stress-adjusted tAI\nin the first 20 codons", x= "",
                      y = "adjusted tAI (2nd-20th codon)", fill="faster translation") +
  #facet_wrap( ~ experiment , ncol =1 ) + 
  scale_fill_manual(values=c("gray60","red")) +
  theme_gul + theme(legend.position="bottom")

ggsave( plot = g, filename = paste0(PATH,"part3/anticodon_enrichment--x.stAI_y.stress_z.faster.translation.png"), 
        dpi=250, width = 6.47, height = 6.47)



# --- attempt to test if the changes in tRNA abundance, reflected by changes in relative adaptiveness are somewhat 
# predictive of the global protein production rate
g <- ggplot(data=subset(SMoPT.data,experiment!="normal"), aes(x=adj.tAI, y = log2( global.protein_synthesis.rate ), 
                                                              group = up_reg )) + 
  geom_point() +
  scale_x_continuous( labels = prettyNum ) +
  stat_smooth(method="lm") +
  labs( x= "stress-adjusted tAI", y = "global protein production\n(log2)") +
  facet_wrap( label ~ up_reg ) + theme_gul


g <- ggplot(data=subset(SMoPT.data,experiment!="normal"), 
       aes(x=adj.tAI, y = log2( ratio_total.time ), group = faster.translation )) + 
  geom_hline( yintercept = 0, color = "gray70", size = 0.5) +
  geom_point(aes(color = faster.translation), size = 1.5, alpha=0.8) +
  scale_x_continuous( labels = prettyNum ) +
  scale_color_manual(values=c("gray40","red")) +
  stat_smooth(method="lm", color = "white") + 
  labs( x= "stress-adjusted tAI", y = "relative translation time\n(log2)") +
  facet_wrap( ~ experiment + time, ncol = 4 ) + theme_gul + theme(legend.position = "none")


ggsave( plot = g, filename = paste0(PATH,"part3/adjusted.tAI--x.stAI--y.total.time.translation_z.stress.time.png"),
        dpi = 250, width= 6.74, height = 4.75)


# [] Increase in initiation frequency correlates with higher protein production rates -----
g <- ggplot(data=subset(SMoPT.data, experiment!="normal"), aes(x= 1/ratio_initiation , y = log2( global.protein_synthesis.rate ) )) + 
  geom_point() +
  stat_smooth(method="gam", formula = y ~ log(x)) +
  scale_x_sqrt() +
  labs( title = "Increase in initiation frequency correlates with higher protein production rates",
         x = expression(paste("frequency of initiation (",tau[i]^s/tau[i]^n,")")),
         y = "log2 protein production rate"
         ) +
  facet_wrap( experiment ~ time, ncol =4  ) + theme_gul

ggsave( plot = g, filename = paste0(PATH,"part3/initiation--x.initiation_freq_y.production_rate_z.stress20min.pdf"),
        dpi = 250, width= 6.74, height = 4.75, useDingbats=F) 


# --------------------------------------------------------------------------- #
g <- ggplot(data= subset(SMoPT.data, !is.na(faster.translation) & experiment !="normal"), aes(x= faster.translation, y = log2(ratio_initiation)  )) + 
  geom_violin(aes(fill=faster.translation)) +
  geom_hline( yintercept= 0 ) +
  scale_fill_manual(values=c("gray60","red")) + 
  labs(x="genes",y="relative waiting time\nfor initiation (log2)") +
  facet_wrap(  experiment ~ time, nrow = 2 ) + theme_gul + 
  theme(legend.position="none", axis.text.x = element_blank() )






#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------# Adaptation to tRNA pool #---------------------------------------------------------------------------------------------#
#--------------------------------# Anticodon Enrichment  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------# FTS genes #---------------------------------------------------------------------------------------------#
#--------------------------------#  s-tAI of FTS genes  #------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

gene.master.melt <- melt(subset(gene.master.table,!is.na(faster.translation)), id.vars = c("ORF", "stress", "experiment", "time", "faster.translation")  )
gene.master.melt$value <- as.numeric(gene.master.melt$value)
selected.features <- c("est.mRNA_abundance", "log2.mRNA_abundance",
                       "length_CDS", "AUGCAI", "IniProb", "global.protein_synthesis.rate")

gene.master.melt$label  <- paste(gene.master.melt$experiment,"\nt=" , gene.master.melt$time, sep="") 
gene.master.melt$label  <- factor( gene.master.melt$label, levels = c(
                                                      "normal\nt=0",
                                                      "diauxic\nt=20",
                                                      "diauxic\nt=120",
                                                      "ox\nt=20",
                                                      "ox\nt=120",
                                                      "osm\nt=20",
                                                      "osm\nt=120",
                                                      "temp\nt=20", 
                                                      "temp\nt=120" )
                            ) 


mannwhitney.tests <- ddply( subset(gene.master.melt, variable %in% selected.features, select = c(stress, faster.translation,variable, value)), .(variable),
      function(x) { 
        # medians
        wilcox.test.batch(x = x)  }
      )
#         stats <- ddply( x, .(stress, faster.translation), summarise, median = median(value), sd = sd(value), n= length(value) )
#         
#         ddply( x, )
#         
#         # wilcoxon tests
#         w <- wilcox.test.batch( x = x )  
#           #pairwise.wilcox.test(x$value, g =  paste( x$stress, x$faster.translation), data = x)
#         
#         # format the whole thing
#         d <- subset( melt(w$p.value), !is.na(value) & as.character(X1) != as.character(X2))  
#         d <- do.call( rbind, apply( d, 1, function(x){ data.frame( t(c(sort(x[1:2]), x[[3]] )), stringsAsFactors=F) } ) )
#         
#         dd <- data.frame( do.call(rbind, strsplit(as.character(d$X1), split = " ")), do.call(rbind, strsplit(as.character(d$X2), split = " ")), 


head( subset(mannwhitney.tests, p.value < 0.05) )

# [] multiplot characterising FTS genes ----
library(grid)
library(gridExtra)

p1 <- ggplot( data = subset(gene.master.melt, variable %in% "est.mRNA_abundance"), 
        aes(x= label, y = ( value ) , fill = faster.translation )) + 
  geom_boxplot(notch=T, outlier.colour = "gray70", outlier.size=1) +
  scale_fill_manual(values=c("gray40","red")) +
  theme_gul + theme(legend.position = "none") +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))
  ) +
  labs( x = "", y = "mRNA abundance") +
  facet_wrap( ~ variable)

p2 <- ggplot( data = subset(gene.master.melt, variable %in% "length_CDS"), 
              aes(x= label, y = ( value/3 ) , fill = faster.translation )) + 
  geom_boxplot(notch=T, outlier.colour = "gray70", outlier.size=1) +
  scale_fill_manual(values=c("gray40","red")) +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))
    ) +
  theme_gul + theme(legend.position = "none") +
  labs( x = "", y = "length (aa)") +
  facet_wrap( ~ variable)


p3 <- ggplot( data = subset(gene.master.melt, variable %in% "IniProb"), 
              aes(x= label, y =  value , fill = faster.translation )) + 
  geom_boxplot(notch=T, outlier.colour = "gray70", outlier.size=1) +
  scale_fill_manual(values=c("gray40","red")) +
  theme_gul + theme(legend.position = "none") +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))
                 )  +
  labs( x = "", y = expression(p[i])) +
  facet_wrap( ~ variable)



p4 <- ggplot( data = subset(gene.master.melt, variable %in% "global.protein_synthesis.rate" & experiment !="normal"), 
              aes(x= label, y =  value , fill = faster.translation )) + 
  geom_boxplot(notch=T, outlier.colour = "gray70", outlier.size=1) +
  scale_fill_manual(values=c("gray40","red")) +
  theme_gul + theme(legend.position = "none") +
  scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))
  )  +
  labs( x = "", y = "") +
  facet_wrap( ~ variable)

require(grid)
require(gridExtra)
g <- grid.arrange(p1, p2, p3, p4, ncol = 2)

png(filename = paste0(PATH,"part3/FST.genes--multiplot.png"), res =250,units = "in", 
    width = 9.6, height = 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


pdf(file = paste0(PATH,"part3/FST.genes--multiplot.pdf"), 
    width = 9.6, height = 8)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()