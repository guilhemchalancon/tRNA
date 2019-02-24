#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#     OBJECTIVES    #--------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# Here I try to map transcriptomes highlighting how translation efficiency is distributed across it, depending on mRNA abundance
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

study.translation.kinetics





stAI.cor <- matrixify(reshape2::dcast( study.translation.kinetics, name ~ experiment + time, value.var = "stAI" ))
stAI.cor <- stAI.cor[complete.cases(stAI.cor),]
cor(stAI.cor)

stAIFC.cor <- matrixify(reshape2::dcast( study.translation.kinetics, name ~ experiment + time, value.var = "stAIFC" ))
stAIFC.cor <- stAIFC.cor[complete.cases(stAIFC.cor),]
cor(stAI.cor)




# Stress response data set
response  <- melt(list(
 diauxic = c("GTT1","GPP2","SIT1","SPI1","YJL045W","PDH1","NQM1","SNZ3","SNO4"),  
 stress = read.table("data/stress responses/response_to_stress_annotations.txt",comment.char = "!", sep="\t", skip =1, head=2 )[,1],
 ox = read.table("data/stress responses/response_to_oxidative_stress_annotations.txt",comment.char = "!", sep="\t", skip =1, head=2 )[,1],
 osm = read.table("data/stress responses/response_to_osmotic_stress_annotations.txt",comment.char = "!", sep="\t", skip =1, head=2 )[,1],
 temp = read.table("data/stress responses/cellular_response_to_heat_annotations.txt",comment.char = "!", sep="\t", skip =1, head=2 )[,1]
))

colnames(response) <- c("name","experiment")



# Map the transcriptome on 2D matrices
cartography <- data.table( subset(gene.master.table, 
                                  select=c(name, experiment, time, 
                                           log2.mRNA_abundance, tAI, stAI, FS.stAI, 
                                           stAIFC, rank.stAI, rank.tAI, gain.stAI, 
                                           mRNA_abundance, mRNA_abundance.normal, 
                                           global.protein_synthesis.rate)
                                  )
                           )
cartography <- subset(cartography, !is.na(mRNA_abundance))
cartography[ , rank := rank(mRNA_abundance), by =list(experiment, time) ]
cartography[ , rank.t0 := rank(mRNA_abundance.normal), by =list(experiment, time) ]
setkey(cartography, name, experiment)
cartography$responding <- 0
cartography[response, responding:=1]
cartography$responding <- factor(cartography$responding, levels=c(0,1), labels=c(F,T))
cartography$label <- paste(cartography$experiment, paste(cartography$time,"min"), sep="\n")
cartography$label <- factor(cartography$label, levels= c("normal\n0 min", paste( rep(c("diauxic","osm","ox","temp"), times=2), rep(paste(c(20,120), "min"),each=4), sep="\n")) )
cartography$delta <- cartography$rank.stAI - cartography$rank.tAI



# JULY 2015:
# abundance of the top-40 genes
ab_40 <- sum(head(arrange(subset(cartography, experiment == "normal"),plyr::desc(rank)), 40)$mRNA_abundance)
ab_all <- sum(subset(cartography, experiment == "normal")$mRNA_abundance)

ab_40/ab_all

# cartography2 <- cartography[, list(rank, rank.t0, FS.stAI, stAI, tAI, stAIFC, copy = seq(1:mRNA_abundance)), by=list(experiment, time, name)]
# cartography2$time <- factor(cartography2$time, levels=c(0,20,120))
# cartography2$time <- as.numeric(as.character(cartography2$time))
# setkey(cartography2, name, experiment, time)
# #cartography2[, matrix( data = stAI[order(rank)] ) , by=list(experiment, time)]
# 
# 
# # for gene-specific maps
# additions.for.cartography3 <- data.table(study.translation.kinetics[,c("name","experiment","time","label","gain.rank","tAI.profile")] ,key = c("name","experiment","time"))
# cartography3 <- cartography2[ additions.for.cartography3 ]#,  by=c("name","experiment","time") )
# 
# cartography3$responding <- 0
# cartography3[response, responding:=1] 


#test <- subset(cartography2, experiment == "normal")
              # 
              # transcriptome.map <- function(test, var = "stAI", dim.factor = 60){
              #   n <- nrow(test)
              #   n.matrix <- ceiling(n/dim.factor)*dim.factor
              #   m <- matrix( data = c(test[,get(var)][order(test$rank)], rep(NA,times=n.matrix-n)), nrow = n.matrix/dim.factor )
              #   return(m)
              # }
              # 
              # #ggplot(test,aes(x=est.mRNA_abundance, y=adj.tAIFC)) + geom_point()
              # #ggplot(test,aes(x=est.mRNA_abundance, y=tRNA_adaptation_index)) + geom_point()
              # 
              # m <- list()
              # m$empty <- matrix(NA)
              # m$normal <- transcriptome.map(subset(cartography2, experiment=="normal" & time == 0), var = "stAI")
              # m$diauxic_20 <- transcriptome.map(subset(cartography2, experiment=="diauxic" & time == 20), var = "stAI")
              # m$diauxic_120 <- transcriptome.map(subset(cartography2, experiment=="diauxic" & time == 120), var = "stAI")
              # m$osm_20 <- transcriptome.map(subset(cartography2, experiment=="osm" & time == 20), var = "stAI")
              # m$osm_120 <- transcriptome.map(subset(cartography2, experiment=="osm" & time == 120), var = "stAI")
              # m$ox_20 <- transcriptome.map(subset(cartography2, experiment=="ox" & time == 20), var = "stAI")
              # m$ox_120 <- transcriptome.map(subset(cartography2, experiment=="ox" & time == 120), var = "stAI")
              # m$temp_20 <- transcriptome.map(subset(cartography2, experiment=="temp" & time == 20), var = "stAI")
              # m$temp_120 <- transcriptome.map(subset(cartography2, experiment=="temp" & time == 120), var = "stAI")
              # 
              # 
              # mFC <- list()
              # mFC$empty <- matrix(NA)
              # mFC$normal <- transcriptome.map(subset(cartography2, experiment=="normal" & time == 0), var = "stAIFC")
              # mFC$diauxic_20 <- transcriptome.map(subset(cartography2, experiment=="diauxic" & time == 20), var = "stAIFC")
              # mFC$diauxic_120 <- transcriptome.map(subset(cartography2, experiment=="diauxic" & time == 120), var = "stAIFC")
              # mFC$osm_20 <- transcriptome.map(subset(cartography2, experiment=="osm" & time == 20), var = "stAIFC")
              # mFC$osm_120 <- transcriptome.map(subset(cartography2, experiment=="osm" & time == 120), var = "stAIFC")
              # mFC$ox_20 <- transcriptome.map(subset(cartography2, experiment=="ox" & time == 20), var = "stAIFC")
              # mFC$ox_120 <- transcriptome.map(subset(cartography2, experiment=="ox" & time == 120), var = "stAIFC")
              # mFC$temp_20 <- transcriptome.map(subset(cartography2, experiment=="temp" & time == 20), var = "stAIFC")
              # mFC$temp_120 <- transcriptome.map(subset(cartography2, experiment=="temp" & time == 120), var = "stAIFC")
              # 


              # require(fields)
              # 
              # image.plot(m$normal, legend.width = 0.5, asp =1, axes=F)
              # image.plot(m$diauxic_20)
              # image.plot(m$diauxic_120)
              # image.plot(m$osm_20)
              # image.plot(m$osm_120)
              # image.plot(m$ox_20)
              # image.plot(m$ox_120)
              # image.plot(m$temp_20, axes=F)
              # 
              

# require(StingRay)
# data(GEN)

        # 
        # 
        # 
        # # ---- Figure tAI 
        # set.panel()
        # # Here is quick but quirky way to add a common legend to several plots. 
        # # The idea is leave some room in the margin and then over plot in this margin
        # par( mai=c(1,0.1,1,0.1), mar=c(1,0.1,1,0.1)) #mar=c(0.1,0.1,0.1,0.1) # margin of 4 spaces width at right hand side #oma=c(0,0,0,0.2),
        # set.panel(9,2) # 2X2 matrix of plots
        # 
        # # now draw all your plots using usual image command
        # image.plot( m$normal, col=tim.colors(), axes=F, asp=1, useRaster=TRUE, zlim=c(0,1)) 
        # image.plot( mFC$normal, col=tim.colors(), axes=F, asp=1, useRaster=TRUE, zlim=c(0,1)) 
        # 
        # frame()
        # for (  k in 2:10){
        #   image.plot( mFC[[k]], col=tim.colors(),   axes=F, asp=1,useRaster=TRUE, zlim=c(0,1), main = names(m)[k]  ) #zlim=c(0,1)
        # #  image.plot( mFC[[k]], col=tim.colors(), axes=F, asp=1,useRaster=TRUE, zlim=c(0,1), main = names(m)[k] ) #zlim=c(0,1)
        # }
        # dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_maps.pdf")
        # #dev.off() #  type='cairo', antialias=NULL
        # set.panel()

# http://stackoverflow.com/questions/3955095/creating-raster-images-using-r


#ggplot(test,aes(x=est.mRNA_abundance, y=adj.tAIFC)) + geom_point()
#ggplot(test,aes(x=est.mRNA_abundance, y=tRNA_adaptation_index)) + geom_point()

              # r <- list()
              # r$diauxic_20 <- transcriptome.map(subset(cartography3, experiment=="diauxic" & time == 20), var = "responding")
              # r$diauxic_120 <- transcriptome.map(subset(cartography3, experiment=="diauxic" & time == 120), var = "responding")
              # r$osm_20 <- transcriptome.map(subset(cartography3, experiment=="osm" & time == 20), var = "responding")
              # r$osm_120 <- transcriptome.map(subset(cartography3, experiment=="osm" & time == 120), var = "responding")
              # r$ox_20 <- transcriptome.map(subset(cartography3, experiment=="ox" & time == 20), var = "responding")
              # r$ox_120 <- transcriptome.map(subset(cartography3, experiment=="ox" & time == 120), var = "responding")
              # r$temp_20 <- transcriptome.map(subset(cartography3, experiment=="temp" & time == 20), var = "responding")
              # r$temp_120 <- transcriptome.map(subset(cartography3, experiment=="temp" & time == 120), var = "responding")
              # 
              # image.plot(r$osm_20, asp =1, axes=F,  col = c("white","red"), border="grey50")
              # image(r$osm_120, asp =1, axes=T,  col = c("white","red"))
              # image(r$ox_20, asp =1, axes=T,  col = c("white","red"))
              # image(r$ox_120, asp =1, axes=T,  col = c("white","red"))
              # image(r$temp_20, asp =1, axes=T,  col = c("white","red"))
              # image(r$temp_120, asp =1, axes=T,  col = c("white","red"))


# VORONOI TREE MAPS
install.packages("treemap")
require(treemap)
data(business)
head(business)

treemap(subset(cartography, experiment=="normal" & time == 0), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,0.8),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.normal_0_FS.pdf")


treemap(subset(cartography, experiment=="temp" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,.95),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.temp_20_FS.pdf")



treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,.95),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.ox_20_FS.pdf")


treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance.normal", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,.95),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.ox_20_FS.mRNAab_fixed.pdf")



treemap(subset(cartography, experiment=="ox" & time == 120), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,.95),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.ox_120_FS.pdf")



treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="delta",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        #range=c(0,1),
        type="value")


treemap(subset(cartography, experiment=="temp" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="FS.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0,1),
        type="manual")



treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance", 
        vColor="responding",
        palette=c("white","skyblue"),
        border.col="gray30",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        #range=c(0.18,0.45),
        type="categorical")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.tAI_transcriptome_treemap.ox_20.resp.pdf")



ggplot(data=cartography, aes(x=log2.mRNA_abundance, y=FS.stAI ) ) + 
  geom_point() +
  facet_wrap( ~ label, nrow=2)
  


treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="mRNA_abundance.normal", 
        vColor="log2.mRNA_abundance",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(-4.5,4.5),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.mRNAab_treemap.ox_20_FS.pdf")



treemap(subset(cartography, experiment=="ox" & time == 20), 
        index=c("name"), 
        vSize="global.protein_synthesis.rate", 
        vColor="gain.stAI",
        palette="-RdYlBu",
        border.col="white",
        lowerbound.cex.labels=0.5,
        force.print.labels=F,
        range=c(0.38,0.59),
        type="manual")
dev.copy2pdf(device = quartz, file = "results/stress tAI/stress.prod_treemap.ox_20_FS.pdf")
