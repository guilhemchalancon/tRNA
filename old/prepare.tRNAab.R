# just need to reformat a bit the data frames to prepare figures showing the reproducibility levels etc. 
# (mainly for the first result section of the tRNA-Chapter of my thesis)
setwd("~/Documents/MRC16")
require(XLConnect)
require(plyr)
require(reshape)
require(reshape2)

# integrate data from XLSX spreadsheet
  wb <- loadWorkbook("data/tRNA abundance/original data/tRNA_steady_state_abundances_update.xlsx") # change in case of update
  tRNAab <- list()
  tRNAab$abundance <- readWorksheet(wb, sheet = getSheets(wb)[2],header=T, startRow=2, endRow=44)
  tRNAab$abundance[,2:ncol(tRNAab$abundance)] <- sapply(tRNAab$abundance[,-1], function(x){ round(as.numeric(x),3) } ) #sapply(tRNAab$abundance[,-1], as.numeric)
  tRNAab$conditions <- readWorksheet(wb, sheet = getSheets(wb)[4],header=T)
  tRNAab$conditions$treatment <- factor( tRNAab$conditions$treatment, labels = c("diauxic","ox","osm","temp") )

  tRNAab$conditions$label <- factor(tRNAab$conditions$treatment, levels=unique(tRNAab$conditions$treatment), labels=c("peroxide","heat_shock","ethanol","sorbitol"))
  tRNAab$conditions$experiment <- factor(tRNAab$conditions$treatment, levels=unique(tRNAab$conditions$treatment), labels=c("ox","temp","diauxic","osm"))
  tRNAab$conditions$stress <- paste(tRNAab$conditions$experiment, tRNAab$conditions$time, sep="_")
  colnames(tRNAab$abundance)[1] <- "anticodon"

#------------------------------------------------------------------------------------------------------------------------#

# updated data table from Marc (September 2014) -- susceptible to not contain mistakes (the previous contained some for CAA, at least)
  wb <- loadWorkbook("data/tRNA abundance/original data/Updated_tRNA_table.xlsx") # change in case of update
  tRNAab2 <- list()
  tRNAab2$abundance <- readWorksheet(wb, sheet = getSheets(wb)[1],header=T, startRow=2, endRow=44)[,-2]
  tRNAab2$abundance[,2:ncol(tRNAab$abundance)] <- sapply(tRNAab2$abundance[,-1], function(x){ round(as.numeric(x),3) } )
  tRNAab2$conditions <- tRNAab$conditions
  colnames(tRNAab2$abundance)[1] <- "anticodon"

test <- merge( merge( melt( tRNAab$abundance  ), tRNAab$conditions[,c(1,4,8)], by.x="variable", by.y="reference" ),
       merge( melt( tRNAab2$abundance  ), tRNAab$conditions[,c(1,4,8)], by.x="variable", by.y="reference" ),
       by = c("variable","anticodon","type","stress") 
      )

subset(test, abs(value.x - value.y) > 0.5)
subset( tRNAab2$abundance, anticodon %in% subset(test, abs(value.x - value.y) > 0.5)$anticodon )

#------------------------------------------------------------------------------------------------------------------------#


  # 1 # Wide format used in the construction of all master.tables
  # only use the average fold change for the wide format table
  tRNAab.foldchange <-  tRNAab2$abundance[ , c("anticodon", subset(tRNAab2$conditions, type == "average")$reference ) ]
  # warnings are normal here
  colnames(tRNAab.foldchange)[-1] <- as.character(factor( colnames(tRNAab.foldchange)[-1],  levels = tRNAab$conditions$reference[seq(1,length(tRNAab$conditions$reference),by=2)], 
                                                                                            labels = tRNAab$conditions$stress[seq(1,length(tRNAab$conditions$reference),by=2)]))


  # ##### scale( tRNAab.foldchange,  ) 
  write.table(tRNAab.foldchange, "data/tRNA abundance/tRNAab_quantified.foldchange.txt", sep="\t", quote=F, row.names=F)



#------------------------------------------------------------------------------------------------------------------------#

  # 2 # Matrix form
  # Transform into a matrix (for pheatmap)
  rawM <- as.matrix(tRNAab2$abundance[,-1]) ; rownames(rawM) <- tRNAab2$abundance$tRNA_anticodon
  M <- rawM[complete.cases(rawM),]
  colnames(M) <- tRNAab$conditions$label
  
  # log2 transformation ( since raw value = abundance(x)/abundance(wt) )
  log2M <- log2(M)
  # scaled (centred on 0 and scaled to 1)
  scaledM <- apply(M, 2, scale); rownames(scaledM) <- rownames(M)
  # scalation of the log2 data
  tM <- apply(log2M, 2, scale); rownames(tM) <- rownames(M)
  save(tM,file="results/tRNA.abundance/tM.Rda")


#------------------------------------------------------------------------------------------------------------------------#


# 3 # Long format with mean fold change and standard deviation
  tRNA.data <- merge( melt(tRNAab2$abundance), tRNAab$conditions[,-5], by.x="variable", by.y="reference" )[,-1]
  colnames(tRNA.data) <- c("anticodon","value","experiment","time","type")
  
  marc.data <- dcast( tRNA.data, formula = anticodon + experiment + time ~ type, value.var = "value")
  write.table(marc.data, file="data/tRNA abundance/tRNAab_quantified.mean_sd.txt",sep="\t", quote=F, row.names=F)


#------------------------------------------------------------------------------------------------------------------------#
