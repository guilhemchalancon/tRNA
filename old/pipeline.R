# pipeline


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

setwd("~/Documents/MRC16")
source("scripts/commands.R")
source("scripts/SMoT.R")


source("scripts/prepare.tRNAab.R") # imports data from XLSX spreadsheets and build tRNAab_quantified.foldchange
source("scripts/rosetta.tables.R"); rm(T, TAB) # builds the rosetta tables
source("scripts/prepare.codon.table.R")
source("scripts/prepare.anticodon.table.R")
source("scripts/prepare.gene.table.R")

source("scripts/tRNA_response_to_stress.analysis.R")
