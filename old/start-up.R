# Installing and playing with GenomicFeatures to later on extract Condon Usage Indexes
# 
source("http://www.bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
#

library("GenomicFeatures")
source("~/Documents/MRC16/scripts/refClass.R") # Problem with initRefFields, likely due to absence of Rcpp
source("~/Documents/Scripts/GenomicFeatures.R")
source("~/Documents/Scripts/IRanges_overview.R")
source("~/Documents/Scripts/IRanges_Tricks.R")