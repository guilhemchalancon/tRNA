require(seqinr)
require(plyr)
require(YET)

setwd("~/Documents/MRC16/")

source("scripts/commands.R")
load("data/Rdata/orf_all.codon.usage_eff.Rda")
load("data/Rdata/orf_all.codon.usage_freq.Rda")
load("data/Rdata/orf_all.codon.usage_rscu.Rda")

load("data/Rdata/gash.all.Rda")
gash.q20 <- quantilize(gash.all,bin.number=20,name="gash")
gash.top20 <- select.bins(gash.q20,c("0-5","5-10","10-15","15-20"))

extract.count <- lapply(gash.top20, function(x) uco.eff[x,])
extract.count <- extract.count [conditions_to_keep]

extract.freq <- lapply(gash.top20, function(x) uco.freq[x,])
extract.freq <- extract.freq [conditions_to_keep]

extract.rscu <- lapply(gash.top20, function(x) uco.rscu[x,])
extract.rscu <- extract.rscu [conditions_to_keep]

clusters <- split(names(clusters.conditions),f=clusters.conditions)

meanNA <- function(x){mean(x,na.rm=T)}
sumNA <- function(x){sum(x,na.rm=T)}

create.codon.df <- function(clusters, data, code=c("severe","mild"), label=c("codon","score"), Function=meanNA){
  
  set1 <- lapply( data[clusters[[1]]], function(x) {colwise(.fun=Function)( as.data.frame( x) )} ); set1 <- do.call(what=rbind,set1)
  set2 <- lapply( data[clusters[[2]]], function(x) {colwise(.fun=Function)( as.data.frame( x) )} ); set2 <- do.call(what=rbind,set2)
  set1 <- as.data.frame(set1); cond1 <- rownames(set1)
  set2 <- as.data.frame(set2); cond2 <- rownames(set2)
  set1 <- melt(set1)
  set2 <- melt(set2)
  set1$condition <- cond1
  set2$condition <- cond2
  set1$type <- code[1]
  set2$type <- code[2]
  
  colnames(set1)[1:2] <- label
  colnames(set2)[1:2] <- label
  return(rbind(set1,set2))
}

COUNT <- create.codon.df(clusters, extract.count, label=c("codon","count"),Function=sum)
RSCU <- create.codon.df(clusters, extract.rscu, label=c("codon","RSCU"))
FREQ <- create.codon.df(clusters, extract.freq, label=c("codon","frequency"))

RSCU$codon <- factor(RSCU$codon, levels=names(CAI[order(CAI)]), ordered= T)
FREQ$codon <- factor(FREQ$codon, levels=names(CAI[order(CAI)]), ordered= T)

RSCU.nonstop <- subset(RSCU, codon !=c("tga","taa","tag"))
FREQ.nonstop <- subset(FREQ, codon !=c("tga","taa","tag"))


require(ggplot2)

p <- ggplot( data = RSCU.nonstop, aes(x= codon, y=RSCU, colour=factor(type), group=1) )
p + geom_smooth(aes(colour=factor(type), group=factor(type)))
p + geom_smooth(aes(colour=factor(type), group=factor(type))) + opts(axis.text.x=theme_text(angle=90))

p <- ggplot( data = FREQ.nonstop, aes(x= codon, y=frequency, colour=factor(type), group=1) )
p + geom_smooth(aes(colour=factor(type), group=factor(type))) + opts(axis.text.x=theme_text(angle=90))

c <- ddply(COUNT, ~ type + codon , .fun=function(x){ mean(x$count) })
c$AA <- apply(c, 1, function(x){ translate(s2c(as.character(x[2]))) })
colnames(c)[3] <- "count"
total.counts <-  ddply(c, ~ type + AA, .fun=function(x){ sum(x$count)} )
colnames(total.counts)[3] <- "total.counts"
c <- merge(c, total.counts, by=c("type","AA"))

f <- ddply(FREQ, ~ type + codon , .fun=function(x){ mean(x$frequency) })
s <- ddply(RSCU, ~ type + codon , .fun=function(x){ mean(x$RSCU) })
ss <- split(s, f=s$type)

# The software seems to take the data from openwetware.org/wiki/
#   
#   After playing a little bit with the software, I found the codon usage table for E. coli, which contains the following information:
#   
#   Codon
# Amino acid (3 letters)
# Amino acid (1 letter)
# Fraction
# Frequency (per 1000)
# Number of codons
# RSCU
# RSCUmax
# Wi
# 
# Basically, given a number of CDSs, they use the total number of codons found in all CDSs (1603901 codons for E. coli) and the number of hits for each codon (e.g 32529 for GCA, 53984 for GCG, 40914 for GCC and 24609 for GCT; 152036 total) then they compute the fraction (e.g. for GCA 53984/152036=0.214) and the frequency per 1000 (e.g for GCA 53984/1603901=20.28 -- I'm not sure about this, I get similar but not identical numbers).
# For each codon they use also the RSCU, RSCU max (the maximum for the given amino acid) and, finally, Wi is the RSCU/RSCUmax.

table01 <- merge(x=s,y=f,by=c("type","codon"))
table02 <- merge(x=table01,y=c,by=c("type","codon"))
colnames(table02)[3:4] <- c("RSCU","frequency")
#table02$AA <- sapply(table02$codon, function(x){ translate(s2c(as.character(x))) })
#table02$aa <- sapply(table02$codon, function(x){ aaa(translate(s2c(as.character(x)))) })
table02$aa <- aaa(table02$AA)
table02$codon  <- toupper(table02$codon)
table02 <- table02[,c("type","codon","aa","AA","frequency","count","RSCU","total.counts")]
table02$max.RSCU <- apply(table02, 1, function(x){ max(table02[ which(table02$AA %in% x[4] & table02$type == x[1] ), ]$RSCU) } )

table02 <- arrange(table02,type,AA) # this is fantastic!!!

table02$frequency <- table02$frequency*1000
table02$Wi <- table02$RSCU/table02$max.RSCU
table02$fraction <- table02$count/table02$total.counts
table02 <- table02[, c("type","codon","aa","AA","fraction","frequency","count","RSCU", "max.RSCU","Wi") ]

table03 <- split(table02, f=table02$type)

table03 <- lapply(table03, function(y){ 
  y <- cbind(y,  t( apply(y, 1, function(x) {range(which(y$aa==x[3] & y$type == x[1])) }) ) );
  colnames(y)[ (ncol(y) -1):ncol(y)]  <-  c("start","end")
  y[, c(5:10) ] <- format(y[, c(5:10) ], digits=3)
  write.table(y[2:ncol(y)], file=paste("results/gfp sequence optimization/profile.table_",unique(y[1]),".txt",sep=""), quote=F, sep=",",row.names=T )
  return(y)
})

save(table03, file= "results/gfp sequence optimization/profile.tables.Rda")
