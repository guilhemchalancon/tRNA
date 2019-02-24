setwd("~/Documents/MRC16/")
require(data.table)

anticodon.usage.table <- fread("results/master tables/anticodon.usage.table.txt",header=T)

go.slim.ORFs <- read.table("data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
colnames(go.slim.ORFs)[2] <- "name"


O <- data.frame(subset(unique(go.slim.ORFs[,c("GO.description","GO.type")]), !GO.description %in% c("not_yet_annotated","biological_process","molecular_function","cellular_component") ), stringsAsFactors=F)
O$GO.description <- as.character(O$GO.description)
O$GO.type <- as.character(O$GO.type)
rownames(O) <- NULL
O <- subset(O, GO.type %in% "P" )# | GO.description %in% names(complexes[complexes>15]) )




# to make a table for Arg-tRNAs
arg.tRNAs <- dlply( subset(anticodon.master.table, aa == "Arg" & experiment == "diauxic", select = c(time, anticodon, total.tRNAab, total.tRNAab.normal ,demand.mRNAab, foldchange)) , .(time))


diauxic.genes <- as.character(subset(go.slim.ORFs, GO.description ==  "invasive growth in response to glucose limitation")$name)


data.arginine.tRNAs <- subset(anticodon.usage.table, anticodon %in% c("ACG","UCU","CCU","CCG") )
data.arginine.tRNAs$diauxic.gene <- with(data=data.arginine.tRNAs, name %in% diauxic.genes )

enrichment.arginine.tRNAs <- arrange(subset( data.arginine.tRNAs, 
                                             select=c(anticodon, name, diauxic.gene, freq.in_CDS, OR.CDS, preference_in.CDS, n.in_CDS )), 
                                     name, anticodon)


ggplot( enrichment.arginine.tRNAs, aes( x = anticodon, y = preference_in.CDS ) ) + geom_boxplot( aes(fill = diauxic.gene), notch=T )

ggplot( enrichment.arginine.tRNAs, aes( x = anticodon, y = freq.in_CDS ) ) + geom_boxplot( aes(fill = diauxic.gene), notch=T )


ggplot( enrichment.arginine.tRNAs, aes( x = anticodon, y = OR.CDS ) ) + geom_boxplot( aes(fill = diauxic.gene), notch=T )


# ggplot( ddply( enrichment.arginine.tRNAs, .(diauxic.gene, anticodon), summarise, preference_in.CDS = median(preference_in.CDS, na.rm=T) ), aes( x = diauxic.gene, y = preference_in.CDS ) ) + geom_histogram( aes(fill = anticodon), stat = "identity", position="fill" )

wc.table <- do.call( rbind, lapply( setNames(c("ACG","CCU","CCG","UCU"), c("ACG","CCU","CCG","UCU")), function(c){ wilcox.test.batch(subset(enrichment.arginine.tRNAs, anticodon == c), grouping.vars = c("diauxic.gene"), value = "preference_in.CDS") }))

ggplot( ddply( enrichment.arginine.tRNAs, .(diauxic.gene, anticodon), summarise, preference_in.CDS = median(preference_in.CDS, na.rm=T) ), 
        aes( x = anticodon, y = preference_in.CDS ) ) + 
  geom_histogram( aes(fill = diauxic.gene), stat = "identity", position="dodge" ) +
  scale_y_continuous(label = percent) + 
  labs(y="preference in CDS\n(compared to other Arg-tRNAs)")


enrichment.arginine.tRNAs.long <- dcast(enrichment.arginine.tRNAs, name + diauxic.gene ~ anticodon, value.var = "OR.CDS" )
subset(enrichment.arginine.tRNAs.long, diauxic.gene == T & (ACG + UCU <= CCG + CCU) & CCG > 1 & CCU > 1  )


# 
subset(enrichment.arginine.tRNAs.long, name %in% c("CAT8","AAP1","SSA1","SSA3","SSA4","TFS1","PGM2","GAC1")  )
 copy2clipboard( subset(enrichment.arginine.tRNAs.long, CCG > 5 | CCU > 5   )$name )

###3


require(scales)
ggplot( ddply( enrichment.arginine.tRNAs, .(diauxic.gene, anticodon), summarise, freq.in_CDS = median(freq.in_CDS, na.rm=T) ), 
        aes( x = anticodon, y = freq.in_CDS ) ) + geom_histogram( aes(fill = diauxic.gene), stat = "identity", position="dodge" ) + 
  scale_y_continuous(labels=percent) 


d <- dcast(data.arginine.tRNAs, name ~ anticodon,  value.var = "n.in_CDS" )

weirdos <- subset( d, CCG + CCU > ACG + UCU  )

bruckmann2008.diff_protab_diauxic <- toupper(as.character(read.table("data/diauxic_shift_diffprotab/Bruckmann2008_diauxic_diffprotab.txt")[,1]))

intersect(bruckmann2008.diff_protab_diauxic,weirdos$name)
intersect(bruckmann2008.diff_protab_diauxic,data.arginine.tRNAs$name)

# 
# total.Arg <- ddply( enrichment.arginine.tRNAs, .(anticodon, diauxic.gene), summarise, n.in_CDS = sum(n.in_CDS) )
# m <- cast( subset(total.Arg, anticodon %in% c("UCU","CCU") ), anticodon ~ diauxic.gene, value = "n.in_CDS" )
# rownames(m) <- m[,1] 
# m <- m[,-1]
# chisq.test(m)  
# 
# m <- cast( subset(total.Arg, anticodon %in% c("UCU","CCG") ), anticodon ~ diauxic.gene, value = "n.in_CDS" )
# rownames(m) <- m[,1] 
# m <- m[,-1]
# chisq.test(m)  
# m
# 
# m <- cast( subset(total.Arg, anticodon %in% c("UCU","ACG") ), anticodon ~ diauxic.gene, value = "n.in_CDS" )
# rownames(m) <- m[,1] 
# m <- m[,-1]
# chisq.test(m)  
enrichment.arginine.tRNAs$is.present <- enrichment.arginine.tRNAs$n.in_CDS > 0
test <- subset(enrichment.arginine.tRNAs, anticodon == "CCG")

with(test, chisq.test(x = is.present, y = diauxic.gene) )


# 
# proportions.Arg <- with(arg.tRNAs[["120"]], data.frame( anticodon = anticodon, 
#                                      normal =  round(total.tRNAab.normal/sum(total.tRNAab.normal),2),
#                                      diauxic.120 = round(total.tRNAab/sum(total.tRNAab),2)
#                                      ) )


proportions.Arg <- subset(anticodon.master.table, time == 120 & experiment == "diauxic" & aa == "Arg", select=c(anticodon, relative.availability.normal, relative.availability))

