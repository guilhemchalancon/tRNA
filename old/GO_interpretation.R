setwd("~/Documents/MRC16/")


go.slim.ORFs <- read.table("data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
colnames(go.slim.ORFs)[2] <- "name"

complexes <- table(subset(go.slim.ORFs, GO.type == "C")$GO.description)
processes <- table(as.character(subset(go.slim.ORFs, GO.type == "P")$GO.description))

summary(processes)

# list of all Gene ontologies (Cell compartment,Process,Function)

O <- data.frame(subset(unique(go.slim.ORFs[,c("GO.description","GO.type")]), !GO.description %in% c("not_yet_annotated","biological_process","molecular_function","cellular_component") ), stringsAsFactors=F)
O$GO.description <- as.character(O$GO.description)
O$GO.type <- as.character(O$GO.type)
rownames(O) <- NULL
O <- subset(O, GO.type %in% "P" )# | GO.description %in% names(complexes[complexes>15]) )

grep("response", O$GO.description, value = T)


GO.results.table <- data.frame(read.csv("results/GO_results_table.txt",sep="\t",header=1, as.is = T, check.names=F ), stringsAsFactors=F)
GO.results.table$p.value <- as.numeric(as.character(GO.results.table$p.value))
GO.results.table$p.hyper <- as.numeric(as.character(GO.results.table$p.hyper))

# significant results ranked by stress (experiment) and time point
results.of.interest <-  subset(GO.results.table, p.hyper <0.01 & diagnostic == "enriched", select= c(experiment, time, .id, y,odds, CI.low, sensitivity,specificity,accuracy, p.value ) )
results.of.interest$p.value <- format.pval( results.of.interest$p.value, digits = 1 )

dlply(results.of.interest, .(.id), arrange, experiment, time )
dlply(results.of.interest, .(experiment), arrange, experiment, time )


dlply( subset(results.of.interest, experiment == "ox"), .(.id), arrange, experiment, time )




# ------ DIAUXIC ----------

"GO:0045471"
nrow(subset(go.slim.ORFs, GO.description=="invasive growth in response to glucose limitation"))
subset(go.slim.ORFs, GO.description=="invasive growth in response to glucose limitation" & name %in% subset(study.translation.kinetics, experiment =="diauxic")$name )

response.to.diauxic.shift <-  arrange( subset(GO.results.table, y=="invasive growth in response to glucose limitation", 
                                               select= c(experiment, time, .id, y, n.TT, n.TF, n.FT, odds, sensitivity,specificity,accuracy, p.value, p.hyper )
), experiment, time, plyr::desc(n.TT/n.FT))
response.to.diauxic.shift$p.value <- format.pval( response.to.diauxic.shift$p.value, digits = 1 )

subset(response.to.diauxic.shift, experiment == "diauxic")

subset(study.translation.kinetics, experiment == "diauxic" & time == 120 & ratio_apparent.translation.rate > 1 & 
         name %in% subset(go.slim.ORFs, GO.description=="invasive growth in response to glucose limitation")$name
         )


arrange(subset(GO.results.table, experiment == "diauxic" & p.hyper < 0.01 & .id == "prod.up", select= c(.id, experiment,time, y, n.TT, sensitivity, adjusted.p)  ), plyr::desc(sensitivity))
arrange(subset(GO.results.table, experiment == "diauxic" & p.hyper < 0.01 & .id == "rank.tAI.up", select= c(.id, experiment,time, y, n.TT, sensitivity, adjusted.p)  ), plyr::desc(sensitivity))


subset(study.translation.kinetics,  experiment == "diauxic" & name == "ICY2", select=c(log2.mRNA_abundance, gain.stAI, tAI.profile, global.protein_synthesis.rate) )
subset(study.translation.kinetics,  experiment == "diauxic" & name == "GSY1", select=c(log2.mRNA_abundance, gain.stAI, tAI.profile, global.protein_synthesis.rate) )



# ------ OX ----------
response.to.oxidative.stress <-  arrange( subset(GO.results.table, y=="response to oxidative stress", 
                                              select= c(experiment, time, .id, y, n.TT, n.TF, n.FT, odds, CI.low, CI.high, sensitivity,specificity,accuracy, p.hyper )
), experiment, time, plyr::desc(n.TT/n.FT))

response.to.oxidative.stress$p.hyper <- format.pval( response.to.oxidative.stress$p.hyper, digits = 1 )

subset(response.to.oxidative.stress, experiment == "ox")

length(subset(go.slim.ORFs, GO.description == "response to oxidative stress")$name)
noquote(sort(as.character(subset(go.slim.ORFs, GO.description == "response to oxidative stress")$name)))


# ------ OSM ----------
response.to.osmotic.stress <-  arrange( subset(GO.results.table, y=="response to osmotic stress", 
                                               select= c(experiment, time, .id, y, n.TT, n.TF, n.FT, odds, sensitivity,specificity,accuracy, p.value )
                                               ), 
                                        experiment, time, plyr::desc(n.TT/n.FT))
response.to.osmotic.stress$p.value <- format.pval( response.to.osmotic.stress$p.value, digits = 1 )

subset(response.to.osmotic.stress, experiment == "osm")

 length(subset(go.slim.ORFs, GO.description == "response to osmotic stress")$name)
 noquote(sort(as.character(subset(go.slim.ORFs, GO.description == "response to osmotic stress")$name)))



# ------ Temp ----------
response.to.temperature.stress <-  arrange( subset(GO.results.table, y=="response to heat", 
                                                 select= c(experiment, time, .id, y, n.TT, n.TF, n.FT, odds, sensitivity,specificity,accuracy, p.value )
), experiment, time, plyr::desc(n.TT/n.FT))


response.to.temperature.stress$p.value <- format.pval( response.to.temperature.stress$p.value, digits = 1 )

length(subset(go.slim.ORFs, GO.description == "response to heat")$name)
noquote(sort(as.character(subset(go.slim.ORFs, GO.description == "response to osmotic stress")$name)))
subset(response.to.temperature.stress, experiment == "temp")



# Quick check, low hope: could I back up the estimations of protein abundance with empirical data?

require(StingRay)
data(XPR)
diff.protab <- with( subset(XPR$protein.abundance_noise_plasticity.Breker2013, gene %in% NAMES$systematic.name), 
      data.frame( fNAME(gene), protein_abundance_H2O2, protein_abundance_ref1, protein_abundance_H2O2/protein_abundance_ref1))
colnames(diff.protab) <- c("name","protab.ox","protab.normal","protab.ratio")

oxidative.stress.data <- merge( subset(gene.master.table, time == 20 & experiment == "ox"), 
                                diff.protab, by="name" )


with(oxidative.stress.data , cor.test( protab.ox, mRNA_abundance, method= "spearman" ) )
with(oxidative.stress.data, cor.test( protab.ratio, global.protein_synthesis.rate, method= "spearman" ) )
with(oxidative.stress.data, cor.test( gain.stAI, protab.ratio, method = "spearman" ) )
with(oxidative.stress.data, cor.test( log2.mRNA_abundance, protab.ratio, method = "spearman" ) )



ggplot(oxidative.stress.data, aes(x= log(global.protein_synthesis.rate), log(protab.ratio) )) + geom_point()

with(subset(oxidative.stress.data, global.protein_synthesis.rate < 5 & protab.ratio < 4 ),
     cor.test( protab.ratio, global.protein_synthesis.rate, method= "spearman" ) )

