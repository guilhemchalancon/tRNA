# Analysis of the stochastic simulation by Marc using Josua Plotkin's source code (Cell theory paper of whole-cell modeling of translation)
setwd("~/Documents/Research/MRC16")
source("scripts/commands.R")
source("scripts/SMoT.R")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

# LOAD DATA 

load("results/SciSignalling 2018 revisions/Datasets/SMoPT.table_May2018.Rd")
study.translation.kinetics <- data.table::fread("results/master tables/study.translation.kinetics.txt", header=T)

go.slim.ORFs <- read.table(file = "data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
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


test <- subset( automate.fisherexact( data = subset(study.translation.kinetics, experiment == "temp" & time == 20), classification.table = go.slim.ORFs,
                              expression = 'gain.rank > 2', keys = O ), 
        p.value < 0.05 & diagnostic %in% c("enriched (significantly)","enriched (perhaps)")  )

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

GO.results <- list()

GO.results$rank.tAI.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'gain.rank > 2', keys = O)  )
})

GO.results$rank.tAI.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'gain.rank < 0.5', keys = O) )
})

# ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), summarise, sum( gain.stAI > 1.25 ))
# 
# GO.results$tAI.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
#   subset( automate.fisherexact( data = x, expression = 'gain.stAI > 1.25', keys = O) )
# })

# 
# GO.results$tAI.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
#   subset( automate.fisherexact( data = x, expression = 'gain.stAI < 0.75', keys = O) )
# })
# 

# ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), summarise, sum( ratio_events > 1.5 ))

GO.results$prod.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'ratio_events > 1.5', keys = O)  )
})

GO.results$prod.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'ratio_events < 0.5', keys = O) )
})



GO.results$ap.translation.rate.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'ratio_apparent.translation.rate > 1', keys = O)  )
})

GO.results$ap.translation.rate.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'ratio_apparent.translation.rate < 1', keys = O)  )
})


#ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), summarise, sum( variation_initiation_frequency > 0.05 ))

GO.results$freq.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'variation_initiation_frequency > 0.05', keys = O)  )
})

GO.results$freq.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'variation_initiation_frequency < -0.25 ', keys = O))
})

# ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), summarise, sum( variation_speed > 0.1 ))
# ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), summarise, sum( variation_speed < -0.25 ))

GO.results$speed.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'variation_speed > 0.1', keys = O)  )
})

GO.results$speed.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'variation_speed < -0.25', keys = O ) )
})

# ---- Additions April 2nd 2015

# load("results/GO_results.Rd")

GO.results$increased.tAI <-  ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'tAI.profile == "increased\nadaptation" ', keys = O ) )
})
GO.results$increased.tAI$x <- factor(GO.results$increased.tAI$x, label = "increased tAI")  


GO.results$decreased.tAI <-  ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'tAI.profile == "decreased\nadaptation" ', keys = O ) )
})
GO.results$decreased.tAI$x <- factor(GO.results$decreased.tAI$x, label = "decreased tAI")  


GO.results$TFS.genes <-  ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'ratio_apparent.translation.rate > 1', keys = O ) )
})


GO.results$mRNA.down <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'log2.mRNA_abundance < -1', keys = O) )
})

GO.results$mRNA.up <- ddply( subset(study.translation.kinetics, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'log2.mRNA_abundance > 1', keys = O) )
})



save(list = "GO.results", file = "results/GO_results.Rd")
GO.results <- lapply(GO.results, unique)
GO.results.table <- ldply( GO.results )
write.table(GO.results.table, file = "results/GO_results_table.txt",row.names=F, quote=T, sep="\t")

# MAY 2018
GO.results.table <- data.table::fread("results/GO_results_table.txt", sep="\t",header=T)

GO_results_summary <- GO.results.table[!is.na(adjusted.p), c(".id","experiment","time","x","y","type","p.hyper","adjusted.p","diagnostic","high.confidence" )]
data.table::setnames(GO_results_summary, old = c(".id","x","y"), new = c("gene_group","criterion","GO_slim"))

write.table(GO_results_summary, file = "results/SciSignalling 2018 revisions/GO_results_summary.csv",sep=",",row.names = FALSE)

GO_results_summary[,list(mean(adjusted.p)), by="gene_group"]
