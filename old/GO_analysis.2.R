# Analysis of the stochastic simulation by Marc using Josua Plotkin's source code (Cell theory paper of whole-cell modeling of translation)
setwd("~/Documents/MRC16")
source("scripts/commands.R")
source("scripts/SMoT.R")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

# LOAD DATA 

load("results/master tables/SMoPT.table.Rd")
study.translation.kinetics <- read.table("results/master tables/study.translation.kinetics.txt", header=T)

go.slim.ORFs <- read.table("data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
colnames(go.slim.ORFs)[2] <- "name"

GO.results.table <- read.table(file = "results/GO_results_table.txt",header=1, sep="\t")

GO.SLIMs <- data.frame(subset(unique(go.slim.ORFs[,c("GO.description","GO.type")]), !GO.description %in% c("not_yet_annotated","biological_process","molecular_function","cellular_component") ), stringsAsFactors=F)
rownames(GO.SLIMs) <- NULL
sort(GO.SLIMs$GO.description)


head(GO.results.table)
p.value < 0.05 & diagnostic %in% c("enriched (significantly)","enriched (perhaps)")

stress.responses <- c("response to chemical","response to osmotic stress", "response to oxidative stress", "response to starvation","response to heat")
subset(GO.results.table, y %in% stress.responses )

as.character(GO.results.table$sensitivity)

ggplot( data = subset(GO.results.table, y %in% stress.responses ), aes( x = depercent(sensitivity), y = depercent(specificity) )  ) + geom_point()


# focus on genes with annotations for response to stress only
ggplot( data = subset(GO.results.table, y %in% stress.responses ), aes( x = depercent(sensitivity), y = -log10(as.numeric(p.value)) )  ) + 
  geom_point( aes(color = experiment, shape = factor(time) ) ) + scale_y_sqrt() +
  geom_hline( yintercept = -log10(0.05), lty = 3 ) +
  facet_wrap(  ~ .id  ) 


ddply( subset(GO.results.table, y %in% stress.responses ), .(experiment, time, .id), nrow)


subset(GO.results.table, y %in% stress.responses & .id == "speed.up" & min.cell >5)


# all genes
ggplot( data = subset(GO.results.table, min.cell>5 ), aes( x = depercent(sensitivity), y = -log10(as.numeric(adjusted.p)) )  ) + 
  geom_point( aes(color = experiment, shape = factor(time) ) ) + scale_y_sqrt() +
  geom_hline( yintercept = -log10(0.05), lty = 3 ) + scale_x_continuous(labels=percent) +
  labs(x="recall rate", y="significance (-log10 p)") +
  facet_wrap(  ~ .id, nrow=2  ) + theme(legend.position="bottom")

subset(GO.results.table, min.cell >5 & .id == "rank.tAI.up" & depercent(sensitivity) > 0.2)
