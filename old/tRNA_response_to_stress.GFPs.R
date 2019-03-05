PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
source("scripts/commands.R"); source("scripts/SMoT.R")
require(Biostrings)
require(seqinr)

# ---- Stress codon adaptiveness -- computed in prepare.codon.table.R ------
stress.codon.adaptiveness  <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1)
anticodon.master.table <- read.table("results/master tables/anticodon.master.table.txt",header=1)
codon.master.table <- read.table("results/master tables/codon.master.table.txt",header=1)

# >> correspondence between codon, anticodons and amino acids
rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1,sep="\t")
rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",header=1,sep="\t")


# ---- Analyse GFP sequences ----
GFPs <- readDNAStringSet(filepath = "~/Documents/MRC16/data/GFP seq/Scer.GFP.fasta") #GFP_optimised.fasta")

# there's a Valine insertion in position 2 in Marc's sequences
lapply(GFPs, function(x){ c2s(translate(s2c(tolower(as.character(x))))[1:10]) } ) 
lapply(GFPs, function(x){ c2s(translate(s2c(tolower(as.character(x))))) } ) 

# ---- Count codons -----
GFP.codon.freq.matrix <- compute.codon.frequency.genes(gene.names = names(GFPs), seq.set = GFPs, l = "all", start = 2)
GFP.codon.count <-  melt( ldply( setNames(names(GFPs),names(GFPs)), 
                        function(x){ get.uco(GFPs[[x]], method="eff") } ), id.vars = ".id" )
colnames(GFP.codon.count)  <- c("GFP.construct","codon","n.codon")

GFP.anticodon.count <- ddply( merge( GFP.codon.count, subset(rosetta.codons, select=c(codon, anticodon)), by = "codon"), .(GFP.construct, anticodon), mutate, 
                              n.anticodon = sum(n.codon))

# compute change in number of anticodon / CDS compared to original GFP
require(reshape2)
M.GFP.anticodon <- reshape2::dcast(unique(GFP.anticodon.count[,c("anticodon","n.anticodon","GFP.construct")]), anticodon  ~ GFP.construct, value.var = "n.anticodon" )
GFP.var.n.anticodon <- melt( data.frame( anticodon = M.GFP.anticodon$anticodon, M.GFP.anticodon[,names(GFPs)[1:6]] - M.GFP.anticodon$eGFP_normal), id.vars = "anticodon")
colnames(GFP.var.n.anticodon) <- c("anticodon","GFP.construct","diff.anticodon")

#GFP.anticodon.count <- dcast(GFP.anticodon.count, formula = anticodon ~ GFP.construct, value.var = "n.anticodon")


M <-  t(GFP.codon.freq.matrix) 
Md.freq  <- data.frame( codon = rownames(M), (M)[,1:6], row.names = NULL )
Md.delta <- data.frame( codon = rownames(M), (M - M[,"eGFP_normal"])[,1:6], row.names = NULL )
Md <- merge(melt(Md.freq, id.vars = "codon"), melt(Md.delta, id.vars = "codon"), by=c("codon","variable"))
colnames(Md) <- c("codon","GFP.construct","codon.freq","diff.codon.freq")

GFP.d <- merge(Md, subset(rosetta.codons,select=c(codon, aa, anticodon, CAI.O)), by = "codon")
GFP.d <- merge(GFP.d, GFP.anticodon.count, by = c("anticodon","codon","GFP.construct"))
GFP.d <- merge(GFP.d, GFP.var.n.anticodon, by = c("anticodon","GFP.construct"))
GFP.d$GFP.construct <- factor( GFP.d$GFP.construct, levels = paste0("eGFP_", c("diauxic","ox","osm","temp","normal","wt")) )

GFP.count.aa <- ddply( unique(subset(GFP.d, GFP.construct == "eGFP_normal", select=c(aa,anticodon,n.anticodon))), .(aa), summarise, n.aa.original = sum(n.anticodon) )
GFP.d               <- merge(GFP.d, GFP.count.aa, by = "aa")
GFP.d$set <- factor( paste0( GFP.d$aa, "-",GFP.d$anticodon, " (/", GFP.d$n.aa.original, ")"), ordered = T ) 

GFP.a <- merge(GFP.d, 
               subset(anticodon.master.table, 
                      select= c(anticodon, experiment, time, aa.with.1.anticodon, demand.mRNAab, adjusted.tGCN, foldchange, log2.foldchange, total.tRNAab.marc, total.tRNAab)
                      ), 
               by = "anticodon" )
GFP.a$experiment    <- factor( GFP.a$experiment, levels = c("diauxic","ox","osm","temp","normal") )


# -------------- codon relative adaptivenes
GFP.a <- merge(GFP.a, subset(codon.master.table, 
                             select=c(codon, anticodon, experiment, time, w , delta_w, av.elongation_time, ratio_elongation) ), 
               by = c("codon","anticodon", "experiment", "time") )


# -------------- layer
layer <- expand.grid( GFP.construct = c("diauxic","ox","osm","temp","normal"), experiment = c("diauxic","ox","osm","temp","normal") )
layer$optimal <- ifelse(layer$GFP.construct == layer$experiment, T, F)
layer$GFP.construct <- paste0("eGFP_", layer$GFP.construct)
layer <- rbind(layer,layer)
layer$time  <- c(rep(20,20), rep(0,5), rep(120,20), rep(0,5))
layer <- unique(layer)

# ignore codons for now
GFP.a.anticodons               <- unique(GFP.a[, setdiff(colnames(GFP.a),c("codon", "codon.freq","diff.codon.freq","n.codon"))])
GFP.a.anticodons$GFP.construct <- factor( GFP.a.anticodons$GFP.construct, levels = paste0("eGFP_", c("diauxic","ox","osm","temp","normal")) )
GFP.a.anticodons$experiment    <- factor( GFP.a.anticodons$experiment, levels = c("diauxic","ox","osm","temp","normal") )


# dlply( subset(GFP.a.anticodons, GFP.construct != "eGFP_normal" & time != 0), .(aa,time), 
#        function(x){ arrange( x[,c("aa","anticodon","GFP.construct","experiment", "time", "n.anticodon", "diff.anticodon")], experiment, time, GFP.construct)
#                     } 
#        )

# 
# ggplot(data=  subset( anticodon.master.table, (time == 120 ) | (experiment=="normal") ), 
#        aes(x=experiment, y=(total.tRNAab.marc - free.tRNAab)/100000,  fill=CAI.O)) +
#   geom_bar(stat="identity", aes(group=anticodon), position=position_dodge() ) + 
#   geom_text( aes(label = anticodon, group=anticodon, ymax = (total.tRNAab.marc - free.tRNAab)/100000), 
#               size=2.25, hjust=-0.2, color="gray30", angle = 90 ) +
#   scale_fill_manual(values=c( "skyblue3","gray40","orange" ) ) +
#   labs (x="\nstress", y="Bound tRNAs\n(x 100,000)\n", 
#         title="Abundance of bound tRNAs 120min after stress induction\n",
#         fill="anticodon") + 
#   facet_wrap( ~ aa) +
#   theme_gul + theme(axis.text.x=element_text(angle=90, hjust=1), legend.position="bottom")


optimals <- subset(layer, optimal== TRUE)[,1:2]


#---- [] GFP ------
g <- ggplot( data = subset(GFP.a.anticodons, GFP.construct != "eGFP_normal" & n.aa.original > 0 & time == 20 ) ) + 
        # background
        geom_bar( stat = "identity", position="dodge", aes(y =  n.aa.original , x = set, group = aa ), fill = "gray80", alpha = 0.5 ) + 
        geom_bar( stat = "identity", position="dodge", aes(y = -n.aa.original , x = set, group = aa ), fill = "gray80", alpha = 0.5 ) + 
        # results
        geom_bar( stat = "identity", position="dodge", aes(y = diff.anticodon , x = set, fill = CAI.O, group = aa ) ) + 
        geom_hline( yintercept = 0, color = "black", size = 0.5) +
        geom_line( data = subset( GFP.a.anticodons, GFP.construct != "eGFP_normal" & n.aa.original > 0 & time == 20 &
                                    paste(GFP.construct, experiment) == paste(optimals$GFP.construct, optimals$experiment)), 
                   aes(x= set, y=log2.foldchange*10, group=GFP.construct ), color = "red" ) + 
        #geom_text( aes(label = diff.anticodon, group = aa, ymax = diff.anticodon ), position=dodge, size = 2.25, hjust=-0.2 ) +
  #data = subset(layer, optimal == T & GFP.construct != "eGFP_normal")
        #geom_rect( aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "gray20", alpha=0, size=1) +
        scale_y_continuous( labels = prettyNum) +
        scale_fill_manual( values = c("skyblue3","orange") ) +
        labs(x = expression(atop("anticodons",atop("(ordered by AA, and from least to most used in the GFP sequence)",""))), 
             y="variation in number of anticodon") +
        facet_wrap( ~GFP.construct, ncol = 1, drop = T ) + 
        theme( axis.text.x = element_text(angle=90, size = 8)) 
ggsave(plot = g, filename = paste0(PATH,"part1/GFPs_x_log2.tRNAab_y.variation.anticodon.pdf"), 
       dpi=250, width=8.1, height=14.1, useDingbats=FALSE )


# -------------- codon relative adaptivenes
g <- ggplot(data = subset(GFP.a, time !=0 & GFP.construct != "eGFP_normal") ) + 
    geom_hline( yintercept = 0, color = "gray") +
    geom_vline( xintercept = 0, color = "gray") +
  geom_point( aes(x = delta_w , y = diff.anticodon, color = CAI.O ) ) + 
    geom_rect( data = subset(layer, optimal == T & time != 0), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "gray20", alpha=0, size=1) +
    scale_y_continuous( labels = prettyNum) +
    scale_color_manual( values = c("skyblue3","orange") ) +
  labs(x = "change in adaptiveness", y="variation in n. of codon occurrences", title = "sequence optimisation and relative adaptiveness (2)") +
  facet_wrap( experiment + time ~ GFP.construct, ncol = 4 )

ggsave(plot = g, filename = paste0(PATH,"part1/GFPs_x_delta.w_y.variation.anticodon.pdf"), 
       dpi=250, width=8.1, height=14.1, useDingbats=FALSE )


# -------------- codon relative adaptivenes
g <- ggplot(data = subset(GFP.a.anticodons, time !=0 & GFP.construct != "eGFP_normal") ) + 
    geom_hline( yintercept = 0, color = "gray") +
    geom_vline( xintercept = 0, color = "gray") +
  geom_point( aes(x = -log2(ratio_elongation) , y = diff.anticodon, color = CAI.O ) ) + 
    geom_rect( data = subset(layer, optimal == T & time != 0), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "gray20", alpha=0, size=1) +
    scale_y_continuous( labels = prettyNum) +
    scale_color_manual( values = c("skyblue3","gray","orange") ) +
  labs(x = "faster elongation?", y="variation in n. of codon occurrences", title = "sequence optimisation and relative adaptiveness (1)") +
  facet_wrap( experiment + time ~ GFP.construct, ncol = 4 )

ggsave(plot = g, filename = paste0(PATH,"part1/GFPs_x_faster.elongation.y_variation.anticodon.pdf"), 
       dpi=250, width=8.1, height=14.1, useDingbats=FALSE )


g <- ggplot(data = subset(GFP.a.anticodons, time !=0 & GFP.construct != "eGFP_normal") ) + 
  geom_hline( yintercept = 0, color = "gray") +
  geom_vline( xintercept = 0, color = "gray") +
  geom_point( aes(x = av.elongation_time , y = diff.anticodon, color = CAI.O ) ) + 
  geom_rect( data = subset(layer, optimal == T & time != 0), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "gray20", alpha=0, size=1) +
  scale_y_continuous( labels = prettyNum) +
  scale_color_manual( values = c("skyblue3","gray","orange") ) +
  labs(x = "av. elongation time", y="variation in n. of codon occurrences", title = "sequence optimisation and relative adaptiveness (1)") +
  facet_wrap( experiment + time ~ GFP.construct, ncol = 4 )


ggplot(data = subset(anticodon.master.table, time !=0 ) ) + 
  geom_hline( yintercept = 0, color = "gray") +
  geom_vline( xintercept = 0, color = "gray") +
  geom_point( aes(x = log2.foldchange , y = diff.anticodon, color = CAI.O ) ) + 
  geom_rect( data = subset(layer, optimal == T & time != 0), aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "gray20", alpha=0, size=1) +
  scale_y_continuous( labels = prettyNum) +
  scale_color_manual( values = c("skyblue3","gray","orange") ) +
  labs(x = "log2 change in tRNA abundance", y="variation in number of anticodon") +
  facet_wrap( experiment + time ~ GFP.construct, ncol = 4 )


# ---- Compute stAI ----
stress.tAI.GFPs <- unique(
  # exclude stop codons and methionine from calculation (As in dosReis2004)
  arrange( ddply(subset(stress.codon.adaptiveness, !codon %in% c("taa","tga","tag") ),  .(experiment, time), 
                 function(d){ 
                   if( unique(as.character(d$experiment)) == "normal" ){
                       data.frame( ORF = rownames(GFP.codon.freq.matrix), 
                                   adj.tAI   = round(get.tai(x=GFP.codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w"),2),
                                   adj.tAIFC = round(get.tai(x=GFP.codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" ),2)
                       )                      
                   } else {
                     data.frame( ORF = rownames(GFP.codon.freq.matrix), 
                                 adj.tAI   = round(get.tai(x=GFP.codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w" ),2),
                                 adj.tAIFC = round(get.tai(x=GFP.codon.freq.matrix[,-c(11,12,15)], w.data = d, score = "w_FC" ),2)
                     )
                   }
                 } 
  ), 
  ORF, experiment, time))



GFP.stAI.comparison <- list(
  m.20  = matrixify( cast( subset(stress.tAI.GFPs, time !=120) , formula = experiment ~ ORF, value = "adj.tAIFC" ) ),
  m.120 = matrixify( cast( subset(stress.tAI.GFPs, time != 20) , formula = experiment ~ ORF, value = "adj.tAIFC" ) ) 
)

GFP.stAI.ggplot <- ldply(GFP.stAI.comparison, function(x){ melt(as.matrix(x))})
GFP.stAI.ggplot$.id <- factor(GFP.stAI.ggplot$.id, levels=c("m.20","m.120"), labels=c("20 min", "120 min"))


GFP.stAI.comparison.correlation     <- ldply( GFP.stAI.comparison, function(x){ melt(cor(x[complete.cases(x),])) })
# GFP.stAI.comparison.correlation$X1  <- factor(GFP.stAI.comparison.correlation$X1, levels = c("normal","diauxic","ox","osm","temp"))
# GFP.stAI.comparison.correlation$X2  <- factor(GFP.stAI.comparison.correlation$X2, levels = c("normal","diauxic","ox","osm","temp"))
# GFP.stAI.comparison.correlation$.id <- factor(GFP.stAI.comparison.correlation$.id, levels = c("m.20","m.120"), labels = c("t=20","t=120"))


g <- ggplot( data = melt(GFP.stAI.comparison.correlation), aes(x=X1, y=X2, fill=value)) + 
  geom_tile() + 
  geom_text( aes(label = round(value,2) ), size = 3, alpha=0.5, color = "white" ) +
  coord_fixed() + labs(x="", y = "", title = "Correlation of stAI") +
  theme_gul +
  theme( axis.text.x = element_text(angle=90)) +
  facet_wrap( ~ .id , drop = T)


ggplot( data = GFP.stAI.ggplot, aes(x=X1, y=X2, fill=value)) + 
  geom_tile() + 
  geom_text( aes(label = round(value,2) ), size = 6, alpha=0.75, color = "white" ) +
  coord_fixed() + labs(x="", y = "", title = "Correlation of stAI") +
  theme_gul +
  labs( fill = "stAI") +
  theme( axis.text.x = element_text(angle=90)) +
  facet_wrap( ~ .id , drop = T)



ggplot( subset(anticodon.master.table, time !=0), aes(x= experiment, y = log2.foldchange) ) + geom_boxplot() + facet_grid( ~ time )

require(reshape2)
# stress.codon.adaptiveness.mat <- dcast(stress.codon.adaptiveness, codon ~ experiment + time ,value.var = "w" )
# stress.codon.adaptiveness.mat[,-1] <- round(stress.codon.adaptiveness.mat[,-1],2)
# stress.codon.adaptiveness.mat <- melt(stress.codon.adaptiveness.mat)

stress.codon.adaptiveness.mat <- stress.codon.adaptiveness
stress.codon.adaptiveness.mat$variable <-  paste(stress.codon.adaptiveness.mat$experiment,stress.codon.adaptiveness.mat$time, sep="_")
stress.codon.adaptiveness.mat[ which(stress.codon.adaptiveness.mat$time == 0), ]$w_FC <- stress.codon.adaptiveness.mat[ which(stress.codon.adaptiveness.mat$time == 0), ]$w

stress.codon.adaptiveness.mat$codon <- factor(stress.codon.adaptiveness.mat$codon, levels = rev(unique(stress.codon.adaptiveness.mat$codon)) )
stress.codon.adaptiveness.mat$variable <- factor(stress.codon.adaptiveness.mat$variable, 
                                                 levels = c("normal_0",paste(rep(c("diauxic","ox","osm","temp"), each=2), rep(c(20,120),times = 4), sep="_")), 
                                                 labels = c("Non-stress",
                                                            "Diauxic\n20min", "Diauxic\n120min",
                                                            "Oxidative\n20min","Oxidative\n120min",
                                                            "Osmotic\n20min","Osmotic\n120min",
                                                            "Temperature\n20min","Temperature\n120min")  )
stress.codon.adaptiveness.mat <- merge(stress.codon.adaptiveness.mat, subset(rosetta.codons, select=c(aa,codon,CAI.O)), by="codon")
stress.codon.adaptiveness.mat$codon <- factor(stress.codon.adaptiveness.mat$codon, levels = as.character(arrange(rosetta.codons, CAI.O, aa, codon)$codon) )
stress.codon.adaptiveness.mat$codon.set <- factor( paste(stress.codon.adaptiveness.mat$codon,stress.codon.adaptiveness.mat$aa,sep="-"), 
                                                   levels = paste(arrange(rosetta.codons, CAI.O, aa, codon)$codon, arrange(rosetta.codons, CAI.O, aa, codon)$aa, sep="-") )


# codon relative adaptiveness across all conditions
g <- ggplot( data = stress.codon.adaptiveness.mat , aes(x=variable, y=codon.set, fill=w_FC) ) + 
  geom_tile( width = 0.95, height= 1) + labs( title = "codon relative adaptiveness", x="", y="")

ggsave(plot = g, filename = paste0(PATH,"part2/adaptiveness--heatmap.pdf"), dpi=250, width=8.1, height=12.1, useDingbats=FALSE )
ggsave(plot = g, filename = paste0(PATH,"part2/adaptiveness--heatmap.png"), dpi=250, width=8.1, height=12.1)


ggplot( data = subset(stress.codon.adaptiveness.mat, time > 0) , aes(x=variable, y=codon.set, fill= w/(w - delta_w) ) ) + 
  geom_tile( width = 0.95, height= 1) + labs( title = "codon relative adaptiveness", x="", y="") + 
  scale_fill_gradient2( low = "red", mid = "white", high = "blue", midpoint = 1)

ggplot( data = subset(gene.master.table, !is.na(up_reg)) , aes(x=experiment, y=adj.tAIFC, fill=up_reg) ) + geom_boxplot() + facet_grid( ~ time )

# --------------------------------------------------------------------------------------------------------------------------------------- #
# Produce tables of s-tAI (or s-tAIFC) for the optimised GFP constructs (with rows and columns ordered to matched with Figure 2)
matrixify(reshape2::dcast( subset(stress.tAI.GFPs, time == 20 & !ORF %in% c("eGFP_normal","eGFP_wt") ), ORF ~ experiment + time, value.var = "adj.tAI" ))[c("eGFP_ox","eGFP_osm","eGFP_diauxic","eGFP_temp"),c("ox_20","osm_20","diauxic_20","temp_20")]
matrixify(reshape2::dcast( subset(stress.tAI.GFPs, time == 20 & !ORF %in% c("eGFP_normal","eGFP_wt") ), ORF ~ experiment + time, value.var = "adj.tAIFC" ))[c("eGFP_ox","eGFP_osm","eGFP_diauxic","eGFP_temp"),c("ox_20","osm_20","diauxic_20","temp_20")]


matrixify(reshape2::dcast( subset(stress.tAI.GFPs, time == 20 & !ORF %in% c("eGFP_wt") ), ORF ~ experiment + time, value.var = "adj.tAIFC" ))[c("eGFP_normal","eGFP_ox","eGFP_osm","eGFP_diauxic","eGFP_temp"),c("ox_20","osm_20","diauxic_20","temp_20")]

# --------------------------------------------------------------------------------------------------------------------------------------- #

matrixify(reshape2::dcast( subset(stress.tAI.GFPs, time == 120 & !ORF %in% c("eGFP_normal","eGFP_wt") ), ORF ~ experiment + time, value.var = "adj.tAI" ))[c("eGFP_ox","eGFP_osm","eGFP_diauxic","eGFP_temp"),c("ox_120","osm_120","diauxic_120","temp_120")]
matrixify(reshape2::dcast( subset(stress.tAI.GFPs, time == 120 & !ORF %in% c("eGFP_normal","eGFP_wt") ), ORF ~ experiment + time, value.var = "adj.tAIFC" ))[c("eGFP_ox","eGFP_osm","eGFP_diauxic","eGFP_temp"),c("ox_120","osm_120","diauxic_120","temp_120")]



ORDER <- rev(c("oxidative","osmotic","temperature","diauxic"))
LABELS.conditions <- rev(c("ox","osm","temp","diauxic"))
LABELS.strains <- rev(c("eGFP_ox","eGFP_osm","eGFP_temp","eGFP_diauxic"))
                         
                         
data.gfp <- read.csv("data/GFP seq/data_gpf_heatmap.csv",header=1)
data.gfp <- ddply(data.gfp, .(optimised.stain), mutate, scaled.value = range01(scale(value)), ranged.value = (value)/max(value) ) 
data.gfp$condition <- factor(data.gfp$condition, levels= ORDER, labels = LABELS.conditions )
data.gfp$optimised.stain <- factor(data.gfp$optimised.stain, levels= ORDER, labels = LABELS.strains )
data.gfp$label.strain <- gsub(d.GFP$strain, pattern="^(.*)_(.*)$", replacement="\\2")


d.GFP.Pi <- subset(data.gfp, select= c(optimised.stain, value, ranged.value, condition, label.strain))
colnames(d.GFP.Pi) <- c("strain","raw.value","prod.rate","condition","label.strain")

d.GFP.stAI <- subset(stress.tAI.GFPs, time == 20 & !ORF %in% c("eGFP_normal","eGFP_wt") )[,c("ORF","adj.tAI","adj.tAIFC","experiment")]
colnames(d.GFP.stAI) <- c("strain","stAI","stAIFC","condition")

d.GFP <- merge(d.GFP.stAI, d.GFP.Pi, by = c("strain","condition"))
  require(tidyr)



ggplot(d.GFP, aes(x=stAIFC, y=raw.value)) + geom_point(aes(colour=condition)) + theme(legend.position = "bottom")
ggplot(d.GFP, aes(x=stAIFC, y=prod.rate)) + geom_point(aes(colour=condition)) + theme(legend.position = "bottom")
ggplot(d.GFP, aes(x=stAI, y=prod.rate)) + geom_point(aes(colour=condition)) + theme(legend.position = "bottom")

g1 <- ggplot(d.GFP, aes(x=stAI, y=raw.value)) + stat_smooth(method="lm", color = "gray20") + 
  geom_point(aes(colour=condition)) + theme(legend.position = "bottom") + labs(x="s-tAI", y= "production rate" )

ggsave(plot=g1, file="~/Documents/MRC16/Paper/draft/july add-ons/GFP_correlation.s-tAI_prod.pdf", dpi=200, width=3.5, height=3.25)

with(d.GFP, cor.test(raw.value, stAI) )
with(d.GFP, cor.test(raw.value, stAIFC) )

with(d.GFP, cor.test(prod.rate, stAI) )
with(d.GFP, cor.test(prod.rate, stAIFC) )

t2 <- ddply(d.GFP, .(condition), function(x){ test <- cor.test(x$raw.value, x$stAI); 
                                                return( data.frame( R = round(test$estimate, 2), p = test$p.value)) } )


t2b <- ddply(d.GFP, .(strain), function(x){ test <- cor.test(x$raw.value, x$stAI); 
return( data.frame( R = round(test$estimate, 2), p = test$p.value)) } )

round(matrixify(reshape2::dcast(d.GFP, strain ~ condition, value.var = "raw.value")),2)[paste("eGFP", c("diauxic","ox","osm","temp"),sep="_"),c("diauxic","ox","osm","temp")]


g2 <- ggplot(d.GFP, aes(x=stAI, y=raw.value)) + stat_smooth(method="lm", color = "gray20", se = F) + 
  geom_point(aes(colour=condition)) + theme(legend.position = "bottom") + labs(x="s-tAI", y= "production rate" ) + 
  geom_text(aes(label= label.strain), size = 4, vjust=1.5) +
  scale_y_continuous(labels=prettyNum, limits=c(0,1)) +
  scale_x_continuous(labels=prettyNum) +
  facet_wrap( ~ condition, scales = "free") + theme_bw() + theme(legend.position="none")

ggsave(plot=g2, file="~/Documents/MRC16/Paper/draft/july add-ons/GFP_correlation.s-tAI_prod_facets.pdf", dpi=200, width=6.5, height = 6, useDingbats=F)


ggplot(d.GFP, aes(x=stAI, y=raw.value)) + 
  stat_smooth(method="lm", se = F, aes(group=condition, color= condition), fullrange = T) + 
  geom_point(aes(colour=condition)) + theme(legend.position = "bottom") + labs(x="s-tAI", y= "production rate" ) + 
  geom_text(aes(label= label.strain), size = 4, vjust=1.5) +
  scale_y_continuous(labels=prettyNum, limits=c(0,1)) +
  scale_x_continuous(labels=prettyNum) + scale_color_manual(values=c("yellow","purple","green","orange")) +
  theme_gul
