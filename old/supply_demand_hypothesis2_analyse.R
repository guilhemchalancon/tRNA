setwd("~/Documents/MRC16")
require(data.table)
require(Biostrings)
require(seqinr)
require(StingRay)
require(plyr)
require(reshape)
source("scripts/commands.R")
source("scripts/SMoT.R")

#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#
# ---- LOAD DATA ----
#-#############-----------##############-----------##############-----------##############-----------##############-----------#############-#

rosetta.anticodons <- read.table(file="data/info/rosetta.anticodons.txt",header=1)

# | processed in supply_demand_hypothesis2.R
anticodon.economy <- read.table(file = "data/anticodon-codon balance/anticodon.economy.txt", header=1)
#anticodon.economy.up_reg <- read.table(file = "data/anticodon-codon balance/anticodon.economy.txt", header=1)
master.table.reduced <- data.table( read.table("results/SMoPT/master.table.reduced.txt", header=1), key=c("ORF","experiment","time") )



#
# Preambule
#
ggplot(rosetta.anticodons, aes(x= aa.with.1.anticodon, y = tGCN)) + geom_point( aes(color=CAI.O))

ggplot(rosetta.anticodons, aes(x= relative.availability, y = tGCN)) + 
  geom_line( data = subset(rosetta.anticodons, relative.availability < 1 ), aes( group = aa), lty=3, color="gray") +
  geom_point( aes(color=CAI.O), size=5) +
  geom_text( aes(label=aa), size=3) + facet_wrap( ~ topology ) + theme_light() 


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Q1 - how often are specific anticodons required (per minute)?
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

data <- subset(anticodon.economy, experiment %in% c("ox","ox","osm") )
data <- anticodon.economy.up_reg
data <- anticodon.economy


ggplot( data = anticodon.economy, aes( x = demand.mRNAab, y = demand.events )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=5,alpha=0.7, aes(color = copy.number) ) + 
  geom_text(aes(label=anticodon), size=2.5) + 
  stat_smooth( method="lm" ) + 
  scale_x_sqrt() + scale_y_sqrt() +
  scale_color_gradient(low="gray90",high="cyan") +
  labs(x="demand\n(number of occurrence throughout the transcriptome)\n",y="demand\n(number of elongation events)\n") +
  theme_bw() + #coord_fixed() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_wrap( experiment ~ time ) # , scales = "free"

# Clearly, anticodon in high copy numbers (and presumbably higher abundance) are read more than once on average


ggplot( data = anticodon.economy, aes( x = demand.events, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=5,alpha=0.7, aes(color = copy.number) ) + 
  geom_text(aes(label=anticodon), size=2.5) + 
  stat_smooth( method="lm" ) + 
  scale_x_sqrt() + scale_y_sqrt() +
  scale_color_gradient(low="gray90",high="cyan") +
  labs(x="demand\n(number of elongation events)\n",y="\nsupply\n(number of free cognate-tRNAs)") +
  theme_bw() + coord_fixed() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_wrap( experiment ~ time, scales = "free" )

# Clearly, anticodon in high copy numbers (and presumbably higher abundance) are read more than once on average


ggplot( data = anticodon.economy, aes( x = copy.number, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=5,alpha=0.7, aes(color = demand.events ) ) + 
  geom_text(aes(label=anticodon), size=2.5) + 
  stat_smooth( method="lm" ) + 
  scale_color_gradient(low="gray90",high="cyan") +
  labs(x="tGCN",y="supply (number of free cognate-tRNAs)",color="demand") +
  theme_bw() + coord_fixed() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_wrap( experiment ~ time, scales = "free" )



ggplot( data = anticodon.economy, aes( x = demand.events, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_line( aes(group=anticodon), lty=1, col="gray" ) +
  geom_point(size=5,alpha=0.7, aes(color=first.pos) ) + 
  geom_text(aes(label=anticodon), size=2.5) + 
  #scale_color_gradient(low="gray90",high="cyan") +
  scale_x_log10() + scale_y_log10() +
  labs(x="usage (elongation events)", y="supply\n(number of free cognate-tRNAs)\n") +
  theme_bw() + 
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_wrap( ~ aa )

ggplot( data = subset(anticodon.economy.up_reg, experiment !="temp"), aes( x = demand.events, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=5,alpha=0.7, aes(color = copy.number) ) + 
  geom_text(aes(label=anticodon), size=2.5) + 
  stat_smooth( method="lm" ) + 
  scale_x_sqrt() + scale_y_sqrt() +
  scale_color_gradient(low="gray90",high="cyan") +
  labs(x="number of elongation events",y="supply (number of free cognate-tRNAs)") +
  theme_bw() + coord_fixed() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_grid( experiment ~ time, scales = "free" )
#facet_wrap( ~ anticodon , scales = "free" )


ggplot( data = master.table.reduced, aes( x = up_reg, y= mean_nTE  )) +
  geom_boxplot(notch=T) +
  theme_bw() + coord_fixed() +
  theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  facet_grid( experiment ~ time, scales = "free" )






ggplot( data = anticodon.economy, aes( x = demand.events, y = supply )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_point(size=4,alpha=0.6, aes(color = copy.number) ) + 
  geom_text(aes(label=anticodon), size=3) + 
  stat_smooth( method="lm" ) + 
  scale_x_sqrt() + scale_y_sqrt() +
  theme_bw() + coord_fixed() +
  facet_grid( experiment ~ time, scales = "free" )




ggplot( data = data, aes( x = log2(supply/demand), y = log2(ratio_tRNAab) )) + 
  geom_hline(yintercept=0, lty=3) +
  geom_vline(xintercept=0,lty=3) +
  geom_point(size=4,alpha=0.4, aes(color = copy.number) ) + 
  geom_text(aes(label=recognised), size=3) + 
  stat_smooth( method="lm" ) + 
  #scale_x_sqrt() + scale_y_sqrt() +
  theme_bw() + coord_fixed() +
  facet_grid( experiment ~ time, scales = "free" )


ggplot( data = data, aes( x = supply, y = log2(ratio_tRNAab) )) + 
  geom_abline(intercept=0, slope=1,lty=3) +
  geom_hline(xintercept=0,lty=3) +
  geom_point(size=4,alpha=0.4, aes(color = copy.number) ) + 
  geom_text(aes(label=recognised), size=3) + 
  stat_smooth( method="lm" ) +
  scale_x_log10() +
  #scale_x_sqrt() + scale_y_sqrt() +
  theme_bw() + coord_fixed() +
  facet_grid( experiment ~ time, scales = "free" )


ggplot( data = data, aes( x = demand.events, y = log2(ratio_tRNAab) )) + 
  geom_hline(yintercept=0, lty=3, color="gray30") +
  stat_smooth( method="lm" ) + 
  geom_point(size=4,alpha=0.4, aes(color = aa) ) + 
  geom_text(aes(label=aa), size=3) + 
  scale_x_log10() +
  theme_bw() + coord_fixed() +
  facet_grid( experiment ~ time, scales = "free" )

# ---------------- CHECK TRPs ------------------ #
classification <- GEN$tandem.repeats[,-c(2:4)]
colnames(classification)[1] <- "ORF"
automate.chi.squared(data =  master.table.reduced, expression = "faster.translation == T", keys = T, classification.table = classification )
# no significant association with TRPs

#------ check delta initiation vs delta mRNA abundance

head(master.table.reduced)



# ggplot( data = master.table.020, aes( x= rand_mRNA*log2.mRNA_abundance - rand_mRNA , y= delta_initiation ) ) + 
#   geom_point() + facet_wrap(~ experiment) + stat_smooth(method="lm") +
#   theme_bw()

# verify that the genes that we defined as up-regulated (top20%) at the transcriptional level do 
# lead to higher protein synthesis rate (ie not how many new proteins are expected to be produced by minute, be it a result of transcription, translation or both)
# ---- 1) with box plots
g <- ggplot( data = subset(master.table.reduced,time>0 & !is.na(up_reg) & !is.na(global.protein_synthesis.rate) ), 
             aes(x=up_reg )  
) +
  geom_boxplot(notch=T, width=0.75, aes(fill=up_reg, y = log10( global.protein_synthesis.rate) )) + 
  #geom_bar(stat="i",width=0.75, aes(fill=up_reg, y = median( global.protein_synthesis.rate) )) + 
  facet_grid( experiment ~ time ) + #
  scale_fill_manual(values=c("gray","cyan")) +
  labs(x="",y="global protein synthesis rate\n", fill="up-regulated transcripts") + 
  theme_bw() + theme( legend.position = "bottom")

ggsave(plot=g, file="results/SMoPT/up_reg.global.protein_synthesis.rate.pdf", dpi=400)

# ---- 2) with bar plots
g <- ggplot( data = ddply( subset(master.table.reduced,time>0 & !is.na(up_reg) & !is.na(global.protein_synthesis.rate) ), .(experiment,time, up_reg), summarize, global.protein_synthesis.rate= mean(global.protein_synthesis.rate) ), 
             aes(x=up_reg, y=global.protein_synthesis.rate )  
) +
  geom_bar(stat="identity", width=0.5, aes(fill=up_reg)) + 
  facet_grid( experiment ~ time ) + #
  scale_fill_manual(values=c("gray","cyan")) +
  labs(x="",y="global protein synthesis rate\n", fill="up-regulated transcripts") + 
  theme_bw() + theme( legend.position = "bottom")

ggsave(plot=g, file="results/SMoPT/up_reg.global.protein_synthesis.rate2.pdf", dpi=400)
