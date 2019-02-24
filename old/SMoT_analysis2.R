setwd("~/Documents/MRC16")
source("scripts/SMoT.R") # commands
require(data.table)
require(reshape)
require(plyr)
#Some references on stress response in yeast:
# Clare E. Simpson and Mark P. Ashe. Adaptation to stress in yeast: to translate or not? Biochem. Society Trans. 2012, 40:794-799
# Audrey P. Gasch et al. Genomic expression programs in the response of yeast cells to environmental changes. Mol. Biol. Cell 2000, 11:4241-4257
# Eulàlia de Nadal et al. Controlling gene expression in response to stress. Nature Rev. Genetics 2011, 12:833-845
# Audrey P. Gasch. The environmental stress response: a common yeast response to diverse environmental stresses (Book Chapter)
# Haruo Saito and Francesc Posas. Response to hyperosmotic stress. Gentics 2012, 192:289-318
# Kevin A. Morano et al. The response to heat shock and oxidative stress in Saccharomyces cerevisiae. Genetics 2011.

# copy2clipboard(unique(arrange(subset(master.table.020, select=c(Gene, rand_mRNA)),Gene))$rand_mRNA)

# fread(input = "results/SMoPT/master.table.reduced.txt", colClasses="character") # bugs... mistakes a character column for a logical one, beh.
master.table.reduced <- data.table( read.table("results/SMoPT/master.table.reduced.txt", header=1), key=c("ORF","experiment","time") )


master.table.reduced$log10.est.mRNA_abundance <- log10(master.table.reduced$est.mRNA_abundance)
master.melt <- melt( subset(master.table.reduced, experiment !="normal" & !is.na(faster.translation) ), id.vars = c("ORF","experiment","time","faster.translation","TRPs","up_reg", "stress", "quant.mRNAab","overal_charge", "delta_initiation"))
variables.of.interest = c("ratio_total.time","global.protein_synthesis.rate", "av.initiation_time","av.initiation_time.normal","ratio_initiation","ratio_elongation","delta_elongation","length_CDS","ratio_elongation", "log10.est.mRNA_abundance")

hlines <- data.frame( variable = variables.of.interest, yint = c(1,1,NA,NA,1,1,0,NA,1,1) )
nrow(master.melt)
head(master.melt)


# generated in supply_demand_hypothesis2.R
expressed.codons_vs_adaptiveness  <-  read.table( file = "results/stress tAI/stress.codon.adaptiveness_vs_expression.txt",sep="\t", header=1)



# Look up ratio of initiation (stress/normal) vs. ratio of elongation (stress/normal) across conditions
ggplot( data = subset(gene.master.table, experiment !="normal"), aes( x= log2(ratio_initiation), y= log2(ratio_elongation) ) ) +
  geom_hline(yintercept=0, lty=2, color="gray60") +
  geom_vline(xintercept=0, lty=2, color="gray60") +
  geom_point(aes(color=faster.translation)) +
  scale_color_manual(values=c("gray40","red")) +
  facet_wrap( experiment ~ time ) +
  theme_bw() + theme(legend.position = "none")


# Look up delta initiation (stress/normal) vs. delta elongation (stress/normal) across conditions
ggplot( data = subset(gene.master.table, experiment !="normal"), aes( x= delta_initiation, y= delta_elongation ) ) +
  geom_hline(yintercept=0, lty=2, color="gray60") +
  geom_vline(xintercept=0, lty=2, color="gray60") +
  geom_point(aes(color=faster.translation)) +
  scale_color_manual(values=c("gray40","red")) +
  facet_wrap( experiment ~ time ) +
  theme_bw() + theme(legend.position = "none")






ggplot( subset( master.melt, experiment == "diauxic" & time == 20 & variable %in% variables.of.interest ), aes( x= as.factor(faster.translation), y = as.numeric(value), fill = faster.translation )  ) + 
  geom_hline( data =  hlines, aes(yintercept= yint), lty=1, color="gray" ) +
  geom_boxplot( notch = T ) + 
#  scale_y_sqrt() +
  scale_fill_manual(values = c("gray70","red") ) +
  facet_wrap( ~ variable, scale="free" ) + labs(x="",y="", title="response to diauxic shock\n") +
  theme_bw() + theme(legend.position = "none")


ggplot( subset( master.melt, experiment == "osm" & time == 20 & variable %in% variables.of.interest ), aes( x= as.factor(faster.translation), y = as.numeric(value), fill = faster.translation )  ) + 
  geom_hline( data =  hlines, aes(yintercept= yint), lty=1, color="gray" ) +
  geom_boxplot( notch = T ) + 
  #  scale_y_sqrt() +
  scale_fill_manual(values = c("gray70","red") ) +
  facet_wrap( ~ variable, scale="free" ) + labs(x="",y="", title="response to osmotic shock\n") +
  theme_bw() + theme(legend.position = "none")


ggplot( subset( master.melt, experiment == "ox" & time == 20 & variable %in% variables.of.interest ), aes( x= as.factor(faster.translation), y = as.numeric(value), fill = faster.translation )  ) + 
  geom_hline( data =  hlines, aes(yintercept= yint), lty=1, color="gray" ) +
  geom_boxplot( notch = T ) + 
  #  scale_y_sqrt() +
  scale_fill_manual(values = c("gray70","red") ) +
  facet_wrap( ~ variable, scale="free" ) + labs(x="",y="", title="response to oxidative stress\n") +
  theme_bw() + theme(legend.position = "none")


ggplot( subset( master.melt, experiment == "temp" & time == 20 & variable %in% variables.of.interest ), aes( x= as.factor(faster.translation), y = as.numeric(value), fill = faster.translation )  ) + 
  geom_hline( data =  hlines, aes(yintercept= yint), lty=1, color="gray" ) +
  geom_boxplot( notch = T ) + 
  #  scale_y_sqrt() +
  scale_fill_manual(values = c("gray70","red") ) +
  facet_wrap( ~ variable, scale="free" ) + labs(x="",y="", title="response to temperature stress\n") +
  theme_bw() + theme(legend.position = "none")


table( master.table.reduced[ , list( faster.elongation = ratio_elongation < 0.8 )  ] )


# ----- EXPRESSED CODONS ----- #

g <- ggplot( data = codon.master.table, aes(x= rank.demand.mRNAab, y=w)) + 
  stat_smooth(method="lm",color="gray40") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\nranked demand\n(codon occurence x mRNA abundance summed over genes)", y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

ggsave(plot = g, filename = "results/stress tAI/codon.adaptiveness_across.conditions_1.pdf")

g <- ggplot( data = codon.master.table, aes(x= rank.demand.mRNAab, y=w)) + 
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_wrap( ~ codon ) + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\nranked demand\n(codon occurence x mRNA abundance summed over genes)", y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom", axis.text.x=element_blank() )

ggsave(plot = g, filename = "results/stress tAI/codon.adaptiveness_across.conditions_2.pdf")




g <- ggplot( data = codon.master.table , aes(x= demand.mRNAab, y=w)) + #expressed.codons_vs_adaptiveness
  stat_smooth(method="lm",color="gray40") +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_wrap(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\ndemand\n(codon occurence x mRNA abundance summed over genes)", y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

ggsave(plot = g, filename = "results/stress tAI/codon.adaptiveness_across.conditions_3.pdf")


ggplot( data = codon.master.table, aes(x= rank.demand.mRNAab, y= rank.demand.mRNAab.up )) + 
  geom_point(alpha=0.8, aes(size = w, color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\ndemand\n(in all mRNAs)", y="demand\n(in top-20% up-regulated mRNAs)\n") +
  theme_bw() + theme(legend.position="bottom")




ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= rank.demand.mRNAab, y= w )) + 
  stat_smooth(method="lm",color="gray40") +
  geom_point(alpha=0.8, aes(color=CAI.O), size=5) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\nranked demand\n(in all mRNAs)", y="w\ncodon relative adaptiveness (to tRNA pool)\n") +
  theme_bw() + theme(legend.position="bottom")

# in top-20% up-regulated genes, blue codons "invade" the right region (increase in relative demand)
ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= rank.demand.mRNAab.up, y= w )) + 
  stat_smooth(method="lm",color="gray40") +
  geom_point(alpha=0.8, aes(color=CAI.O), size=5) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\nranked demand\n(in top-20%-up regulated mRNAs)", y="w\ncodon relative adaptiveness (to tRNA pool)\n") +
  theme_bw() + theme(legend.position="bottom")


# gain in ranks for codon demand (mRNA abundance)
ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= rank.demand.mRNAab.up - rank.demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x="\ndifferential ranked demand\n(codon occurence x mRNA abundance summed over genes)\nup-regulated mRNAs - all mRNAs",
       y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

# same, based on actual demand (not rank)
ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= demand.mRNAab.up - demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x=expression(Delta[demand]), y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

# relative change (compared to all mRNAs)
ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (demand.mRNAab.up - demand.mRNAab)/demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

# now showing the relative change (compared to all mRNAs) in relative demand (compared to other codons)
g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

ggsave(plot=g, filename = "results/stress tAI/codon.adaptiveness_across.conditions_upregulated.genes_vs_all.genes.pdf",dpi=400)



g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y=delta_w)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="w\nchange in relative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")

ggsave(plot=g, filename = "results/stress tAI/codon.adaptiveness.change_across.conditions_upregulated.genes_vs_all.genes.pdf",dpi=400)


g <- ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= (relative.demand.mRNAab.up - relative.demand.mRNAab)/relative.demand.mRNAab, y= av.elongation_time )) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x='relative difference in demand\n(top-20% vs all mRNAs)', y="elongation time (s)") +
  theme_bw() + theme(legend.position="bottom")



# same. based on actual demand from elongation events
ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= demand.events.up - demand.events, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment ~ time, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x=expression(Delta[events]), y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")





ggplot( data = subset(codon.master.table, experiment !="normal"), aes(x= demand.events.up - demand.events, y=w)) + 
  geom_vline(xintercept=0, lty=2) +
  geom_point(alpha=0.8, size=4.25, aes(color=CAI.O)) + 
  geom_text(aes(label=codon), size=2.25, color="gray20") +
  facet_grid(experiment + time ~ aa, scale="free") + 
  scale_color_manual(values=c("skyblue3","orange")) +
  labs(x=expression(Delta[events]), y="w\nrelative adaptiveness\n",color="optimality\n(normal cond.)") +
  theme_bw() + theme(legend.position="bottom")


