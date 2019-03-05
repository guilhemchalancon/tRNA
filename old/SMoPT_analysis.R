setwd("~/Documents/MRC16/")
source("scripts/commands.R")

PATH  <-  "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"

# l'affaire est ketchup
anticodon.usage.table <- read.table("results/master tables/anticodon.usage.table.txt",header=1)


require(scales)
SMoPT.data <- read.table(file = "results/master tables/SMoPT.table.txt",sep="\t", header=1)
 SMoPT.data$label <- factor(SMoPT.data$label, 
                           levels= c(paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"),"normal\n0 min")
                           )
head(SMoPT.data)


COLFAST <- "#ff610d"

# remove entries for non-stress conditions, and ignore potential problematic entries that have null mRNA abundance 
# (March2015: I corrected these by attributing them their continuous mRNA abundance estimations, but I think Marc gave them a value of 1 in the simulation)
study.translation.kinetics <- subset(SMoPT.data, experiment !="normal" & mRNA_abundance >0)
# % variation in elongation speed
study.translation.kinetics$variation_speed <- (study.translation.kinetics$elongation.speed - study.translation.kinetics$elongation.speed.normal) / study.translation.kinetics$elongation.speed.normal
# % variation in tAI 
study.translation.kinetics$variation_tAI <- (study.translation.kinetics$stAI - study.translation.kinetics$tAI)/study.translation.kinetics$tAI
# % variaiton in initiation frequency
study.translation.kinetics$variation_initiation_frequency <- (1/study.translation.kinetics$av.initiation_time - 1/study.translation.kinetics$av.initiation_time.normal)*study.translation.kinetics$av.initiation_time.normal
# % variation in translation rate (alternative to log2 of ratio)
study.translation.kinetics$variation_translation_rate <- (study.translation.kinetics$expected.translation_rate - study.translation.kinetics$expected.translation_rate.normal)/study.translation.kinetics$expected.translation_rate.normal
# scaled variation values in each condition
study.translation.kinetics <- ddply(study.translation.kinetics, .(experiment,time), mutate, 
                                    scaled.change.tAI = scale(variation_tAI), 
                                    scaled.variation_speed = scale(variation_speed),
                                    scaled.variation_initiation_frequency = scale(variation_initiation_frequency),
                                    scaled.variation_rate = scale(variation_translation_rate)
                                    )

write.table(study.translation.kinetics, file="results/master tables/study.translation.kinetics.txt",sep="\t",row.names=F)
study.translation.kinetics <- read.table("results/master tables/study.translation.kinetics.txt",sep="\t",header=1)
study.translation.kinetics$label <- factor(study.translation.kinetics$label, 
                                     levels= paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n"))
#anticodon.usage.table <- read.table("results/master tables/anticodon.usage.table.txt",sep="\t",header=1)


go.slim.ORFs <- read.table("data/annotations/go_slim_ORFs.txt",header=1,sep="\t")
colnames(go.slim.ORFs)[2] <- "name"
  O <- data.frame(subset(unique(go.slim.ORFs[,c("GO.description","GO.type")]), !GO.description %in% c("not_yet_annotated","biological_process","molecular_function","cellular_component") ), stringsAsFactors=F)
  O$GO.description <- as.character(O$GO.description)
  O$GO.type <- as.character(O$GO.type)
  rownames(O) <- NULL
  O <- subset(O, GO.type %in% "P" )# | GO.description %in% names(complexes[complexes>15]) )


# Rapid analysis of GO slim enrichment in ttg-rich mRNAs
ttg.rich.genes <- subset(anticodon.usage.table, anticodon == "CAA")
ttg.rich.genes$ttg.rich <- with(data=ttg.rich.genes, OR.CDS > 2 & OR.CDS_CI.low > 1 )
ttg.GO <- automate.fisherexact( data = ttg.rich.genes, expression = 'ttg.rich == TRUE', keys = O) 

View(subset(ttg.GO, high.confidence == "yes", select = c(y, sensitivity, specificity, n.TT, adjusted.p)))



# Rapid analysis of GO slim enrichment in ttg-rich mRNAs
gcc.rich.genes <- subset(anticodon.usage.table, anticodon == "CCG")
gcc.rich.genes$gcc.rich <- with(data=gcc.rich.genes, OR.CDS > 1 & OR.CDS_CI.low > 1 )
gcc.GO <- automate.fisherexact( data = gcc.rich.genes, expression = 'gcc.rich == TRUE', keys = O) 

View(subset(gcc.GO, high.confidence == "yes", select = c(y, sensitivity, specificity, n.TT, adjusted.p)))

# CAA (Leu) enrichment
CAA.table <- arrange( merge(study.translation.kinetics , subset(anticodon.usage.table, anticodon == "CAA"), by = "name", all.x=T), plyr::desc(freq.in_CDS) )
CAA.table$ttg.rich <-  with(data=CAA.table,  OR.CDS > 2 & OR.CDS_CI.low > 1)

colnames(CAA.table)

# Alternative definitions of translation rates have little impact on the results above. They highly correlate:
with( data = study.translation.kinetics, expr = cor(1/ratio_total.time,ratio_expected.translation.rate)) # good. 0.86
with( data = study.translation.kinetics, expr = cor(1/ratio_total.time,ratio_apparent.translation.rate)) # good. 0.74
with( data = study.translation.kinetics, expr = cor(expected.translation_rate,apparent.translation.rate)) # good. 0.87


# LINEAR REGRESSIONS NEEDED IN SUPPORT OF PLOTS 

regressions.speed.tAI <- ddply( study.translation.kinetics, .(label), function(x){ batch.regressions(x, x = "stAI",y = "elongation.speed") })
regressions.var.speed.var.tAI <- ddply( study.translation.kinetics, .(label), function(x){ batch.regressions(x, x = "variation_tAI",y = "variation_speed") })
regressions.speed <- ddply( study.translation.kinetics, .(label), function(x){ batch.regressions(x, x = "elongation.speed.normal",y = "elongation.speed") })
regressions.tAI.initiation_freq <- ddply( study.translation.kinetics, .(label), function(x){ batch.regressions(x, x = "stAI",y = "variation_initiation_frequency") })
regressions.CAA.elongation <- ddply( CAA.table, .(label), function(x){ batch.regressions(x, x = "freq.in_CDS",y = "elongation.speed") })


# Are genes that have faster translation have a significantly increased intiation frequency during stress?
wc.initiation.faster.translation <- ddply( study.translation.kinetics, .(experiment, time), function(x){ 
  wilcox.test.batch(x =x , grouping.vars = c("faster.translation"), 
                    value = "variation_initiation_frequency" ) } )





ggplot(data= study.translation.kinetics, 
       aes( x = delta_initiation, y = delta_elongation, 
            #color = tAI.profile
            color = faster.translation
            #color = delta_initiation < - delta_elongation
       )  ) + 
  geom_point()  +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  geom_abline(intercept=0, slope=-1) +
  labs(x=expression(paste(Delta["i"]~"\n(change in waiting time for initiation)")), 
       y=expression(paste(Delta["e"]~"\n(change in average elongation time)"))
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_color_manual(values = c("gray70",COLFAST)) +
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")





###### FIGURE ELONGATION KINETICS ########




g1.elongation <- ggplot(data= study.translation.kinetics, 
       aes( x = elongation.speed.normal, y = elongation.speed, 
            #color = tAI.profile
            #color = faster.translation
            #color = delta_initiation < - delta_elongation
       )  ) + 
  geom_point(pch=21, fill="white", size=1.25)  +
  #geom_density2d() +
  geom_text( data = regressions.speed, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =4 ) +
  geom_abline(intercept=0, slope=1, lty=3) +
  labs(x="elongation speed\nin normal conditions (codon/s)", 
       y="elongation speed\nduring stress (codon/s)"
  ) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  #scale_color_manual(values = c("gray20","red")) +
  stat_smooth(se = T,method = "gam", color="skyblue3") +
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")



g2.elongation <- ggplot(data= study.translation.kinetics, 
            aes( x = stAI, y = elongation.speed, 
                 #color = tAI.profile
                 #color = faster.translation
                 #color = delta_initiation < - delta_elongation
            )  ) + 
  geom_point(pch=21, fill="white", size=1.25)  +
  #geom_density2d() +
  #geom_hline(yintercept=0, lty=3) +
 # geom_vline(xintercept=0, lty=3) +
  geom_text( data = regressions.speed.tAI, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =4 ) +
  #geom_abline(intercept=0, slope=1, lty=3) +
  labs(x="stress-adjusted tAI", 
       y="elongation speed\nduring stress (codon/s)"
  ) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  #scale_color_manual(values = c("gray20","red")) +
  stat_smooth(se = T,method = "gam", color="red") +
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

require(scales)
g3.elongation <- ggplot(data= study.translation.kinetics, 
       aes( x = variation_tAI, y = variation_speed, 
            #color = tAI.profile
            #color = faster.translation
            #color = delta_initiation < - delta_elongation
       )  ) + 
  geom_point(pch=21, fill="white", size=1.25)  +
  #geom_density2d() +
  geom_hline(yintercept=0, lty=3) +
  geom_vline(xintercept=0, lty=3) +
  geom_text( data = regressions.var.speed.var.tAI, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =4 ) +
  #geom_abline(intercept=0, slope=1, lty=3) +
    labs(x="variation in tAI\n(%)", 
         y="variation in average elongation speed\n(%)"
    ) +
  scale_x_continuous( labels = percent ) +
  scale_y_continuous( labels = percent ) +
  stat_smooth(se = T,method = "gam", color="orange") +
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

##### PLOT
#pdf(paste0(PATH,"fig_results_elongation_kinetics.pdf"), width = 10,  height = 15,  useDingbats=F)
png(paste0(PATH,"fig_results_elongation_kinetics.png"), width = 10,  height = 15, res=150, units="in")
require(gridExtra)
grid.draw(rbind_gtable_max(ggplotGrob(g2.elongation),
                           ggplotGrob(g1.elongation),
                           ggplotGrob(g3.elongation)
                           ))
dev.off()

#####




# 
# 
# ggplot(data= study.translation.kinetics, 
#        aes( x = apparent.translation.rate, y = 1/av.initiation_time*1000, 
#             #color = tAI.profile
#             #color = faster.translation
#             #color = delta_initiation < - delta_elongation
#        )  ) + 
#   geom_point(pch=21, fill="white", size=1.25)  +
#   geom_abline(intercept=0, slope=1, lty=1, color="blue") +
#   stat_smooth(se = T,method = "gam", color="orange") +
#   facet_wrap( ~ label, nrow=2) + 
#   theme_gul +
#   labs(x= "apparent translation rate\n(n.events/mRNA)",y="predicted initiation frequency") +
#   theme(legend.position="none")
# 


# how much elongation time explains translation rate? 
# ggplot(data= study.translation.kinetics, 
#        aes( x = apparent.translation.rate, y = 1/av.elongation_time*1000, 
#             #color = tAI.profile
#             #color = faster.translation
#             #color = delta_initiation < - delta_elongation
#        )  ) + 
#   geom_point(pch=21, fill="white", size=1.25)  +
#   geom_abline(intercept=0, slope=1, lty=1, color="blue") +
#   stat_smooth(se = T,method = "gam", color="orange") +
#   facet_wrap( ~ label, nrow=2) + 
#   theme_gul +
#   theme(legend.position="none")

# 
# ggplot(data= study.translation.kinetics, 
#        aes( x = log2(ratio_apparent.translation.rate), y = variation_speed, 
#             #color = tAI.profile
#             #color = faster.translation
#             #color = delta_initiation < - delta_elongation
#        )  ) + 
#   geom_point(pch=21, fill="white", size=1.25)  +
#   geom_abline(intercept=0, slope=1, lty=1, color="blue") +
#   stat_smooth(se = T,method = "gam", color="orange") +
#   facet_wrap( ~ label, nrow=2) + 
#   theme_gul +
#   theme(legend.position="none")

# from the point where translation rate increase, so do initiation frequency

require(scales)
g.data <- subset(study.translation.kinetics, experiment =="diauxic" & time == 20)
gA <- ggplot( data= subset(g.data, ratio_initiation > 0 ), 
  #ggplot( data=subset(study.translation.kinetics, experiment !="normal" & ratio_initiation > 0 ), 
             aes( x = variation_speed, y = log2(ratio_total.time) )  ) + 
  stat_smooth( data = subset(g.data, experiment !="normal" & (tAI >=0.5 | name %in% c("DPM1","PUB1"))),
               method = "loess", color ="gray30",size=1.5, alpha=1) +
  stat_smooth( data = subset(g.data, experiment !="normal" & (tAI >=0.5 | name %in% c("DPM1","PUB1"))),
               method = "loess", color ="white",size=1, alpha=1) +
  geom_point(  data = subset(g.data, experiment !="normal" & tAI < 0.5),
    aes(fill = faster.translation, shape = tRNA_adaptation_index >=0.5),
    pch=21
  )  +
  geom_point(  data = subset(g.data, experiment !="normal" & tAI >=0.5), 
    fill="white", shape=24) +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  labs( x="variation in elongation speed (%)", y="Changes in\ntranslation time (log2)" ) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  scale_shape_manual(values = c(16,24) ) +
  scale_x_continuous(labels=percent) +
  facet_wrap( ~ label, nrow = 2) + 
  theme_gul + 
  theme(legend.position="none", plot.margin=unit(c(2,2,8,2),"mm"))
#log2(1/ratio_initiation)

gB <- ggplot( data=subset(g.data, ratio_initiation > 0 ),  
  #ggplot( data=subset(study.translation.kinetics, experiment !="normal" & ratio_initiation > 0 ), 
        aes( x = log2(1/ratio_initiation), y = log2(ratio_apparent.translation.rate) )  ) + 
  geom_point(pch=21, aes(fill = faster.translation))  +
  stat_smooth( method = "lm", color ="gray30",size=1.75, alpha=1) +
  stat_smooth( method = "lm", color ="white",size=1, alpha=1) +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  labs( x="Changes in initiation frequency (log2)", y="Changes in\ntranslation time (log2)" ) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap( ~ label, nrow = 2) + 
  theme_gul + 
  theme(legend.position="none", plot.margin=unit(c(2,2,8,2),"mm"))
# 
# ggsave( plot = g, filename = paste0(PATH,"fig_results_elongation.riboprot.png"), 
#         dpi = 250, width = 9.75, height = 6)

pdf( file = paste0(PATH,"fig_results_elong_init_dependency.pdf"), width=7.5, height=3.5, useDingbats = F)
require(gridExtra)
grid.draw(cbind_gtable_max(ggplotGrob(gA),ggplotGrob(gB)))
dev.off()

 





g <- ggplot( data=subset(study.translation.kinetics, experiment !="normal" & ratio_initiation > 0 ), 
             aes( x = (tAI >=0.5) , y= mRNA_abundance)
           ) + 
  geom_boxplot(notch=T) +
  facet_wrap( ~ label, nrow = 2) + 
  theme_gul + 
  scale_y_log10() +
  theme(legend.position="none")


ggplot( data=subset(study.translation.kinetics, experiment !="normal" & ratio_initiation > 0 ), 
        aes( x = log2(gain.rank) , y= mRNA_abundance)
) + 
  geom_point(aes(color=tAI.profile)) +
  geom_vline(xintercept=0, lty=3, color="gray") +
  facet_wrap( ~ label, nrow = 2) + 
  theme_gul + 
  scale_color_manual(values=c("red","orange", "gray30")) +
  scale_y_log10(breaks=c(1,10,100,1000)) +
  theme(legend.position="none")


######## CHARACTERISING FTS GENES ########

# DELTA FIGURE
ggplot(data= study.translation.kinetics, 
       aes( x = delta_initiation, y = delta_elongation, 
            #color = tAI.profile
            color = faster.translation
            #color = delta_initiation < - delta_elongation
       )  ) + 
  geom_point()  +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  geom_abline(intercept=0, slope=-1) +
  labs(x=expression(paste(Delta["i"]~"\n(change in waiting time for initiation)")), 
       y=expression(paste(Delta["e"]~"\n(change in average elongation time)"))
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_color_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")


# 
# ggplot(data= study.translation.kinetics, 
#        aes( x = delta_initiation, y = delta_elongation, 
#             #color = tAI.profile
#             color = tAI.profile
#             #color = delta_initiation < - delta_elongation
#        )  ) + 
#   geom_point()  +
#   geom_vline(xintercept=0, lty=2, colour = "gray50") + 
#   geom_hline(yintercept=0, lty=2, colour = "gray50") + 
#   geom_abline(intercept=0, slope=-1) +
#   labs(x=expression(paste(Delta["i"]~"\n(change in waiting time for initiation)")), 
#        y=expression(paste(Delta["e"]~"\n(change in average elongation time)"))
#   ) +
#   #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
#   scale_x_continuous( labels = prettyNum ) +
#   scale_y_continuous( labels = prettyNum ) +
#   scale_color_manual(values = c("red","orange","gray30")) +
#   facet_wrap(  ~ label, nrow=2) + 
#   theme_gul #+
#  theme(legend.position="none")



# genes that are translated faster tend to have an increased initiation frequency
ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 30), 
       aes( x = faster.translation, y = (variation_initiation_frequency)
       )  ) + 
  geom_boxplot(notch=T, aes(fill=faster.translation)) +
  scale_fill_manual(values=c("gray70",COLFAST)) +
  scale_y_continuous(labels=percent)+
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")



d.FTS <- ddply(study.translation.kinetics, .(label), summarise, 
      increased.freq = percent(round(sum(variation_initiation_frequency > 0)/length(variation_initiation_frequency),2)),
      faster.elong   = round(sum(variation_speed > 0)/length(variation_speed),2),
      faster.translation = percent(round(sum(faster.translation==T)/length(faster.translation),2))
)
d.FTS$faster.elong <- percent(d.FTS$faster.elong) # really can't get why it works this way
d.FTS$faster.elong <- with(d.FTS, ifelse(faster.elong == "0%", "<0.1%", faster.elong) )

require(scales)
# NEW DELTA FIGURE
g1.FTS <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10), 
       aes( x = variation_initiation_frequency, y = variation_speed, 
       )  ) + 
  geom_point(aes(fill = faster.translation), pch=21, alpha=0.5)  +
  geom_vline(xintercept=0, lty=2, colour = "gray50") + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  geom_text( data = subset(study.translation.kinetics, name %in% c("PUB1","DPM1","INH1")), 
             aes(label=name), size=2.25, hjust=-0.25) +
  geom_text( data = subset(study.translation.kinetics, name %in% c("DSS1","RDS1")  & variation_initiation_frequency < 10), 
             aes(label=name), size=2.25, vjust=1.75) +
  geom_text( data = d.FTS, aes( label = (faster.elong) ), x=-Inf, y=+Inf,    hjust = -0.1, vjust = 1.5,  size =3 ) +
  geom_text( data = d.FTS, aes( label = (increased.freq)), x=Inf, y=-Inf,    hjust = 1.5,  vjust =-0.75, size =3 ) +
  geom_text( data = d.FTS, aes( label = (faster.translation)), x=Inf, y=Inf, hjust = 1.5,  vjust = 1.5,  size =3, color=COLFAST ) +
  # geom_abline(intercept=0, slope=-1) +
  labs(x=expression(i^{"%"}*" (percent change in initiation frequency)"), 
       y=expression(e^{"%"}*" (percent change in elongation speed)")
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_continuous( labels = percent ) +
  scale_y_continuous( labels = percent ) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none", plot.margin = unit(c(2,2,8,2),"mm"))



ggsave(plot=g1.FTS, filename =  paste0(PATH,"part3/g1.FTS.pdf"), useDingbats=F, width=9, height=5.45)
ggsave(plot=g1.FTS, filename =  paste0(PATH,"fig_results_delta.png"), dpi=190, width=9.5, height=5.45)

PUB1.OR <- arrange( subset(anticodon.usage.table, name %in% c("PUB1","DPM1") & OR.CDS_valid.interval == 1,
                select=c(anticodon, name, preference_in.CDS, preference_ramp, n.in_CDS, OR.CDS, OR.ramp) ),
         name, plyr::desc(OR.CDS), plyr::desc(OR.ramp))
PUB1.OR[,-c(1:2)] <- sapply(PUB1.OR[,-c(1:2)], round, 2)

anticodon.data <- subset(anticodon.master.table, select=c(anticodon, aa, experiment, time, recognised, log2.foldchange, CAI.O))


PUB1.trans <- subset(study.translation.kinetics, name %in% c("PUB1","DPM1"), select=c(name,experiment, time, n.events, 
                                                                        elongation.speed, elongation.speed.normal, FS.gain.tAI) )
PUB1.trans[,-c(1:4)] <- sapply(PUB1.trans[,-c(1:4)], round, 2)

PUB1.data  <- merge( merge(PUB1.OR, anticodon.data, by=c("anticodon")), PUB1.trans, by= c("experiment","time","name"))

View(arrange( PUB1.data, name, plyr::desc(OR.CDS), plyr::desc(OR.ramp)))

arrange( subset( PUB1.data, log2.foldchange > 0), name, experiment, time, plyr::desc(log2.foldchange) )

get.sequence.info <- function(
  name = "PUB1"
  ){
  require(StingRay)
  require(seqinr)
  data(SEQ)
  seq <- tolower(get.codons(SEQ$ORF_CDS[[fORF(name)]]))
  data <- arrange( merge ( data.frame( pos=1:length(seq), codon = seq ), subset(rosetta.codons, select=c(codon, anticodon, aa, CAI.O)), by = "codon"), pos)
  return(data)
}


seq.data <- get.sequence.info("DPM1")
seq.map <- arrange(merge(seq.data, codon.master.table, by=c("codon","aa","anticodon","CAI.O")), pos)
seq.map$label <- paste(seq.map$experiment, seq.map$time, sep="\n")
seq.map$label <- factor(seq.map$label, levels=  paste( rep(c("diauxic","ox","osm","temp"), each = 2 ), rep(c(20,120), time = 4 ), sep="\n"),
                               labels = paste( rep(c("diauxic","oxidative","osmotic","temperature"), each = 2 ), rep(c("20min","120min"), time = 4 ), sep="\n") )

ddply(seq.map, .(label), summarise, 
      faster.codon.elongation =(round(sum(ratio_elongation < 1)/length(ratio_elongation),2)),
      total.elongation = sum(delta_elongation)
      )


ggplot(data = subset(seq.map, time > 0), aes(x=pos, y = factor(1), fill=log2.foldchange, label=codon)) + geom_tile() +
  coord_fixed(ratio = 10 ) +
  theme_minimal() + theme(legend.position="bottom", axis.text.y = element_blank()) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high="#2166AC", midpoint = 0) + 
  facet_wrap( ~ label, ncol =2 )

ggplot(data = subset(seq.map, time>0), aes(x=pos, y = factor(1), fill= log2(ratio_elongation), label=codon)) + geom_tile() +
  coord_fixed(ratio = 10 ) +
  theme_minimal() + theme(legend.position="bottom", axis.text.y = element_blank()) + 
  labs(x="",y="") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high="#2166AC", midpoint = 0) +
  facet_wrap( ~ label, ncol =2 )

subset(study.translation.kinetics, name == "PUB1")

ggplot(data = subset(seq.map, time == 0), aes(x=pos, y = factor(1), fill=CAI.O, label=codon)) + geom_tile() +
  coord_fixed(ratio = 10 ) +
  theme_minimal() + theme(legend.position="bottom", axis.text.y = element_blank()) + 
  scale_fill_manual(values=c("skyblue","orange"))

# 
# ggplot(data = subset(seq.map, time>0), aes(x=log2(ratio_elongation), y = log2.foldchange, label=codon, fill = CAI.O)) + 
#   geom_point(pch=21) +
#   geom_text(size=2.5, vjust=1.4)+
#   theme(legend.position="bottom", axis.text.y = element_blank()) + 
#   scale_fill_manual(values=c("skyblue","orange")) +
#   facet_wrap( ~ label, nrow =2 )



#### FIGURE INITIATION KINETICS ####

# extreme outliers (var initiation freq >20) have been removed (to explicitly mention in figure caption)
data <- subset(study.translation.kinetics, variation_initiation_frequency < 10)
g1.initiation <- ggplot(data= data, 
       aes( x = log2(ratio_apparent.translation.rate), y = variation_initiation_frequency, 
            #color = tAI.profile
            #color = faster.translation
            #color = delta_initiation < - delta_elongation
       )  ) + 
  geom_point(pch=21, fill="white", size=1.25, alpha=0.5)  +
  geom_text(data=subset(data, variation_initiation_frequency>10), aes(label=name), size=2.25)  +
  geom_abline(intercept=0, slope=1, lty=1, color="blue") +
  geom_vline(xintercept=0,lty=3)+
  scale_y_continuous(labels=percent) +
  stat_smooth(se = T,method = "loess", color="orange", n = 50) +
  facet_wrap( ~ label, nrow=2) + 
  theme_gul +
  labs(x="change in translation rate\n(log2)", y="variation in initiation frequency (%)") +
  theme(legend.position="none")

ggsave(plot=g1.initiation, filename =  paste0(PATH,"part3/g1.initiation.pdf"), useDingbats=F, width=9, height=5.45)


data <- subset(study.translation.kinetics, variation_initiation_frequency < 10)
g2.medians <- ddply(data, .(label,faster.translation), summarise, 
                 av.initiation_time = median(av.initiation_time), 
                 av.elongation_time = median(av.elongation_time),
                 elongation.speed= median(elongation.speed))

g2.initiation <- ggplot(data= data, 
       aes( x = av.initiation_time/60, y = av.elongation_time/60, 
       )  ) + 
  geom_point(aes(fill = faster.translation), pch=21, alpha=0.35)  +
  geom_rug(data=subset(data,faster.translation==F), sides="br", color="gray70", alpha=0.5) +
  geom_rug(data=subset(data,faster.translation==T), sides="tl", color=COLFAST, alpha=0.5) +
  geom_rug(data= ddply(subset(data,faster.translation==T), .(label), summarise, 
                       av.initiation_time = median(av.initiation_time), 
                       av.elongation_time= median(av.elongation_time) ), sides="tl", color="blue") +
  geom_rug(data= ddply(subset(data,faster.translation==F), .(label), summarise, 
                       av.initiation_time = median(av.initiation_time), 
                       av.elongation_time= median(av.elongation_time) ), sides="br", color="blue") +
  geom_vline(data = g2.medians, aes(xintercept=av.initiation_time/60), lty=3, colour = "gray30") + 
  geom_hline(data = g2.medians, aes(yintercept=av.elongation_time/60), lty=3, colour = "gray30") + 
  geom_point( data= g2.medians, shape=23,size=3.5, fill="white", aes(color=faster.translation) ) +
  labs(x="waiting time between initiation events (min)", 
       y="duration of elongation (min)"
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_sqrt( labels = prettyNum, breaks=c(1,5,10,20,40,60) ) +
  scale_y_sqrt( labels = prettyNum, breaks=c(1,5,10,15,20,25) ) +
  scale_color_manual(values = c("gray70",COLFAST)) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")



ggsave(plot=g2.initiation, filename =  paste0(PATH,"fig_results_TFS_map.pdf"), useDingbats=F, width=11, height=6.54)
ggsave(plot=g2.initiation, filename =  paste0(PATH,"fig_results_TFS_map.png"),  width=11, height=6.54, dpi=150)




#################################################################################

id.vars = c("ORF", "name", "valid","experiment", "time", "faster.translation","tAI.profile","label") 
features <- c("mRNA_abundance", "log2.mRNA_abundance", "length_CDS", "AUGCAI", "IniProb", "global.protein_synthesis.rate",
              "stAI.20codons","stAIFC.20codons", "n.codons",
              "gain.stAI.20codons","gain.stAIFC.20codons",
              "gain.stAI.20codons.scaled", "gain.stAIFC.20codons.scaled")

gene.master.melt <- melt( study.translation.kinetics[, c(id.vars, features)], 
                         id.vars = id.vars  )
gene.master.melt$value <- as.numeric(gene.master.melt$value)
gene.master.melt$faster.translation <- ifelse(gene.master.melt$faster.translation == T, "TFS", "other")

mannwhitneys <- list(
  diauxic = ddply( subset(gene.master.melt, experiment == "diauxic" & time == 20 & variable %in% features, select = c(label, faster.translation,variable, value)), .(variable),
                            function(x) { 
                              wilcox.test.batch(x = x, grouping.vars = c("label","faster.translation"), value = "value")  }
  ),
  oxidative = ddply( subset(gene.master.melt, experiment == "ox" & time == 20  & variable %in% features, select = c(label, faster.translation,variable, value)), .(variable),
                    function(x) { 
                      wilcox.test.batch(x = x, grouping.vars = c("label","faster.translation"), value = "value")  }
  ),
  osmotic = ddply( subset(gene.master.melt, experiment == "osm" & time == 20  & variable %in% features, select = c(label, faster.translation,variable, value)), .(variable),
                    function(x) { 
                      wilcox.test.batch(x = x, grouping.vars = c("label","faster.translation"), value = "value")  }
  ),
  temperature = ddply( subset(gene.master.melt, experiment == "temp" & time == 20  & variable %in% features, select = c(label, faster.translation,variable, value)), .(variable),
                    function(x) { 
                      wilcox.test.batch(x = x, grouping.vars = c("label","faster.translation"), value = "value")  }
  )
)
  
View(mannwhitneys$diauxic)



p1 <- ggplot( data = subset(study.translation.kinetics, experiment == "diauxic" & time == 20) ) + 
  geom_boxplot( aes(x= faster.translation, y = mRNA_abundance, fill=faster.translation), notch=T ) +
  scale_y_log10( breaks=c(1,10,100,1000)) +
  geom_text( data = subset(mannwhitneys$diauxic, variable == "mRNA_abundance"), aes(label= paste(paste("r=",(r.wendt),sep=""),paste("p<",p.adjusted), sep="\n") ), y=2.8, x=2, hjust=0.5, vjust=0.4, size =3 ) +
  scale_fill_manual(values=c("gray70", COLFAST)) +
  labs( x = "", y="mRNA abundance") + theme(legend.position = "none")

p2 <- ggplot( data = subset(study.translation.kinetics, experiment == "diauxic" & time == 20) ) + 
  geom_boxplot( aes(x= faster.translation, y = n.codons, fill=faster.translation), notch=T ) +
  scale_y_log10(breaks=c(50,100,250,500,1000,2000, 4000) ) +
  geom_text( data = subset(mannwhitneys$diauxic, variable == "n.codons"), aes(label= paste(paste("r=",(r.wendt),sep=""),paste("p<",p.adjusted), sep="\n") ), y=2, x=2, hjust=0.5, vjust=0.4, size =3 ) +
  scale_fill_manual(values=c("gray70", COLFAST)) +
labs( x = "", y="CDS length") + theme(legend.position = "none")

p3 <- ggplot( data = subset(study.translation.kinetics, experiment == "diauxic" & time == 20) ) + 
  geom_boxplot( aes(x= faster.translation, y = AUGCAI, fill=faster.translation), notch=T ) +
  #scale_y_log10() +
  geom_text( data = subset(mannwhitneys$diauxic, variable == "AUGCAI"), aes(label= paste(paste("r=",(r.wendt),sep=""),paste("p<",p.adjusted), sep="\n") ), y=0.9, x=2, hjust=0.5, vjust=0.4, size =3 ) +
  scale_fill_manual(values=c("gray70", COLFAST)) +
  labs( x = "", y="Start codon adaptation\n(AUG-CAI)") + theme(legend.position = "none")

p4 <- ggplot( data = subset(study.translation.kinetics, experiment == "diauxic" & time == 20 ) ) + 
  geom_boxplot( aes(x= faster.translation, y = (gain.stAI.20codons), fill=faster.translation), notch=T ) +
  #scale_y_log10() +
  geom_text( data = subset(mannwhitneys$diauxic, variable == "gain.stAI.20codons"), aes(label= paste(paste("r=",(r.wendt),sep=""),paste("p<",p.adjusted), sep="\n") ), y=1, x=2, hjust=0.5, vjust=0.4, size =3 ) +
  scale_fill_manual(values=c("gray70", COLFAST)) +
  labs( x = "", y="Change in\ncodon ramp adaptation (stAI)") + theme(legend.position = "none")



# 
# ggplot( data = subset(study.translation.kinetics, experiment == "diauxic" & time == 20 ) ) + 
#   geom_point( aes(x =  variation_initiation_frequency, y = (gain.stAI.20codons), fill=faster.translation), pch=21 ) +
#   #scale_y_log10() +
#   scale_fill_manual(values=c("gray70", COLFAST)) 

require(gridExtra)
png(filename = paste0(PATH,"fig_results_TFS_multiplot.png"), res =250,units = "in", 
    width = 9.6, height = 11.6)
  grid.arrange( ggplotGrob(g2.initiation),
    arrangeGrob(empty, rbind_gtable_max( cbind_gtable_max(ggplotGrob(p1),ggplotGrob(p2)), 
                                         cbind_gtable_max(ggplotGrob(p3),ggplotGrob(p4))
                                        ),
                empty,
                widths=c(1/4,1/2,1/4), ncol=3, main = "TFS genes in early response to diauxic shift" 
                ), 
    heights=c(5/10,5/10))
dev.off()


pdf(file = paste0(PATH,"fig_results_TFS_multiplot.pdf"), useDingbats=F,
    width = 9.6, height = 11.6)
grid.arrange( ggplotGrob(g2.initiation),
              arrangeGrob(empty, rbind_gtable_max( cbind_gtable_max(ggplotGrob(p1),ggplotGrob(p2)), 
                                                   cbind_gtable_max(ggplotGrob(p3),ggplotGrob(p4))
              ),
              empty,
              widths=c(1/4,1/2,1/4), ncol=3, main = "TFS genes in early response to diauxic shift" 
              ), 
              heights=c(5/10,5/10))
dev.off()


#################################################################################
            g3.initiation <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10), 
                   aes( x = faster.translation, y = elongation.speed, 
                   )  ) + 
              geom_hline(yintercept= mean(subset(SMoPT.data, time ==0)$elongation.speed), lty=3, colour = "#758316") + 
              geom_boxplot(notch=T,alpha=0.65, outlier.size = 0, aes(fill=faster.translation)) +
              labs(x="faster translation", 
                   y="elongation speed\n(codons/s)"
              ) +
              #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
              #scale_x_continuous( labels = prettyNum ) +
              scale_y_continuous( labels = prettyNum ) +
              scale_fill_manual(values = c("gray70",COLFAST)) +
              facet_wrap(  ~ label, nrow=2) + 
              theme_gul +
              theme(legend.position="none")
            
            g4.initiation <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10 & 1/av.initiation_time < 0.25), 
                         aes( x = faster.translation, y = 1/av.initiation_time)  ) + 
              stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                           geom = "crossbar", width = 0.65) +
              geom_hline(yintercept= 1/mean(subset(SMoPT.data, time ==0)$av.initiation_time), lty=3, colour = "#758316") + 
              geom_boxplot(notch=T,alpha=0.65, outlier.size = 0, aes(fill=faster.translation)) +
              #   geom_text( data = subset(study.translation.kinetics, name %in% c("PUB1","DPM1","INH1")), 
              #              aes(label=name), size=2.5, hjust=-0.25) +
              #   geom_text( data = subset(study.translation.kinetics, name %in% c("DSS1","RDS1")  & variation_initiation_frequency < 10), 
              #              aes(label=name), size=2.5, vjust=1.75) +
              # geom_abline(intercept=0, slope=-1) +
              labs(x="faster translation", 
                   y=expression("initiation frequency "*s^{-1})
              ) +
              #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
              #scale_x_continuous( labels = prettyNum ) +
              scale_y_continuous( labels = prettyNum ) +
              scale_fill_manual(values = c("gray70",COLFAST)) +
              facet_wrap(  ~ label, nrow=2) + 
              theme_gul +
              theme(legend.position="none")


data <- subset(study.translation.kinetics, variation_initiation_frequency < 10)
g5.medians <- ddply(data, .(label,faster.translation), summarise, 
                 av.initiation_time.normal = median(av.initiation_time.normal), 
                 variation_initiation_frequency= median(variation_initiation_frequency))

g5.initiation <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10), 
       aes( x = 1/av.initiation_time.normal, y = variation_initiation_frequency)  ) + 
  geom_vline(data = g5.medians, aes(xintercept=1/av.initiation_time.normal, color=faster.translation), lty=6) + 
  geom_hline(data = g5.medians, aes(yintercept=variation_initiation_frequency, color=faster.translation), lty=6) + 
  geom_point(aes(fill = faster.translation), pch=21,alpha=.5)  +
  geom_hline(yintercept=0, lty=2, colour = "gray30") + 
#   geom_text( data = subset(study.translation.kinetics, name %in% c("PUB1","DPM1","INH1")), 
#              aes(label=name), size=2.5, hjust=-0.25) +
  geom_text( data = subset(study.translation.kinetics, name %in% c("DSS1","RDS1")  & variation_initiation_frequency < 10), 
             aes(label=name), size=2.5, vjust=1.75) +
  geom_point( data= g5.medians, shape=23,size=3.5, fill="white", aes(color=faster.translation) ) +
  # geom_abline(intercept=0, slope=-1) +
  labs(x="initiation frequency (normal conditions)", 
       y="%i\nchanges in initiation frequency (%)"))
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_sqrt( labels = prettyNum ) +
  scale_y_continuous( labels = percent ) +
  scale_color_manual(values = c("gray50","orange")) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

ggsave(plot=g5.initiation, filename =  paste0(PATH,"part3/g5.initiation.pdf"), useDingbats=F, width=11, height=6.54)


g6.medians <- ddply(data, .(label,faster.translation), summarise, 
                    elongation.speed.normal = median(elongation.speed.normal), 
                    variation_speed= median(variation_speed))

g6.speed <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10), 
             aes( x = elongation.speed.normal, y = variation_speed)  ) + 
  geom_vline(data = g6.medians, aes(xintercept=elongation.speed.normal, color=faster.translation), lty=6) + 
  geom_hline(data = g6.medians, aes(yintercept=variation_speed, color=faster.translation), lty=6) + 
  geom_point(aes(fill = faster.translation), pch=21,alpha=.5)  +
  geom_hline(yintercept=0, lty=2, colour = "gray30") + 
  #   geom_text( data = subset(study.translation.kinetics, name %in% c("PUB1","DPM1","INH1")), 
  #              aes(label=name), size=2.5, hjust=-0.25) +
#   geom_text( data = subset(study.translation.kinetics, name %in% c("DSS1","RDS1")  & variation_initiation_frequency < 10), 
#              aes(label=name), size=2.5, vjust=1.75) +
   geom_text( data = subset(study.translation.kinetics, name %in% c("PUB1","DPM1","INH1")), 
             aes(label=name), size=2.25, hjust=-0.25) +
  geom_point( data= g6.medians, shape=23,size=3.5, fill="white", aes(color=faster.translation) ) +
  # geom_abline(intercept=0, slope=-1) +
  labs(x="elongation speed in normal conditions (codons/s)", 
       y="%e\nchanges in elongation speed"
  ) +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_sqrt( labels = prettyNum ) +
  scale_y_continuous( labels = percent ) +
  scale_color_manual(values = c("gray50","orange")) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

ggsave(plot=g6.speed, filename =  paste0(PATH,"part3/g6.speed.pdf"), useDingbats=F, width=11, height=6.54)


require(gridExtra)
pdf( file = paste0(PATH,"fig_results_TFS_kinetics.pdf"), useDingbats=F, width=11, height=12)
grid.draw(rbind_gtable_max(ggplotGrob(g5.initiation),ggplotGrob(g6.speed)))
dev.off()

png( file = paste0(PATH,"fig_results_TFS_kinetics.png"), width=11, height=12, units = "in", res=150 )
grid.draw(rbind_gtable_max(ggplotGrob(g5.initiation),ggplotGrob(g6.speed)))
dev.off()



g7.medians <- ddply(data, .(label,faster.translation), summarise, 
                    apparent.translation.rate = median(apparent.translation.rate), 
                    global.protein_synthesis.rate= median(global.protein_synthesis.rate))
   
g7.production <- ggplot(data= subset(study.translation.kinetics, variation_initiation_frequency < 10), 
       aes( x = apparent.translation.rate, y = global.protein_synthesis.rate )  ) + 
  geom_hline(yintercept = 1, color="gray40", lty=3) + 
  geom_vline(data = g7.medians, aes(xintercept=apparent.translation.rate,     color=faster.translation), lty=6) + 
  geom_hline(data = g7.medians, aes(yintercept=global.protein_synthesis.rate, color=faster.translation), lty=6) + 
  geom_point(aes(fill = faster.translation), pch=21,alpha=.5)  +
  geom_text( data = subset(study.translation.kinetics, name %in% c("DSS1","RDS1")  & variation_initiation_frequency < 10), 
             aes(label=name), size=2.5, vjust=1.75) +
  geom_point( data= g7.medians, shape=23,size=3.5, fill="white", aes(color=faster.translation) ) +
  # geom_abline(intercept=0, slope=-1) +
  labs(y=expression(Pi~"global protein production rate"), 
       x="translation rate") +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_sqrt( labels = prettyNum ) +
  scale_y_sqrt( labels = prettyNum, breaks=c(0,1,5,10,20,30) ) +
  scale_color_manual(values = c("gray50","orange")) +
  scale_fill_manual(values = c("gray70",COLFAST)) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

ggsave(plot=g7.production, filename =  paste0(PATH,"part3/g7.production.pdf"), useDingbats=F, width=11, height=6.54)

require(gridExtra)


# bit too early to bring "fast.translation", since in the dissertation i look at it from the perspective of whether CAA influences elongation
  # g1.CAA <- ggplot(data= subset(CAA.table, elongation.speed < 10), 
  #        aes( x = paste( faster.translation,  OR.CDS > 2), y = elongation.speed
  #        )  ) + 
  #   geom_hline(yintercept= mean(subset(SMoPT.data, time ==0)$elongation.speed), lty=3, colour = "#758316") + 
  #   geom_boxplot(notch=T,alpha=0.65, outlier.size = 0, aes(fill= OR.CDS > 2)) +
  #   labs(x="faster translation", 
  #        y="elongation speed\n(codons/s)"
  #   ) +
  #   #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  #   #scale_x_continuous( labels = prettyNum ) +
  #   scale_y_continuous( labels = prettyNum ) +
  #   scale_fill_manual(values = c("gray70",COLFAST)) +
  #   facet_wrap(  ~ label, nrow=2) + 
  #   theme_gul +
  #   theme(legend.position="none")

wc.g2.CAA <-  ddply( subset(CAA.table, elongation.speed < 15), .(label), function(x){ 
    wilcox.test.batch(x =x , grouping.vars = c("ttg.rich"), 
                      value = "elongation.speed" ) } )
require(scales)
g1.CAA <- ggplot(data= subset(CAA.table, elongation.speed < 15), 
                 aes( x = freq.in_CDS, y = elongation.speed )  ) + 
  geom_point(aes(fill = ttg.rich), pch=21,alpha=.5)  +
  stat_smooth(se = T,method = "loess", color="orange") +
  geom_text( data = regressions.CAA.elongation, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =4 ) +
  labs(y="elongation speed\n(codons/s)", 
       x="CAA demand in CDS\n(% of ttg codons)") +
  #geom_boxplot( notch=T, outlier.size=0, alpha=0.7) + 
  scale_x_continuous( labels = percent ) +
  scale_y_sqrt( labels = prettyNum ) +
  #scale_color_manual(values = c("gray50","orange")) +
  scale_fill_manual(values = c("gray70","#0899FF")) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")

#ggsave(plot=g1.CAA, filename =  paste0(PATH,"part3/g1.CAA.pdf"), useDingbats=F, width=11, height=6.54)
ggsave(plot=g1.CAA, filename =  paste0(PATH,"fig_results_CAA_scatterplot.png"), width=11, height=6.54, dpi=150)



#####

empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )


g2a.CAA <- ggplot(data=subset(CAA.table, elongation.speed <15 & time == 20 & experiment == "diauxic")) + 
  geom_histogram(aes(x=freq.in_CDS, fill=ttg.rich), binwidth= range(CAA.table$freq.in_CDS)[2]/50 ) + 
  labs(x="frequency of ttg codon in CDS") +
  scale_x_continuous(label=percent) +
  theme(legend.position="none") + scale_fill_manual(values=c("gray40","#0899FF"))

g2b.CAA <- ggplot(data= subset(CAA.table, elongation.speed < 15), 
       aes( x = ttg.rich, y = elongation.speed
       )  ) + 
  geom_hline(yintercept= mean(subset(SMoPT.data, time ==0)$elongation.speed), lty=3, colour = "#758316") + 
  geom_boxplot(notch=T,alpha=0.65, outlier.size = 0, aes(fill= ttg.rich)) +
  labs(x="Genes enriched in ttg codon\n(decoded by Leu-CAA tRNA)", 
       y="elongation speed\n(codons/s)"
  ) +
  #geom_text( data = wc.g1.CAA, x=1, y=Inf, aes(label=(r)), parse=T, vjust=1.5 ) + 
  scale_y_continuous( labels = prettyNum ) +
  scale_fill_manual(values = c("gray40","#0899FF")) +
  facet_wrap(  ~ label, nrow=2) + 
  theme_gul +
  theme(legend.position="none")


  png(paste0(PATH,"fig_results_CAA_boxplot.png"), width = 10,  height = 7, res=150, units="in")
  require(gridExtra)
  grid.arrange( arrangeGrob(g2a.CAA, empty, ncol=1) , g2b.CAA, ncol=2, heights=c(1/3,2/3), widths=c(1/3,2/3))
  dev.off()

  pdf(paste0(PATH,"fig_results_CAA_boxplot.pdf"), width = 10,  height = 7, usingDingbats=F, units="in")
  require(gridExtra)
  grid.arrange( arrangeGrob(g2a.CAA, empty, ncol=1) , g2b.CAA, ncol=2, heights=c(1/3,2/3), widths=c(1/3,2/3))
  dev.off()



data.t0 <- merge( subset(SMoPT.data, experiment =="normal" & mRNA_abundance >0), subset(anticodon.usage.table, anticodon == "CAA"), by = "name")
data.t0$ttg.rich <-  with(data=data.t0,  OR.CDS > 2 & OR.CDS_CI.low > 1)

ggplot(data= data.t0, aes(x=ttg.rich,y=elongation.speed) )+ geom_boxplot()



#------ APRIL 2015 ---- Alternative approach

# data.frame( experiment = c("diauxic", "diauxic", "ox", "ox", "osm", "osm", "temp", "temp"), time=rep(c(20,120), times = 4),
#             alpha = c() )
CAA.balance <- subset(anticodon.master.table, anticodon == "CAA", select = c(experiment, time, total.tRNAab, demand.mRNAab, foldchange) )
CAA.balance$SD_ratio <- with(CAA.balance, total.tRNAab / demand.mRNAab)


kinetics.ttg.rich <- ddply(CAA.table, .(experiment, time, ttg.rich), summarise, 
      median.speed = median(elongation.speed, na.rm=T), 
      median.sc.var.speed = median(scaled.variation_speed, na.rm=T),
      median.variation.speed = median( variation_speed , na.rm = T),
      median.log2.change.speed = median( log2(1/ratio_elongation)  ),
      median.ratio_apparent.translation.rate = median( ratio_apparent.translation.rate, na.rm=T)
        )

data <- merge(CAA.balance, kinetics.ttg.rich, by = c("experiment","time"))

ggplot(data = data, aes(x=SD_ratio, y = median.ratio_apparent.translation.rate, color = ttg.rich )) + geom_line()



# 



