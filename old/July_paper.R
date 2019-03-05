# ESR and mRNA abundance (Steady-state) ---------------------------------------------------------------------------------------------------------------------------------
ESR.table$tAI.profile <- factor(ESR.table$tAI.profile, levels=c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
ESR.table$cluster <- factor(ESR.table$cluster, levels=c("ESR.down","other","ESR.up"))
ESR.table$is.ESR.up <- factor( with( ESR.table, ifelse( cluster == "ESR.up", "ESR.up", "other.genes") ), levels = c("other.genes","ESR.up"))



DATA <- subset(ESR.table, time == "20 min" & stAI > 0)
gESR_mRNA <- ggplot(DATA, aes( x = as.factor(experiment),  y = mRNA_abundance )) +
  geom_boxplot(aes(fill=cluster), notch=T) + scale_y_log10() + scale_fill_manual(values=c("red","gray","skyblue")) + labs(x="stress",y="mRNA abundance")


ggsave(plot=gESR_mRNA, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_mRNA.pdf",dpi=250, useDingbats=F, width=6, height=4.5)

tESR_mRNAab <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "cluster"), value = "mRNA_abundance"  )
}
))

write.table(tESR_mRNAab, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_mRNAab.wc.txt", quote=F, sep="\t", row.names=F)

# ESR and log2 mRNA abundance fold change ---------------------------------------------------------------------------------------------------------------------------------

DATA <- subset(ESR.table, time == "20 min" & stAI > 0)
gESR_log2.mRNA <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2.mRNA_abundance )) +
  geom_boxplot(aes(fill=cluster), notch=T) + scale_fill_manual(values=c("red","gray","skyblue")) + labs(x="stress",y="log2 mRNA abundance")

ggsave(plot=gESR_log2.mRNA, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_log2.mRNA.pdf",dpi=250, useDingbats=F, width=6, height=4.5)


tESR_log2.mRNAab <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "cluster"), value = "log2.mRNA_abundance"  )
}
))

write.table(tESR_log2.mRNAab, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_mRNAab.wc.txt", quote=F, sep="\t", row.names=F)


# ESR and tAI ---------------------------------------------------------------------------------------------------------------------------------

DATA <- subset(ESR.table, time == "20 min" & stAI > 0)
gESR_01 <- ggplot(DATA, aes( x = as.factor(experiment),  y = stAI )) +
  geom_boxplot(aes(fill=cluster), notch=T) + scale_fill_manual(values=c("red","gray","skyblue")) + labs(x="stress",y="s-tAI")

ggsave(plot=gESR_01, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_tAI.pdf",dpi=250, useDingbats=F, width=6, height=4.5)

tESR_01 <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "cluster"), value = "stAI"  )
}
))

write.table(tESR_01, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_tAI.wc.txt", quote=F, sep="\t", row.names=F)






# ESR and change in tAI -----------------------------------------------------------------------------------------------------------------------

DATA <- subset(ESR.table, time == "20 min" & stAI > 0 & gain.stAI.scaled< 80)
gESR_gain.tAI <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(gain.FS.stAI) )) +
  geom_boxplot(aes(fill=cluster), notch=T) + scale_fill_manual(values=c("red","gray","skyblue")) + labs(x="stress",y="gain s-tAI") + ylim(-2,2)

ggsave(plot=gESR_gain.tAI, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_gain.tAI.pdf",dpi=250, useDingbats=F, width=6, height=4.5)


ggplot(DATA, aes( x = as.factor(experiment),  y = (stAI) )) + facet_wrap( ~ experiment, scales = "free" ) +
  geom_boxplot(aes(fill=is.ESR.up), notch=T) + scale_fill_manual(values=c("gray","skyblue")) + labs(x="stress",y="s-tAI") 

ggplot(DATA, aes( x = as.factor(experiment),  y = (stAI) )) + facet_wrap( ~ experiment, scales = "free" ) +
  geom_boxplot(aes(fill=tAI.profile), notch=T) + scale_fill_manual(values=c("red","gray","orange")) + labs(x="stress",y="s-tAI") 


gProfile_tAI <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) + facet_wrap( ~ experiment, scales = "free" ) +
  geom_boxplot(aes(fill=tAI.profile), notch=T, outlier.size = 0) + 
  scale_fill_manual(values=c("red","gray","orange")) + 
  labs(x="stress",y="protein production rate fold change") + 
  theme(legend.position="none") + ylim(-4,2)



gProfile_tAI <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) + facet_wrap( ~ experiment, scales = "free" ) +
  geom_boxplot(aes(fill=tAI.profile), notch=T, outlier.size = 0) + 
  scale_fill_manual(values=c("red","gray","orange")) + 
  labs(x="stress",y="protein production rate fold change") + 
  theme(legend.position="none") + ylim(-4,2)


ggsave(plot=gProfile_tAI, filename = "~/Documents/MRC16/Paper/draft/july add-ons/Profile_tAI.pdf",dpi=250, useDingbats=F, width=3, height=4.5)


tProfile_tAI <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "tAI.profile"), value = "global.protein_synthesis.rate"  )
}
))




tESR_tAI <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "cluster"), value = "gain.stAI"  )
}
))

write.table(tESR_tAI, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_gain.tAI.wc.txt", quote=F, sep="\t", row.names=F)








# ESR and protein production rates ------------------------------------------------------------------------------------------------------------

DATA <- subset(ESR.table, time == "20 min" & global.protein_synthesis.rate > 0)
gESR_Pi <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) +
  geom_boxplot(aes(fill=cluster), notch=T) + scale_fill_manual(values=c("red","gray","skyblue")) + labs(x="stress",y=expression(Pi))

ggsave(plot=gESR_Pi, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_Pi.pdf",dpi=250, useDingbats=F, width=6, height=4.5)



    
  tESR_Pi <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "cluster"), value = "global.protein_synthesis.rate"  )
}
))

write.table(tESR_Pi, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR_Pi.wc.txt", quote=F, sep="\t", row.names=F)


gESR.up_Pi <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) +
  geom_boxplot(aes(fill=is.ESR.up), notch=T, outlier.size = 0 ) + ylim(-3.5,3.5) + scale_fill_manual(values=c("gray","skyblue")) + labs(x="stress",y=expression(Pi) ) + theme(legend.position="none")

ggsave(plot=gESR.up_Pi, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR.up_Pi.pdf",dpi=250, useDingbats=F, width=3, height=2.5)

tESR.up_Pi <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "is.ESR.up"), value = "global.protein_synthesis.rate"  )
}
))

write.table(tESR.up_Pi, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR.up_Pi.wc.txt", quote=F, sep="\t", row.names=F)


gESR.up_stAI <- ggplot(DATA, aes( x = as.factor(experiment),  y = FS.stAI )) +
  geom_boxplot(aes(fill=is.ESR.up), notch=T, outlier.size = 0 ) + scale_fill_manual(values=c("gray","skyblue")) + labs(x="stress",y="s-tAI (scaled)" ) + theme(legend.position="none") + ylim(0,0.8)

ggsave(plot=gESR.up_stAI, filename = "~/Documents/MRC16/Paper/draft/july add-ons/ESR.up_stAI.pdf",dpi=250, useDingbats=F, width=3, height=2.5)


tESR.up_stAI <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "is.ESR.up"), value = "FS.stAI"  )
}
))

write.table(tESR.up_stAI, file = "~/Documents/MRC16/Paper/draft/july add-ons/ESR.up_stAI.wc.txt", quote=F, sep="\t", row.names=F)


# protein production fold change at 20min ----------------------------------------------------------------------------------------------------
DATA <- subset(ESR.table, time == "20 min" & global.protein_synthesis.rate > 0)

gtAI_Pi20 <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) +
  geom_boxplot(aes(fill=tAI.profile), notch = T) + scale_fill_manual(values=c("red","gray","orange")) + labs(x="stress",y=expression(Pi))

    ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
      wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "tAI.profile"), value = "global.protein_synthesis.rate"  )
    }
    ))

ggsave(plot=gtAI_Pi20, filename = "~/Documents/MRC16/Paper/draft/july add-ons/tAI_Pi20.pdf",dpi=250, useDingbats=F, width=10, height=4.5)
    
ttAI_Pi20 <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "tAI.profile"), value = "global.protein_synthesis.rate"  )
}
))

write.table(ttAI_Pi20, file = "~/Documents/MRC16/Paper/draft/july add-ons/tAI_Pi20.wc.txt", quote=F, sep="\t", row.names=F)










    
# protein production fold change at 120min ---------------------------------------------------------------------------------------------------
DATA <- subset(ESR.table, time == "120 min" & global.protein_synthesis.rate > 0)

gtAI_Pi120 <- ggplot(DATA, aes( x = as.factor(experiment),  y = log2(global.protein_synthesis.rate) )) +
  geom_boxplot(aes(fill=tAI.profile), notch = T) + scale_fill_manual(values=c("red","gray","orange")) + labs(x="stress",y=expression(Pi))

ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "tAI.profile"), value = "global.protein_synthesis.rate"  )
}
))

ggsave(plot=gtAI_Pi120, filename = "~/Documents/MRC16/Paper/draft/july add-ons/tAI_Pi120.pdf",dpi=250, useDingbats=F, width=6, height=4.5)

ttAI_Pi120 <- ldply(lapply( setNames(c("diauxic","ox","osm","temp"),c("diauxic","ox","osm","temp")), function(e){
  wilcox.test.batch(  subset(DATA, experiment == e), grouping.vars = c( "tAI.profile"), value = "global.protein_synthesis.rate"  )
}
))

write.table(ttAI_Pi120, file = "~/Documents/MRC16/Paper/draft/july add-ons/tAI_Pi120.wc.txt", quote=F, sep="\t", row.names=F)







# Gene table ----------------------------------------------------------------------------------------------------------------------------------

GSH.table <- subset( SMoPT.data, name %in% c("GSH1","GSH2"), select = 
                       c(name, experiment, time, faster.translation, mRNA_abundance, mRNA_abundance.normal, log2.mRNA_abundance, stAI, FS.stAI, FS.tAI, gain.stAI, tAI.profile  ) )


arrange(GSH.table, name, time, log2.mRNA_abundance)




# Gene enrichment sets:
  
subset(GO.results$decreased.tAI,high.confidence=="yes")
subset(GO.results$increased.tAI,high.confidence=="yes")


ggplot(data=subset(SMoPT.data, time !=0), aes(x=FS.tAI, y=stAIFC)) + geom_point(aes(color=tAI.profile1)) + facet_wrap( ~ time + experiment, nrow = 2)
ggplot(data=subset(ESR.table, time !=0), aes(x=FS.tAI, y=stAIFC)) + geom_point(aes(color=cluster)) + facet_wrap( ~ time + experiment, nrow = 2)

# tAI analysis: genes with decreasing ranks

GO_decreased.tAI <-  ddply( subset(SMoPT.data, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'tAI.profile1 == "decreased\nadaptation" ', keys = O ) )
})
GO_decreased.tAI$x <- factor(GO_decreased.tAI$x, label = "decreased tAI")  

results_GO_decreased.tAI <- subset(GO_decreased.tAI, high.confidence=="yes", select = c(experiment, time, x, y, diagnostic, high.confidence,  p.value, adjusted.p, sensitivity, accuracy))
subset(results_GO_decreased.tAI, diagnostic == "enriched")

write.table(  results_GO_decreased.tAI, path.expand("~/Documents/MRC16/Paper/draft/july add-ons/GO_decreased-tAI.txt"), row.names = F, sep="\t" )

# tAI analysis: genes with decreasing ranks

GO_increased.tAI <-  ddply( subset(SMoPT.data, experiment!="normal"), .(experiment, time), function(x) {
  subset( automate.fisherexact( data = x, expression = 'tAI.profile1 == "increased\nadaptation" ', keys = O ) )
})
GO_increased.tAI$x <- factor(GO_increased.tAI$x, label = "increased tAI")  


results_GO_increased.tAI <- subset(GO_increased.tAI, high.confidence=="yes", select = c(experiment, time, x, y, diagnostic, high.confidence,  p.value, adjusted.p, sensitivity, accuracy))
subset(results_GO_increased.tAI, diagnostic == "enriched")


write.table(  results_GO_increased.tAI, path.expand("~/Documents/MRC16/Paper/draft/july add-ons/GO_increased-tAI.txt"), row.names = F, sep = "\t" )





# Gene table for Madan ------------------------------------------------------------------------------------------------------------------------------------------------------------------
variables.to.use <- c("name","ORF","experiment","time","mRNA_abundance","mRNA_abundance.normal", "log2.mRNA_abundance","tAI","stAIFC","tAI.profile1",
                      "faster.translation","ratio_events","av.initiation_time","av.initiation_time.normal","av.elongation_time","av.elongation_time.normal","ratio_elongation","global.protein_synthesis.rate") 
to.round <- c("tAI","stAIFC","ratio_events","av.initiation_time","av.initiation_time.normal","av.elongation_time","av.elongation_time.normal","ratio_elongation","global.protein_synthesis.rate") 

input_madan_tables <- subset(SMoPT.data, select= eval(variables.to.use))
input_madan_tables[,to.round] <- sapply( input_madan_tables[,to.round], round, 2)
input_madan_tables$tAI.profile1 <- factor(input_madan_tables$tAI.profile1, levels=c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"), labels=c("decreased adaptation","similar adaptation","increased adaptation"))

colnames(input_madan_tables) <- c("gene name","ORF","stress","time","mRNA abundance","mRNA abundance (normal)", "log2 mRNA abundance change","n-tAI","s-tAI","tAI-group", "is TFS", "change in translation rate","av. waiting time initiation", "av. waiting time initiation (normal)","av. time to elongate", "av. time to elongate (normal)", "ratio elongation time", "Pi (protein production fold change)")

madan_tables <- dlply(input_madan_tables, .(stress, time) )

library(WriteXLS)
WriteXLS("madan_tables", path.expand("~/Documents/MRC16/Paper/draft/july add-ons/gene_table.xlsx"), col.names = T, row.names = F, BoldHeaderRow = T)


##########

l_up <- dlply(subset(results_GO_increased.tAI, time == 20 & diagnostic == "enriched"), .(experiment), function(x) as.character(x$y))
l_down <- dlply(subset(results_GO_decreased.tAI, time == 20 & diagnostic == "enriched"), .(experiment), function(x) as.character(x$y))

Reduce(intersect, l_up)
Reduce(intersect, l_down)

intersect( l_up$diauxic, l_up$ox)
intersect( l_up$temp, intersect( l_up$diauxic, l_up$ox))
intersect( l_up$diauxic, l_up$osm)


intersect( l_down$diauxic, l_down$ox)
intersect( l_down$temp, intersect( l_down$diauxic, l_down$ox))
intersect( l_down$diauxic, l_down$osm)
