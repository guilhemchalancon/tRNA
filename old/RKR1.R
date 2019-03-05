# RKR1 
# RING domain E3 ubiquitin ligase; involved in ubiquitin-mediated degradation of non-stop proteins; c
# omponent of ribosome-bound RQC (ribosome quality control) complex required for degradation of polypeptides arising from stalled translation; 
# degrades products of mRNAs lacking a termination codon regardless of a poly(A) tail; functional connections to chromatin modification; homolog of mouse Listerin, mutations in which reported to cause neurodegeneration


subset(gene.master.table, name == "RKR1")
subset(stAI.data, name == "RKR1")


study.translation.kinetics$quant.initiation <- quantile.it( study.translation.kinetics$variation_initiation_frequency )
study.translation.kinetics$quant.speed <- quantile.it( study.translation.kinetics$variation_speed )
study.translation.kinetics$quant.prod <- quantile.it( study.translation.kinetics$ratio_events )
study.translation.kinetics$quant.tAI <- quantile.it( study.translation.kinetics$variation_tAI )
#ggplot(study.translation.kinetics, aes(x=quant.initiation, y= log2(ratio_initiation) )) + geom_boxplot(notch=T)
#ggplot(study.translation.kinetics, aes(x=quant.elongation, y= variation_speed )) + geom_boxplot(notch=T)



# ------------------------------ Interplay between changes in mRNA abundance and changes in tAI --------------------------------------------




# Genes that go up or down in transcription, and gain or lose a lot in tAI
Profiles <- subset(study.translation.kinetics, quant.mRNAab %in% c("0-20","80-100") & tAI.profile %in% c("decreased\nadaptation","increased\nadaptation") )
Profiles$tAI.profile4 <- "other" 
Profiles$tAI.profile4[ Profiles$tAI.profile == "increased\nadaptation" & Profiles$FS.stAI > mean(Profiles$FS.stAI)  ] <- "now_well_adapted" #+ sd(Profiles$FS.tAI)
Profiles$tAI.profile4[ Profiles$tAI.profile == "decreased\nadaptation" & Profiles$FS.stAI < mean(Profiles$FS.stAI)  ] <- "now_poorly_adapted" #- sd(Profiles$FS.tAI)
Profiles <- subset(Profiles, tAI.profile4 !="other")

Profiles$profile <- factor( paste(Profiles$quant.mRNAab, Profiles$tAI.profile4), 
                            levels = sort(unique(paste(Profiles$quant.mRNAab, Profiles$tAI.profile4))), #"0-20 now_poorly_adapted"   "0-20 now_well_adapted"     "80-100 now_poorly_adapted" "80-100 now_well_adapted" 
                            labels = c("incons.down","consis.up","consis.down","incons.up")
                           )

Profiles$profile <- factor(Profiles$profile, levels= c ("consis.down","incons.up","incons.down","consis.up") )

dlply ( Profiles, .(experiment, time, profile ), function(x) {as.character(x$name)})
dlply ( Profiles, .(experiment, time, profile ), function(x) {x$global.protein_synthesis.rate})




# Wilcoxon tests Translation Efficiency ~ stAI groups
wc.table <- wilcox.test.batch(Profiles, grouping.vars = "profile", value = "global.protein_synthesis.rate")
wc.table


ggplot( data = Profiles, aes(x=profile, y=global.protein_synthesis.rate)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=T, outlier.size = 0) + 
  labs(y=expression(Pi), x="", title = "Protein synthesis rates are more influenced\nby transcription than by elongation")+
  scale_y_continuous(limits=c(0,3))

g1 <- ggplot( data = Profiles, aes(x= FS.stAI, y = log2.mRNA_abundance  )) + 
  geom_point( aes(fill=profile), pch=21, color="gray20" ) +
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  geom_hline( yintercept = 0, lty=3) + geom_vline( data = data.frame( mean = mean(Profiles$FS.stAI)) ,aes(xintercept = mean), lty=3) +
  labs( x = "stAI (scaled)", y = "fold change mRNA (log2)", title="Aggregated\ndata" ) + theme(legend.position="none")

g2 <- ggplot( data = Profiles, aes(x=profile, y=ratio_events)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=T, outlier.size = 0, aes(fill=profile)) + 
  labs(y=expression(Pi), x="", title = "Protein production rates are more influenced\nby transcription than by elongation")+
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  scale_y_continuous(limits=c(0,3)) + theme(legend.position="none")

require(gridExtra)
pdf(file = paste0(PATH,"part3/SMoPT_tAI.transcription.protein.synthesis.pdf"), width=9.24, height=5.03, useDingbats=F)
grid.arrange(g1, g2, ncol=2)
dev.off()
  #scale_y_sqrt() 


# ------------------------------ Interplay between changes in mRNA abundance and changes in initiation frequency --------------------------------------------


Profiles.init <- subset(study.translation.kinetics, quant.mRNAab %in% c("0-20","80-100") & quant.initiation %in% c("0-20","80-100"))


Profiles.init$profile <- factor( paste(Profiles.init$quant.mRNAab, Profiles.init$quant.initiation), 
                            levels = sort(unique(paste(Profiles.init$quant.mRNAab, Profiles.init$quant.initiation))),
                            labels = c("consis.up","incons.down","incons.up","consis.down")
)

Profiles.init$profile <- factor(Profiles.init$profile, levels= c ("consis.down","incons.up","incons.down","consis.up") )




Profiles.elong <- subset(study.translation.kinetics, quant.mRNAab %in% c("0-20","80-100") & quant.speed %in% c("0-20","80-100"))


Profiles.elong$profile <- factor( paste(Profiles.elong$quant.mRNAab, Profiles.elong$quant.speed), 
                                 levels = sort(unique(paste(Profiles.elong$quant.mRNAab, Profiles.elong$quant.speed))),
                                 labels = c("consis.up","incons.down","incons.up","consis.down")
)
Profiles.elong$profile <- factor(Profiles.elong$profile, levels= c ("consis.down","incons.up","incons.down","consis.up") )



g1 <- ggplot( data = Profiles, aes(x= FS.stAI, y = log2.mRNA_abundance  )) + 
  geom_point( aes(fill=profile), pch=21, color="gray20" ) +
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  geom_hline( yintercept = 0, lty=3) + geom_vline( data = data.frame( mean = mean(Profiles$FS.stAI)) ,aes(xintercept = mean), lty=3) +
  labs( x = "stAI (scaled)", y = "fold change mRNA (log2)", title="Aggregated\ndata" ) + theme(legend.position="none")

g2 <- ggplot( data = Profiles, aes(x=profile, y=ratio_events)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=T, outlier.size = 0, aes(fill=profile)) + 
  labs(y=expression(Pi), x="", title = "Protein production rates are more influenced\nby transcription than by elongation")+
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  scale_y_continuous(limits=c(0,3)) + theme(legend.position="none")

g3 <- ggplot( data = Profiles.init, aes(x= log2(1/ratio_initiation), y = log2.mRNA_abundance  )) + 
  geom_point( aes(fill=profile), pch=21, color="gray20" ) +
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  geom_hline( yintercept = 0, lty=3) + geom_vline( xintercept = 0, lty=3) +
  labs( x = "fold change initiation frequency (log2)", y = "fold change mRNA (log2)", title="Aggregated\ndata" ) + theme(legend.position="none")


g4 <- ggplot( data = Profiles.init, aes(x=profile, y=ratio_events)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=T, outlier.size = 0, aes(fill=profile)) + 
  labs(y=expression(Pi), x="", title = "Increase in initiation frequency strongly\ninfluences protein production rates")+
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  scale_y_continuous(limits=c(0,3)) + theme(legend.position="none")

ggplot( data = Profiles.elong, aes(x=profile, y=ratio_events)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=T, outlier.size = 0, aes(fill=profile)) + 
  labs(y=expression(Pi), x="", title = "Protein production rates are more influenced\nby transcription than by elongation")+
  scale_fill_manual(values=c("white","gray30","gray30","white")) +
  scale_y_continuous(limits=c(0,3)) + theme(legend.position="none")



require(gridExtra)
pdf(file = paste0(PATH,"part3/SMoPT_tAI.transcription.protein.synthesis.pdf"), width=9.24, height=5.03, useDingbats=F)
grid.arrange(g1, g2, g3, g4, ncol=2)
dev.off()




