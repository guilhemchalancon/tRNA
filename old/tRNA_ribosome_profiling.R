# table1 <- read.csv("ribosome_profiling_gerashchenko2012/sd01.txt",header=2, sep="\t", skip=1)#[,1:2]
# colnames(table1) <- c("gene","05", "gene", "30")
# 
# table <- subset( melt(table1, id.vars = "gene"), !is.na(value))
# 
# head(table)
# 
# colnames(table) <- c("gene","time","log2_TE.change")
# write.table(table, "ribosome_profiling_gerashchenko2012/sd01_melt.txt",sep="\t", quote= F)

########################################################################################################################################################
#### INPUT ####
PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
setwd("~/Documents/MRC16/")
source("scripts/commands.R"); source("scripts/SMoT.R")

require(StingRay)
data(XPR)

# Define 
tRNAab.ox <- subset( tRNAab_quantified.foldchange_long, stress == "ox_20" )

# Ribosome profiling data in oxidative stress from Gerashchenko2012: Translation Efficiency
translation.efficiency.ox <- XPR$ribosome.profiling.Gerashchenko2012$sd01
colnames(translation.efficiency.ox)[1] <- "ORF"
translation.efficiency.ox$name  <- fNAME(translation.efficiency.ox$ORF)

  # Translational efficiency (TE) is a measure of how well translated a particular gene is relative to its mRNA abundance. TE can be defined as the number of footprints
  # divided by the number of mRNA-seq reads normalized to gene length and total number of reads, i.e., footprint in reads per kilobase per million mapped reads (rpkm)/mRNA rpkm. A higher TE value represents greater potency of mRNA for translation. TE was
  # used to examine translationally regulated genes. If a gene had a log2 (TE change) above 1.5 or below 1.5, it was considered up- or down-regulated, respectively.

write.table(translation.efficiency.ox, "results/master tables/translation.efficiency.ox.txt", sep="\t", quote=F, row.names=F)


# Ribosome profiling data in oxidative stress from Gerashchenko2012: RPKM measurements
ribosome.profiling.ox <- XPR$ribosome.profiling.Gerashchenko2012$sd02[1:5508,] # after 5508, clusters of gene IDs (e.g. YPL036W_YGL008C)
colnames(ribosome.profiling.ox)[1] <- "ORF"
ribosome.profiling.ox$name <- fNAME(ribosome.profiling.ox$ORF)
ribosome.profiling.ox <- subset(ribosome.profiling.ox, (Fp1_rpkm_30min+Fp2_rpkm_30min) > 10 & (mRNA1_rpkm_30min+mRNA2_rpkm_30min) > 10 ) 

data(PRO)
  lengths <- subset(PRO$properties, select=c(gene,protein.length))
  colnames(lengths)
  lengths$name <- fNAME(lengths$gene)
  lengths$normalised.length <- range01(lengths$protein.length)

ribosome.profiling.ox <- merge( ribosome.profiling.ox[,-1],  lengths[,-1], by="name" )
ribosome.profiling.ox <- within(ribosome.profiling.ox, TE.30 <- (Fp1_rpkm_30min+Fp2_rpkm_30min)/2 / (mRNA1_rpkm_30min + mRNA2_rpkm_30min)/2 )
ribosome.profiling.ox <- within(ribosome.profiling.ox, TE.00 <- (Initial_Fp1_rpkm+Initial_Fp2_rpkm)/2 / (Initial_mRNA1_rpkm) )


# Merge with stress-adjuted tRNA Adaptation Index data
#load("results/master tables/stAI.data.R")
load("results/master tables/SMoPT.table.Rd")
stAI.data.ox <- merge( subset(SMoPT.data, experiment == "ox" ), 
                       subset(translation.efficiency.ox, time == 30)[,c("name", "log2_TE.change")], 
                       by = "name")
stAI.data.ox <- merge( stAI.data.ox, ribosome.profiling.ox, by = "name")
stAI.data.ox$tAI.profile1 <- factor(stAI.data.ox$tAI.profile1, levels=c("decreased\nadaptation","similar\nadaptation","increased\nadaptation"))
save(stAI.data.ox, file = "results/master tables/stAI.data.ox.R")



# Enrichment data for all anticodons in all the gene set
anticodon.enrichment.per.gene.ox <- fread("results/master tables/anticodon.enrichment.per.gene.ox.txt", sep="\t", header = T, skip = 0)

# Anticodon data for ribosomal density analysis
anticodons <- list()
  anticodons$n.in.CDS    <-  cast( subset(anticodon.enrichment.per.gene.ox, time == "20"), ORF ~ anticodon, value = "n.in_CDS" )
  anticodons$freq.in.CDS <-  cast( subset(anticodon.enrichment.per.gene.ox, time == "20"), ORF ~ anticodon, value = "freq.in_CDS" )
  anticodons$OR.CDS      <-  cast( subset(anticodon.enrichment.per.gene.ox, time == "20"), ORF ~ anticodon, value = "OR.CDS" )
  anticodons$stAI        <-  subset(stAI.data.ox, experiment == "ox" )

  # add gene names
  anticodons[1:3] <- lapply(1:3, function(x) { anticodons[[x]]$name <- fNAME(anticodons[[x]]$ORF); return(anticodons[[x]]) })
  

# build a tables combining ribosome profiling and computed measures of anticodon frequency and stress-adjusted tRNA adaptation index 
profiling <- anticodons
profiling[1:3]  <- lapply(1:3, function(x){ 
  merge( subset(translation.efficiency.ox, time == "30"),
         anticodons[[x]],
         by = c("name","ORF") 
  ) }
)
profiling$stAI <- merge( subset(translation.efficiency.ox, time == "30"), anticodons$stAI, by = "name" )


########################################################################################################################################################
# tRNAs matching with codons that gained adaptability in Oxidative stress: CCG, UCC and UAA
improved <- as.character( subset(codon.master.table, experiment == "ox" & time == 20 & delta_w_FC > 0 )$anticodon )

# Do genes with the highest frequency of CCG|UCC|UAA anticodons show a higher translation efficiency?
improved.anticodons.freq <- merge( data.frame( 
  name = profiling$n.in.CDS$name,
  log2_TE.change = profiling$n.in.CDS$log2_TE.change, 
  f.CCG.UCC.UAA  = rowSums( profiling$n.in.CDS[,improved, drop=F] ) / rowSums( profiling$n.in.CDS[,  5:45]  )
  ),
  anticodons$stAI, by = "name"
)


ggplot( data = improved.anticodons.freq, aes( x =  f.CCG.UCC.UAA, y = log2_TE.change) ) + geom_point()

ggplot( data = improved.anticodons.freq, aes( x =  f.CCG.UCC.UAA, y = stAI) ) + geom_point()
ggplot( data = improved.anticodons.freq, aes( x =  f.CCG.UCC.UAA > 0.05, y = gain.rank) ) + geom_boxplot(notch=T)



# ------ # # ------ #
a <- melt(profiling$freq.in.CDS, id.vars = c("name","ORF","time","log2_TE.change"))

ggplot( data = a, aes( x = value, y = log2_TE.change )  ) + geom_point(alpha=0.5) + facet_wrap( ~ variable ) + theme_gul

# ------ # # ------ #

ggplot( data = profiling$stAI, aes( x = stAI, y = log2_TE.change )  ) + geom_point(alpha=0.5) + stat_smooth(method="lm") +  theme_gul
ggplot( data = profiling$stAI, aes( x = quant.gain, y = log2_TE.change )  ) + geom_boxplot(alpha=0.5) + theme_gul
ggplot( data = profiling$stAI, aes( x = quant.tAI, y = log2_TE.change )  ) + geom_boxplot(alpha=0.5) + theme_gul


rm(a)
rm(b)

########################################################################################################################################################
### [] Variations in s-tAI ####

g <- ggplot( subset(SMoPT.data,experiment!="normal"), aes( x = FS.tAI, y = FS.stAI ) ) + 
  geom_point(aes(color=tAI.profile)) + 
  geom_abline( intercept=0, slope =1, lty=2, color="gray20", alpha=0.95) +
  stat_smooth( method = "lm", color = "white", alpha=0.5) + coord_fixed() + scale_color_manual(values=c("red","gray40","orange")) +
  facet_wrap(  ~ experiment) + ggtitle("Variations in tRNA Adaptation Index during stress")

ggsave(plot=g, filename = paste0(PATH,"part3/adjusted.tAIFC--scatterplot.pdf"), width=8.1, height=8.1, useDingbats=FALSE )

########################################################################################################################################################
### [] Variations in ribosome profiling data ####

load("results/master tables/stAI.data.ox.R") # 
# [1] "name"                          "experiment"                    "tAI"                           "rank.tAI"                      "stAI"                         
# [6] "rank.stAI"                     "gain.tAI"                      "gain.rank"                     "quant.gain"                    "quant.tAI"                    
# [11] "FS.stAI"                       "FS.tAI"                        "FS.gain.tAI"                   "est.mRNA_abundance"            "log2.mRNA_abundance"          
# [16] "global.protein_synthesis.rate" "tAI.profile"                   "tAI.profile2"                  "log2_TE.change"               


      # Partial correlations
      require(ppcor)
      d <- stAI.data.ox[,c("log2_TE.change","log2.mRNA_abundance", "gain.stAI", "global.protein_synthesis.rate")]
      pcor( d[complete.cases(d),]  )
      
      # Wilcoxon tests Translation Efficiency ~ stAI groups
      wc.table <- wilcox.test.batch(stAI.data.ox, grouping.vars = "tAI.profile", value = "log2_TE.change")
     
    
      
      # Box plots Translation Efficiency ~ stAI groups
      g <- ggplot( data = stAI.data.ox, aes(x=tAI.profile1, y=log2_TE.change, fill=tAI.profile1)  ) + 
#          geom_hline(y= 1.4)+
#          geom_hline(y=-1.4)+
        scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
        labs(y="variation in translation efficiency (log2)\n(Gerashchenko et al., 2012)", x="", fill="") +
        theme_minimal() +
        theme(axis.text.x=element_blank())

        
      write.table(subset(stAI.data.ox, select = c(name, experiment, time, tAI.profile1, log2_TE.change, TE.00, TE.30 )), "~/Dropbox/Madan-Guilhem/update 0901/translation_efficiency_vs_tAIprofile.txt",sep="\t",quote=F, row.names=F )

ggplot( data = stAI.data.ox, aes(x=quantile.it(x = stAI, N=3), y=TE.30)  ) + 
  #          geom_hline(y= 1.4)+
  #          geom_hline(y=-1.4)+
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
#  labs(y="variation in translation efficiency (log2)\n(Gerashchenko et al., 2012)", x="", fill="") +
  theme_minimal() +
  theme(axis.text.x=element_blank())




# Gerashchenko et al. 2012 estimated changes in translation efficiency (log2-TE), based on measurements of ribosome occupancy made during response to oxidative stress. 
# We used these data to test if the genes that we predict to have "increased adaptation" based on the changes in tRNA abundance we measured during oxidative stress show an expected
# gain in translation efficiency during this stress (log2-TE>0). Indeed, we found that genes with increased adaptation had a median log2-TE of 0.12 
# (corresponding to a +8% increase in efficiency), and that their log2-TE distribution was significantly higher compared to 
# that of genes whose adaptation was predicted to remain similar during stress.
# (Wilcoxson Rank Sum test, P<4.28e-10, r=0.11), 
# Further, we found that genes whose adaptation was estimated to decrease consistantly showed a negative median log2-TE (-0.19, indicating a decrease of 12% in translation efficiency),
# and that this shift was significantly lower than genes whose adaptation was predicted to remain similar during stress (P<5.7e-05, r=0.07).
# Consistently, we found that differences in log2-TE where strongest between genes with increased and decreased adaptation (P<1.12e-14, r=0.21).
# Overall, these results indicate that in the case of oxidative stress, estimations of tRNA adaptation index are consistent with measured changes in ribosome occupancy on transcripts.

require(gridExtra)
      # Figure
      pdf(paste0(PATH,"part3/adjusted.tAIFC--translation.efficiency.pdf"), width=8.1, height=8.1, useDingbats=FALSE )
      grid.arrange( g %+% theme(legend.position="none"), arrangeGrob( tableGrob(wc.table[,c("X1","X2","Z","p.value","r","median.1","median.2","n.1","n.2")], 
                                                                                widths=unit(1,"null"), 
                                                                                heights=unit(1/(nrow(wc.table)),"npc"), 
                                                                                gpar.coltext=gpar(fontsize=8),
                                                                                gpar.coretext=gpar(fontsize=8), show.rownames =F ), # %+% theme_minimal()
                                                                      g_legend(g)
      ),
      nrow=1,
      widths = c(1/4, 3/4), main = "Ribosome profiling data shows changes in translation efficiency\n that are consistent with changes in tAI (oxidative stress)"  
      )
      dev.off()

########################################################################################################################################################

# Box plots Ribosome occupancy ~ stAI groups

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=protein.length, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  scale_y_log10() +
  theme(axis.text.x=element_blank())

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=Fp1_rpkm_30min, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  #labs(y="ribosome occupancy, t=30 (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  scale_y_log10() +
  theme(axis.text.x=element_blank())

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=Initial_Fp1_rpkm, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  #labs(y="ribosome occupancy, t=0 (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  scale_y_log10() +
  theme(axis.text.x=element_blank())

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=Initial_mRNA1_rpkm, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  #labs(y="mRNA abundance, t=0 (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  scale_y_log10() +
  theme(axis.text.x=element_blank())

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=Initial_Fp1_rpkm/Initial_mRNA1_rpkm, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  #labs(y="mRNA abundance, t=0 (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  scale_y_log10() +
  theme(axis.text.x=element_blank())

ggplot( data = stAI.data.ox, aes(x=tAI.profile, y=Fp1_rpkm_30min/mRNA1_rpkm_30min, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) + 
  #labs(y="mRNA abundance, t=0 (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  scale_y_log10() +
  theme(axis.text.x=element_blank())


ggplot( data = stAI.data.ox, aes(x=tAI.profile, y= (Fp1_rpkm_30min/Initial_Fp1_rpkm * Initial_mRNA1_rpkm/mRNA1_rpkm_30min), fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) +
  scale_y_log10() +
  #labs(y="ribosome occupancy (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  theme(axis.text.x=element_blank())

# not sure here.
ggplot( data = stAI.data.ox, aes(x=tAI.profile, y= TE.30/TE.00, fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) +
  scale_y_log10() +
  labs(y="ribosome occupancy (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  theme(axis.text.x=element_blank())


ggplot( data = stAI.data.ox, aes(x=tAI.profile, y= TE.30 * normalised.length , fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) +
  scale_y_log10() +
  labs(y="ribosome occupancy (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  theme(axis.text.x=element_blank())


ggplot( data = stAI.data.ox, aes(x=tAI.profile, y= TE.00 * normalised.length , fill=tAI.profile)  ) + 
  scale_fill_manual(values=c("red","gray40","orange")) + geom_boxplot(notch=T) +
  scale_y_log10() +
  labs(y="ribosome occupancy (RPKM)\n(Gerashchenko et al., 2012) - log scale", x="", fill="") +
  theme(axis.text.x=element_blank())
