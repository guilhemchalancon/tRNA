setwd("~/Documents/MRC16/")
source("scripts/commands.R")
#     NEW APPROACH
# >> anticodon enrichment ----- 
# l'affaire est ketchup
anticodon.usage.table <- read.table("results/master tables/anticodon.usage.table_10.txt",header=1)

anticodon.enrichment.per.gene <- merge( anticodon.usage.table, 
                                        ddply(subset(study.translation.kinetics, select = c(experiment, time, stress, ORF, ratio_initiation, ratio_total.time, variation_initiation_frequency, global.protein_synthesis.rate) ) , .(experiment, time, stress), function(x){  
                                          x$ratio_initiation.q   <- as.numeric(quantile.it( x$ratio_initiation, N = 7 ));
                                          x$variation_initiation.cat <- ifelse( x$variation_initiation_frequency > 0.5, "faster", ifelse( abs(x$variation_initiation_frequency) <= 0.5, "stable", "slower" ))
                                          x$ratio_total.time.q   <- as.numeric(quantile.it( x$ratio_total.time, N = 7 )); 
                                          x$gppr.q <- as.numeric(quantile.it( x$global.protein_synthesis.rate, N = 7 ))
                                          return(x) }), by = "ORF"
)

anticodon.enrichment.per.gene <- subset(anticodon.enrichment.per.gene, !is.na(ratio_initiation) )
anticodon.enrichment.per.gene$conditions <- paste(anticodon.enrichment.per.gene$experiment, counts.anticodon_in.ramp$time, sep="\n")
anticodon.enrichment.per.gene$conditions <- factor(anticodon.enrichment.per.gene$conditions, 
                                              levels=  paste( rep(c("diauxic","ox","osm","temp"), time = 2 ), rep(c(20,120), each = 4 ), sep="\n"),
                                              labels = paste( rep(c("diauxic shift","oxidative stress","osmotic stress","temperature stress"), time = 2 ), rep(c("20min","120min"), each = 4 ), sep="\n") )


CCG.enrichment.table <- subset(anticodon.enrichment.per.gene, anticodon == "CCG")

# this tables gives the aggregated frequency of each codon in the ramp a
anticodon.enrichment.aggregated <- ddply(subset(anticodon.enrichment.per.gene,time != 0), 
                                         .(experiment, time, stress, variation_initiation.cat, anticodon), summarise, 
                                         f.ramp = sum(n.in_ramp) / sum(L.ramp), a = sum(n.in_ramp) , b= sum(L.ramp) - sum(n.in_ramp) )


write.table(subset(anticodon.enrichment.per.gene, experiment=="ox"), "results/master tables/anticodon.enrichment.per.gene.ox.txt", sep="\t", quote=F, row.names=F)


#anticodon.enrichment.aggregated <- ddply(anticodon.enrichment.aggregated, .(experiment, time, stress, ratio_initiation.q), summarise, sum( f.ramp ))

head(anticodon.enrichment.aggregated)
require(reshape2)


a <- dcast( anticodon.enrichment.aggregated, formula = experiment + time + stress + anticodon ~  variation_initiation.cat, value.var = "a" )
colnames(a)[5:ncol(a)] <- paste( "a", colnames(a)[5:ncol(a)], sep=".")
b <- dcast( anticodon.enrichment.aggregated, formula = experiment + time + stress + anticodon ~  variation_initiation.cat, value.var = "b" )
colnames(b)[5:ncol(b)] <- paste( "b", colnames(b)[5:ncol(b)], sep=".")
counts.anticodon_in.ramp <- merge(a,b, by=c("experiment","time","stress","anticodon"))

head(counts.anticodon_in.ramp)
z <- qnorm(1-0.05/2)

# compare the top bin to the bottom bin
#     counts.anticodon_in.ramp$RR.1_7 <- counts.anticodon_in.ramp$a1 / counts.anticodon_in.ramp$a7 
#     counts.anticodon_in.ramp$OR.1_7 <- counts.anticodon_in.ramp$a1 / counts.anticodon_in.ramp$a7 * counts.anticodon_in.ramp$b7 / counts.anticodon_in.ramp$b1
#     counts.anticodon_in.ramp$SE.1_7 <- round( sqrt( 1 / counts.anticodon_in.ramp$a1 + 1 / counts.anticodon_in.ramp$a7 + 1 / counts.anticodon_in.ramp$b1 + 1 / counts.anticodon_in.ramp$b7 ), 3)
#     counts.anticodon_in.ramp$CIlow.1_7 <-  round( exp( log( counts.anticodon_in.ramp$OR.1_7 ) - 1.96*counts.anticodon_in.ramp$SE.1_7 ), 3)
#     counts.anticodon_in.ramp$CIhigh.1_7 <- round( exp( log( counts.anticodon_in.ramp$OR.1_7 ) + 1.96*counts.anticodon_in.ramp$SE.1_7 ), 3)
#     counts.anticodon_in.ramp$consistency.1_7  <- sign( log(counts.anticodon_in.ramp$CIlow.1_7) * log(counts.anticodon_in.ramp$CIhigh.1_7)  ) 
#     
#     # compare the top bin to all other bins (meh?)
#     counts.anticodon_in.ramp$OR.1_rest <- counts.anticodon_in.ramp$a1 / ( rowSums(counts.anticodon_in.ramp[,paste0("a",2:7)]) ) * rowSums(counts.anticodon_in.ramp[,paste0("b",2:7)]) / counts.anticodon_in.ramp$b1
#     counts.anticodon_in.ramp$SE.1_rest <- round( sqrt( 1 / counts.anticodon_in.ramp$a1 + 1 / rowSums(counts.anticodon_in.ramp[,paste0("a",2:7)]) + 1 / counts.anticodon_in.ramp$b1 + 1 / rowSums(counts.anticodon_in.ramp[,paste0("b",2:7)]) ), 3)
#     counts.anticodon_in.ramp$CIlow.1_rest <-  round( exp( log( counts.anticodon_in.ramp$OR.1_rest ) - 1.96*counts.anticodon_in.ramp$SE.1_rest ), 3)
#     counts.anticodon_in.ramp$CIhigh.1_rest <- round( exp( log( counts.anticodon_in.ramp$OR.1_rest ) + 1.96*counts.anticodon_in.ramp$SE.1_rest ), 3)
#     counts.anticodon_in.ramp$consistency.1_rest  <- sign( log(counts.anticodon_in.ramp$CIlow.1_rest) * log(counts.anticodon_in.ramp$CIhigh.1_rest)  ) 

# rm(a,b, z)

# compare the top bin to the bottom bin
    counts.anticodon_in.ramp$RR.faster_slower <- counts.anticodon_in.ramp$a.faster / counts.anticodon_in.ramp$a.slower 
    counts.anticodon_in.ramp$OR.faster_slower <- counts.anticodon_in.ramp$a.faster / counts.anticodon_in.ramp$a.slower * counts.anticodon_in.ramp$b.slower / counts.anticodon_in.ramp$b.faster
    counts.anticodon_in.ramp$SE.faster_slower <- round( sqrt( 1 / counts.anticodon_in.ramp$a.faster + 1 / counts.anticodon_in.ramp$a.slower + 1 / counts.anticodon_in.ramp$b.faster + 1 / counts.anticodon_in.ramp$b.slower ), 3)
    counts.anticodon_in.ramp$CIlow.faster_slower <-  round( exp( log( counts.anticodon_in.ramp$OR.faster_slower ) - 1.96*counts.anticodon_in.ramp$SE.faster_slower ), 3)
    counts.anticodon_in.ramp$CIhigh.faster_slower <- round( exp( log( counts.anticodon_in.ramp$OR.faster_slower ) + 1.96*counts.anticodon_in.ramp$SE.faster_slower ), 3)
    counts.anticodon_in.ramp$consistency.faster_slower  <- sign( log(counts.anticodon_in.ramp$CIlow.faster_slower) * log(counts.anticodon_in.ramp$CIhigh.faster_slower )  ) 
    
    # compare the top bin to all other bins (meh?)
    counts.anticodon_in.ramp$OR.faster_rest <- counts.anticodon_in.ramp$a.faster / ( rowSums(counts.anticodon_in.ramp[,c("a.stable","a.slower")]) ) * rowSums(counts.anticodon_in.ramp[,c("b.stable","b.slower")]) / counts.anticodon_in.ramp$b.faster
    counts.anticodon_in.ramp$SE.faster_rest <- round( sqrt( 1 / counts.anticodon_in.ramp$a.faster + 1 / rowSums(counts.anticodon_in.ramp[,c("a.stable","a.slower")]) + 1 / counts.anticodon_in.ramp$b.faster + 1 / rowSums(counts.anticodon_in.ramp[,c("b.stable","b.slower")]) ), 3)
    counts.anticodon_in.ramp$CIlow.faster_rest <-  round( exp( log( counts.anticodon_in.ramp$OR.faster_rest ) - 1.96*counts.anticodon_in.ramp$SE.faster_rest ), 3)
    counts.anticodon_in.ramp$CIhigh.faster_rest <- round( exp( log( counts.anticodon_in.ramp$OR.faster_rest ) + 1.96*counts.anticodon_in.ramp$SE.faster_rest ), 3)
    counts.anticodon_in.ramp$consistency.faster_rest  <- sign( log(counts.anticodon_in.ramp$CIlow.faster_rest) * log(counts.anticodon_in.ramp$CIhigh.faster_rest)  ) 



counts.anticodon_in.ramp <- merge( counts.anticodon_in.ramp, anticodon.master.table, by = c("anticodon","experiment","time","stress") )

counts.anticodon_in.ramp$conditions <- paste(counts.anticodon_in.ramp$experiment, counts.anticodon_in.ramp$time, sep="\n")
counts.anticodon_in.ramp$conditions <- factor(counts.anticodon_in.ramp$conditions, 
                        levels=  paste( rep(c("diauxic","ox","osm","temp"), time = 2 ), rep(c(20,120), each = 4 ), sep="\n"),
                        labels = paste( rep(c("diauxic shift","oxidative stress","osmotic stress","temperature stress"), time = 2 ), rep(c("20min","120min"), each = 4 ), sep="\n") )


# focus on the special case of anticodons that seem most consistently enriched / depleted
#counts.anticodon_in.ramp <- subset(counts.anticodon_in.ramp, select = c(experiment, time, stress, anticodon, OR.1_7, CIlow.1_7, CIhigh.1_7, consistency) )
# anticodon.difference.ramp.FTSvsSTS <- subset(counts.anticodon_in.ramp, consistency.1_7 == 1)
# anticodon.difference.ramp.FTSvsSTS <- ddply(anticodon.difference.ramp.FTSvsSTS , .(stress), transform, n.tRNAs = length(unique(anticodon)) )

anticodon.enrichment.ramp_faster.vs.rest <- subset(counts.anticodon_in.ramp, consistency.faster_rest == 1)
anticodon.enrichment.ramp_faster.vs.rest <- ddply(anticodon.enrichment.ramp_faster.vs.rest , .(stress), transform, n.tRNAs = length(unique(anticodon)) )

#CCG.enrichment.table <- ddply(CCG.enrichment.table, .(experiment, time, stress), function(x){  
#  x$ratio_initiation.q <- as.numeric(quantile.it( x$ratio_initiation, N = 7 )); 
#  return(x) })





# [] consistent depletion/enrichment of anticodons in the ramp of FTS and STS genes -----
g <- ggplot( data= anticodon.enrichment.ramp_faster.vs.rest, aes( x = reorder(anticodon, OR.faster_rest), y = OR.faster_rest - 1 )) + 
  geom_bar( stat="identity", origin = 1, aes(width = 0.9*n.tRNAs/max(n.tRNAs), fill = CAI.O ) ) +
  geom_hline( yintercept = 0 ) +
  geom_errorbar( aes(ymin = CIlow.faster_rest -1, ymax= CIhigh.faster_rest -1), 
                 width=0.1, size=0.3, color="gray40" ) + 
  facet_wrap( ~ conditions, scales = "free", nrow = 2 ) + scale_y_continuous( labels = function(x){x+1} ) +
  scale_fill_manual(values = c("skyblue3","gray","orange") ) + 
  labs(x="",y = "OR",
       title="depletion/enrichment of anticodons in the ramp of genes initiated faster") +
  theme_gul + theme(legend.position="none")

ggsave(plot = g, filename = paste0(PATH,"part3/fig_results_ramp.enrichment.pdf"), 
       dpi=250, width = 12.7, height = 6.5, useDingbats=F )



## too big figure to plot frequently
# g <- ggplot( data = subset(anticodon.enrichment, time == 20 & !is.na(preference_in.CDS)), aes( x = faster.translation, y = preference_in.CDS, fill = faster.translation )  ) + 
#   geom_boxplot( notch = T) + 
#   scale_y_continuous( labels = prettyNum ) +
#   scale_fill_manual( values = c("gray40","red")) +
#   facet_wrap( anticodon ~ stress )
# 
# ggplot(plot=g, filename=paste0(PATH,"part3/anticodon_enrichment--x.faster.translation_y.enrichment_z.stress20min_t.anticodon.png"), 
#        dpi=250, width=12, height=12)
# 
# # [] distribution of anticodon frequency based on ratio initiation ----
# g <- ggplot( data = subset(anticodon.enrichment, time==20 & !is.na(ratio_initiation) & !is.na(preference_in.CDS) ), 
#              aes( x= log2(ratio_initiation), y = preference_in.CDS ) ) + 
#   geom_point( aes(color = faster.translation ), size = 1 ) +
#   geom_vline(xintercept=0, size=0.25) +
#   scale_x_continuous( labels = prettyNum ) +
#   scale_y_continuous( labels = prettyNum ) + 
#   facet_wrap( anticodon ~ experiment ) + 
#   scale_color_manual(values = c("gray40","red")) + 
#   theme(legend.position="none")
# 
# ggsave(plot=g, filename = paste0(PATH,"part3/anticodon_enrichment--x.ratio_initiation_y.anticodon.frequency_z.stress_t.anticodon.png"),
#         dpi=250, width = 11.7*1.10, height=9.39*1.10)

#--------------------------------------------------------------#
require(grid)
require(scales)
# grutte de grulle
# [] tRNA abundance vs anticodon enrichment in the first 20 codons -----
g <- ggplot( data = counts.anticodon_in.ramp, aes( x = total.tRNAab, y = OR.faster_rest ) ) +
  geom_hline( yintercept = 1, lty= 3 ) +
  geom_point( data = subset(counts.anticodon_in.ramp, consistency.faster_rest == 1), size = 8, color="gray40", fill="white", pch=21) + 
  geom_point( aes(color = CAI.O), alpha = 0.7, size = 5 ) +
  geom_text(  aes(label= anticodon), size = 2.25) +
  scale_y_continuous( labels = prettyNum ) +
  scale_x_log10( labels = scientific_10 ) +
#   scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x),
#                  labels = trans_format("log10", math_format(10^.x))
#   ) +
  scale_size_continuous( range = c(6,9)) +
  facet_wrap(  ~ conditions, nrow = 2 ) + 
  labs( y= "Odd ratios\nof finding a given anticodon in the ramp\nof genes with faster initiation compared to other genes",
        x= "tRNA abundance", 
        title = ""
  ) +
  scale_colour_manual( values = c("skyblue3","gray60","orange") ) + theme(legend.position="none")

ggsave( plot = g, dpi = 250, width = 12.7, height = 7.5,useDingbats=F,
        filename = paste0(PATH,"part3/fig_results_CCG_tRNAab.pdf") )




g <- ggplot( data = counts.anticodon_in.ramp, aes( x = log2.foldchange, y = OR.faster_rest ) ) +
  geom_hline( yintercept = 1, lty = 3 ) +
  geom_vline( xintercept = 0, lty = 3 ) +
  geom_point( data = subset(counts.anticodon_in.ramp, consistency.faster_rest == 1), size = 8, color="gray40", fill="white", pch=21) + 
  geom_point( aes(color = CAI.O), alpha = 0.7, size = 5 ) +
  geom_text(  aes(label= anticodon), size = 2.25) +
  scale_y_continuous( labels = prettyNum ) +
  scale_x_continuous( labels = prettyNum ) +
  #   scale_x_log10( breaks = trans_breaks("log10", function(x) 10^x),
  #                  labels = trans_format("log10", math_format(10^.x))
  #   ) +
  scale_size_continuous( range = c(6,9)) +
  facet_wrap(  ~ conditions, nrow = 2 ) + 
  labs( y= "Odd ratios\nof finding a given anticodon in the ramp\nof genes with faster initiation compared to other genes",
        x= "tRNA abundance fold change (log2)", 
        title = ""
  ) +
  scale_colour_manual( values = c("skyblue3","gray60","orange") ) + theme(legend.position="none")

ggsave( plot = g, dpi = 250, width = 12.7, height = 7.5,useDingbats=F,
        filename = paste0(PATH,"part3/fig_results_CCG_log2.tRNAab.pdf") )












# [] relative tRNA availability vs anticodon enrichment in the first 20 codons -----
g <- ggplot( data = counts.anticodon_in.ramp, aes( x = OR.faster_rest, y = relative.availability ) ) +
  geom_vline( xintercept = 1 ) +
  geom_point( data = subset(counts.anticodon_in.ramp, consistency.faster_rest == 1), size = 8, color="gray40", fill="white", pch=21) + 
  geom_point( aes(color = CAI.O), alpha = 0.7, size = 5 ) +
  geom_text(  aes(label= anticodon), size = 2.25) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_size_continuous( range = c(6,9)) +
  facet_wrap(  ~ conditions, nrow = 2 ) + 
  labs( x= expression(atop("Odd ratio",atop("of finding anticodon j in the ramp of FTS genes compared to STS genes",""))),
        y= expression(atop("Relative availability of tRNA",atop("compared to synonymous tRNAs",""))), 
        title = ""
  ) +
  scale_colour_manual( values = c("skyblue3","gray60","orange") ) + theme(legend.position="none")

ggsave( plot = g, dpi = 250, width = 12.7, height = 7.5,useDingbats=F,
        filename = paste0(PATH,"part3/anticodon_enrichment--x.enrichment_y.relative.availability_z.stress.pdf") )




# [] change in relative tRNA availability vs anticodon enrichment in the first 20 codons -----
g <- ggplot( data = subset(counts.anticodon_in.ramp,relative.availability.normal !=1), aes( x = OR.faster_rest, y = relative.availability / relative.availability.normal ) ) +
  geom_vline( xintercept = 1, color = "gray40" ) + geom_hline( yintercept = 1, color = "gray40" ) +
  geom_point( data = subset(counts.anticodon_in.ramp, relative.availability.normal !=1 & relative.availability / relative.availability.normal > 1 & OR.faster_rest > 1 ), 
              size = 8, color="gray40", fill="white", pch=21) + 
  geom_point( aes(color = CAI.O), alpha = 0.7, size = 5 ) +
  geom_text(  aes(label= anticodon), size = 2.25) +
  scale_x_continuous( labels = prettyNum ) +
  scale_y_continuous( labels = prettyNum ) +
  scale_size_continuous( range = c(6,9)) +
  facet_wrap(  ~ conditions, nrow = 2 ) + 
  labs( x= expression(atop("Odd ratio",atop("of finding anticodon j in the ramp of FTS genes compared to other genes",""))),
        y= expression(atop("Change in relative availability of tRNA",atop("compared to synonymous tRNAs",""))), 
        title = ""
  ) +
  scale_colour_manual( values = c("skyblue3","gray60","orange") ) + theme(legend.position="none")

ggsave( plot = g, dpi = 250, width = 12.7, height = 7.5,useDingbats=F,
        filename = paste0(PATH,"part3/anticodon_enrichment--x.enrichment_y.chagne.relative.availability_z.stress.pdf") )







# [] does CGG enrichment associates with initiation kinetics? ----- 
g <- ggplot( data = subset(anticodon.enrichment.per.gene, anticodon == "CCG" & time !=0 & !is.na(ratio_initiation)), 
             aes( x = ifelse(variation_initiation.cat == "faster","faster","other") ,
                  y = preference_in.CDS, 
                  fill = factor(ifelse(variation_initiation.cat == "faster","faster","other"))) ) + 
  geom_boxplot(notch=T) +
  facet_wrap( ~ conditions,  ncol = 4 ) +
  labs(title= "CCG enrichment in genes with various changes in initiation kinetics", 
       x = "relative decrease in total translation time",
       y = "%CCG among Arg position / gene") +
  #scale_x_discrete( labels = c("faster","slowed\ndown","intermediate")) +
  scale_fill_brewer( "clarity" ) +
  theme_gul + theme(legend.position="none")

ggsave( plot = g, filename = paste0(PATH,"part3/fig_results_CCG_boxplot_prefCDS.pdf"),
        dpi = 250, width= 12.8, height = 7.6, useDingbats=F)



# [] OD ratios for having CCG in the CDS in group of genes with increasing gain in translation speed during stress -----
g <- ggplot( data = subset(anticodon.enrichment.per.gene, anticodon == "CCG" & time !=0 & !is.na(ratio_initiation)), 
             aes( x = ifelse(variation_initiation.cat == "faster","faster","other"),
                  y = OR.CDS, 
                  fill = ifelse(variation_initiation.cat == "faster","faster","other") )) + 
  geom_boxplot( notch =T ) +
  facet_wrap( ~ conditions,  ncol = 4 ) +
  labs(title= "CCG enrichment in genes with various changes in initiation kinetics", 
       x = "relative decrease in total translation time",
       y = "Odd-ratio (CCG in CDS)/ gene") +
  scale_x_discrete( labels = c("slowed\ndown","","","intermediate","","","speed\nup")) +
  scale_fill_brewer( "clarity" ) +
  theme_gul + theme(legend.position="none")

ggsave( plot = g, filename = paste0(PATH,"part3/fig_results_CCG_boxplot_OR.pdf"),
        dpi = 250, width= 12.8, height = 7.6)

# ---------------------------------------------------------------------------------------------------------------- # 
# [*] CCG influence in the codon ramp ####


#anticodon.master.table  <- read.table("results/master tables/anticodon.master.table.txt", header=1) # especially adjusted tGCN, anticodon supply, relative availability
codon.master.table      <- read.table("results/master tables/codon.master.table.txt", header=1) # especially relative adaptiveness and anticodon deman
# codon.freq.matrix <- as.matrix(read.table( file="results/stress tAI/codon.freq.mat",sep="\t",header=1))
# codon.freq.matrix.ramp <- as.matrix(read.table( file="results/stress tAI/codon.freq.matrix_first20codons.mat",sep="\t",header=1))
# codon.freq.matrix.outside.ramp <- as.matrix(read.table( file="results/stress tAI/codon.freq_outside.ramp.mat",sep="\t",header=1))
codon.eff.matrix <- as.matrix(read.table( file="results/stress tAI/codon.eff.mat",sep="\t",header=1))
codon.eff.matrix.ramp <- as.matrix(read.table( file="results/stress tAI/codon.eff.matrix_first10codons.mat",sep="\t",header=1))
codon.eff.matrix.outside.ramp <- as.matrix(read.table( file="results/stress tAI/codon.eff_outside.ramp.mat",sep="\t",header=1))

#     rownames(codon.eff.matrix.20codons) == rownames(codon.eff.matrix.outside.ramp) 
#     rownames(codon.eff.matrix.20codons) == rownames(codon.eff.matrix)
    all( rowSums(codon.eff.matrix.ramp) + rowSums(codon.eff.matrix.outside.ramp) - rowSums(codon.eff.matrix) == 0 )



all( colnames(codon.eff.matrix) == colnames(codon.eff.matrix.ramp) )
to.ignore <- which( colnames(codon.eff.matrix) %in% c("cgg","tga","taa","tag") )
stop.codons <- which( colnames(codon.eff.matrix) %in% c("tga","taa","tag") )
target.codon <- which( colnames(codon.eff.matrix) == "cgg" )

# computed in prepare.codon.table.R
stress.codon.adaptiveness  <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1)


cgg.table <- subset(codon.master.table, codon == "cgg", select = c(codon, anticodon, experiment, time, demand.mRNAab, demand.events,
                                                                   w.2, av.elongation_time, num_events,
                                                                   foldchange, log2.foldchange, total.tRNAab
) )

 
stress.tAI_noCCG.genes <- unique(
  # exclude stop codons and methionine from calculation (As in dosReis2004)
  arrange( ddply(subset(stress.codon.adaptiveness, !codon %in% c("taa","tga","tag") ),  .(experiment, time), 
                 function(d){ 
                   data.frame( ORF = rownames(codon.eff.matrix), 
                               L = rowSums(codon.eff.matrix),
                               stAI.2       =   get.tai(x=codon.eff.matrix[,-c(11,12,15)], w.data = d, score = "w.2" ),
                               stAI.ramp           =   get.tai(x=codon.eff.matrix.ramp[,-stop.codons],   w.data = d, score = "w.2" ),
                               stAI.not.ramp       =   get.tai(x=codon.eff.matrix.outside.ramp[,-stop.codons],   w.data = d, score = "w.2" ),
                               stAI.not.cgg.ramp   =   get.tai(x=codon.eff.matrix.ramp[,-to.ignore],   w.data = d, score = "w.2" ),
                               f.ramp.cgg   =  codon.eff.matrix.ramp[,target.codon]/9,
                               f.cgg   =   codon.eff.matrix[,target.codon]/rowSums(codon.eff.matrix),
                               w.cgg   =   unique(subset(d, codon=="cgg")$w.2)
                   ) 
                 }), 
           ORF, experiment, time)
)

# get.tai works equally well with codon frequencies or absolute numbers (see command's script)
# all( abs(get.tai(x=codon.eff.matrix[,-stop.codons],   w.data = d, score = "w.2" ) - get.tai(x=codon.freq.matrix[rownames(codon.eff.matrix),-stop.codons],   w.data = d, score = "w.2" )) < 10^{-5} )

#stress.tAI_noCAA.genes$hypothetical.stAI <- with(stress.tAI_noCAA.genes, w.ttg^f.ttg * stAI.ttg^{1-f.ttg} )
stress.tAI_noCCG.genes$stAI.ramp_f           <- with(stress.tAI_noCCG.genes, stAI.ramp^{9/L} )
stress.tAI_noCCG.genes$stAI.not.ramp_1_f     <- with(stress.tAI_noCCG.genes, stAI.not.ramp^{1- 9/L} )
stress.tAI_noCCG.genes$stAI.not.cgg.ramp_1_f <- with(stress.tAI_noCCG.genes, stAI.not.cgg.ramp^{1-f.ramp.cgg} )
stress.tAI_noCCG.genes$cgg.impact.ramp       <- with(stress.tAI_noCCG.genes, w.cgg^f.ramp.cgg )
# indicator variable:
stress.tAI_noCCG.genes$cgg.in.ramp           <- with(stress.tAI_noCCG.genes, f.ramp.cgg > 0 )



# prepare data set
CCG.table <- arrange( merge(study.translation.kinetics , subset(anticodon.usage.table, anticodon == "CCG"), by = c("name","ORF"), all.x=T), plyr::desc(freq.in_CDS) )
CCG.table <- ddply(CCG.table, .(label), mutate, scaled.initiation_frequency = scale( 1/av.initiation_time ) )


CCG.influence <- merge(CCG.table, stress.tAI_noCCG.genes, by = c("ORF","experiment","time"))
CCG.influence$initiation_freq <- 1/CCG.influence$av.initiation_time * 60
  write.table(CCG.influence, file="results/master tables/CCG.influence.txt",sep="\t",quote=F, row.names=F)
  save(CCG.influence, file="results/master tables/CCG.influence.Rd")

      head(with( CCG.influence, data.frame(stAI, stAI.ramp, stAI.not.ramp, close = abs(stAI - stAI.ramp) < 10^{-5} )) )




# verify calculations are correct
table( with(CCG.influence, abs(stAI.ramp_f * stAI.not.ramp_1_f - stAI) < 10^{-5} ) ) # not ok at all (but explained: I changed the calculation of tGCN_s which should also modify w^s and thus stAI )
table( with(CCG.influence, abs(stAI.ramp_f * stAI.not.ramp_1_f - stAI.2) < 10^{-5} ) ) 
with(CCG.influence, cor.test(stAI, stAI.2)) # but this is not a big deal, since overall stAI with round(tGCN_s) or round(tGCN_s,1) correlate at >0.994) 
  all( with(CCG.influence, abs(cgg.impact.ramp * stAI.not.cgg.ramp_1_f - stAI.ramp) < 10^{-5} ) ) # ok 
table( with(CCG.influence, abs(stAI.ramp - stAI) < 10^{-5} ) )  # ok-ish (should be all false)
  all( with(CAA.influence, abs(ttg.impact * stAI.ttg ^ {1-f.ttg} - stAI) < 10^{-5} ) ) # ok



  CCG.influence


#####-----------------------------------------------------------------------------------------------------------------------######

load("results/master tables/CCG.influence.Rd")

notCCG.cor <- matrixify(reshape2::dcast( CCG.influence, name ~ experiment + time, value.var = "stAI.not.cgg.ramp" ))
notCCG.cor <- notCCG.cor[complete.cases(notCCG.cor),]
cor(notCCG.cor)

summary(CCG.influence$cgg.impact.ramp)
subset(stress.codon.adaptiveness, codon == "cgg")

# select genes with measurements in all 8 conditions (exluded genes will be those for which the stochastic simulation occcasionally returned NaN for at least one of their conditions)
genes.to.keep.ramp <- as.character( subset( ddply(CCG.influence, .(name), summarise, n.measures = length(time)), n.measures == 8)$name )
length(genes.to.keep.ramp)

# select genes with measurements in all 8 conditions (exluded genes will be those for which the stochastic simulation occcasionally returned NaN for at least one of their conditions)
genes.to.keep <- as.character( subset( ddply(subset(CCG.influence,f.ramp.cgg>0), .(name), summarise, n.measures = length(time)), n.measures == 8)$name )
length(genes.to.keep)




# [] Apply multiple regressions for all genes, and measure the relative importance of the  #####
require(yhat)

# --- influence of the ramp vs rest of CDS on initiation frequency ----
all.genes.regressions.ramp_vs_rest <- lapply( setNames(genes.to.keep.ramp, genes.to.keep.ramp), function(g){
  #print(g)
  data <- subset(CCG.influence, name == g)
  #data$initiation_freq_fc <- 1/data$ratio_initiation
  ## Regression
  lm.out<-lm( initiation_freq ~ stAI.ramp_f + stAI.not.ramp_1_f, data= data) #  
  regrOut<-calc.yhat(lm.out)
}
)


ramp.relative.importance <- do.call(rbind, lapply( names( all.genes.regressions.ramp_vs_rest ), function(x){
  pred <- data.frame(all.genes.regressions.ramp_vs_rest[[x]]$PredictorMetrics)
  pred$variable  <- rownames(pred)
  
  unique <- reshape2::dcast(data = pred, x ~ variable, value.var = "Unique")[,c("x",rownames(pred))]
  colnames(unique) <- c("name","U.ramp","U.out.ramp","U.total")
  
  common <- reshape2::dcast(data = pred, x ~ variable, value.var = "Common")[,c("x",rownames(pred))]
  colnames(common) <- c("name","C.ramp","C.out.ramp","C.total")
  
  r <- reshape2::dcast(data = pred, x ~ variable, value.var = "r")[,c("x",rownames(pred))]
  colnames(r) <- c("name","r.ramp","r.out.ramp","r.total")
  
  rs <- reshape2::dcast(data = pred, x ~ variable, value.var = "rs")[,c("x",rownames(pred))]
  colnames(rs) <- c("name","rs.ramp","rs.out.ramp","rs.total")
  
  rs2 <- reshape2::dcast(data = pred, x ~ variable, value.var = "rs2")[,c("x",rownames(pred))]
  colnames(rs2) <- c("name","rs2.ramp","rs2.out.ramp","rs2.total")
  
  # output
  data.frame(  unique, 
               common[,-1],
               r[,-1],
               rs[,-1],
               rs2[,-1], 
               'ramp.dominanted' = with( rs2, rs2.ramp > rs2.out.ramp ),
               check.names = F
  )
}
))

dim(ramp.relative.importance)


# --- adaptiveness of cgg in ramp on initiation frequency ----
require(yhat)
all.genes.regressions.cgg_ramp <- lapply( setNames(genes.to.keep, genes.to.keep), function(g){
  #print(g)
  data <- subset(CCG.influence, name == g)
  #data$initiation_freq_fc <- 1/data$ratio_initiation # already there
  ## Regression
  lm.out<-lm( initiation_freq ~ cgg.impact.ramp + stAI.not.cgg.ramp_1_f  + stAI.not.ramp, data= data) #  
  regrOut<-calc.yhat(lm.out)
}
)

    # Retrieve specific data from the list 
    CCG.relative.importance <- do.call(rbind, lapply( names( all.genes.regressions.cgg_ramp ), function(x){
      pred <- data.frame(all.genes.regressions.cgg_ramp[[x]]$PredictorMetrics)
      pred$variable  <- rownames(pred)
      
      unique <- reshape2::dcast(data = pred, x ~ variable, value.var = "Unique")[,c("x",rownames(pred))]
      colnames(unique) <- c("name","U.cgg","U.stAI.cgg","U.not.ramp","U.total")
      
      common <- reshape2::dcast(data = pred, x ~ variable, value.var = "Common")[,c("x",rownames(pred))]
      colnames(common) <- c("name","C.cgg","C.stAI.cgg","C.not.ramp","C.total")
      
      r <- reshape2::dcast(data = pred, x ~ variable, value.var = "r")[,c("x",rownames(pred))]
      colnames(r) <- c("name","r.cgg","r.stAI.cgg","r.not.ramp","r.total")
      
      rs <- reshape2::dcast(data = pred, x ~ variable, value.var = "rs")[,c("x",rownames(pred))]
      colnames(rs) <- c("name","rs.cgg","rs.stAI.cgg","rs.not.ramp","rs.total")
      
      rs2 <- reshape2::dcast(data = pred, x ~ variable, value.var = "rs2")[,c("x",rownames(pred))]
      colnames(rs2) <- c("name","rs2.cgg","rs2.stAI.cgg","rs2.not.ramp","rs2.total")
      
      # output
      data.frame(  unique, 
                   common[,-1],
                   r[,-1],
                   rs[,-1],
                   rs2[,-1], 
                   'cgg.dominanted' = with( rs2, rs2.cgg > rs2.stAI.cgg & rs2.cgg > rs2.not.ramp),
                   check.names = F
      )
    }
    ))
   
    CCG.relative.importance.data <- merge(CCG.relative.importance, unique(subset(CCG.influence, select = c(name, n.in_ramp, n.in_CDS, L.CDS, freq.in_CDS, OR.CDS, OR.CDS_CI.low, cgg.in.ramp  )) ), by="name")

# ------ Quick look-up ----

  subset(CCG.relative.importance.data, cgg.dominanted == T)
  subset(ramp.relative.importance, ramp.dominanted == T)


#-------- ANALYSIS OF THE INFLUENCE OF Arg-CCG ON INITIATION FREQUENCY ------


    # Do I find an association between the presence of a cgg codon in the ramp and an increase in initiation frequency?
    chisqr.CGG <- ddply( CCG.influence, .(experiment,time), function(x){
      #x <- subset(CCG.influence, experiment == "diauxic" & time == 20)
      test <- with(x, chisq.test( x = n.in_ramp > 0, y = variation_initiation_frequency > 0   ))
      data.frame( experiment = unique(x$experiment), time = unique(x$time), chi2= test$statistic, n= sum(test$observed), phi = sqrt(test$statistic/sum(test$observed)), df = test$parameter, p.value = test$p.value, p.sign = ifelse( test$p.value < 0.05, T,F) )  
    } )
    chisqr.CGG$adjusted.p <- noquote( format.pval( p.adjust(chisqr.CGG$p.value, method = "fdr"),digits = 3 ))
    
    data.chisqr.CGG <- merge( chisqr.CGG, subset(anticodon.master.table, anticodon == "CCG", select=c(experiment, time, log2.foldchange, foldchange)), by = c("experiment","time"))   
      
    
    ggplot(data.chisqr.CGG, aes(x= phi, y= log2.foldchange )) + geom_point(pch=21, fill="white") + 
      labs(y="tRNA abundance fold change (log2)", x=expression(phi)) +
      stat_smooth(method="lm", color = "gray30", se=F)
#     ggplot(data.chisqr.CGG, aes(x= -log10(p.value), y= log2.foldchange )) + geom_point()
#     ggplot(data.chisqr.CGG, aes(x= -log10(p.value), y= phi )) + geom_point()





g.cgg1 <- ggplot( CCG.relative.importance.data, aes(x=U.stAI.cgg,  y = U.cgg )) + 
#  geom_point(data = subset(CAA.relative.importance.data, ttg.rich == F), pch=21, fill="gray70") + 
  geom_point(data = subset(CCG.relative.importance.data), aes(fill =  U.stAI.cgg < U.cgg ), pch=21) + 
  geom_abline(intercept=0, slope=1, lty=3) + 
  geom_text( data = subset(CCG.relative.importance, U.stAI.cgg < U.cgg & U.cgg > 0.2), aes(label = name), size = 3, vjust = 1.9  ) +
  theme(legend.position="none") +
  scale_x_continuous(labels=prettyNum) + scale_y_continuous(labels=prettyNum) + scale_fill_manual(values =c("gray","#0899FF"))



g.cgg2 <- ggplot( CCG.relative.importance.data, aes(x=U.not.ramp,  y = U.cgg )) + 
  #  geom_point(data = subset(CAA.relative.importance.data, ttg.rich == F), pch=21, fill="gray70") + 
  geom_point(data = subset(CCG.relative.importance.data), aes(fill =  U.not.ramp < U.cgg ), pch=21) + 
  geom_abline(intercept=0, slope=1, lty=3) +
  geom_text( data = subset(CCG.relative.importance, U.not.ramp < U.cgg & U.cgg > 0.2), aes(label = name), size = 3, vjust = 1.9  ) +
  theme(legend.position="none") +
  scale_x_continuous(labels=prettyNum) + scale_y_continuous(labels=prettyNum) + scale_fill_manual(values =c("gray","#0899FF"))

pdf(paste0(PATH,"part3/fig_results_cgg.ramp.pdf"), width = 8, height = 2.5, useDingbats=F)
grid.draw( cbind_gtable_max( ggplotGrob(g.cgg1), ggplotGrob(g.cgg2) ) )
dev.off()


g1 <- ggplot( ramp.relative.importance, aes(x=rs2.out.ramp,  y = rs2.ramp )) + 
  #  geom_point(data = subset(CAA.relative.importance.data, ttg.rich == F), pch=21, fill="gray70") + 
  geom_point(data = subset(ramp.relative.importance), aes(fill =  rs2.out.ramp < rs2.ramp), pch=21) + 
  geom_abline(intercept=0, slope=1, lty=3) + coord_equal() +
  #geom_text( data = subset(ramp.relative.importance, rs2.out.ramp < rs2.ramp), aes(label = name), size = 3, vjust = 1.9  ) +
  labs( x = expression(r[s]^2*"(outside ramp)"), y = expression(r[s]^2*"(codon ramp)") ) +
  theme(legend.position="none") +
  scale_x_continuous(labels=prettyNum) + scale_y_continuous(labels=prettyNum) + scale_fill_manual(values =c("gray","#0899FF"))


g2 <- ggplot( ramp.relative.importance, aes(x=U.out.ramp,  y = U.ramp )) + 
  #  geom_point(data = subset(CAA.relative.importance.data, ttg.rich == F), pch=21, fill="gray70") + 
  geom_point(data = subset(ramp.relative.importance), aes(fill =  U.out.ramp < U.ramp), pch=21) + 
  geom_abline(intercept=0, slope=1, lty=3) + coord_equal() +
  geom_text( data = subset(ramp.relative.importance, U.out.ramp < U.ramp & U.ramp > 0.4), aes(label = name), size = 3, vjust = 1.9  ) +
  theme(legend.position="none") +
  scale_x_continuous(labels=prettyNum) + scale_y_continuous(labels=prettyNum) + scale_fill_manual(values =c("gray","#0899FF"))


require(gridExtra)
pdf( paste0(PATH,"part3/fig_results_ramp.pdf"), width = 10, height=9, useDingbats=F )
grid.draw( cbind_gtable_max( ggplotGrob(g1), ggplotGrob(g2) ) )
dev.off()

TFS.genes.diauxic_ox <- unique(as.character(subset(SMoPT.data, faster.translation == T & experiment %in% c("diauxic","ox") )$name))
faster.initiation.diauxic_ox <- unique(as.character(subset(SMoPT.data, ratio_initiation < 1 & experiment %in% c("diauxic","ox") )$name))
ramp.relative.importance$TFS <- ramp.relative.importance$name %in% TFS.genes.diauxic_ox
ramp.relative.importance$i   <- ramp.relative.importance$name %in% faster.initiation.diauxic_ox
ramp.relative.importance$ramp.driven <- with(ramp.relative.importance, rs2.ramp > rs2.out.ramp )
ramp.driven.genes <- ramp.relative.importance$name[ with(ramp.relative.importance, rs2.ramp > rs2.out.ramp) ]

SMoPT.data$ramp.driven  <- SMoPT.data$name %in% ramp.driven.genes

ggplot(data = subset(study.translation.kinetics, name %in% as.character(ramp.relative.importance$name) & variation_initiation_frequency < 4),
                     aes(x= name %in% ramp.driven.genes, y = log2(1/ratio_initiation)  ) ) + geom_boxplot(notch=T) + facet_wrap( ~ label)

ggplot(data = subset(SMoPT.data, name %in% as.character(ramp.relative.importance$name) & time == 0),
       aes(x= name %in% ramp.driven.genes, y = tAI  ) ) + geom_boxplot(notch=T) + facet_wrap( ~ label)

ggplot(data = subset(SMoPT.data, name %in% as.character(ramp.relative.importance$name) & time == 0),
       aes(x= name %in% ramp.driven.genes, y = log10(IniProb)  ) ) + geom_boxplot(notch=T) + facet_wrap( ~ label)


wilcox.test.batch(x = subset(SMoPT.data, name %in% as.character(ramp.relative.importance$name) & time == 0), grouping.vars = "ramp.driven", value = "IniProb" )
wilcox.test.batch(x = subset(SMoPT.data, name %in% as.character(ramp.relative.importance$name) & time == 0), grouping.vars = "ramp.driven", value = "tAI" )
wilcox.test.batch(x = subset(SMoPT.data, name %in% as.character(ramp.relative.importance$name) & time == 0), grouping.vars = "ramp.driven", value = "AUGCAI" )

head(ramp.relative.importance)
 head(study.translation.kinetics)
table(ramp.relative.importance[,c("TFS","ramp.driven")] )
table(ramp.relative.importance[,c("i","ramp.driven")] )

ggplot(data = subset(study.translation.kinetics, experiment %in% c("diauxic","ox") ),
       aes(x= name %in% ramp.driven.genes, y = log2(1/ratio_initiation)  ) ) + geom_boxplot(notch=T) + facet_wrap( ~ label)

enrichment.ramp.driven.and.TFS.genes <- with(ramp.relative.importance, chisq.test( TFS, ramp.driven ))
enrichment.ramp.driven.and.TFS.genes$effect.size <-  sqrt( enrichment.ramp.driven.and.TFS.genes$statistic  / sum(enrichment.ramp.driven.and.TFS.genes$observed) )

enrichment.ramp.driven.and.faster.initiation <- with(ramp.relative.importance, chisq.test( i, ramp.driven ))
enrichment.ramp.driven.and.faster.initiation$effect.size <- sqrt( enrichment.ramp.driven.and.faster.initiation$statistic  /  sum(enrichment.ramp.driven.and.faster.initiation$observed) )



# PLOT 1x
ggplot(data=cgg.table, aes( x = total.tRNAab, y= 1/av.elongation_time ) ) + 
  #stat_smooth( method="gam", formula = (y) ~ x*(1-exp(x)), se = F, color = "gray30")
  stat_smooth(method = "lm", level=0, color = "gray30") +
  geom_point(pch=21,fill="white", size=3) + 
  labs( y = "codon elongation speed (codon/s)", x="anticodon supply (Arg-CCG)" ) +
  scale_x_continuous(labels= scientific_10)


# PLOT 2
g1 <- ggplot( data = subset(CCG.influence, name %in% c("MRS2") ), aes( x = stAI, y = elongation.speed ) ) + 
  stat_smooth(method="lm", color="gray30") + 
  geom_line( lty=3 ) + 
  geom_point( pch = 21, fill = "white", size =3 ) +
  labs( title = "MRS2", y = "elongation speed (aa/s)") +
  scale_x_continuous(labels=prettyNum) +  scale_y_continuous(labels=prettyNum)

g2 <- ggplot( data = subset(CCG.influence, name %in% c("MRS2") ), aes( x = stAI, y = initiation_freq ) ) + 
  stat_smooth(method="lm", color="gray30") + 
  geom_line( lty=3 ) + 
  geom_point( pch = 21, fill = "white", size =3 ) +
  labs( title = "MRS2", y = "initiation frequency (min-1)") +
  scale_x_continuous(labels=prettyNum) +  scale_y_continuous(labels=prettyNum)

g3 <- ggplot( data = subset(CCG.influence, name %in% c("MRS2") ), aes( x = cgg.impact.ramp, y = elongation.speed ) ) + 
  stat_smooth(method="lm", color="gray30") + 
  geom_line( lty=3 ) + 
  geom_point( pch = 21, fill = "white", size =3 ) +
  labs( title = "MRS2", y = "elongation speed (aa/s)") +
  scale_x_continuous(labels=prettyNum) +  scale_y_continuous(labels=prettyNum)

g4 <- ggplot( data = subset(CCG.influence, name %in% c("MRS2") ), aes( x = cgg.impact.ramp, y = initiation_freq ) ) + 
  stat_smooth(method="lm", color="gray30") + 
  geom_line( lty=3 ) + 
  geom_point( pch = 21, fill = "white", size =3 ) +
  labs( title = "MRS2", y = "initiation frequency (min-1)") +
  scale_x_continuous(labels=prettyNum) +  scale_y_continuous(labels=prettyNum)

require(gridExtra)
grid.draw( rbind_gtable_max( cbind_gtable_max( ggplotGrob(g1), ggplotGrob(g2) ), cbind_gtable_max( ggplotGrob(g3), ggplotGrob(g4) ) ) )



# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# mulitplot figure

d <- with(CCG.relative.importance.data, data.frame( name = name, cgg.in.ramp = cgg.in.ramp,
                                                    r.cgg=r.cgg, q.cgg = quantile.it(U.cgg/U.not.ramp), ratio.U = U.cgg/U.not.ramp, freq.cgg = freq.in_CDS, U.stAI.cgg = U.stAI.cgg, U.cgg = U.cgg   ) )

require(scales)

g_reg1 <- ggplot( d, aes(x=freq.cgg, y=r.cgg) ) + geom_point(aes(color=cgg.in.ramp)) + 
  geom_text(data=subset(d, r.cgg > 0.6),aes(label=name), face="bold", color = "skyblue4", size = 3, vjust=1.5) +
  scale_x_continuous(label=percent) + scale_y_continuous(label=percent) + labs( y = expression(r[cgg]), x = expression(f[cgg])) +
  scale_color_manual(values=c("gray25","#0899FF")) + theme(legend.position="none")

g_reg2 <- ggplot( d, aes(x=freq.cgg, y=U.cgg) ) + geom_point(aes(color=cgg.in.ramp)) + geom_text(data=subset(d, U.cgg > 0.1),aes(label=name), face="bold", 
                                                                                              color = "skyblue4", size = 3, vjust=1.5) +
  scale_x_continuous(label=percent) + scale_y_continuous(label=percent) + labs( y = "Unique contribution of cgg", x = expression(f[cgg]) ) +
  scale_color_manual(values=c("gray25","#0899FF"))+ theme(legend.position="none")




require(gridExtra)



