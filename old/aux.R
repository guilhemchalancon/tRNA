ddply(marc.data, .(experiment, time), function(x){ cor(x$average, x$std, method = "pearson")})

ddply(marc.data, .(experiment, time), function(x){ summary(lm( std ~ average, data=x))$r.squared })

range.20 <- summary( 2^(subset(anticodon.master.table, time == 20)$log2.foldchange) )
range.60 <- summary( 2^(subset(anticodon.master.table, time == 60)$log2.foldchange) )
range.120 <- summary( 2^(subset(anticodon.master.table, time == 120)$log2.foldchange) )


max(range.20)/min(range.20)
max(range.60)/min(range.60)
max(range.120)/min(range.120)


DNR <- data.frame( ddply( marc.data, .(time), summarise, max = max(average), min = min(average), DNR = max(average)/min(average), dispersion = mean(std^2/average) ), 
            cor = sapply(COR, function(x){diag(x) <- NA; mean(x,na.rm=T);}))

ggsave( plot = ggplot(data=DNR, aes(x=cor, y=DNR, fill=time)) + geom_line() + geom_point(aes(fill=time), size=5, pch=21) + 
  geom_text(aes(label=paste(time,"min")) , hjust=-0.4, size = 4, angle=90) + ylim(50, 350) + 
  labs(x="Pearson's r", y="Dynamic range") + theme(legend.position="none"), filename="~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part1/dynamic.range.pdf", useDingbats=F, width=3.32, height=2.79) 


set  <- subset(anticodon.master.table, time == 0)
ordered.anticodons <- as.character( arrange(set, plyr::desc(total.tRNAab))$anticodon )
ordered.aa <- as.character( arrange(set, plyr::desc(total.tRNAab))$aa )

set$anticodon <- factor( set$anticodon, levels = ordered.anticodons ) 
set$label  <- factor( set$anticodon, levels = ordered.anticodons, labels = paste( ordered.aa, ordered.anticodons, sep="-") )

g1 <- ggplot(set, aes(x=factor(0), y=factor(0))) + geom_tile(aes(fill=total.tRNAab)) +   
  facet_wrap( ~ label) + 
  scale_fill_continuous(low = "#DEEBF7",high = "#3182BD") + coord_fixed(ratio = 0.3) +
  labs(x="",y="") + theme(legend.position="none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

set2  <- subset(anticodon.master.table, time == 20 & experiment == "osm")
set2$anticodon <- factor( set2$anticodon, levels = ordered.anticodons ) 
set2$label  <- factor( set2$anticodon, levels = ordered.anticodons, labels = paste( ordered.aa, ordered.anticodons, sep="-") )

g2 <- ggplot(set2, aes(x=factor(0), y=factor(0))) + geom_tile(aes(fill=total.tRNAab)) +   
  facet_wrap( ~ label) + 
  scale_fill_continuous(low = "#DEEBF7",high = "#3182BD") + coord_fixed(ratio = 0.3) +
  labs(x="",y="") + theme(legend.position="none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())


require(gridExtra)
pdf("~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part1/variations.tRNA_illustration.pdf",  width=33.2/3, height=27.9/4)
  grid.arrange(g1,g2, nrow=1)
dev.off()

rosetta.codons$w <- 0.9
rosetta.codons <- ddply(rosetta.codons, .(aa), mutate, w = w*length(w)/7 )
#rosetta.codons$w[df$TIME == 5] <- 0.9 * 3/4

requires(scales)
g <- ggplot(data=rosetta.codons, aes(x=codon, y=genomic.freq )) + geom_bar(aes(width = w, fill=CAI.O), stat="identity") + 
  facet_grid( ~ aa, scales="free_x") + 
  scale_fill_manual(values=c("skyblue","orange")) + 
  #scale_fill_manual(values=c("gray30","gray30")) + 
  labs(x="codons", y="genomic\frequency") +
  theme(legend.position="none",axis.text.x=element_text(angle=90, vjust=0.5)) + scale_y_continuous(labels=percent)
ggsave(plot=g, filename = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part1/genomic.frequency.codons.pdf", width=12.7, height=1.85)



ddply(anticodon.master.table, .(experiment, time), 
      function(x){ 42 - rank(x$total.tRNAab)[ which(x$anticodon == "CAA") ] } 
      )


ggplot( marc.data, aes( x= quantcut(average, q = seq(0,1,by=0.1) ), y = std/average )) + 
  geom_boxplot(notch=T) + 
  labs(x="average fold change",y="CV", title="Stationary CV in function of mean\n(expected from log-normal distribution)
       ")

### ------------------------------------------------------------ 21 April 2014 -------------------------------------------------- ###
setwd("~/Documents/MRC16")
data.gfp <- read.csv("data/GFP seq/data_gpf_heatmap.csv",header=1)
data.gfp <- ddply(data.gfp, .(optimised.stain), mutate, scaled.value = range01(scale(value)), ranged.value = (value)/max(value) ) 
data.gfp$condition <- factor(data.gfp$condition, levels= rev(c("oxidative","osmotic","temperature","diauxic")) )
data.gfp$optimised.stain <- factor(data.gfp$optimised.stain, levels= c("oxidative","osmotic","temperature","diauxic"), 
                                   labels =c("eGFP_ox","eGFP_osm","eGFP_temp","eGFP_diauxic") )

g <- ggplot(data=data.gfp, aes(x=optimised.stain,y=condition, fill=ranged.value)) + geom_tile(height=0.9) + 
  geom_text( aes(label=round(ranged.value,2)),color="white", size =5) +
  scale_fill_gradient(low = "gray90", high="gray30") +
  #scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2878B8", midpoint = 0.5) + 
  coord_equal() +
  labs(x="optimised strain", y="stress condition", fill="scaled\nfold\nchange") + theme(legend.position="none")

ggsave(plot=g, filename =  "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/fig_results_gfp.2.pdf", useDingbats=F,width=5.39 , height=3.76)
require(RColorBrewer)
rev(colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))[1:80]

g <- ggplot( data = subset(study.translation.kinetics, time == 120 ), aes(x=experiment, y= log2( global.protein_synthesis.rate) )) + geom_boxplot( notch = T, width=0.75 ) + geom_hline( yintercept= 0, lty=3)
ggsave(plot=g, filename =  "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/fig_results_Pi.pdf", useDingbats=F,width=6.28 , height=3.68  )


### MADAN PAPER MARC 

# g1 <- ggplot( data = subset(study.translation.kinetics, time == 20  & valid == T & variation_initiation_frequency < 4), aes(x= label , y= variation_initiation_frequency) ) + 
#   geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T) +  theme_minimal() + labs(x="",y= "difference in initiation frequency (%)") + scale_y_continuous(labels=percent)
# 
# g2 <- ggplot( data = subset(study.translation.kinetics, time == 20  & valid == T), aes(x= experiment , y= variation_speed) ) + 
#    geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T) + theme_minimal() + labs(x="",y= "difference in elongation speed (%)") + scale_y_continuous(labels=percent)
study.translation.kinetics$label <- factor(study.translation.kinetics$lab, levels= levels(study.translation.kinetics$label)[c(1,5,2,6,3,7,4,8)] )

g1a <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_initiation_frequency < 4), aes(x= label , y= variation_initiation_frequency) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) +  theme_minimal() + labs(x="",y= "difference in initiation frequency (%)") + scale_y_continuous(labels=percent) 

g2a <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 ), aes(x= label , y= variation_speed) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_minimal() + labs(x="",y= "difference in elongation speed (%)") + scale_y_continuous(labels=percent)


g3a <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 ), aes(x= label , y= log2(global.protein_synthesis.rate)) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_minimal() + labs(x="",y= "change in protein production (log2)") 


require(gridExtra)
#pdf("~/Dropbox/tRNA_paper/Paper_figure3B_boxplots_raw.pdf")
pdf("~/Desktop/tRNA_boxplots_a.pdf", height = 5, width = 8)
grid.draw( cbind_gtable_max( ggplotGrob(g1a), ggplotGrob(g2a) ) )
dev.off()


pdf("~/Desktop/tRNA_boxplots_protein.prod_a.pdf", height = 5, width = 4)
g3a
dev.off()



g1b <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_initiation_frequency < 4 & time !=120), aes(x= label , y= variation_initiation_frequency) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "difference in initiation frequency (%)") + scale_y_continuous(labels=percent) 

g2b <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time !=120), aes(x= label , y= variation_speed) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "difference in elongation speed (%)") + scale_y_continuous(labels=percent)

g3b <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time !=120), aes(x= label , y= log2(global.protein_synthesis.rate)) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "change in protein production (log2)") 


require(gridExtra)
#pdf("~/Dropbox/tRNA_paper/Paper_figure3B_boxplots_raw.pdf")
pdf("~/Desktop/figs_4_marc/tRNA_boxplots_b.pdf", height = 7.5, width = 9)
grid.draw( cbind_gtable_max( ggplotGrob(g1b), ggplotGrob(g2b) ) )
dev.off()


pdf("~/Desktop/figs_4_marc/tRNA_boxplots_protein.prod_b.pdf", height = 7.5, width = 6)
g3b
dev.off()

# MADAN PAPER MARC 27 MAY 2015

g1c <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_initiation_frequency < 4 & time !=120), aes(x= label , y= variation_initiation_frequency) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) + labs(x="",y= "difference in initiation frequency (%)") + 
  scale_y_continuous(labels=percent, limits=c(-1,1)) + theme_gul

g2c <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time !=120), aes(x= label , y= variation_speed) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) +  labs(x="",y= "difference in elongation speed (%)") + 
  scale_y_continuous(labels=percent, limits=c(-1,1)) + theme_gul

g3c <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time !=120), aes(x= label , y= log2(global.protein_synthesis.rate)) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) + labs(x="",y= "change in protein production (log2)") + ylim(-3,3) + theme_gul


require(gridExtra)
#pdf("~/Dropbox/tRNA_paper/Paper_figure3B_boxplots_raw.pdf")
pdf("~/Desktop/figs_4_marc/tRNA_boxplots_c.pdf", height = 7.5, width = 9)
grid.draw( cbind_gtable_max( ggplotGrob(g1c), ggplotGrob(g2c) ) )
dev.off()


pdf("~/Desktop/figs_4_marc/tRNA_boxplots_protein.prod_c.pdf", height = 7.5, width = 6)
g3c
dev.off()



# 28 May NOW for 120min

require(scales)

g1d <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_initiation_frequency < 4 & time ==120), aes(x= label , y= variation_initiation_frequency) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "difference in initiation frequency (%)") + scale_y_continuous(labels=percent) 

g2d <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time ==120), aes(x= label , y= variation_speed) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "difference in elongation speed (%)") + scale_y_continuous(labels=percent)

g3d <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time ==120), aes(x= label , y= log2(global.protein_synthesis.rate)) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 1) + theme_gul + labs(x="",y= "change in protein production (log2)") 


require(gridExtra)
#pdf("~/Dropbox/tRNA_paper/Paper_figure3B_boxplots_raw.pdf")
pdf("~/Desktop/figs_4_marc/tRNA_boxplots_d.pdf", height = 7.5, width = 9)
grid.draw( cbind_gtable_max( ggplotGrob(g1d), ggplotGrob(g2d) ) )
dev.off()


pdf("~/Desktop/figs_4_marc/tRNA_boxplots_protein.prod_d.pdf", height = 7.5, width = 6)
g3d
dev.off()

# MADAN PAPER MARC 27 MAY 2015

g1e <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_initiation_frequency < 4 & time ==120), aes(x= label , y= variation_initiation_frequency) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) + labs(x="",y= "difference in initiation frequency (%)") + 
  scale_y_continuous(labels=percent, limits=c(-1,1)) + theme_gul

g2e <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time ==120), aes(x= label , y= variation_speed) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) +  labs(x="",y= "difference in elongation speed (%)") + 
  scale_y_continuous(labels=percent, limits=c(-1,1)) + theme_gul

g3e <- ggplot( data = subset(study.translation.kinetics, valid == T & variation_speed < 1 & time ==120), aes(x= label , y= log2(global.protein_synthesis.rate)) ) + 
  geom_hline(yintercept=0, lty=3) + geom_boxplot(notch=T, outlier.size = 0) + labs(x="",y= "change in protein production (log2)") + ylim(-3,3) + theme_gul


require(gridExtra)
#pdf("~/Dropbox/tRNA_paper/Paper_figure3B_boxplots_raw.pdf")
pdf("~/Desktop/figs_4_marc/tRNA_boxplots_e.pdf", height = 7.5, width = 9)
grid.draw( cbind_gtable_max( ggplotGrob(g1e), ggplotGrob(g2e) ) )
dev.off()


pdf("~/Desktop/figs_4_marc/tRNA_boxplots_protein.prod_e.pdf", height = 7.5, width = 6)
g3e
dev.off()



# repeatability 
View(head(marc.data))

# percentage of measurements within 2SDs
1 - length(with(marc.data, which( abs(average - mean(average)) > 2* sd(average) ) )  )/nrow(marc.data)

marc.data$sem <- with( marc.data, sqrt( std^2/3  ))

hist(log(marc.data$sem) )

####### table anticodons
require(stargazer)

stargazer( arrange( subset(rosetta.anticodons, select=c(anticodon, aa, codon, tGCN, CAI.O)), anticodon ), summary = F,
           out = "~/Dropbox/PhD/thesis/TeX/items/tables/appendix/table_anticodons.draft.tex"
           )


arrange( subset(rosetta.codons, select=c(codon, anticodon, aa)), codon )

# third appendix table
codon.master.table      <- read.table("results/master tables/codon.master.table.txt", header=1) # especially relative adaptiveness and anticodon deman

table <-  dcast(subset(codon.master.table, select = c(codon, anticodon, experiment,time, w.2)), codon + anticodon ~ experiment + time, value.var = "w.2" )

stargazer(table[,c("codon","anticodon","normal_0","diauxic_20","diauxic_120","ox_20","ox_120","osm_20","osm_120","temp_20","temp_120")], 
           summary = F,digits = 2,
           out  "~/Dropbox/PhD/thesis/TeX/items/tables/appendix/table_adaptiveness.draft.tex"
)

lm <- with(codon.master.table, cor.test( delta_elongation, delta_w, method="spearman"  ) )
summary(lm)


#### JULY 2015
load("~/Documents/MRC09/data/PLS-PM/df.PLSPM_bivariate.outliers.removed.Rda")
df <- subset(df.PLSPM.nTnS.o,select=c(gene, PARS_average_5UTR, PARS_average_3UTR, PARS_average_CDS, length_3UTR, length_5UTR,length_polyA_tail, PARS_mad_3UTR, PARS_mad_5UTR, structuredness, EFE_20up23down, freq_U.last_30nt ))
data_wc <- merge( gene.master.table, df, by.x="name", by.y="gene", all.x=T)

tAI.TFS <- wilcox.test.batch( x = subset(gene.master.table, stress!="normal_0" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )


wilcox.test.batch( x = subset(data_wc, stress=="ox_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="diauxic_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="diauxic_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="osm_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="osm_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="temp_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )
wilcox.test.batch( x = subset(data_wc, stress=="temp_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "tAI" )

wilcox.test.batch( x = subset(data_wc, stress=="ox_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = "length_polyA_tail" )

to.try <- c( "PARS_average_5UTR", "PARS_average_3UTR", "PARS_average_CDS", "length_3UTR", "length_5UTR","length_polyA_tail", "PARS_mad_3UTR", "PARS_mad_5UTR", "structuredness", "EFE_20up23down", "freq_U.last_30nt")

painfully_laborious <- list(
  diauxic_20 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="diauxic_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = x ) )
  }) ),
  diauxic_120 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="diauxic_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = x ) )
  }) ) ,
  ox_20 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="ox_20" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = x ) )
  }) ),
  ox_120 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="ox_120" & !is.na(faster.translation) ), grouping.vars = c("stress","faster.translation"), value = x ) )
  }) ),
  osm_20 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="osm_20" & !is.na(faster.translation) ), grouping.vars = c("faster.translation"), value = x ) )
  }) ),
  osm_120 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="osm_120" & !is.na(faster.translation) ), grouping.vars = c("faster.translation"), value = x ) )
  }) ),
  temp_20 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="temp_20" & !is.na(faster.translation) ), grouping.vars = c("faster.translation"), value = x ) )
  }) ),
  temp_120 = do.call(rbind, lapply( to.try, function(x){ 
    data.frame( var = x,  wilcox.test.batch( x = subset(data_wc, stress=="temp_120" & !is.na(faster.translation) ), grouping.vars = c("faster.translation"), value = x ) )
  }) )
)


TFS_extra <- arrange( ldply(painfully_laborious), var, abs(r.wendt) )
View(TFS_extra)
