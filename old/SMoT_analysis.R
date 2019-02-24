setwd("~/Documents/MRC16")
source("scripts/SMoT.R") # commands

#Some references on stress response in yeast:
# Clare E. Simpson and Mark P. Ashe. Adaptation to stress in yeast: to translate or not? Biochem. Society Trans. 2012, 40:794-799
# Audrey P. Gasch et al. Genomic expression programs in the response of yeast cells to environmental changes. Mol. Biol. Cell 2000, 11:4241-4257
# Eulàlia de Nadal et al. Controlling gene expression in response to stress. Nature Rev. Genetics 2011, 12:833-845
# Audrey P. Gasch. The environmental stress response: a common yeast response to diverse environmental stresses (Book Chapter)
# Haruo Saito and Francesc Posas. Response to hyperosmotic stress. Gentics 2012, 192:289-318
# Kevin A. Morano et al. The response to heat shock and oxidative stress in Saccharomyces cerevisiae. Genetics 2011.

# copy2clipboard(unique(arrange(subset(master.table.020, select=c(Gene, rand_mRNA)),Gene))$rand_mRNA)


# DATA ######
require(data.table)
go.slim <- fread("data/annotations/go_slim_mapping.txt",header=F,sep="\t")
setnames(go.slim, colnames(go.slim), c("ORF","Name", "sgdid", "GO.type","GO.description","GO.ID","molecule"))
go.slim.ORFs <- subset(go.slim, ORF %in% go.slim$ORF[ grep( "ORF", go.slim$molecule ) ] )
go.slim.ORFs$Name <- fNAME(go.slim.ORFs$ORF)

go.slim.ORFs <- arrange(subset(go.slim.ORFs, GO.description !="other", select=c(ORF,Name,GO.type,GO.description, GO.ID)), Name, GO.type)
setnames(x = go.slim.ORFs, "Name", "name")
write.table(go.slim.ORFs, file = "data/annotations/GO.slim_ORFs.txt", quote=F, sep="\t", row.names=F)

unique(subset(go.slim.ORFs, GO.type == "P")$GO.description)

go.slim.interest <- read.table("data/annotations/go_slim_of_interest.txt",header=1, sep="\t", stringsAsFactors = F)

GOs <- as.character(go.slim.interest$GO.description)
GO.slim <-   melt( lapply( setNames(GOs,GOs), function(x){ subset( go.slim.ORFs, GO.description %in% x )$ORF }  ) )
colnames(GO.slim) <- c("ORF","GO.description")
require(StingRay)
GO.slim$name <- fNAME(GO.slim$ORF)

write.table(GO.slim, "data/annotations/GO.slim_SMoT.txt",sep="\t",row.names=F)

# build master tables: compile data on outcomes of the stochastic simulation of translation in different conditions and time points, together with variables from df.PLSPM
master.table.020 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
master.table.120 <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
master.table.controls <- compile.master.table( path = "data/SMoPT/stochastic.simulation/unformated/controls/")

write.table( master.table.020, file = "data/SMoPT/stochastic.simulation/master.table.020min.txt", quote=F, row.names=F, sep="\t")
write.table( master.table.120, file = "data/SMoPT/stochastic.simulation/master.table.120min.txt", quote=F, row.names=F, sep="\t")
write.table( master.table.controls, file = "data/SMoPT/stochastic.simulation/master.table.controls.txt", quote=F, row.names=F, sep="\t")

# build codon tables: check the elongation time of individual codons in different conditions and time points
codon.table.020 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/20 minutes/")
codon.table.120 <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/120 minutes/")
codon.table.controls <- compile.codon.elongation.matrix( path = "data/SMoPT/stochastic.simulation/unformated/controls/")




  # summary table enable time-course analysis of codon elongation dynamics upon stress
  codon.table.adaptation <- melt( data.frame( codon = codon.table.020$codon, 
                                              experiment = codon.table.020$experiment, 
                                              tAI = codon.table.020$tAI, 
                                              av.time.normal=codon.table.020$av.elongation_time.normal,
                                              't0'=    1, 
                                              't20' = codon.table.020$ratio_elongation, 
                                              't120'=codon.table.120$ratio_elongation ), 
                                  id.vars = c("codon","experiment","tAI","av.time.normal") )
  colnames(codon.table.adaptation)[5] <- "time"
  codon.table.adaptation$high_tAI  <- ifelse( codon.table.adaptation$tAI > 0.2, T, F )
  codon.table.adaptation$av.time <- codon.table.adaptation$av.time.normal * codon.table.adaptation$value
  codon.table.adaptation$experiment <- factor(codon.table.adaptation$experiment, levels=levels(codon.table.adaptation$experiment), labels=c("diauxic","osmolarity","oxidative","temperature"))
  codon.table.adaptation$time <- factor( codon.table.adaptation$time, labels=c(0,20,120) )

write.table( codon.table.020, file = "data/SMoPT/stochastic.simulation/codon.table.020min.txt", quote=F, row.names=F, sep="\t")
write.table( codon.table.120, file = "data/SMoPT/stochastic.simulation/codon.table.120min.txt", quote=F, row.names=F, sep="\t")
write.table( codon.table.controls, file = "data/SMoPT/stochastic.simulation/codon.table.controls.txt", quote=F, row.names=F, sep="\t")
write.table( codon.table.adaptation, file = "data/SMoPT/stochastic.simulation/codon.table.adaptation.txt", quote=F, row.names=F, sep="\t")


# Venn diagrams showing the overlap between gene that get translated faster in different stress conditions and time points
venn_020 <- extract.venn.data(dataset = master.table.020)
venn_120 <- extract.venn.data(dataset = master.table.120)

plot(venn_020$venn)
plot(venn_120$venn)



# compare translation dynamics in short and long term
keys <- c("experiment","Name")
var.of.interest <- c("ORF","time","Name","av.initiation_time","av.elongation_time","ratio_initiation","ratio_elongation","ratio_total.time","est.mRNA_abundance","log2.mRNA_abundance","global.protein_synthesis.rate")
not.desired <- c("Gene","ORF")

df.master.dynamics <- merge( subset(master.table.020, select= c(keys,var.of.interest)), 
                        subset(master.table.120,select= setdiff(colnames(master.table.120),not.desired) ),
                        by= keys  )

# beware that there are 331 genes with NAs in either case (because they had +Inf in av.initiation time I reckon)
df.master.dynamics <- subset(df.master.dynamics, !(is.na(ratio_total.time.x) | is.na(ratio_total.time.y)))
df.master.dynamics$ratio_mRNA.x  <- df.master.dynamics$est.mRNA_abundance.x / df.master.dynamics$rand_mRNA
df.master.dynamics$ratio_mRNA.y  <- df.master.dynamics$est.mRNA_abundance.y / df.master.dynamics$rand_mRNA
df.master.dynamics$transcription.dynamics <- paste( ifelse(df.master.dynamics$log2.mRNA_abundance.x > 0, "up", "down") ,       ifelse(df.master.dynamics$log2.mRNA_abundance.y > 0, "up", "down"), sep="->"  )
df.master.dynamics$translation.dynamics <- paste( ifelse(df.master.dynamics$ratio_total.time.x > 1, "slow", "fast") , ifelse(df.master.dynamics$ratio_total.time.y > 1, "slow", "fast"), sep="->"  )





# df.master.dynamics$fast_20  <-  ifelse(df.master.dynamics$ratio_total.time.x > 1, 0, 1) 
# df.master.dynamics$fast_120 <- ifelse(df.master.dynamics$ratio_total.time.y > 1,  0, 1)
write.table(x = df.master.dynamics, file="data/SMoPT/stochastic.simulation/master.table.dynamics.txt",sep="\t", quote=F, row.names=F)



             


map.dynamics(data = df.master.dynamics, GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="consistency")$GO.description ))
map.dynamics(data = df.master.dynamics, GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="translational response")$GO.description ))
map.dynamics(data = df.master.dynamics, GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="metabolic response")$GO.description ))
map.dynamics(data = df.master.dynamics, GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="cell cycle")$GO.description ))

map.dynamics(data = df.master.dynamics, what = "transcription.dynamics", GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="consistency")$GO.description ))
map.dynamics(data = df.master.dynamics, what = "transcription.dynamics", GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="translational response")$GO.description ))
map.dynamics(data = df.master.dynamics, what = "transcription.dynamics", GO.slim = subset(GO.slim, GO.description %in% subset(go.slim.interest, motivation=="cell cycle")$GO.description ))


data.dynamics <- melt( table(df.master.dynamics[,c("transcription.dynamics","translation.dynamics","experiment")])   )




ggplot( data = data.dynamics, aes(x=transcription.dynamics, y=translation.dynamics, fill=value) ) + 
  geom_tile( data= subset(data.dynamics, transcription.dynamics %in% c("up->down","up->up") & translation.dynamics %in% c("fast->slow","fast->fast") ), fill = "orange1") + 
  geom_tile( data= subset(data.dynamics, transcription.dynamics == "up->down" & translation.dynamics == "fast->slow"), fill = "orange3") + 
  geom_tile( width=0.65, height=0.65) + 
  geom_text( aes(label=value), color="white" ) +
  facet_wrap( ~ experiment) +
  coord_equal() + theme_bw() + labs(fill="n")


# dynamics.venn.diagrams <- dlply( subset(df.master.dynamics, select=c(experiment, Name, fast_20, fast_120)) , .(experiment), function(x){ y <- as.matrix(x[,3:4]); rownames(y) <- x[,2]; venneuler(y) } )
# plot(dynamics.venn.diagrams$diauxic)

head(df.master.dynamics)

transcription.dynamics.genes <- dlply( df.master.dynamics, .(experiment), function(x) dlply( x, .(transcription.dynamics), function(y) as.character( y$ORF ) ))
translation.dynamics.genes   <- dlply( df.master.dynamics, .(experiment), function(x) dlply( x, .(translation.dynamics), function(y) as.character( y$ORF ) ))



# Get the up-regulated genes in Gasch

load("data/Rdata/gash.all.Rda")
colnames(gash.all)

"Heat_Shock_20min_hs-1","Heat_Shock_80min_hs-1"
"1M_sorbitol_15_min","1M_sorbitol_120_min"



# Gene ontologies to check
            # 
            # copy2clipboard( Reduce(union, lapply(1:4, function(x){ dynamics.genes[[x]]$'slow->fast' } ) ) )
            # # slow then fast
            # reproduction
            # nuclear division
            # meoitic nuclear division
            # cell cycle
            # cell division
            # tRNA gene clustering
            # rDNA condensation
            # cellular response to stimulus
            # response to division
            # 
            # 
            # copy2clipboard( Reduce(union, lapply(1:4, function(x){ dynamics.genes[[x]]$'fast->fast' } ) ) )
            # # fast then fast
            # response to DNA damage stimulus
            # DNA repair
            # cell cycle
            # "GO:0007049"
            # 
            # DNA metabolic process
            # cellular aromatic compound metabolic proecess
            # cellular response to stimulus
            # cellular response to stress
            # response to stimulus
            # 
            # copy2clipboard( Reduce(union, lapply(1:3, function(x){ dynamics.genes[[x]]$'slow->fast' } ) ) )
            # # fast then slow, except temperature stress
            # rDNA condensation
            # "GO:0070550"
            # cell cycle
            # tRNA gene clustering
            # "GO:0070058"
            # 



ggplot( data = master.table.020, aes(x=global.protein_synthesis.rate, y=ratio_initiation) ) + 
  geom_vline( xintercept = 1) + geom_hline( yintercept= 1) +
  geom_point() +
  facet_wrap( ~ experiment ) + scale_x_log10() + scale_y_log10() 
 

(1500/60)/mean(subset(master.table.020, experiment=="temp")$av.elongation_time, na.rm=T) * 2e5




automate.chi.squared( data = subset(data, experiment == "ox"), expression = 'transcription.dynamics == "up->down" ', keys = as.character(go.slim.interest$GO.description) )



chisquare.tests <- list()

chisquare.tests$global.protein_synthesis.020 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression="global.protein_synthesis.rate.x > 1"  ) )
chisquare.tests$global.protein_synthesis.120 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression="global.protein_synthesis.rate.y > 1"   ) )

chisquare.tests$fast.translation.020 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression = "ratio_total.time.x < 1"  ) )
chisquare.tests$fast.translation.120 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression = "ratio_total.time.y < 1"  ) )

chisquare.tests$increased.transcription.020 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression = "log2.mRNA_abundance.x > 0"  ) )
chisquare.tests$increased.transcription.120 <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression = "log2.mRNA_abundance.y > 0"  ) )

chisquare.tests$transient.guys <- ldply( dlply( master.table.dynamics, .(experiment), automate.chi.squared, keys = as.character(go.slim.interest$GO.description), expression = 'transcription.dynamics == "up->down" & translation.dynamics == "fast->slow" '  ) )

save(chisquare.tests, file="results/SMoPT/chisquare_tests.Rda")

# Enriched
View( data.frame(subset( chisquare.tests$fast.translation.020, log2_xy.TT > 0 & p.value < 0.025 ), row.names=NULL) )
# Depleted
View( data.frame(subset( chisquare.tests$fast.translation.020, log2_xy.TT < 0 & p.value < 0.025 ), row.names=NULL) )

# Enriched
View( data.frame(subset( chisquare.tests$global.protein_synthesis.020, log2_xy.TT > 0 & p.value < 0.025 ), row.names=NULL) )
# Depleted
View( data.frame(subset( chisquare.tests$global.protein_synthesis.020, log2_xy.TT < 0 & p.value < 0.025 ), row.names=NULL) )

# Enriched
View( data.frame(subset( chisquare.tests$increased.transcription.020, log2_xy.TT > 0 & p.value < 0.025 ), row.names=NULL) )
# Depleted
View( data.frame(subset( chisquare.tests$increased.transcription.020, log2_xy.TT < 0 & p.value < 0.025 ), row.names=NULL) )

# No-thing for translation :(
subset( chisquare.tests_fast.translation, log2_xy.TT > 0 & p.value < 0.025 )
subset( chisquare.tests_fast.translation, log2_xy.TT < 0 & p.value < 0.025 )



# --------------------------------------------------------------------------- #
# Analysis

translate.better.under.stress.2  <- subset(master.table.2, faster.translation == T)$Name
copy2clipboard(translate.better.under.stress.2)
# cell cycle [GO:0007049]  0.027211  47

ggplot( data = master.table.2, aes(x=expected.translation_time) ) + geom_histogram(binwidth=1) + facet_grid( experiment ~ . ) + labs( x = "exp. translation time (min)")

ggplot( data = subset(master.table.2, experiment=="genome_diauxic_120"), aes(x=expected.translation_time.normal, y = protein_abundance) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "exp. translation time (min)", y = "protein abundance (PAXDB)")


ggplot( data = subset(master.table.2, experiment=="genome_diauxic_120"), aes(x=expected.translation_time.normal, y = ribosome_density_YPD) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "exp. translation time (min)", y = "ribosome density (YPD)")


ggplot( data = subset(master.table.2, experiment=="genome_diauxic_120"), aes(x=expected.translation_time.normal, y = length_CDS/3) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "exp. translation time (min)", y = "protein length")

ggplot( data = subset(master.table.2, experiment=="genome_diauxic_120"), aes(x=length_CDS/3, y = ribosome_density_YPD) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "protein length", y = "ribosome density (YPD)")


ggplot( data = subset(master.table.2, experiment=="genome_diauxic_120"), aes(x=length_CDS/3, y = residual_cost) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() +  
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "protein length", y = "metabolic cost")


ggplot( data = master.table.2, aes(x=expected.translation_time, y = residual_cost) ) + 
  geom_point() + stat_smooth(method="lm") +
  scale_x_log10() +  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2)) +
  labs( x = "expected.translation_time", y = "metabolic cost")






# - 1 -  Dependency on gene length

# One would expect that the initiation time doesn't depend on gene length
ggplot( data = master.table.2, aes(x = length_CDS, y=av.initiation_time) ) + 
  geom_point() + 
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

# In contrast, the elongation time is linearly dependent on gene length
ggplot( data = master.table.2, aes(x = length_CDS, y=av.elongation_time) ) + 
  geom_point() +
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

# The sum of initiation time and elongation time also (marginally) depend on gene length
ggplot( data = master.table.2, aes(x = length_CDS, y=expected.translation_time) ) + 
  geom_point() + stat_smooth(method="lm") +
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

# By taking the residuals of the expected time ~ gene length regression, the length-dependency is controlled for
ggplot( data = master.table.2, aes(x=expected.translation_time, y = residual.translation.time) ) + 
  geom_point() + stat_smooth(method="lm") +
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

# Inded this has now vanished
ggplot( data = master.table.2, aes(x=length_CDS, y = residual.translation.time) ) + 
  geom_point() + stat_smooth(method="lm") +
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))


ggplot( data = master.table.2, aes(x=expected.translation_time.normal, y = residual.translation.time) ) + 
  geom_point() + stat_smooth(method="lm") +
  geom_hline(yintercept=0)+
  scale_x_log10() + facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))


# check dependency of gene lenght and delta elongation

ggplot( data = master.table.2, aes(x=length_CDS, y = delta_elongation) ) + 
  geom_point() + stat_smooth(method="lm") +
  geom_hline(yintercept=0)+
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

ggplot( data = master.table.2, aes(x=length_CDS, y = residual.delta_elongation) ) + 
  geom_point() + stat_smooth(method="lm") +
  geom_hline(yintercept=0)+
  facet_grid( experiment ~ . ) + 
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2))

# now we can talk 
ggplot( data = master.table.2, aes(x=delta_initiation, y = residual.delta_elongation, color = is.outliers(delta_initiation, h = 1.5) & is.outliers(residual.delta_elongation, h = 1.5) ) ) + 
  geom_point() +
  geom_hline(yintercept=0) +   geom_vline(xintercept=0)+
  facet_grid( experiment ~ . ) +
  scale_color_manual( values =  c("gray","orange2")) +
  labs( colour = "outliers" ) +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")


# now we can talk 
ggplot( data = master.table.2, aes(x=delta_initiation, y = residual.delta_elongation, color = (abs(zscore.delta_initiation) > 1.65) | (abs(zscore.residual.delta_elongation) > 1.65) ) ) + 
  geom_point() +
  geom_hline(yintercept=0) +   geom_vline(xintercept=0)+
  facet_grid( experiment ~ . ) +
  scale_color_manual( values =  c("gray","orange2")) +
  labs( colour = "outliers" ) +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")

translate.better.under.stress <- subset( master.table.2, ( (abs(zscore.delta_initiation) > 1.65) | (abs(zscore.residual.delta_elongation) > 1.65) ) & residual.delta_elongation < 0 & delta_initiation < 0 ) 


dlply( translate.better.under.stress, .(experiment) )

genes.of.interest <- unique(translate.better.under.stress[,c("ORF","Name","experiment")] )
genes.of.interest$Name <- as.character(genes.of.interest$Name)
table(genes.of.interest$Name)

copy2clipboard( unique(genes.of.interest$Name) )

ggplot( data = master.table.2, aes( x = Name %in% genes.of.interest$Name, y = tRNA_adaptation_index )) + geom_boxplot(notch=T) +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")

ggplot( data = master.table.2, aes( x = Name %in% genes.of.interest$Name, y = IniProb )) + geom_boxplot(notch=T) + scale_y_log10() +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")

ggplot( data = master.table.2, aes( x = Name %in% genes.of.interest$Name, y = length_CDS )) + geom_boxplot(notch=T) + scale_y_log10() +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")

ggplot( data = master.table.2, aes( x = Name %in% genes.of.interest$Name, y = av.elongation_time.normal )) + geom_boxplot(notch=T) + scale_y_log10() +
  theme_bw() + theme(axis.title.x=element_text(hjust=0.5, vjust=0.2), axis.title.y=element_text(hjust=0.5, vjust=0.2), legend.position= "bottom")




  ggplot( subset( master.table.020, !is.na(ratio_initiation) ), aes( x = est.mRNA_abundance/rand_mRNA , y = delta_initiation )  ) + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  scale_x_log10() + 
  geom_point() +
  facet_grid( . ~ experiment) + theme_bw()


  ggplot( subset( master.table.120, !is.na(ratio_initiation) ) , aes( x = est.mRNA_abundance/rand_mRNA , y = ratio_initiation )  ) + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  scale_x_log10() + scale_y_log10() + 
  geom_point() +
  facet_grid( . ~ experiment) + theme_bw()

  ggplot( subset( master.table.020, !is.na(ratio_initiation) & av.initiation_time > 0 ), aes( x = est.mRNA_abundance , y = av.initiation_time )  ) + 
  geom_hline(yintercept=0, lty=2, colour = "gray50") + 
  scale_x_log10() +   scale_y_log10() + 
  geom_point() + stat_smooth(method="lm") +
  facet_grid( . ~ experiment) + theme_bw()
