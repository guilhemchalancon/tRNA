require(tidyr)
setwd("~/Documents/MRC16/")
s = c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68)
names(s) <- c("A:t","G:c","U:a","C:g","G:t","I:c","I:a","U:g")


load("data/Lister/orf_all.codon.usage_eff.Rda")
load("data/Lister/orf_all.codon.usage_freq.Rda")
load("data/Lister/orf_all.codon.usage_rscu.Rda")
EFF.transcriptome  <- ddply( arrange(melt(uco.eff),  gene), .(variable), summarise, mean.eff = mean(value,na.rm=T), SD.eff = sd(value,na.rm=T) )
FREQ.transcriptome <- ddply( arrange(melt(uco.freq), gene), .(variable), summarise, mean.freq = mean(value,na.rm=T), SD.freq = sd(value,na.rm=T) )
RSCU.transcriptome <- ddply( arrange(melt(uco.rscu), gene), .(variable), summarise, mean.RSCU = mean(value,na.rm=T), SD.RSCU = sd(value,na.rm=T) )
colnames(EFF.transcriptome)[1] <- c("codon")
colnames(FREQ.transcriptome)[1] <- c("codon")
colnames(RSCU.transcriptome)[1] <- c("codon")

selective.constrains <- data.frame( pos.1_anticodon = c("A","G","U","C","G","A","A","U"), pos.3_codon = c("t","c","a","g","t","c","a","g"), s = s )
stress.codon.adaptiveness <- read.table("results/stress tAI/stress.codon.adaptiveness.txt",header=1)

data <- subset(rosetta.codons, select=c(anticodon,codon,pos.1_anticodon, pos.3_codon, CAI, tAI, aa))
data <- merge(data, selective.constrains, by=c("pos.1_anticodon", "pos.3_codon") )
data <- merge(data, subset(rosetta.anticodons, select=c(anticodon, tGCN, relative.availability.normal, nTE.O)))
data <- merge(data, subset(stress.codon.adaptiveness, time==0, select=c(codon, w, w_FC)), by = "codon" )
data <- merge(data, EFF.transcriptome,by= "codon" )
data <- merge(data, FREQ.transcriptome,by= "codon" )
data <- merge(data, RSCU.transcriptome,by= "codon" )
data$wobbling <- ifelse(data$s>0, "wobbling","cognate")
data$wobbling.group  <- factor( data$s, labels=c("Watson-Crick","I:c","G:t","U:g","I:a") )

isoacceptor.tRNAs <- as.character(subset(ddply(data, .(anticodon), nrow), V1>1)$anticodon)
aa.with.syn.tRNAs <- as.character(subset(data, aa %in% as.character(subset(ddply(data, .(aa), nrow), V1>1)$aa))$anticodon)


g <- ggplot( subset(data, anticodon %in% isoacceptor.tRNAs), aes(x=wobbling.group, y = mean.RSCU) ) + 
  geom_hline(yintercept=1, lty=3) +
  geom_boxplot(notch=F) +
  geom_point(aes(fill= wobbling, size = tGCN ), pch=21) +
  scale_fill_manual( values=c("gray30","gray65") ) +
  labs(x="Anticodon-codon pairing\n(ordered by selective constraint s)", y="Codon usage (RSCU)", fill="tRNA")

ggsave(plot=g, filename = paste0(PATH,"part1/codon_usage_selective_constraint.pdf"), useDingbats=F, width=10.8, height=6.22)


ggplot( data=subset(data, anticodon %in% aa.with.syn.tRNAs), aes(x = mean.RSCU, y=(1-s)) ) + 
  #geom_boxplot(notch=F) +
  geom_vline(xintercept=1, lty=3)+
  geom_point(aes(shape= wobbling, color= aa), size = 3) +
  geom_text(aes(label=codon),vjust=1.5, size=2.5) +
  scale_color_manual( values = findPalette(aa.with.syn.tRNAs) ) + 
  #stat_smooth(method="lm", color = "gray30") +
  facet_wrap( aa ~ anticodon ) +
  scale_y_continuous(labels=prettyNum) +
labs(x="Codon usage (RSCU)", y="1-s", fill="near-cognate tRNA?")


