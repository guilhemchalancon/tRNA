PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
source("scripts/commands.R"); source("scripts/SMoT.R")
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #

require(XLConnect)
require(plyr)
require(reshape)
require(reshape2)

# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
#                                                                               DATA
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------------------------------------- #

# >> orgininal quantification (raw data)
# fold change data with standard errors for the triplicates (aim is to evaluate the quality of the quantification)
marc.data <- read.table("data/tRNA abundance/tRNAab_quantified.mean_sd.txt",header=1)
marc.data <- merge(marc.data, rosetta.anticodons, by="anticodon", all.x=T)
marc.data$anticodon <- factor( marc.data$anticodon, levels = as.character(unique(arrange(marc.data, aa)$anticodon)) )

# integrate data from XLSX spreadsheet
wb <- loadWorkbook("data/tRNA abundance/Pang2014/tRNAab_Pang2014.xlsx") # change in case of update
pang.data <- readWorksheet(wb, sheet = getSheets(wb)[1],header=T, startRow=1 )


# integrate data from XLSX spreadsheet recomputed manually from SI table 3
wb <- loadWorkbook("data/tRNA abundance/Pang2014/raw_data_parsed.xlsx") # change in case of update
pang.data.2 <- readWorksheet(wb, sheet = getSheets(wb)[3],header=T, startRow=1 )
# pull tRNA isoforms together to quantify anticodon supply 
pang.data.2 <- ddply(pang.data.2, .(anticodon), summarise, mean_normal = sum(mean_normal), mean_ox = sum(mean_ox), foldchange_ox = sum(mean_ox)/sum(mean_normal)  )

pang.data.2$scaled.foldchange <- scale(pang.data.2$mean_ox) / scale(pang.data.2$mean_normal)
pang.data.2$FS.foldchange  <- (pang.data.2$scaled.foldchange - min(pang.data.2$scaled.foldchange))/(max(pang.data.2$scaled.foldchange)-min(pang.data.2$scaled.foldchange))

verification <- merge( pang.data.2[,-2], merge( pang.data, subset(marc.data, experiment == "ox" & time == 60), by = "anticodon" ), by="anticodon")
verification$foldchange.hyp <- ifelse( verification$tRNAab.foldchange_ox > 0, verification$tRNAab.foldchange_ox, - 1/(verification$tRNAab.foldchange_ox) ) 



# our data and theirs based on what they present as "mean fold change"
g1 <- ggplot( data = verification, aes(x= average, y =  tRNAab.foldchange_ox ) ) + 
  geom_vline(xintercept=1, color="gray", lty=3) +
  geom_hline(yintercept=1, color="gray", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray") +
  # stat_smooth(data = subset(verification, anticodon != "UAG"), method="lm") +
  geom_point(aes(color = CAI.O, size = std)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tRNA fold change (Marc)", y = "tRNA fold change (Pang)" ) +
  theme(legend.position="top", axis.title.x = element_text(vjust=0.15))


# are there quantification reported in table S3 relate to mean fold change rations computable from the XLS supplementary file?
g2 <- ggplot( data = verification, aes(x= foldchange_ox, y =  tRNAab.foldchange_ox ) ) + 
  stat_smooth( method="lm", color = "white") +
  geom_point(aes(color = CAI.O), size = 5) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tRNA fold change (recomputed)", y = "tRNA fold change (Table s3)", 
        title = "the measurements Pang et al. show in the paper vs. the ones\nthat are deduced from their raw data" ) +
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))

# are there quantification reported in table S3 relate to mean fold change rations computable from the XLS supplementary file?
g3 <- ggplot( data = verification, aes(x= foldchange_ox, y =  foldchange.hyp ) ) + 
  geom_abline(slope=1, intercept=0, color="gray") +
  stat_smooth( method="lm", color = "white") +
  geom_point(aes(color = CAI.O), size = 5) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tRNA fold change (recomputed)", y = "tRNA fold change (Table s3)",
        title = "recalibrating data from Table S3 gets us closer\n to what can be computed from the raw data" ) +
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))




# our data and theirs based on the interpretation of actual fold changes that they likely heavily modified to increase contrasts between up and down regulated tRNAs
g4 <- ggplot( data = verification, aes(x= average, y =  foldchange_ox ) ) + 
  geom_vline(xintercept=1, color="gray", lty=3) +
  geom_hline(yintercept=1, color="gray", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray") +
  # stat_smooth(data = subset(verification, anticodon != "UAG"), method="lm") +
  geom_point(aes(color = CAI.O, size = std)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tRNA fold change (Marc)", y = "tRNA fold change (Pang)",
        title = "Using actual fold changes from Pang et al." ) + 
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))


# and now ranking: can we salvage something by normalising our results?
g5 <- ggplot( data = verification, aes(x= rank(average), y = rank(foldchange_ox) ) ) + 
  geom_vline(xintercept=1, color="gray", lty=3) +
  geom_hline(yintercept=1, color="gray", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray") +
  # stat_smooth(data = subset(verification, anticodon != "UAG"), method="lm") +
  geom_point(aes(color = CAI.O, size = std)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "ranks in tRNA fold change (Marc)", y = "ranks in tRNA fold change (Pang)",
        title = "Are the ranks similar?") +
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))


# are there quantification of tRNA in normal condition proportional to tGCN?
g6 <- ggplot( data = verification, aes(x= tGCN, y =  mean_normal ) ) + 
  stat_smooth( method="lm") +
  geom_point(aes(color = CAI.O, size = SD_ox)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_log10(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tGCN", y = "tRNA RPKM (normal)", title = "Do RPKM measurements in non-stressed conditions\nshow proportionality with tGCN?" ) +
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))



ggplot( data = verification, aes(x= foldchange_ox, y =  scaled.foldchange ) ) + 
  geom_vline(xintercept=1, color="gray", lty=3) +
  geom_hline(yintercept=1, color="gray", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray") +
  # stat_smooth(data = subset(verification, anticodon != "UAG"), method="lm") +
  geom_point(aes(color = CAI.O, size = std)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "fold change (recomputed)", y = "scaled fold change",
        title = "Using actual fold changes from Pang et al." ) + 
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))


ggplot( data = verification, aes(x= scaled.foldchange, y =  FS.foldchange ) ) + 
  geom_vline(xintercept=1, color="gray", lty=3) +
  geom_hline(yintercept=1, color="gray", lty=3) +
  geom_abline(slope=1, intercept=0, color="gray") +
  # stat_smooth(data = subset(verification, anticodon != "UAG"), method="lm") +
  geom_point(aes(color = CAI.O, size = std)) +
  geom_text(aes(label= anticodon), size =3 ) +
  scale_x_continuous(labels = prettyNum) +
  scale_y_continuous(labels = prettyNum) +
  scale_color_manual(values = c("skyblue3","gray", "orange")) +
  labs( x = "tRNA fold change (Marc)", y = "tRNA fold change (Pang)",
        title = "Using actual fold changes from Pang et al." ) + 
  theme(legend.position="none", axis.title.x = element_text(vjust=0.15))


require(gridExtra)

pdf(file = "results/tRNA.abundance/Pang2014_analysis.pdf", width = 12, height = 17, useDingbats = F )
grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 2)
dev.off()
