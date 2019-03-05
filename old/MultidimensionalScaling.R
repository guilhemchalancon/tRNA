
require(FactoMineR)
require(StingRay)
data(TAB)


rosetta.anticodons$GC <- gc.content(rosetta.anticodons$anticodon)
rosetta.anticodons$essential <- ifelse( rosetta.anticodons$anticodon %in% c("CCU","CCG","CUG","GAG","CGU"), "essential", "non.essential" )

explanatory.variables <- merge(rosetta.anticodons, TAB$aminoacid.properties, by = "aa") 
explanatory.variables$heatmap.label <- paste(explanatory.variables$anticodon, explanatory.variables$aa, sep="-")
explanatory.variables <- merge(explanatory.variables, subset(anticodon.master.table, time==0, select=c(anticodon, demand.mRNAab.normal), by="anticodon"))
explanatory.variables$tGCN <- ifelse( explanatory.variables$tGCN < 3, "low", ifelse(explanatory.variables$tGCN > 9, "high", "medium") ) 
explanatory.variables <- data.frame(matrixify(subset(explanatory.variables, select=c(heatmap.label, CAI.O, tGCN, hydrophobicity_category, essential, aa.with.1.anticodon,
                                                                                     metab.cost.resp, metab.cost.ferm, GC, demand.mRNAab.normal))), 
                                    stringsAsFactors = F)
head(explanatory.variables)

cor(condition.data[,1:4])


################################################
################################################
################################################
################################################


# Construct data tables for PCA
anticodon.master.table$label  <- paste(anticodon.master.table$anticodon,anticodon.master.table$aa,sep="-") # labels

# Recode per time (each tRNA has 4 data points, for each condition)
### via anticodon.master.table (no time 60)
coded.by.time      <- reshape2::dcast( subset(anticodon.master.table, time >0), formula = label + experiment ~ time, value.var = "foldchange")

### include t60 
coded.by.time.2 <- melt(2^M.tRNAab$all.together)
require(tidyr)
coded.by.time.2 <- separate(coded.by.time.2, col = "X2", into = c("experiment","time"), sep="_" )
coded.by.time.2 <- reshape2::dcast(coded.by.time.2, formula = X1 + experiment ~ time)
colnames(coded.by.time.2)[1] <- "label"


# Recode per codnition (each tRNA has 3 data points, for each time point)
coded.by.condition <- reshape2::dcast( subset(anticodon.master.table, time >0), formula = label + time ~ experiment, value.var = "foldchange")

### include t60
coded.by.condition.2 <- melt(2^M.tRNAab$all.together)
require(tidyr)
coded.by.condition.2 <- separate(coded.by.condition.2, col = "X2", into = c("experiment","time"), sep="_" )
coded.by.condition.2 <- reshape2::dcast(coded.by.condition.2, formula = X1 + time ~ experiment)
colnames(coded.by.condition.2)[1] <- "label"

# Prepare input tables for PCA (including explanatory variables)

## includes t60
condition.data <- data.frame(matrixify( merge(coded.by.condition.2,explanatory.variables, by.x="label", by.y=0), ref.n = c(1,2) ))
condition.data$time <- as.numeric(gsub(rownames(condition.data), pattern="^(.*)_(.*)$", replacement="\\2"))
condition.data[,c(1:4,10:13)] <- sapply(condition.data[,c(1:4,10:13)], function(x){as.numeric(as.character(x))})

## includes t60
     time.data <- data.frame(matrixify( merge(coded.by.time.2,explanatory.variables, by.x="label", by.y=0), ref.n = c(1,2) ))
     time.data$condition <- as.factor(gsub(rownames(time.data), pattern="^(.*)_(.*)$", replacement="\\2"))
     time.data[,c(1:3, 9:12)] <- sapply(time.data[,c(1:3,9:12)], function(x){as.numeric(as.character(x))})

require(FactoMineR)
# includes quanti explanatory var. in the model (ie they are active)
PCA.data <- subset(condition.data, time !=60)
PCA.condition <- PCA(X = PCA.data, 
                     ind.sup = c(grep("CAA",rownames(PCA.data)), grep("UAG",rownames(PCA.data))),
                     quali.sup = c(5:9,14), 
                     scale.unit = T)
     #PCA.time <- PCA(X = time.data, quali.sup = c(3:8,14), scale.unit = T)

plot(PCA.condition,  cex = 0.75, axis = c(1,2), select = "cos2 0.7", xlim=c(-5,7),
     habillage=5,col.hab=c("skyblue3","gray40","orange")  )
plot(PCA.condition, choix="var", cex = 0.75, axis = c(1,2), select = "cos2 0.3")
plotellipses(PCA.condition,    keepvar = c(5:9), axis = c(1,2))
plotellipses(PCA.condition,    keepvar = c(5:9), axis = c(2,3))


FA.condition <- factanal( subset(condition.data, time !=60)[,c(1:4,10:12)], factors = 3, scores = "regression", rotation = "promax")
FA.condition

# correlation of all active variables
require(pheatmap)
pheatmap(cor(condition.data[,c(1:4,10:13)]), cluster_rows = F, cluster_cols = F)

# correlation of explicative variables
cor(condition.data[,10:13])
summary(condition.data)
summary(PCA.condition)
# PCA CONDITION MODEL
head(condition.data)
estim_ncp(X = PCA.data[,-c(5:9,14)])
estim_ncp(X = condition.data[,-c(5:9,14)])

dimdesc(PCA.condition)
PCA.condition$var$cos2
PCA.condition$var$coord



setwd("~/Documents/MRC16/")
summary(PCA.condition,   nbelements=Inf,file="results/PCA.txt")
PCA.condition$var$cos2
PCA.condition$ind$contrib
dimPCA.condition  <- dimdesc(PCA.condition)

dimPCA.condition$Dim.1$quali
dimPCA.condition$Dim.1$quanti

colnames(condition.data)
plot(PCA.condition,choix="ind", cex=0.6, habillage = 5, select="cos2 0.75", unselect=0, col.quali = "blue", col.hab=c("skyblue","gray","orange")  )
par(mfrow=c(1,2))
plot(PCA.condition,choix="var",new.plot=FALSE, select="contrib 5")

plot(PCA.time,choix="var",new.plot=FALSE, select="cos2 0.6")
plot(PCA.time,choix="ind",new.plot=FALSE, habillage = 3,select="cos2 0.6", col.hab=c("skyblue","gray","orange"))

summary(PCA.condition)
PCA.condition$var
cor(condition.data[,-c(5:10,16)])


# PCA ALL MODEL (t=20 and t=120)
quantitative.data <- cbind(2^M.tRNAab$all.together)
all.data <- data.frame(matrixify(merge(quantitative.data,explanatory.variables, by=0)))
#all.data <- all.data[,-10] # remove AAs

quanti.vars <- c(colnames(all.data)[1:ncol(quantitative.data)], colnames(explanatory.variables)[-c(1:5)]) 
all.data[,quanti.vars] <- sapply(all.data[,quanti.vars], function(x){as.numeric(as.character(x))})


cor(all.data[,1:ncol(quantitative.data)])

require(FactoMineR)
PCA.all <- PCA(X = all.data, scale.unit = T, quali.sup = 13:17 ) # ,quanti.sup=19:23 -4
plotellipses(PCA.all, keepvar = (c(13,15:18) - 4))

dimdesc(PCA.all)

plot(PCA.all, habillage=9, col.hab=c("skyblue3","gray40","orange"), select="cos2 0.6", unselect=0, cex=0.8)
plot(PCA.all,choix="var",new.plot=FALSE, select="cos2 0.6")

PCA.all$var$cos2
PCA.all$var$coord
dimdesc(PCA.all)



plot(PCA.condition, habillage=5, col.hab=c("skyblue3","gray40","orange"), select="cos2 0.6", unselect=0, cex=0.8)
plot(PCA.condition,choix="var",new.plot=FALSE, select="cos2 0.6")

plotellipses(PCA.all, cex=1, axis = c(2,3))
plotellipses(PCA.condition, cex=0.85, magnify=1.62, select= "cos2 0.6", level=0.99)
dev.copy2pdf(device = quartz, file = paste0(PATH,"part1/PCA_ellipses.pdf"), useDingbats=F )




##############################################################

require("tsne")
data(TAB)

# - variables coded by condition (time points as individuals)
indesirables <- c(grep("Thr", rownames(condition.data)), grep("CAA", rownames(condition.data)), grep("UAG", rownames(condition.data)), grep("CAU", rownames(condition.data)) )
M.condition  <- scale( subset(condition.data, time > 0 )[ - indesirables , -c(5:9,14)],center = T,scale = T)

      head(M.condition)
      summary(M.condition)
      cor(M.condition)

M.time  <- scale( time.data[ - indesirables , -c(5:9,14)],center = T,scale = T)

      head(M.time)
      summary(M.time)
      cor(M.time)

indesirables <- c( grep("CAA", rownames(all.data)), grep("UAG", rownames(all.data)), grep("CAU", rownames(all.data)) )
M.all  <- scale( all.data[ - indesirables , -c(13:17,20)],center = T,scale = T)

# RUN tSNE
  M <- M.all[,1:12]
  summary(M)
  data <- all.data
  tRNAab.tsne <- data.frame(tsne( dist(M)/100, perplexity = 5, k = 3, max_iter = 2500 ))
  colnames(tRNAab.tsne)[1:3] <- c("Dim1","Dim2","Dim3")
  tRNAab.tsne$ind <- rownames(M)
  tRNAab.tsne.data <- merge(tRNAab.tsne, data,by.x="ind", by.y=0)
  #tRNAab.tsne.data <- separate(tRNAab.tsne.data, col = "ind", sep = "_", into = c("label","time"))

tRNAab.tsne.data$profile  <- paste( tRNAab.tsne.data$essential, tRNAab.tsne.data$aa.with.1.anticodon, tRNAab.tsne.data$tGCN )
tRNAab.tsne.data$profile <- factor( tRNAab.tsne.data$profile,
                                    level= sort(unique(tRNAab.tsne.data$profile)),
                                    label= c("essential-isoacceptor-rare","viable-unique-abundant",
                                             "viable-unique-interm","viable-isoacceptor-abundant",
                                             "viable-isoacceptor-rare","viable-isoacceptor-mixed")
                                  )

#write.table(tRNAab.tsne.data, file = "results/tSNE.txt", quote=F)

tRNAab.tsne.data <- read.table("results/tSNE.txt", header=1)
require(FactoMineR)
 
# res.pca.tsne <- PCA(X = all.data[ - indesirables,], scale.unit = T, quali.sup=c(13:17,20), quanti.sup = 18:21)
# plot(res.pca.tsne, habillage = 16)
# plot(res.pca.tsne, habillage = 16, axes=c(1,3))
# plotellipses(res.pca.tsne)
# estim_ncp(X = all.data[ - indesirables,-c(13:21)])
# dimdesc(res.pca.tsne)



Group1 <- as.character(tRNAab.tsne.data$ind[tRNAab.tsne.data$Dim2 > 0 & tRNAab.tsne.data$Dim1<0])
group1 <- gsub(Group1, pattern="^(.*)\\-(.*)$", replacement="\\1")
Group3 <- as.character(tRNAab.tsne.data$ind[tRNAab.tsne.data$Dim2 < 0 & tRNAab.tsne.data$Dim1>100])
group3 <- gsub(Group3, pattern="^(.*)\\-(.*)$", replacement="\\1")

tRNAab.tsne.data$is.G1 <- ifelse( tRNAab.tsne.data$ind %in% Group1, "G1", "G2")
tRNAab.tsne.data$group <- ifelse( tRNAab.tsne.data$ind %in% Group3, "G3", ifelse(tRNAab.tsne.data$ind %in% Group1, "G1", "G2")  )


# minimal dataset for subsequent analysis of Gr.1-tRNAs in other scripts
d <- arrange(data.frame( 
            anticodon = as.character(rosetta.anticodons$anticodon), 
            aa = as.character(rosetta.anticodons$aa), 
            group = with( rosetta.anticodons, ifelse( anticodon %in%  group1, "Gr.1", ifelse(anticodon %in% group3, "Gr.3", "Gr.2") ) ),
            coded.by.essential.gene = with( rosetta.anticodons, ifelse( anticodon %in%  c("CCU","CCG","CUG","GAG","CGU"), "essential.tRNA.gene", "multiple.tRNA.genes")  ), 
            stringsAsFactors = F
          ), group, anticodon, aa)

d[ with(d, anticodon %in% c("CAA","UAG","CAU")), "group"] <- "outgroup"
write.table(d, "results/tSNE_groups.txt", sep="\t", quote=F, row.names=F)
 



### PLOTS ###
g.tSNE <- ggplot(data=tRNAab.tsne.data, aes(x=Dim1, y=Dim2)) + 
  geom_vline(xintercept=0, lty=3, color="gray40") + 
  geom_hline(yintercept=0, lty=3, color="gray40") +
  geom_point( aes(fill= paste(essential, aa.with.1.anticodon), color= group), pch=21, size=5.5) + 
  geom_text( aes(label = ind), size = 3 ) +
  scale_fill_manual(values=c("red","skyblue3","purple")) +
  #scale_shape_manual(values = c(15,16)) +
  #theme(legend.position = "bottom") +
  labs( color = "profile", color = "profile", title = "t-SNE of differential tRNA expression") 

ggsave(plot=g.tSNE, filename = paste0(PATH,"part1/tSNE.pdf"), useDingbats=F, width=8.51, height=5.51)

# --------------------------------- # START OF REVISION MADE IN 2018 # --------------------------------- #
      # AMMENDUM 2018: Science Signaling review
      require(magrittr) 
      
      # first I reloaded tRNAab.tsne.data
      tRNAab.tsne.data <- read.table("results/tSNE.txt", header=1)
      
      # then took the tSNE dimensions' scores and created a numeric matrix to conform to required inputs into kmeans()
      tRNAab.tsne.m <- as.matrix(tRNAab.tsne.data[,c("Dim1","Dim2","Dim3")])
      rownames(tRNAab.tsne.m) <- tRNAab.tsne.data$ind
      
      # perform k-means clustering on outputs of the t-SNE
      tRNAab.tsne.k <- tRNAab.tsne.m %>% kmeans( centers = 3) 
      
      # column-bind cluster results to input data
      tRNAab.tsne.d <- data.frame( tRNAab.tsne.data, cluster = tRNAab.tsne.k$cluster )
      
      # plot clustering along Dim1 and Dim2
      ggplot(tRNAab.tsne.d, aes(x=Dim1, y=Dim2)) +
        geom_vline(xintercept=0, lty=3, color="gray40") + 
        geom_hline(yintercept=0, lty=3, color="gray40") +
        geom_point( aes(fill= factor(cluster)), color = "white", pch=21, size=5.5) + 
        geom_text( aes(label = ind), size = 2 ) + labs( fill = "k-means")
      
      # save cluster data as csv file 
      write.csv( tRNAab.tsne.d[,c("ind","tGCN","essential","Dim1","Dim2","Dim3","cluster")], 
                 file ="results/SciSignalling 2018 revisions/tRNAab.tse.clusters.csv",
                 quote = FALSE)
      
      # save results of k-means as serialised data file
      save(tRNAab.tsne.k,file = "results/SciSignalling 2018 revisions/tRNAab.tse.kmeans.Rda")

# --------------------------------- # END OF REVISION MADE IN 2018 # --------------------------------- #

# late additions to anticodon.master.table
anticodon.tSNE <- subset(anticodon.master.table, time>0)
anticodon.tSNE$essential <- ifelse( anticodon.tSNE$anticodon %in% c("CCU","CCG","CUG","GAG","CGU"), T, F )
anticodon.tSNE$is.G1 <- ifelse(anticodon.tSNE$anticodon %in% group1, T, F )
anticodon.tSNE$group <- ifelse(anticodon.tSNE$anticodon %in% group1, "G1", ifelse(anticodon.tSNE$anticodon %in% group3, "G3", "G2")  )
anticodon.tSNE$set <- factor(paste(anticodon.tSNE$time, anticodon.tSNE$group), levels = c("20 G1","20 G2", "20 G3","120 G1", "120 G2","120 G3") )
#anticodon.tSNE$set <- factor(paste(anticodon.tSNE$time, anticodon.tSNE$group), levels = c("20 TRUE","20 FALSE", "120 TRUE", "120 FALSE") )
# anticodon.tSNE$set.triad <- factor(paste(anticodon.tSNE$time, anticodon.tSNE$triad), 
#                                      levels = c("20 G1","20 G2", "20 G3","120 G1", "120 G2","120 G3") )
anticodon.tSNE$category <- ifelse( anticodon.tSNE$aa.with.1.anticodon == T, "unique", ifelse( anticodon.tSNE$essential == T, "essential", "other" ) )

g <- ggplot(data= anticodon.tSNE, aes(x=set, 
                                                        y=log2.foldchange) ) +
  geom_hline(yintercept=log2(0.5), lty=3) +
  geom_hline(yintercept=0, lty=3) +
  geom_hline(yintercept=log2(1.5), lty=3) +
  geom_dotplot( aes(fill=category), binaxis = "y", stackdir = "center", dotsize=1.2, stackratio=1.15) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.65) +
  facet_wrap( ~ experiment, ncol=1) +
  scale_fill_manual(values=c("red","skyblue3","purple")) +
  labs( x = "", y = "tRNA fold change\n(log2)" ) + theme(legend.position="none")


wc.diff <- ddply( anticodon.tSNE, .(experiment, time), function(x){ wilcox.test.batch(x =x , grouping.vars = c("is.G1"), value = "log2.foldchange" ) } )

wc.groups <- ddply( anticodon.tSNE, .(experiment, time), function(x){ wilcox.test.batch(x =x , grouping.vars = c("group"), value = "log2.foldchange" ) } )
wc.groups$p.value <- as.numeric(as.character(wc.groups$p.value))
wc.groups$significance <- ifelse( wc.groups$p.value < 0.001, "***", ifelse( wc.groups$p.value < 0.01, "**", ifelse( wc.groups$p.value < 0.05, "*", "")))

dlply( subset(wc.groups, select=c(experiment,time,X1,X2, p.value, r, median.1, median.2, significance)), .(experiment, time), head)

with( data = subset(anticodon.tSNE, !anticodon %in% c("CAU","CAA","UAG") & experiment == "ox" & time >20 ), expr = table(group, CAI.O )  )

rosetta.codons$group  <- with( data = rosetta.codons, expr = ifelse( anticodon %in% group1, "G1", ifelse( anticodon %in% group3, "G3", "G2") ) )

with( subset(rosetta.codons, !anticodon %in% c("CAU","CAA","UAG") ), table( group, CAI.O) )


g <- ggplot(data= subset(anticodon.tSNE, !anticodon %in% c("UAG","CAU","CAA")), 
       aes(x=triad, y=log2.foldchange) ) +
  geom_hline(yintercept=log2(0.5), lty=3) +
  geom_hline(yintercept=0, lty=1, color="gray35") +
  geom_hline(yintercept=log2(1.5), lty=3) +
  geom_dotplot( aes(fill=category), binaxis = "y", stackdir = "center",  stackratio=1.15) + #dotsize=1.2,
  #geom_text( data = subset(anticodon.tSNE, !anticodon %in% c("UAG","CAU","CAA") & category=="unique"), size=4, aes(label=anticodon)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.65) +
  geom_hline(yintercept=0, lty=3) +
  #geom_line(aes(group=anticodon)) +
  facet_wrap(  ~ experiment +time, nrow=1, scales = "free_x") +
  scale_fill_manual(values=c("red","skyblue3","purple")) +
  labs( x = "", y = "tRNA fold change\n(log2)" ) + theme(legend.position="none")

ggsave(plot=g, filename= paste0(PATH,"part1/tSNE_dotplot.pdf"), useDingbats=F, width=10.1, height=2.56 )




# plot tSNE scalign
ggplot(data=tRNAab.tsne.data, aes(x=Dim1, y=Dim2)) + 
  geom_vline(xintercept=0, lty=3, color="gray40") + 
  geom_hline(yintercept=0, lty=3, color="gray40") +
  #geom_point( aes(color= paste(essential, aa.with.1.anticodon), shape=CAI.O), size=5.5) + 
  geom_point( aes(color= paste(essential, aa.with.1.anticodon), shape=CAI.O), size=5.5) + 
  geom_text( aes(label = ind), size = 3 ) +
  scale_color_manual(values=c("red","skyblue3","purple")) +
  #scale_shape_manual(values = c(15,16)) +
  theme(legend.position = "bottom") +
  labs( color = "profile", color = "profile", title = "t-SNE of differential tRNA expression") 




ddply( anticodon.tSNE, .(experiment, time), function(x){ wilcox.test.batch(x =x , grouping.vars = c("group"), value = "demand.mRNAab" ) } )

ggplot(data=tRNAab.tsne.data, aes(x=Dim1, y=Dim2)) + geom_point( aes(size=as.numeric(demand.mRNAab.normal))) + coord_fixed()
ggplot(data=tRNAab.tsne.data, aes(x=Dim1, y=Dim3)) + geom_point( aes(color=CAI.O)) + coord_fixed()










#------------------------------------------------------------------------------------------------------------------------------#

