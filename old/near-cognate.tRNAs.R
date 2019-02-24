setwd("~/Documents/MRC16")
source("scripts/commands.R")
# near-cognate tRNAs have at most one single mismatch with the cognate motif, and are coding for the wrong amino acid
require(Biostrings)
require(seqinr)

PATH  <- "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"


rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1,sep="\t")
rosetta.anticodons <- read.table("data/info/rosetta.anticodons.txt",header=1,sep="\t")

    codons <- as.character(rosetta.codons$codon)
    anticodons <- as.character(rosetta.anticodons$anticodon)
    
    codons.reverse.complement <- lapply(setNames(codons,codons), function(x) gsub( "T","U", toupper(c2s(comp(rev(s2c(x))))) ) )
    
#     pair.able <- data.frame( on.codon = c("a","c","t","g","t","c","a","g"), 
#                              on.anticodon = c("U","G","A","C","G","A","A","U"), 
#                              complement = c("U","G","A","C","A","G","U","C"), 
#                              type = c(rep("watson.crick",4), rep("wobble",4)) )
#     
old.approach = F
    
bindable.tRNAs <- do.call( rbind, lapply( setNames(codons, codons), function(x){
      
    if(old.approach == T) {  
      # find those tRNAs whose anticodon as at best one mismatch (substitution) with the reverse-complement of the codon
      bindable.tRNA  <- anticodons[ agrep(pattern = codons.reverse.complement[[x]] , x=anticodons, max.distance = list(insertions=0, substitutions=1) ) ]
      
      # correct (cognate) pairing
      correct <- as.character(subset(rosetta.codons, codons == x)$anticodon)
      
      
      if( correct == codons.reverse.complement[[x]]) correct.one <- "cognate" else correct.one  <- "wobble"
    
      result <- data.frame( codon = x, anticodon = bindable.tRNA, 
                  pairing = ifelse( bindable.tRNA %in% correct, correct.one, "near-cognate"),
                  anticodon.prefix = as.character(lapply( bindable.tRNA, function(x) c2s(s2c(x)[2:3]) ))
                )
      result$risk <- ifelse( result$anticodon.prefix == subset(result, !pairing %in% "near-cognate")$anticodon.prefix & result$pairing %in% "near-cognate", "high", ifelse( result$pairing !="near-cognate","no","medium") )
      return(result)
    } else {
      
      # alternative approach APRIL 2015
      all.possible.pairings <- ldply(anticodons, function(u){ 
        d <- data.frame( codon = x, anticodon = u, position.on.codon = 1:3, on.codon = s2c(x), on.anticodon = rev(s2c(u)), anticodon.prefix = c2s(s2c(u)[2:3]) ,
                         complement = rev(s2c(codons.reverse.complement[[x]])))
        d$pairing  <- type.pairing( codon.nt = d$on.codon, anticodon.nt = d$on.anticodon )
        
        # if the putative pairing is I:a elsewhere than in 3rd codon position, then I guess it should count as "bad" since this pairing only happens when A34 is modified, which I think does not happen on other anticodon positions
        impossible.I_a <- with( d, which( position.on.codon < 3 & on.codon == "a" & on.anticodon == "A") )
        
        if( length(impossible.I_a) > 0 ){
          #print(d[impossible.I_a,] )
          #cat("Assuming an I:a pairing in the 1st two codon positions. Corrected to 'bad'\n")
          d[ with( d, which( as.numeric(position.on.codon) < 3 & on.codon == "a" & on.anticodon == "A") ), "pairing" ] <- "bad"
        }
        
        d$position.on.codon <- paste("pos.",1:3,sep="")
        dd <- reshape2::dcast( d, codon + anticodon + anticodon.prefix ~ position.on.codon, value.var = "pairing" )
        return(dd)
      }
      )
      
      all.possible.pairings$pairing <- with( all.possible.pairings, ifelse( pos.1 == "watson.crick" & pos.2 == "watson.crick"  & pos.3 !="bad" , ifelse( pos.3 == "watson.crick", "cognate_wc", "cognate_wobbling"),
                                                                               ifelse( pos.1 != "bad" & pos.2 != "bad" & pos.3 != "bad", "possible.helix" , "unlikely" )
      ) ) 
#       all.possible.pairings$risk <-  with(all.possible.pairings, ifelse( anticodon.prefix == unique(subset(all.possible.pairings, pairing %in% c("cognate_wc", "cognate_wobbling") )$anticodon.prefix) & pairing %in% "possible.helix" , "high",
#                                             ifelse( pairing !="possible.helix","no","medium")
#                                             )
#              )

      all.possible.pairings$propensity <- apply( all.possible.pairings[,c("pos.1","pos.2","pos.3")], 1, function(x) sum(x!="bad") )
      all.possible.pairings$risk <- with(all.possible.pairings, ifelse( propensity < 2, "no", ifelse(propensity > 2, "high", "medium" ) ) )

      
      return(all.possible.pairings)
    }
        
} # end else 
) )
rownames(bindable.tRNAs) <- NULL
    
    bindable.tRNAs <- arrange( merge( bindable.tRNAs, subset(rosetta.codons, select=c(codon,aa, anticodon)), by = "codon"), codon)
    bindable.tRNAs <- bindable.tRNAs[,c("codon","aa","anticodon.y","anticodon.x","pairing","anticodon.prefix","propensity","risk")]
    colnames(bindable.tRNAs)[2:4] <- c("aa","anticodon","incorporated.anticodon")
    bindable.tRNAs <-  arrange( merge( bindable.tRNAs, unique(subset(rosetta.codons, select=c(anticodon,aa))), by.x = "incorporated.anticodon", by.y="anticodon"), codon)
    bindable.tRNAs <- bindable.tRNAs[,c("codon","anticodon","aa.x","incorporated.anticodon","aa.y","pairing","anticodon.prefix","propensity","risk")]
    colnames(bindable.tRNAs)[c(3,5)] <- c("aa","incorporated.aa")
    bindable.tRNAs$misincorporation <- bindable.tRNAs$incorporated.aa != bindable.tRNAs$aa
    bindable.tRNAs$type  <- with( bindable.tRNAs, ifelse( anticodon == incorporated.anticodon, "cognate",
                                                          ifelse( anticodon != incorporated.anticodon & pairing %in% c("cognate_wc","cognate_wobbling"), "near-cognate", 
                                                                  ifelse( pairing == "possible.helix", "possible.near-cognate", "non-cognate" )
                                                                 )
                                                          ) 
                                ) # with
  

# Save table
  write.table(bindable.tRNAs, file = "data/Genetic code/near-cognate.tRNAs.txt", quote=F, sep="\t")

#

# bindable.tRNAs$value <- as.numeric(as.character(factor( paste(bindable.tRNAs$pairing, bindable.tRNAs$risk), 
#                                                         levels=c("cognate no", "near-cognate high","near-cognate medium","wobble no"), 
#                                                         labels=c(1, 0.55,0.35,0.85))))
bindable.tRNAs$value <- as.numeric(
  as.character(
    factor( paste(bindable.tRNAs$type, bindable.tRNAs$risk), 
            levels=c("cognate high", "near-cognate high", "possible.near-cognate high", "non-cognate medium","non-cognate no" ),
            labels=c(1, 0.75,0.5,0.25,0.01)
                                                       # levels=c("cognate no", "near-cognate high","near-cognate medium","wobble no"), 
                                                      #  labels=c(1, 0.55,0.35,0.85))))
    )
  )
)

near.cognate.matrix <- matrixify(reshape2::dcast( bindable.tRNAs, incorporated.anticodon ~ codon, value.var = "value" ))
near.cognate.matrix[is.na(near.cognate.matrix)] <- 0
require(pheatmap)

require(RColorBrewer)
annotations_colors= list( CAI.O = c(N="skyblue",N.O="gray", O="orange"),
                          aa =  colorRampPalette(brewer.pal(n = 8, name = "Set1"))(20))
names(annotations_colors$aa) <- levels(rosetta.codons$aa)




columns    <- unique(as.character(arrange(bindable.tRNAs,  aa, anticodon, incorporated.anticodon, codon )$codon))
rows       <- unique(as.character(arrange(bindable.tRNAs,  aa, anticodon,  incorporated.anticodon, codon )$anticodon))

pdf(paste0(PATH,"part1/genetic.code.near_cognate.tRNAs_new.pdf"), width = 11, height = 11)
pt <- pheatmap( t(near.cognate.matrix[rows,columns]), treeheight_col=15, 
         treeheight_row=15, cellwidth = 10, cellheight = 10,
         cluster_row = F, cluster_col = F,
       #  annotation = data.frame(matrixify(rosetta.codons[,c("codon","aa","CAI.O")])),
         annotation = data.frame(matrixify(rosetta.anticodons[,c("anticodon","aa","CAI.O")])),
         annotation_colors = annotations_colors,
         clustering_method = "ward", 
         color = colorRampPalette(brewer.pal(n = 4, name = "Blues"))(10) )

pheatmap( near.cognate.matrix[rows,columns], treeheight_col=15, 
          treeheight_row=15, cellwidth = 10, cellheight = 10,
          cluster_row = F, cluster_col = F,
          #  annotation = data.frame(matrixify(rosetta.codons[,c("codon","aa","CAI.O")])),
          annotation = data.frame(matrixify(rosetta.codons[,c("codon","aa","CAI.O")])),
          annotation_colors = annotations_colors,
          clustering_method = "ward", 
          color = colorRampPalette(brewer.pal(n = 4, name = "Blues"))(10) )

dev.off()

# ptd <- data.frame( ID=1:41, anticodon = pt$tree_col$labels[pt$tree_col$order])
## or ## 
ptd <- data.frame( ID=1:41, anticodon = rows)

ptd <- arrange(merge(subset(rosetta.anticodons, select=c(anticodon,aa)), ptd, by = "anticodon"), plyr::desc(ID))



copy2clipboard(paste(ptd$anticodon, ptd$aa,sep="-"))

# 
# the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", 
# "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

# ------------------------------------------------------------------------------------------- #
bindable.tRNAs  <- read.table("data/Genetic code/near-cognate.tRNAs.txt", header=T)
#upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)
#bindable.tRNAs <- merge(subset(bindable.tRNAs, type !="non-cognate"), upstress_codons, by = "codon")


df.near.cognate <- arrange( merge(bindable.tRNAs, subset(anticodon.master.table, select=c(anticodon,recognised, experiment, time, adjusted.tGCN, log2.foldchange,relative.availability, relative.availability.normal, total.tRNAab, total.tRNAab.normal)), 
                       by.x = "incorporated.anticodon", by.y="anticodon"), codon, incorporated.aa, experiment, time)
df.near.cognate <- data.table(merge(df.near.cognate, subset(codon.master.table, select=c(codon, experiment, time, demand.mRNAab, demand.mRNAab.normal)), 
                                    by = c("codon","experiment","time") ))

f1 <- function(x){ correct <- x[x$type %in% c("cognate"), ] 
                   others  <- x[x$type %in% c("near-cognate"),] 
                   competition = sum(others[, "total.tRNAab"])/correct[, "total.tRNAab"]
                   competition.normal = sum(others[, "total.tRNAab.normal"])/correct[, "total.tRNAab.normal"]
                   competition.change = ifelse(competition.normal!=0, competition / competition.normal, NA)
                   data.frame(competition  = competition,
                              competition.normal = competition.normal,
                              competition.change = competition.change
                               )
}


codon.competition <- merge( codon.master.table, ddply(df.near.cognate, .(experiment, time, codon), f1), by=c("experiment","time","codon"))
codon.competition$log2.competition <- with(codon.competition, ifelse( competition >0, log2(codon.competition$competition), NA))
#upstress_codons <- read.table("results/lister/clusters_codons.txt",header=T)
#codon.competition <- merge(upstress_codons, codon.competition, by="codon")
write.table(codon.competition, file = "results/master tables/codon.competition.txt",sep="\t",row.names=F)



d=subset( df.near.cognate, experiment == "diauxic" & risk %in% c("high","no") )
 
d.0 <- unique( d[,!colnames(d)%in%c("time","log2.foldchange","relative.availability","total.tRNAab", "demand.mRNAab"), with=F])
d <- rbind(d, data.frame( d.0, time = 0, 
                          log2.foldchange=0, 
                          relative.availability = d.0$relative.availability.normal, 
                          total.tRNAab = d.0$total.tRNAab.normal,
                          demand.mRNAab = d.0$demand.mRNAab.normal) )
d$risk <- factor(as.character(d$risk), levels=c("high","no"), labels = c("near-cognate","cognate/wobble"))


ggplot(d, aes( x=demand.mRNAab, y = total.tRNAab, color = risk)) + 
  #geom_hline(yintercept=0, lty=3) +
  geom_abline(intercept=0, slope=1, lty=3) +
  geom_line( aes(group = incorporated.anticodon)) +
  geom_point(aes(shape=risk) ) +
  scale_color_manual(values=c("red", "skyblue3") ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap( ~ aa + up_in + codon )

ggplot(d, aes( x=time, y = relative.availability/relative.availability.normal, color = misincorporation)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_line( aes(group = paste(risk, incorporated.anticodon))) +
  geom_point(aes(shape= risk) ) +
  scale_color_manual(values=c("skyblue3","red") ) +
  facet_wrap( up_in ~ aa + codon)

ggplot(d, aes( x=time, y = relative.availability, color = misincorporation)) + 
  geom_hline(yintercept=1, lty=3) +
  geom_line( aes(group = paste(risk, incorporated.anticodon))) +
  geom_point(aes(shape= risk) ) +
  scale_color_manual(values=c("skyblue3","red") ) +
  facet_wrap( up_in ~ aa + codon)

ggplot(d, aes( x=time, y = total.tRNAab, color = misincorporation)) + 
  geom_hline(yintercept=0, lty=3) +
  geom_line( aes(group = paste(risk, incorporated.anticodon))) +
  geom_point(aes(shape=risk) ) +
  scale_color_manual(values=c("skyblue3","red") ) +
  facet_wrap( ~ aa + up_in + codon )


ggplot(data= subset(codon.competition, time>0), aes( x = factor(time), y = as.character(codon))) + 
  geom_tile(aes(fill= log2(competition.change) )) +
  facet_grid( up_in.y ~ experiment, drop = T, scales = "free_y") + scale_fill_gradient2(low = "skyblue",mid="black", high="red",midpoint = 0)

ggplot(data= subset(codon.competition, time>0), aes( x = factor(time), y = as.character(codon))) + 
  geom_tile(aes(fill= log2.competition )) +
  facet_grid( up_in.y ~ experiment, drop = T, scales = "free_y") + scale_fill_gradient2(low = "skyblue",mid="black",high="red",midpoint = 0)



require(pheatmap)
d2 <- matrixify(reshape2::dcast( codon.competition, codon ~ experiment + time, value.var = "log2.competition"   ))
d2 <- d2[,c("normal_0","diauxic_20","ox_20","osm_20","temp_20", "diauxic_120","ox_120","osm_120","temp_120")]
pheatmap(d2, cluster_cols = F, cellwidth=20, cellheight=12, main = "cognate vs. near-cognate tRNA competition")

