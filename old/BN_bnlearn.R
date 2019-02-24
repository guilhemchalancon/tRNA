setwd("~/Documents/MRC16/")
#----------------------------------------------------------------------------------
require(bnlearn)

PATH  <-  "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/"
DOWN <- "#185665"
UP <- "#A36A23"


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                         DATA PREPARATION                                                                                               #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


# 1 # data frame
study.translation.kinetics <- read.table("results/master tables/study.translation.kinetics.txt", header=T)



# I am going to pool together data from various experimental conditions for the ABN analysis: I must first scale data 


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                               SCALING                                                                                                  #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


source.data <- subset(study.translation.kinetics, time >0)
# no not sure I should!
source.data <- ddply(source.data, .(experiment, time), mutate, 
                     s_log2.mRNA_abundance = scale(log2.mRNA_abundance),
                     s_scaled.change.tAI   = scale(scaled.change.tAI),
                     s_variation_speed     = scale(variation_speed),
                     s_variation_initiation_frequency = scale(log2(1/ratio_initiation)),
                     s_protein_production_rate = scale(global.protein_synthesis.rate)
)
# scaling
source.data$quant.mRNAab <- quantile.it(source.data$s_log2.mRNA_abundance, N = 3 ) # already there      
source.data$quant.tAI   <- quantile.it( source.data$s_scaled.change.tAI, N = 3 )
source.data$quant.speed <- quantile.it( source.data$s_variation_speed, N = 3 )
source.data$quant.freq <- quantile.it( source.data$s_variation_initiation_frequency, N = 3 )
source.data$quant.prod <- quantile.it( source.data$s_protein_production_rate, N = 3 )


source.data$quant.mRNAab.10 <- quantile.it(source.data$s_log2.mRNA_abundance, N = 10 ) # already there      
source.data$quant.tAI.10   <- quantile.it( source.data$s_scaled.change.tAI, N = 10 )
source.data$quant.speed.10 <- quantile.it( source.data$s_variation_speed, N = 10 )
source.data$quant.freq.10 <- quantile.it( source.data$s_variation_initiation_frequency, N = 10 )
source.data$quant.prod.10 <- quantile.it( source.data$s_protein_production_rate, N = 10 )

source.data$quant.mRNAab.5 <- quantile.it(source.data$s_log2.mRNA_abundance, N = 9 ) # already there      
source.data$quant.tAI.5   <- quantile.it( source.data$s_scaled.change.tAI, N = 9 )
source.data$quant.speed.5 <- quantile.it( source.data$s_variation_speed, N = 9 )
source.data$quant.freq.5 <- quantile.it( source.data$s_variation_initiation_frequency, N = 9 )
source.data$quant.prod.5 <- quantile.it( source.data$s_protein_production_rate, N = 9 )


source.data$quant.mRNAab.unsc.10<- quantile.it(source.data$log2.mRNA_abundance, N = 10 ) # already there      
source.data$quant.tAI.unsc.10  <- quantile.it( source.data$s_scaled.change.tAI, N = 10 )
source.data$quant.speed.unsc.10 <- quantile.it( source.data$variation_speed, N = 10 )
source.data$quant.freq.unsc.10 <- quantile.it( log2(1/source.data$ratio_initiation), N = 10 )
source.data$quant.prod.unsc.10 <- quantile.it( source.data$global.protein_synthesis.rate, N = 10 )


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                           DISCRETISATION                                                                                               #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


source.data <- with(source.data, expr = data.frame(source.data,
                                                   quant.mRNAab.pooled = quantile.it( log2.mRNA_abundance, N = 3 ),
                                                   quant.tAI.pooled   = quantile.it( variation_tAI, N = 3 ),
                                                   quant.speed.pooled = quantile.it( variation_speed, N = 3 ),
                                                   quant.freq.pooled = quantile.it( log2(1/ratio_initiation), N = 3 ),
                                                   prod = scale(global.protein_synthesis.rate)
)
)


thres.down <- -0.2
thres.up <- 0.2
source.data <- with(source.data, expr = data.frame(source.data,
                                                   thresh.freq.pooled   = cut(variation_initiation_frequency, 
                                                                              breaks = c( min(variation_initiation_frequency), thres.down, thres.up, max(variation_initiation_frequency)), 
                                                                              labels = c("down","interm","up"), ordered_result = T, include.lowest = T),
                                                   thresh.mRNAab.pooled = cut(log2.mRNA_abundance,
                                                                              breaks = c( min(log2.mRNA_abundance), thres.down, thres.up, max(log2.mRNA_abundance)), 
                                                                              labels = c("down","interm","up"), ordered_result = T, include.lowest = T),
                                                   thresh.tAI.pooled    = cut(variation_tAI,
                                                                              breaks = c( min(variation_tAI), thres.down, thres.up, max(variation_tAI)), 
                                                                              labels = c("down","interm","up"), ordered_result = T, include.lowest = T),
                                                   thresh.speed.pooled  = cut(variation_speed,
                                                                              breaks = c( min(variation_speed), thres.down, thres.up, max(variation_speed)), 
                                                                              labels = c("down","interm","up"), ordered_result = T, include.lowest = T),
                                                   thresh.prod.pooled = cut(global.protein_synthesis.rate,
                                                                              breaks = c( min(global.protein_synthesis.rate), (1+thres.down), (1+thres.up), max(global.protein_synthesis.rate)), 
                                                                              labels = c("down","interm","up"), ordered_result = T, include.lowest = T)
                                                   )
)


summary(source.data)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                             SUBSETTING                                                                                                 #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

# extract data subset with variables of interest
data.qt <- subset( source.data, select = c(name, label, quant.mRNAab, quant.tAI, quant.freq, quant.speed, quant.prod, prod )) 
data.qt.10 <- subset( source.data, select = c(name, label, quant.mRNAab.10, quant.tAI.10, quant.freq.10, quant.speed.10, quant.prod.10, prod )) 
data.qt.5  <- subset( source.data, select = c(name, label, quant.mRNAab.5, quant.tAI.5, quant.freq.5, quant.speed.5, quant.prod.5, prod )) 
data.th <- subset( source.data, select = c(name, label, thresh.mRNAab.pooled, thresh.tAI.pooled, thresh.freq.pooled, thresh.speed.pooled, thresh.prod.pooled, prod )) 
data.unsc.10 <- subset( source.data, select = c(name, label, quant.mRNAab.unsc.10, quant.tAI.unsc.10, quant.freq.unsc.10, quant.speed.unsc.10, quant.prod.unsc.10, prod )) 

# select variables that will be used for the BN learning
model.var <- c("a","tAI","e","i","P")

# SEL  <- c("0-20","80-100") # more conservative
#SEL  <- c("0-20","20-40","60-80","80-100")
SEL  <- c("0-33","67-100")

# Bayesian Network Input Data (restrict to only those observations that are in the right space)
#up.data <- subset( data,  quant.mRNAab %in% SEL & quant.tAI %in% SEL & quant.freq %in% SEL & quant.speed %in% SEL & quant.prod %in% SEL) 
colnames(data.qt)[grep("quant",colnames(data.qt))]  <- c("a","tAI","i","e","P") # check the order is fine!
colnames(data.qt.10)[grep("quant",colnames(data.qt.10))]  <- c("a","tAI","i","e","P") # check the order is fine!
colnames(data.qt.5)[grep("quant",colnames(data.qt.5))]  <- c("a","tAI","i","e","P") # check the order is fine!
colnames(data.th)[grep("thresh",colnames(data.th))] <- c("a","tAI","i","e","P") # check the order is fine!
colnames(data.unsc.10)[grep("unsc",colnames(data.unsc.10))] <- c("a","tAI","i","e","P")


data.qt$label <- as.character(data.qt$label)
data.qt$ID <- 1:nrow(data.qt)

data.qt.10$label <- as.character(data.qt.10$label)
data.qt.10$ID <- 1:nrow(data.qt.10)

data.qt.5$label <- as.character(data.qt.5$label)
data.qt.5$ID <- 1:nrow(data.qt.5)


data.th$label <- as.character(data.th$label)
data.th$ID <- 1:nrow(data.th)


data.unsc.10$label <- as.character(data.unsc.10$label)
data.unsc.10$ID <- 1:nrow(data.unsc.10)

# same as data.th, except that all entries with at least one "interm" have been removed
data.th.2 <- subset(data.th, a %in% c("down","up") & tAI %in% c("down","up") & i %in% c("down","up") & e %in% c("down","up") & P %in% c("down","up") )
data.th.2[,model.var] <- data.frame(sapply(data.th.2[,model.var], as.character))
data.th.2$label <- as.character(data.th.2$label)
data.th.2$ID <- 1:nrow(data.th.2)

# data.th recoded to only have 2 categories, up and not up
data.th.up <- data.th
data.th.up[,model.var] <- data.frame( sapply( data.th.up[,model.var], function(x){ ifelse( as.character(x) == "up","up","not_up") }  ) )


data.unsc.up <- data.unsc.10
data.unsc.up[,model.var] <- data.frame( sapply( data.unsc.up[,model.var], function(x){ ifelse( x %in% c("0-10","10-20"),"up","not_up") }  ) )


data.th.down <- data.th
data.th.down[,model.var] <- data.frame( sapply( data.th.down[,model.var], function(x){ ifelse( as.character(x) == "down","down","not_down") }  ) )



# Observations: wow
          # subset(source.data, n.events > 140000)
          # summary(data)
          # #summary(down.data)
          # nrow(data) # 29,817 observations across the data
          # 

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                              LEARNING                                                                                                  #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
# directd arcs that are not permitted (nonsensical)
blacklist <- data.frame(from= c(rep("P",times=4),"tAI","a","e","i","i","e"), to=c("i","e","a","tAI","a","tAI","a","a","tAI","tAI"))
# direct arcs that should be included
whitelist <- data.frame(from="tAI", to="e")

# two types of methods: cosntraint-based algorithms, and score-based algorithms.
# the first learn the structure by analysing the probabilistic relations entailed by the Markov property of Bayesian networks with conditional independence tests
# and then construct a graph which satisfies the corresponding d-separation statements.
# Score-based algorithms assign a score to each candidate BN and try to maximise the score (typically likelihood, or mutual information, or AIC?)
# require(prob)
# 
# prob::prob(x = up.data[,model.var])


# Learning by Constraint-based methods
# 1 - Growth Shrink
  # categories based on tertiles, ordinal 
  bn.gs.qt     <- gs(     data.qt[,model.var], blacklist = blacklist, whitelist = whitelist )
  bn.gs.qt     <- set.arc(bn.gs.qt, "i", "e")
  # deciles categories based on scaled data, ordinal 
  bn.gs.qt.10  <- gs(     data.qt.10[,model.var], blacklist = blacklist, whitelist = whitelist )
  bn.gs.qt.5   <- gs(     data.qt.5[,model.var], blacklist = blacklist, whitelist = whitelist )
  # deciles categories based on scaled data, ordinal  
  bn.gs.unsc.10  <- gs(     data.unsc.10[,model.var], blacklist = blacklist, whitelist = whitelist )
  bn.gs.unsc.up  <- gs(     data.unsc.up[,model.var], blacklist = blacklist, whitelist = whitelist )
  #bn.gs.qt.10     <- set.arc(bn.gs.qt.10, "i", "e")
  # categories based on threshold, ordinal
  bn.gs.th     <- gs(     data.th[,model.var], blacklist = blacklist, whitelist = whitelist )
  bn.gs.th     <- set.arc(bn.gs.th, "i", "e")
  # categories based on threshold, binary
  bn.gs.th2    <- gs(data.th.2[,model.var], blacklist = blacklist, whitelist = whitelist )
  # up vs. not up
  bn.gs.up     <- gs(data.th.up[,model.var], blacklist = blacklist, whitelist = whitelist )
  # down vs. not down
  bn.gs.down   <- gs(data.th.down[,model.var], blacklist = blacklist, whitelist = whitelist )
  # Learning by Score-based methods # 3 - Hill-climbing
  bn.hc <- hc(data.qt.10[,model.var], score = "aic", blacklist = blacklist, whitelist = whitelist )


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                          VISUALISE GRAPHS                                                                                              #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


require(Rgraphviz)
highlight.opts <- list(nodes = c("P"), fill = "gray30", col="black", textCol="white")

par(mfrow=c(2,2))
graphviz.plot(bn.gs.qt,  highlight = highlight.opts)
graphviz.plot(bn.gs.qt.10,  highlight = highlight.opts)
graphviz.plot(iamb(     data.qt.10[,model.var], blacklist = blacklist, whitelist = whitelist ))
graphviz.plot(bn.gs.qt.5,  highlight = highlight.opts)
graphviz.plot(bn.hc, highlight = highlight.opts)
graphviz.plot(gs(     data.qt.10[,model.var] )  )
graphviz.plot(gs(     data.qt.10[,model.var], blacklist = blacklist, whitelist = whitelist )  )

graphviz.plot(gs( data.qt[,model.var] ),   highlight = highlight.opts)
graphviz.plot(gs( data.qt[,c("a","i","e","tAI")], blacklist = blacklist[7:9,] ) )
graphviz.plot(gs( data.th[,c("a","i","e","tAI","P")] ) )
graphviz.plot(bn.gs.th,   highlight = highlight.opts)
graphviz.plot(bn.gs.th2,  highlight = highlight.opts)
graphviz.plot(bn.gs.up,   highlight = highlight.opts)
graphviz.plot(bn.gs.down, highlight = highlight.opts)


graphviz.plot( gs(     data.qt[,model.var] ) )
graphviz.plot( gs(     data.th[,model.var], blacklist = blacklist[1:6,]  ))
graphviz.plot( gs(     data.th.2[,model.var], blacklist = blacklist  ))
graphviz.plot( gs(     data.th.2[,model.var] ) )



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                   CONDITIONAL PROBABILITIES                                                                                              #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


# Fit BN structures to obtain conditional probability tables
BN.GS.qt   <- bn.fit(x = bn.gs.qt, data = data.qt[,model.var]) #, ordinal=T)
BN.GS.qt.10   <- bn.fit(x = bn.gs.qt.10, data = data.qt.10[,model.var]) #, ordinal=T)
BN.GS.th   <- bn.fit(x = bn.gs.th, data = data.th[,model.var]) #, ordinal=T)



probas.BN.GS.qt   <- lapply( BN.GS.qt,   function(x){ y <- as.data.frame(x$prob); y[,"Freq"] <- round(y[,"Freq"],2); y <- y[complete.cases(y),]; return(y)} ) 
probas.BN.GS.qt.10   <- lapply( BN.GS.qt.10,   function(x){ y <- as.data.frame(x$prob); y[,"Freq"] <- round(y[,"Freq"],2); y <- y[complete.cases(y),]; return(y)} ) 
probas.BN.GS.th   <- lapply( BN.GS.th,   function(x){ y <- as.data.frame(x$prob); y[,"Freq"] <- round(y[,"Freq"],2); y <- y[complete.cases(y),]; return(y)} ) 

# bn.gs.qt.10$learning
# bn.gs.qt.10$nodes
# 
# nparams(bn.gs, data = data[,model.var],debug = T)

# Test the conditional independence of X and Y given Z
    # ci.test( "i", "a", "tAI", data = data.th[,model.var] )
    # ci.test( "i", "tAI", "a", data = data.th[,model.var] )
    # ci.test( "i", "tAI", "a", data = data.th[,model.var] )
    # 
    # ci.test( "a", "tAI",  data = data.th[,model.var] )
    # ci.test( "a", "i",  data = data.th[,model.var] )
    # ci.test( "a", "P",  data = data.th[,model.var] )
    # ci.test( "a", "e",  data = data.th[,model.var] )



# Systematically compute conditional independence tests
d <- do.call(rbind, apply( data.frame(t(combn(model.var, 2))), 1, function(x){ 
  X <- as.character(x[[1]])
  Y <- as.character(x[[2]])
  Zs <- setdiff( model.var, c(X,Y) )
  data.frame( X = X, Y = Y, Z = Zs)
 }
))


posterior.links <- with(data.frame(bn.gs.qt.10$arcs), paste(from, to, sep="->"))

# ------------------------------------------------------------------------------------------------ #
#                        Verify validity of conditional independence statements                    #
# ------------------------------------------------------------------------------------------------ #

# ---- model based on binary thresholds (excluding the "interm" category)
arcs <- data.frame(bn.gs.th2$arcs)
verify.ci.th2 <-   dlply(arcs,.(to), function(x){ 
  all.nodes <- unique(c(as.character(arcs$from), as.character(arcs$to)))
  child <- unique(as.character(x$to))
  parents <- as.character(x$from) 
  descendants  <- as.character(subset(arcs, from == child )$to)
  non.descendants <- setdiff( all.nodes, c(child,parents,descendants)  ) 
  if(length(non.descendants)>0){
    lapply(non.descendants, function(N) { 
      return( ci.test( x = child, y = N, z=parents, data = data.th.2[, model.var] ) )
    })
  } else { return( "not applicable" ) }
})


#-------

# ---- model based on deciles (excluding the "interm" category)
arcs <- data.frame(bn.gs.qt.10$arcs)
verify.ci.qt10 <- do.call( rbind, dlply(.data = arcs, .variables = .(to), .inform=T , .fun = function(x){ 
  all.nodes <- unique(c(as.character(arcs$from), as.character(arcs$to)))
  child <- unique(as.character(x$to))
  parents <- as.character(x$from) 
  descendants  <- as.character(subset(arcs, from == child )$to)
  non.descendants <- setdiff( all.nodes, c(child,parents,descendants)  ) 
  if(length(non.descendants)>0){
    do.call(rbind, lapply(non.descendants, function(N) { 
      t <- ci.test( x = child, y = N, z=parents, data = data.qt.10[, model.var], test = "jt" )
        return( data.frame(
          x = child, y = N, z = paste(parents,collapse="+"),
          JT = round(t$statistic,2), p.value = t$p.value, 
          method = t$method, 
          label = t$data.name, 
          alternative = t$alternative
          )  )
    }))
  } else { return( "not applicable" ) }
}))



# independence tests

d.ind <- arrange(do.call(rbind, (apply(data.frame(t(combn(model.var, 2))), 1, 
                                       function(x){ t <- ci.test( x[[1]], x[[2]] , data = data.qt.10[,model.var] );
                                                           data.frame( rel1 = paste(x[[2]], x[[1]], sep="_|_"),
                                                                       stat = round(t$statistic,2), 
                                                                       p.value = t$p.value, method = t$method, 
                                                                       label = t$data.name, alternative = t$alternative )
}))), p.value)


# conditional independence tests
# Each variable is conditionally independent of all its nondescendants
# in the graph given the value of all its parents.
d.cind <- arrange(do.call(rbind, (apply(d, 1, function(x){ t <- ci.test( x[[1]], x[[2]], x[[3]] , test = "jt", 
                                                                     data = data.qt.10[,model.var] );
                         data.frame( rel1 = paste(x[[3]], x[[1]], sep="->"), rel2 = paste(x[[3]], x[[2]], sep="->"), 
                                     stat = round(t$statistic,2), p.value = t$p.value, method = t$method, label = t$data.name, alternative = t$alternative )
                         }))), p.value)

d.cind$adjusted.p.val <- format.pval(p.adjust(d.cind$p.value, method = "fdr"), digits=3)
#subset(dd, p.value < 0.05)
d.cind$p.value <- format.pval(d.cind$p.value,digits = 3)


# nodes that are not connected should be cond.indep (P>0.05)
subset(d.cind, (!rel1 %in% posterior.links) | (!rel2 %in% posterior.links)  )

subset(d.cind, !rel1 %in% posterior.links | !rel2 %in% posterior.links  )

# nodes that are connected may or may not be cond dependent (P>0.05)
subset(d.cind, (rel1 %in% posterior.links & rel2 %in% posterior.links) )



# ------------------------------------------------------------------------------------------------ #
#                                     Format tables for LaTeX                                      #
# ------------------------------------------------------------------------------------------------ #

tex.dd <- verify.ci.qt10[,c("label", "JT","p.value")]
colnames(tex.dd) <- c("Condition","JT", "p-value")

require(stargazer)
#rownames(dd) <- NULL
dd.tex <- stargazer(tex.dd, summary = F, digits = 2, perl = T, 
                    title = "Conditional independence of determinants of protein production rate fold-change P",
                    out = "~/Dropbox/PhD/thesis/TeX/items/tables/chapter5/table_results_conditional_independence_draft.tex"  )


# ------------------------------------------------------------------------------------------------ #
#                                          Priors (thresholds)                                     #
# ------------------------------------------------------------------------------------------------ #

priors <- p.dist <- ddply(melt.bn.data.th, .(variable), summarise, state = names(table(value)), p = table(value)/length(value))

g=ggplot(priors, aes(x=state, y = p, fill=state)) + 
  geom_bar(width=0.75, stat = "identity" ) +
  facet_wrap( ~ variable, ncol=1) + 
  theme_minimal() + scale_fill_manual(values=c(DOWN,"gray60",UP)) +
  geom_text(aes(label=round(p,2), colour = ifelse(p > .57, "white","black") ), size=3.75, y = 0.5) +
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "none" ) + 
  scale_y_continuous(labels=prettyNum, breaks=c(0,0.25,0.5,.75)) +
  scale_color_identity()
ggsave(g, filename=paste0(PATH,"part3/BN_distributions.th.pdf"),width = 2.75, height=5.65,  useDingbats=F)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                   ARC STRENGTHS                                                                                              #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

a.strengths <- arc.strength(x= bn.gs.qt.10, data = data.qt[,model.var] )

b.strengths <- boot.strength(data.qt[,model.var], algorithm = "hc")
b.strengths_bn.gs.qt.10 <- merge(b.strenghts, bn.gs.qt.10$arcs, by = c("from","to"))



get.magnitude.impact <- function( BN.fit, x0 = '50-60', from = "e", to = "P", RND = 2){
  
  probas <- melt(BN.fit[[ to ]]$prob)
  
  query <- paste0(from, '=="', x0, '"')
  converse.query <- paste0(from, '!="', x0, '"')
  
  reference.set <- subset( probas, eval(parse(text=query)) )
  comparison.set <- subset( probas, eval(parse(text=converse.query)) )
  
  other.variables <- setdiff( colnames(probas), c(from,"value") )
  
  comparison.table <- merge(reference.set, comparison.set, by = other.variables)
  return( data.frame( from = from, to = to, MI = round(with( comparison.table,  max(abs(value.x-value.y))),RND ))  )
}

get.magnitude.impact(BN.fit = BN.GS.qt.10, from="a", to = "P")
get.magnitude.impact(BN.fit = BN.GS.qt.10, from="i", to = "P")
get.magnitude.impact(BN.fit = BN.GS.qt.10, from="e", to = "P")

get.magnitude.impact(BN.fit = BN.GS.qt.10, from="tAI", to = "i")
get.magnitude.impact(BN.fit = BN.GS.qt.10, from="tAI", to = "e")



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                   VISUALISE OBSERVATIONS                                                                                               #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

# approach based on scaling (tertiles)
melt.bn.data.qt <- melt(data.qt, id.vars=c("name","label","prod","ID"))
g=ggplot(data=melt.bn.data.qt, aes(x=factor(variable, levels=rev(levels(melt.bn.data.qt$variable))),
                                   y=ID)) + geom_tile(aes(fill=value), width=0.95) + 
  scale_fill_manual(values=c(DOWN, "gray70", UP)) + 
  theme_minimal() + coord_flip() +
  theme(legend.position="bottom", axis.text.x=element_text(size=6), axis.title.x=element_text(size=9)) + labs(x="",y="observations", fill="variation")

ggsave(plot=g, filename = paste0(PATH,"part3/BN_data.qt.pdf"),  width = 6, height=4,useDingbats=F  )


# approach based on scaling (deciles)
melt.bn.data.qt.10 <- melt(data.qt.10, id.vars=c("name","label","prod","ID"))
g=ggplot(data=melt.bn.data.qt.10, aes(x=factor(variable, levels=rev(levels(melt.bn.data.qt.10$variable))),
                                      y=ID)) + 
  geom_raster(aes(fill=as.numeric(value)), width=0.95) + 
  #scale_fill_manual(values=c(DOWN, "gray70", UP)) + 
  scale_fill_gradient2(low=DOWN, mid = "white", high=UP, midpoint = 5.5) +
  theme_minimal() + coord_flip() +
  theme(legend.position="none", 
        axis.text.x=element_text(size=6), 
        axis.title.x=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) + 
  labs(x="",y="observations")

ggsave(plot=g, filename = paste0(PATH,"part3/BN_data.qt10.pdf"),  width = 6, height=4,useDingbats=F  )


# Threshold based approach (no scaling, it's all based on the actual observed values -- a gene with neg.var in tAI has a decreasing tAI, and if p(tAI) is biased towards decreases, then so be it.)
melt.bn.data.th <- melt(data.th, id.vars=c("name","label","prod","ID"))
melt.bn.data.th$label <- factor(melt.bn.data.th$label, levels = 
                                  paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n")
)
g=ggplot(data=melt.bn.data.th, aes(y=factor(variable, levels=rev(levels(melt.bn.data.th$variable))),
                                   x=ID)) + geom_raster(aes(fill=value), width=0.95) + 
  scale_fill_manual(values=c(DOWN, "gray70", UP)) + 
  theme_minimal() + facet_wrap( ~ label, scales = "free_x", nrow=1) +
  theme(legend.position="bottom", axis.text.x=element_text(size=6), axis.title.x=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  labs(x="",y="observations", fill="variation")



ggsave(plot=g, filename = paste0(PATH,"part3/BN_data.th.pdf"),  width = 6, height=4,useDingbats=F  )





melt.bn.data.th.2 <- melt(data.th.2, id.vars=c("name","label","prod","ID"))
melt.bn.data.th.2$label <- factor(melt.bn.data.th.2$label, levels = 
                                    paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n")
)
g=ggplot(data=melt.bn.data.th.2, aes(y=factor(variable, levels=rev(levels(melt.bn.data.th.2$variable))),
                                     x=ID)) + geom_raster(aes(fill=value), width=0.95) + 
  scale_fill_manual(values=c(DOWN, UP)) + 
  theme_minimal() + facet_wrap( ~ label, scales = "free_x", nrow=1) +
  theme(legend.position="bottom", axis.text.x=element_text(size=6), axis.title.x=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  labs(x="",y="observations", fill="variation")


# ----- Up vs. not up

melt.bn.data.th.up <- melt(data.th.up, id.vars=c("name","label","prod","ID"))
melt.bn.data.th.up$label <- factor(melt.bn.data.th.up$label, levels = 
                                     paste( rep(c("diauxic","oxidative","osmotic","temperature"), times = 2 ), paste(rep(c(20,120), each = 4 ),"min",sep=" "), sep="\n")
)
g=ggplot(data=melt.bn.data.th.up, aes(y=factor(variable, levels=rev(levels(melt.bn.data.th.up$variable))),
                                      x=ID)) + geom_raster(aes(fill=value), width=0.95) + 
  scale_fill_manual(values=c(DOWN, UP)) + 
  theme_minimal() + facet_wrap( ~ label, scales = "free_x", nrow=1) +
  theme(legend.position="bottom", axis.text.x=element_text(size=6), axis.title.x=element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  labs(x="",y="observations", fill="variation")


# ------------------------------------------------------------------------------------------------ #
#                                  CPTs shown as tiles (thresholds)                                #
# ------------------------------------------------------------------------------------------------ #

cpt <- probas.BN.GS.th$P
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
g <- ggplot( data=cpt, aes_string(x=colnames(cpt)[[1]], fill="Freq", y= colnames(cpt)[[2]]) ) + 
  geom_tile() + facet_wrap( as.formula(paste("~",given.that,sep=" "))  ) +
  geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  scale_color_identity() + theme(legend.position="none") + coord_equal() 

ggsave(g, filename=paste0(PATH,"part3/BN_P_CPT.th.pdf"), width=8.2, height=7.75,  useDingbats=F)


cpt <- probas.BN.GS.th$e
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
ggplot( data=cpt, aes_string(x=colnames(cpt)[[1]], fill="Freq", y= colnames(cpt)[[2]]) ) + 
  geom_tile() + facet_wrap( as.formula(paste("~",given.that,sep=" ")), ncol=3  ) +
  geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  scale_color_identity() + theme(legend.position="none") + coord_equal() 

cpt <- probas.BN.GS.th$i
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
g <- ggplot( data=cpt, aes_string(x=colnames(cpt)[[1]], fill="Freq", y= colnames(cpt)[[2]]) ) + 
  geom_tile() + facet_wrap( as.formula(paste("~",given.that,sep=" "))  ) +
  geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  scale_color_identity() + theme(legend.position="none") + coord_equal() 

ggsave(g, filename=paste0(PATH,"part3/BN_i_CPT.th.pdf"), width=8.3, height=3,  useDingbats=F)




# ------------------------------------------------------------------------------------------------ #
#                                  CPTs shown as tiles (deciles)                                   #
# ------------------------------------------------------------------------------------------------ #

cpt <- probas.BN.GS.qt.10$P
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
g <- ggplot( data=cpt, aes_string(x=colnames(cpt)[[2]], fill="Freq", y= colnames(cpt)[[1]]) ) + 
  geom_raster() + 
  geom_abline(slope=1, intercept=0, color="gray30") +
  facet_wrap( as.formula(paste("~",given.that,sep=" "))  ) +
  #geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  labs(y=expression(Pi)) +
  scale_fill_gradient(low = "white",high = "#026c8f") +
  scale_color_identity() + 
#   theme_grey() +
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "bottom",
          axis.text.y = element_blank(), axis.text.x = element_blank()
          ) +
  coord_equal() 

ggsave(g, filename=paste0(PATH,"part3/BN_P_CPT.qt10.pdf"), width=12.6, height=12.6,  useDingbats=F)


cpt <- probas.BN.GS.qt.10$e
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
ggplot( data=cpt, aes_string(x=colnames(cpt)[[1]], fill="Freq", y= colnames(cpt)[[2]]) ) + 
  geom_tile() +
  geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  scale_color_identity() + 
  #theme(legend.position="none", axis.text.y = element_blank(), axis.text.x = element_blank() ) +
  coord_equal() 



cpt <- probas.BN.GS.qt.10$i
given.that <- paste( setdiff( colnames(cpt), c("Freq", colnames(cpt)[1:2]) ), collapse=" + ")
g <- ggplot( data=cpt, aes_string(x=colnames(cpt)[[1]], fill="Freq", y= colnames(cpt)[[2]]) ) + 
  geom_tile() + 
  geom_text( aes(label = Freq, colour = ifelse(Freq>0.7, "black", "white"), size = 3 )) +
  labs(title=given.that) +
  scale_color_identity() + theme(legend.position="none") + coord_equal() 

ggsave(g, filename=paste0(PATH,"part3/BN_i_CPT.th.pdf"), width=8.3, height=3,  useDingbats=F)



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
#                                                                                         DISTRIBUTIONS                                                                                                  #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #


g1 <- ggplot( data = source.data, aes(x=s_log2.mRNA_abundance) ) + 
  geom_histogram( aes(fill = ifelse( quant.mRNAab %in% "0-33", UP, ifelse( quant.mRNAab %in% "67-100", DOWN, "gray60" ) )), 
                  binwidth= range(source.data$log2.mRNA_abundance)[2]/90 ) +  
  scale_fill_identity() +
  scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +  scale_y_continuous( breaks=c(0,500,1000)) +
  labs(x="mRNA abundance fold change\n(log2, scaled across conditions)", y="") + 
  theme_bw() 

g2 <- ggplot( data = source.data, aes(x=s_scaled.change.tAI) ) + 
  geom_histogram( aes(fill = ifelse( quant.tAI %in% "0-33", UP, ifelse( quant.tAI %in% "67-100", DOWN, "gray60" ) )), 
                  binwidth= range(source.data$s_scaled.change.tAI, na.rm = T)[2]/150 ) +  
  scale_fill_identity() +
  scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +  scale_y_continuous( breaks=c(0,250,500)) +
  labs(x="variation in tAI\n(scaled across conditions)", y="") +
  theme_bw() 

g3 <- ggplot( data = source.data, aes(x=s_variation_speed) ) + 
  geom_histogram( aes(fill = ifelse( quant.speed %in% "0-33", UP, ifelse( quant.speed %in% "67-100", DOWN, "gray60" ) )), 
                  binwidth= range(source.data$s_variation_speed)[2]/350 ) +  
  scale_fill_identity() +
  scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +  scale_y_continuous( breaks=c(0,500,1000)) +
  labs(x="variation in elongation speed\n(scaled across conditions)", y="") +
  theme_bw() 

g4 <- ggplot( data = source.data, aes(x=s_variation_initiation_frequency) ) + 
  geom_histogram( aes(fill = ifelse( quant.freq %in% "0-33", UP, ifelse( quant.freq %in% "67-100", DOWN, "gray60" ) )), 
                  binwidth= range(source.data$s_variation_initiation_frequency)[2]/250 ) +  
  scale_fill_identity() +
  scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +  scale_y_continuous( breaks=c(0,500,1000)) +
  labs(x="initiation frequency fold change\n(log2, scaled across conditions)", y="") +
  theme_bw() 

g5 <- ggplot( data = source.data, aes(x=s_protein_production_rate) ) + 
  geom_histogram( aes(fill = ifelse( quant.prod %in% "0-33", UP, ifelse( quant.prod %in% "67-100", DOWN, "gray60" ) )),
                  binwidth= range(source.data$s_protein_production_rate)[2]/290 ) +  
  scale_fill_identity() +
  scale_x_continuous(labels=prettyNum, limits=c(-5,5)) + scale_y_continuous( breaks=c(0,500,1500,2500)) +
  labs(x="protein production rate fold change\n(scaled across conditions)", y="") +
  theme_bw() 


pdf(file = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part3/BN_scaling_categorisation.pooled.pdf", useDingbats = F, height=10, width=2.5, compress = T )
require(gridExtra)
grid.draw( rbind_gtable_max(ggplotGrob(g1), ggplotGrob(g2) , ggplotGrob(g3), ggplotGrob(g4), ggplotGrob(g5)) )
dev.off()






# -------------------------------------------------------------------------------------------------------


# read the BN graph as a string
modelstring(bn.hc)

# list undirected edges
undirected.arcs(up.bn.gs)
undirected.arcs(down.bn.gs)

# get the nodes of the graph
nodes(bn.hc)

#bn.hc.AE <- set.arc(bn.hc, "A", "E")



dsep
require(ggm)

dSep( amat =  get.adjacency(graph.data.frame(bn.gs$arcs), sparse = F), first = "a", second = "tAI", cond = c("i","e"))
dSep( amat =  get.adjacency(graph.data.frame(bn.gs$arcs), sparse = F), first = "P", second = "e", cond = c("tAI"))

dSep(DAG(y ~ x, x ~ z), first="y", second="z", cond = "x")

