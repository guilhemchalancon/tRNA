# ----------------- Additive Bayesian Networks --------
setwd("~/Documents/MRC16/")
source("scripts/commands.R")
source("scripts/SMoT.R")

require(abn) # Creation and analysis of Additive Bayesian Networks
require(igraph)
 
# # # ----------------- Tutorial example
#                                           mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
#                                           
#                                           ## setup distribution list for each node
#                                           mydists<-list(b1="binomial",
#                                                         b2="binomial",
#                                                         b3="binomial",
#                                                         g1="gaussian",
#                                                         b4="binomial",
#                                                         p2="poisson",
#                                                         p4="poisson"
#                                           );
#                                           
#                                           ## define model
#                                           mydag<-matrix(data=c(
#                                             0,0,1,0,0,0,0, # b1<-b3
#                                             1,0,0,0,0,0,0, # b2<-b1
#                                             0,0,0,0,0,0,0, #
#                                             0,0,0,0,1,0,0, # g1<-b4
#                                             0,0,0,0,0,0,0, #
#                                             0,0,0,0,0,0,0, #
#                                             0,0,0,0,0,0,0 #
#                                           ), byrow=TRUE,ncol=7);
#                                           colnames(mydag)<-rownames(mydag)<-names(mydat);
#                                           
#                                           myres.c<-fitabn(dag.m=mydag,data.df=mydat,data.dists=mydists, create.graph = T);
#                                           
#                                           
#                                           myres.c
#                                           print(myres.c$mlik);
#                                           plot(myres.c$graph)
# 

# ----------------- Specification of model and data

# 1 # data frame
      study.translation.kinetics <- read.table("results/master tables/study.translation.kinetics.txt", header=T)



      # I am going to pool together data from various experimental conditions for the ABN analysis: I must first scale data 

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
      source.data$quant.mRNAab  <- quantile.it(source.data$s_log2.mRNA_abundance ) # already there      
      source.data$quant.tAI   <- quantile.it( source.data$s_scaled.change.tAI )
      source.data$quant.speed <- quantile.it( source.data$s_variation_speed )
      source.data$quant.freq <- quantile.it( source.data$s_variation_initiation_frequency )
      source.data$quant.prod <- quantile.it( study.translation.kinetics$ratio_events )

# no really not!
# source.data <- ddply(source.data, .(experiment, time), mutate, 
#                      quant.mRNAab = quantile.it( scale(log2.mRNA_abundance) ),
#                      quant.tAI   = quantile.it( scale(variation_tAI) ),
#                      quant.speed     = quantile.it( scale(variation_speed) ),
#                      quant.freq = quantile.it( scale(log2(1/ratio_initiation)) ),
#                      s_protein_production_rate = scale(global.protein_synthesis.rate)
# )

      source.data <- with(source.data, expr = data.frame(source.data,
                           quant.mRNAab.pooled = quantile.it( log2.mRNA_abundance ),
                           quant.tAI.pooled   = quantile.it( variation_tAI ),
                           quant.speed.pooled = quantile.it( variation_speed ),
                           quant.freq.pooled = quantile.it( log2(1/ratio_initiation) ),
                           prod = scale(global.protein_synthesis.rate)
        )
      )

# wow
subset(source.data, n.events > 140000)
      # check scaling worked
#        ggplot(source.data, aes(log2.mRNA_abundance, y =s_log2.mRNA_abundance ) ) + geom_point() 
#        ggplot(source.data, aes(scaled.change.tAI, y =s_scaled.change.tAI ) ) + geom_point() 
#        ggplot(source.data, aes(variation_speed, y =s_variation_speed ) ) + geom_point() 
#        ggplot(source.data, aes(variation_initiation_frequency, y =s_variation_initiation_frequency ) ) + geom_point() 



# extract data subset with variables of interest
data <- subset( source.data, select = c(name, label, quant.mRNAab, quant.tAI, quant.freq, quant.speed, quant.prod, prod )) 


#model.var <- c("v.mRNAab","v.tAI","v.speed","v.freq","prod")
# SEL  <- c("0-20","80-100") # more conservative

# select variables that will be used for the BN learning
model.var <- c("v.mRNAab","v.tAI","v.speed","v.freq","v.prod")
SEL  <- c("0-20","20-40","60-80","80-100")

# Bayesian Network Input Data (restrict to only those observations that are in the right space)
abn.data <- subset( data,  quant.mRNAab %in% SEL & quant.tAI %in% SEL & quant.freq %in% SEL & quant.speed %in% SEL & quant.prod %in% SEL) 
colnames(abn.data)[grep("quant",colnames(abn.data))] <- c("v.mRNAab", "v.tAI","v.freq","v.speed","v.prod")


nrow(abn.data) # 10,201 observations across the data




# encode quantile categories in "up" and "down"
#abn.data[,c("v.mRNAab","v.tAI","v.speed","v.freq")] <- data.frame(sapply( abn.data[,c("v.mRNAab","v.tAI","v.speed","v.freq")], function(x){ ifelse(x == "0-20", "up", "down" )} ))
abn.data[,c("v.mRNAab","v.tAI","v.speed","v.freq","v.prod")] <- data.frame(sapply( abn.data[,c("v.mRNAab","v.tAI","v.speed","v.freq","v.prod")], 
                                                                                   function(x){ ifelse(x %in% c("0-20","20-40"), "up", "down" )} ))

# Nodes in the Bayesian network
abn.data$label <- as.character(abn.data$label)
abn.data$ID <- 1:nrow(abn.data)
write.table(abn.data, "results/ABN.data.txt",quote=F, row.names=F)


table(abn.data$label) # not in all conditions, by the way
table(abn.data$v.tAI) # not in all conditions, by the way

# before selecting data, first plot their distributions

DOWN <- "#185665"
UP <- "#A36A23"

# CONTROL DISTRIBUTIONS
      g1 <- ggplot(data=source.data, aes(x=label, y=s_log2.mRNA_abundance)) + geom_jitter(aes(color= ifelse( quant.mRNAab %in% "0-20", UP, ifelse( quant.mRNAab %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation mRNA abundance")  +theme(axis.text.x=element_blank())
      
      g2 <- ggplot(data=source.data, aes(x=label, y=s_scaled.change.tAI)) + geom_jitter(aes(color= ifelse( quant.tAI %in% "0-20", UP, ifelse( quant.tAI %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation tAI") # + theme(axis.text.x=element_blank())
      
      g3 <- ggplot(data=subset(source.data ) , aes(x=label, y=s_variation_speed)) + geom_jitter(aes(color= ifelse( quant.speed %in% "0-20", UP, ifelse( quant.speed %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation elongation speed") + theme(axis.text.x=element_blank())
      
      g4 <- ggplot(data=source.data, aes(x=label, y=s_variation_initiation_frequency)) + geom_jitter(aes(color= ifelse( quant.freq %in% "0-20", UP, ifelse( quant.freq %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation initiation frequency") + theme(axis.text.x=element_blank())
      
      
            pdf(file = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part3/BN_scaling_categorisation.pooled.pdf", useDingbats = F, height=10, width=10, compress = T )
            require(gridExtra)
            grid.draw( rbind_gtable_max(cbind_gtable_max( ggplotGrob(g1), ggplotGrob(g2)) , cbind_gtable_max(ggplotGrob(g3), ggplotGrob(g4))  ) )
            dev.off


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

      
      g1 <- ggplot(data=source.data, aes(x=label, y= s_log2.mRNA_abundance)) + geom_jitter(aes(color= ifelse( quant.mRNAab.pooled %in% "0-20", UP, ifelse( quant.mRNAab %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation mRNA abundance")  +theme(axis.text.x=element_blank())
      
      g2 <- ggplot(data=source.data, aes(x=label, y= s_variation_tAI)) + geom_jitter(aes(color= ifelse( quant.tAI.pooled %in% "0-20", UP, ifelse( quant.tAI %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation tAI") # + theme(axis.text.x=element_blank())
      
      g3 <- ggplot(data=subset(source.data ) , aes(x=label, y=variation_speed)) + geom_jitter(aes(color= ifelse( quant.speed.pooled %in% "0-20", UP, ifelse( quant.speed %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation elongation speed") + theme(axis.text.x=element_blank())
      
      g4 <- ggplot(data=source.data, aes(x=label, y=s_variation_initiation_frequency)) + geom_jitter(aes(color= ifelse( quant.freq %in% "0-20", UP, ifelse( quant.freq %in% "80-100", DOWN, "gray30" ) ) )) + 
        geom_boxplot(notch=T, outlier.size = F, alpha=0.8) + scale_color_identity() + labs(x="",y="variation initiation frequency") + theme(axis.text.x=element_blank())
      

      pdf(file = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part3/BN_scaling_categorisation.pdf", useDingbats = F, height=10, width=10, compress = T )
      require(gridExtra)
      grid.draw( rbind_gtable_max(cbind_gtable_max( ggplotGrob(g1), ggplotGrob(g2)) , cbind_gtable_max(ggplotGrob(g3), ggplotGrob(g4))  ) )
      dev.off


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

      
      g1 <- ggplot( data = source.data, aes(x=s_log2.mRNA_abundance) ) + 
        geom_histogram( aes(fill = ifelse( quant.mRNAab %in% "0-20", UP, ifelse( quant.mRNAab %in% "80-100", DOWN, "gray60" ) )), 
                        binwidth= range(source.data$log2.mRNA_abundance)[2]/90 ) +  
        scale_fill_identity() +
        scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +
        labs(x="mRNA abundance fold change\n(log2, scaled across conditions)", y="") + 
        theme_bw() 
      
      g2 <- ggplot( data = source.data, aes(x=s_scaled.change.tAI) ) + 
        geom_histogram( aes(fill = ifelse( quant.tAI %in% "0-20", UP, ifelse( quant.tAI %in% "80-100", DOWN, "gray60" ) )), 
                        binwidth= range(source.data$s_scaled.change.tAI, na.rm = T)[2]/150 ) +  
        scale_fill_identity() +
        scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +
        labs(x="variation in tAI\n(scaled across conditions)", y="") +
        theme_bw() 
      
      g3 <- ggplot( data = source.data, aes(x=s_variation_speed) ) + 
        geom_histogram( aes(fill = ifelse( quant.speed %in% "0-20", UP, ifelse( quant.speed %in% "80-100", DOWN, "gray60" ) )), 
                        binwidth= range(source.data$s_variation_speed)[2]/350 ) +  
        scale_fill_identity() +
        scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +
        labs(x="variation in elongation speed\n(scaled across conditions)", y="") +
        theme_bw() 
      
      g4 <- ggplot( data = source.data, aes(x=s_variation_initiation_frequency) ) + 
        geom_histogram( aes(fill = ifelse( quant.freq %in% "0-20", UP, ifelse( quant.freq %in% "80-100", DOWN, "gray60" ) )), 
                        binwidth= range(source.data$s_variation_initiation_frequency)[2]/250 ) +  
        scale_fill_identity() +
        scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +
        labs(x="initiation frequency fold change\n(log2, scaled across conditions)", y="") +
        theme_bw() 
      
      g5 <- ggplot( data = source.data, aes(x=s_protein_production_rate) ) + 
        geom_histogram( binwidth= range(source.data$s_protein_production_rate)[2]/250 ) +  
        scale_fill_identity() +
        scale_x_continuous(labels=prettyNum, limits=c(-5,5)) +
        labs(x="protein production rate fold change\n(scaled across conditions)", y="") +
        theme_bw() 


pdf(file = "~/Dropbox/PhD/thesis/TeX/items/figures/chapter5/part3/BN_distributions.pdf", useDingbats = F )
  require(gridExtra)
  grid.draw( rbind_gtable_max( ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), ggplotGrob(g4), ggplotGrob(g5)  ) )
dev.off()




      
      
    


# 2 # specification of the Directed Acyclic Graph
tAI.dag<-matrix(data=c(
  0,0,0,0,0, # v.mRNAab exogenous
  0,0,0,0,0, # v.tAI exogenous
  0,1,0,0,0, # v.speed <- v.tAI
  0,0,0,0,0, # v.freq (so far) exogenous
  1,0,1,1,0  # v.prod < - v.speed, v.prod <- v.freq
), byrow=TRUE,ncol=length(model.var));
colnames(tAI.dag)<-rownames(tAI.dag)<-model.var;

# 3 # prior distributions of all the model variables
distributions<-list(
                v.mRNAab="binomial",
                v.tAI="binomial",
                v.speed="binomial",
                v.freq="binomial",
                prod="gaussian"
);


# 4 # defined banned edges and retained edges in the DAG search
ban <-matrix(data=c(
  0,1,1,1,1, # v.mRNAab exogenous
  1,0,1,1,1, # v.tAI exogenous
  0,0,0,0,1, # v.speed <- v.tAI
  0,0,0,0,1, # v.freq (so far) exogenous
  0,0,0,0,0  # v.prod < - v.speed, v.prod <- v.freq
), byrow=TRUE,ncol=length(model.var));
colnames(ban)<-rownames(ban) <- model.var


retain <-matrix(data=c(
  0,0,0,0,0, # v.mRNAab exogenous
  0,0,0,0,0, # v.tAI exogenous
  0,1,0,0,0, # v.speed <- v.tAI
  0,0,0,0,0, # v.freq (so far) exogenous
  0,0,0,0,0  # v.prod < - v.speed, v.prod <- v.freq
), byrow=TRUE,ncol=length(model.var));
colnames(retain)<-rownames(retain) <- model.var

# also define each node maximum number of permitted parents 
max.parents<-list("v.mRNAab"=0, "v.tAI"=0,"v.speed"=3,"v.freq"=3,"prod"=4);



# 5 # fit additive bayesian network 
require(abn)
tAI.abn.results <- fitabn( dag.m = tAI.dag, data.df = abn.data[,model.var], data.dists = distributions[model.var], create.graph = T )
 #tAI.bn.results <- fitbn( dag.m = tAI.dag, data.df = abn.data[-c(1:2,8)], data.dists = distributions, create.graph = T )
plot(tAI.abn.results$graph)

plot( graph.adjacency( tAI.abn.results$graph@adjMat ) )




# 6 # in search of a good DAG
mycache<-buildscorecache(data.df=abn.data[,model.var],data.dists=distributions[model.var],
                         dag.banned=ban, dag.retained=retain, max.parents =max.parents);

# Most probable dag (based on maximum likelihood)
mp.dag<-mostprobable(score.cache=mycache);

max.likelihood.mostprobable  <- fitabn(dag.m=mp.dag,data.df=abn.data[,model.var],data.dists=distributions[model.var])$mlik

tAI.abn.results$mlik         # max. likelihood prior BM 
max.likelihood.mostprobable  # most probable BN's max. likelihood

myres<-fitabn(dag.m=mp.dag,data.df=abn.data[,model.var],data.dists=distributions[model.var],create.graph=TRUE);
plot(myres$graph)

# Heuristic search
heur.res2<-search.hillclimber(score.cache=mycache,num.searches=1000,seed=0,verbose=FALSE, dag.retained = mp.dag,
                              trace=F,timing.on=FALSE);
dev.off()

plot(graph.adjacency(t(heur.res2$consensus)))



