
aa.data <- ddply(anticodon.master.table, .(experiment, time, aa), 
                   summarise,
                   total.demand = sum(demand.mRNAab),
                   total.supply = sum(total.tRNAab)
                   )
aa.data.no.CAA <- ddply(subset(anticodon.master.table, anticodon!="CAA"), .(experiment, time, aa), 
                 summarise,
                 total.demand = sum(demand.mRNAab),
                 total.supply = sum(total.tRNAab)
)

subset(aa.data, time == 0)

lm_aa = function(df=aa.data, df2=aa.data.no.CAA, a="total.supply", b="total.demand"){
  formula   <- as.formula(paste(a," ~ ",b,sep = ""))
  formula.2 <- as.formula(paste("~",a,"+",b,sep=" "))
  m1 = lm( formula, df);
  m2 = lm( formula, df2);
  pears = cor.test( formula = formula.2, data = df2, method = "pearson",exact = F);
  line1 <- substitute(italic(r)~"="~pears~" | "~alpha~"="~coeff, 
                      list( pears = format(pears$estimate, digits = 2),
                            coeff = format(summary(m2)$coefficients[[2]], digits = 2)
                      )
  )
  line2 <- substitute(italic(R)^2~"="~r2~" | "~italic(R)[+CAA]^2~"="~r2.caa,
                      list( r2.caa     = format(summary(m1)$r.squared, digits = 2), 
                            r2         = format(summary(m2)$r.squared, digits = 2)
                      )
  )
  as.character( as.expression( paste( c('atop(', line1, ',', line2,')'), collapse=" ") ) )
}

R2.aa <- ddply( aa.data, .(experiment, time), function(x){ lm_aa(x, df2 = subset(aa.data.no.CAA, experiment == unique(x$experiment) & time == unique(x$time)) ) })

ggplot( data = aa.data, aes(x=total.demand/10^5, y=total.supply/10^5)
         ) + 
  stat_smooth(method="lm") +
  geom_point( size = 3.25, color="gray80", alpha=0.9 ) + 
  geom_text( size = 1.5, aes(label=aa) ) +
  geom_text( data = R2.aa, aes_string( label = expression(V1) ), parse =T, x=-Inf, y=+Inf, hjust = -0.1, vjust = 1.5, size =3 ) +
  scale_x_sqrt( labels = prettyNum ) + coord_fixed() +
  scale_y_sqrt( labels = prettyNum ) +
  labs( x = expression(atop("aa demand",atop("x100,000",""))), 
        y = expression(atop("aa supply",atop("tRNA abundance, x100,000",""))),
        title="Balance between amino acid demand and supply") +
  facet_wrap( time ~ experiment, scale="free_x", nrow =2 ) + 
  theme_gul + theme(legend.position="bottom") +
  facet_wrap( time ~ experiment, nrow=2)


#ketchup
# figure out whether isoacceptor tRNA correlate or not
get.isoacceptor.correlations <- function(M=M.tRNAab$all.together){
aa.with.isoacceptor.tRNAs <- names(which(table(subset(rosetta.anticodons, select=c(aa,anticodon))$aa) > 1))
m <- cor(t(M))
require(gdata)
upperTriangle(m, diag = T) <- NA # remove duplicates, as well as diagonal (Pearson = 1 by definition here)
correlations.per.aa <- melt( m )
colnames(correlations.per.aa) <- c("tRNA.1","tRNA.2", "Pearson")
correlations.per.aa <- subset(correlations.per.aa, !is.na(Pearson))
require(tidyr)
correlations.per.aa <- separate(correlations.per.aa, col = "tRNA.1", into = c("anticodon.1","aa.1"), sep="-" )
correlations.per.aa <- separate(correlations.per.aa, col = "tRNA.2", into = c("anticodon.2","aa.2"), sep="-" )
correlations.per.aa$ID <- apply(correlations.per.aa, 1, function(x) paste(sort(x[c("anticodon.1","anticodon.2")]),collapse="-") )
correlations.per.aa <- correlations.per.aa[, !colnames(correlations.per.aa) %in% c("anticodon.1","anticodon.2") ]
correlations.per.aa$isoacceptor <- ifelse(correlations.per.aa$aa.1 == correlations.per.aa$aa.2, T,F)
correlations.per.aa[,c("aa.1","aa.2")] <- t(apply(correlations.per.aa[,c("aa.1","aa.2")], 1, function(x) sort(x)))

duplications <- correlations.per.aa[,c("aa.2","aa.1","Pearson","ID","isoacceptor")]
colnames(duplications) <- colnames(correlations.per.aa)
correlations.per.aa <- unique(rbind( correlations.per.aa,duplications))
correlations.per.aa <- subset(correlations.per.aa, 
                              aa.1 %in% aa.with.isoacceptor.tRNAs & 
                                aa.2 %in% aa.with.isoacceptor.tRNAs )
return(correlations.per.aa)
}

get.isoacceptor.plot <- function(data=correlations.per.aa){
  medians <- ddply(subset(data, isoacceptor==T), .(aa.1, isoacceptor), summarise, Pearson = median(Pearson))
 g <-  ggplot(data= data, aes(x=isoacceptor,y=Pearson) ) + 
    geom_hline( yintercept=1, lty=1, color="gray" ) +
    geom_hline( yintercept=0.75, lty=3 ) +
    geom_hline( yintercept=0, lty=3 ) +
    geom_hline( yintercept=-0.75, lty=3 ) +
    geom_hline( yintercept=-1, lty=1, color="gray" ) +
    geom_dotplot(aes(color=aa.1, fill=aa.1),binaxis = "y", stackdir = "center", binwidth=0.035) +
    geom_boxplot( data = subset(data, isoacceptor==F), notch=T, width=0.5,alpha=0.6, outlier.size = 0) +
    geom_point( data = medians,
                shape = 5, size =4
    ) + ylim(-1,1) +
    facet_wrap( ~ aa.1) + 
    labs(y="Pearson's correlation coefficient", title="Correlation of tRNA expression patterns") + theme(legend.position="none")
  return(g)
}

correlations.per.aa <- lapply( setNames(names(M.tRNAab[-c(2,7)]),names(M.tRNAab[-c(2,7)])), function(x){ 
                          data <- get.isoacceptor.correlations(M = M.tRNAab[[x]] )
                          get.isoacceptor.plot(data = data) %+% labs(title=x)
  }
)

d.correlations <- get.isoacceptor.correlations(M = M.tRNAab[[1]] )
g <- correlations.per.aa$all.together
correlations.per.aa$diauxic  
correlations.per.aa$ox
correlations.per.aa$osm
correlations.per.aa$temp

ggsave(plot = g, filename = paste0(PATH,"part1/aa_correlation.iscoacceptors.pdf"), useDingbats=F, width=6.93, height=9.39)

subset(get.isoacceptor.correlations(M = M.tRNAab$all.together), isoacceptor==T & aa.1 == "Gln")
d <- as.data.frame(t( M.tRNAab$all.together[c("CUG-Gln","UUG-Gln"),]))
colnames(d) <- c("CUG","UUG")
medians <- ddply(subset(get.isoacceptor.correlations(M = M.tRNAab$all.together), isoacceptor==T), .(aa.1, isoacceptor), summarise, Pearson = round(median(Pearson),2))

g <- ggplot( data = d, aes_string(x="CUG", y="UUG") ) + 
  geom_point(size = 3) + coord_fixed() +
  stat_smooth(method="lm")
ggsave(plot = g, filename = paste0(PATH,"part1/aa_correlation.iscoacceptors_medaillon.pdf"), useDingbats=F,
       width=2.4, height=2.4
       )


# check if the level of correlation of isoacceptor tRNAs is biased by amino acid demand
test.correlation <- merge( subset(d.correlations, isoacceptor == T), subset(aa.data, time == 0), by.x="aa.1", by.y="aa")
test.correlation <- subset(test.correlation, select=c(aa.1, Pearson, ID, total.demand,total.supply))

ggplot( data = test.correlation, aes(x= rank(total.demand), y = Pearson)) + geom_point(aes(fill=Pearson > 0.7), pch =21, size =4 ) + geom_text(aes(label=aa.1), size=3, vjust=2)
