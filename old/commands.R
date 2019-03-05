# ---- Graphs settings ----

theme_gul <- theme_update(
  # axis
  axis.text =  element_text( color = "gray25", size = rel(0.75)),
  axis.ticks = element_line( color = "gray70"),
  axis.title.y=element_text(hjust=0.5, vjust=1, angle=90),
  axis.title.x=element_text(hjust=0.5, vjust=-0.05),
  # title
  plot.title = element_text(hjust=0.5, vjust=1.5),
  # background
  panel.background = element_rect(fill="gray94", color="gray25", size=0),
  panel.grid.major = element_line(color="gray85", size = rel(0.6)),#,
  panel.grid.minor = element_line(color="gray88", size = rel(0.4)),
  # strips
  strip.background = element_rect(fill="gray55", color="gray25",size=0.2),
  strip.text       = element_text(color="white", size = rel(0.85) )
)


copy2clipboard <- function( dat, col.names=T, row.names= F, ... ){
  
  if(class(dat)=="data.frame"){
   file <- pipe("pbcopy")
          write.table(dat, file = file, row.names = row.names, col.names = col.names, sep = "\t", quote = F)
   #close.connection(file)  
  }
  
  if(class(dat)=="character"){
   file <- pipe("pbcopy")
          cat(dat, file = file, sep = "\n")
   close(file)    
  }
}

# https://github.com/hadley/plyr/blob/master/R/round-any.r#L28
gul_percent <- function (x) {
  require(scales)
  if(x!=0){ percent(x)} else {
    paste0(comma(0), "%")  
  }
}

depercent <- function(x){
  as.numeric(gsub(as.character(x), pattern="^(.*)\\%$", replacement = "\\1"))/100
}

range01 <- function (x) {(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}

# consider including in YET
quantilize <- function(data, name, bin.number=5){
  require(gtools)
  bin.number <- round(bin.number)
  LABELS <- paste( round((1 -seq(1:bin.number)/bin.number)*100 ), round(100/bin.number) + round((1 -seq(1:bin.number)/bin.number)*100 ),sep="-" )
  data.quantiles <- data.frame(apply( data, 2, function(x) quantcut(x, q=seq(0,1,by=1/bin.number),na.rm=T, labels=LABELS, ordered_result=T) ))
  rownames(data.quantiles) <- rownames(data)
#  save(data.quantiles,file=paste("data/Rdata/",name,".quantiles_",bin.number,".Rda",sep=""),safe=T)
  return(data.quantiles)
}

# Quartile dispertion coefficient
qdc <- function(x){
  q <- quantile(x, na.rm = T)
  return( as.numeric((q[4] - q[2]) / (q[2] + q[4])) ) # (Q3-Q1)/(Q1+Q3)
}

# Automated
select.bins <- function(data,selection){
  output.list <- apply(data, 2, function(x) names(x[which(x %in% selection)]) )
  names(output.list) <- colnames(data)
  return(output.list)
}

# Can't figure out which R command already does this
split.by.bin <- function(data){
  output.list <- list()
  for(i in names(table(data))){
   output.list[[i]] <- names(data[which(data %in% i)])
  }
   return(output.list)
}

show.binning <- function(qdata,odata){
  output <- list()
  binning <- names(table(qdata[,1]))
  for(i in binning ){
    output[[i]] <- odata[which(qdata[,1] %in% i),1]  
  }
  boxplot(output,main=paste("Binning in ",length(binning),"-quantiles",sep=""),las=2,ylab="log transformed expression levels")
}


is.outliers <- function (x,na.rm=TRUE, h=3,...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- h * IQR(x, na.rm = na.rm)
  outliers <- ( x < (qnt[1] - H) | x > (qnt[2] + H) )
  return(outliers)
}



get.sym.matrix <- function(df=primer.alignments, value.var="identity.all", identifiers=c("id1","id2")){
  require(reshape2)
  df2 <- df[,c(rev(identifiers),value.var)]
  colnames(df2)[1:2] <- identifiers
  df <- rbind(df[,c(identifiers,value.var)],df2) 
  ordering  <- sort(unique(as.character(df[,identifiers[1]])))
  formula  <- paste(identifiers[1],identifiers[2],sep = " ~ ")
  identity.matrix <- (matrixify(dcast(df, as.formula(formula), value.var = value.var )  ))
  identity.matrix <- identity.matrix[ ordering , ordering  ]
  return(identity.matrix)
}


conf.int <- function( mu=3.006, sigma=0.707, n = 3, conf = 0.05 ){
  SE  <- sigma/sqrt(n) # std error 
  E <- qt(p = 1 - conf/2, df=n-1)*SE
}

gulqqFitPlot <- function (df, y="average", method="I", distribution="lognormal", ref = "label", 
                          xlab = "Predicted fold change", ylab = "Observed fold change", 
                          fat = FALSE, ...) 
{
  require(extremevalues)
  L <- getOutliers(y = df[,y], method = method, distribution=distribution, ...)
  if (L$method == "Method I") {
    X <- getOutliersII(df[,y], FLim = c(L$Fmin, L$Fmax), distribution = L$distribution)
    L$residuals <- X$residuals
    L$yMin <- X$yMin
    L$yMax <- X$yMax
  }
  df$yHat <- switch(L$distribution, normal = df[,y] - L$residuals, 
                 exponential = df[,y] - L$residuals, exp(log(df[,y]) - L$residuals))
  lgAxis <- switch(L$distribution, normal = "", exponential = "", "xy")
  vLeft  <- head( df$yHat[df[,y] == L$yMin], 1)
  vRight <- head(df$yHat[df[,y] == L$yMax], 1)
  if (is.na(xlab)) 
    xlab <- "Predicted"
  if (is.na(ylab)) 
    ylab <- "Observed"

  df$outlier <- ifelse( df[,y] > L$limit[2] | df[,y] < L$limit[1], T, F )
  
 
 myTitle = bquote( 
   atop("QQ plot " *R^2*" = "*.(round( L$R2, 2)), "("*.(L$distribution)*") " )
 ) 

 mag <- 0.2
 prettyLims <- c( min(df[,c(y,"yHat")]), max(df[,c(y,"yHat")]) ) * c(1-mag, 1+mag)
 
  ggplot(data=df, aes_string(x="yHat", y= y) ) + 
    geom_hline(yintercept=L$limit[1], lty=2, color="gray30") +
    geom_hline(yintercept=L$limit[2], lty=2, color="gray30") +
    geom_vline(xintercept=vLeft, lty=2, color="gray30") +
    geom_vline(xintercept=vRight, lty=2, color="gray30") +
    geom_abline(slope=1,intercept=0)+
    #geom_point(color="black",size=2.6) +
    geom_point(data = subset(df,outlier==F),pch=21,size=3, fill="white") +
    geom_point(data = subset(df,outlier==T),pch=21,size=3, fill="red") +
    geom_text( data = data.frame(label=paste( "Left outliers:", L$nOut[[1]], "\nRight outliers:", L$nOut[[2]] )), 
             aes(label=label), x = log10(prettyLims[1]), y= log10(prettyLims[2]), hjust=-0.05, vjust=1.05,size=4) +
   geom_text( data = subset(df, eval(parse(text=y)) < L$limit[[1]]  ), aes_string(label=ref), hjust=-0.6, size =4 ) +
   geom_text( data = subset(df, eval(parse(text=y)) > L$limit[[2]]  ), aes_string(label=ref), hjust=1.4, size =4 ) +
    scale_y_log10(label=prettyNum) + scale_x_log10(label=prettyNum) + coord_fixed(xlim = prettyLims, ylim=prettyLims ) +
    #scale_fill_manual(values = c("white","red")) +
    labs(x=xlab, y=ylab) + ggtitle(myTitle) + theme(legend.position="none")
  
}

generate.matrices <- function(data){
  require(seqinr)
  # uco.freq and uco.freq should be loaded in memory, otherwise:
    if(!exists("uco.freq")){load("data/Rdata/orf_all.codon.usage_freq.Rda")}
    if(!exists("uco.rscu")){load("data/Rdata/orf_all.codon.usage_rscu.Rda")}
  # Codon usage based on frequency
  data.uco_freq <- sapply( data, function (x) uco.freq[ x, ] )
        #  save(m3d.uco_freq,file="data/Rdata/m3d.top10.codon_usage_freq.Rda")
  # Codon usage based on RSCU
  data.uco_rscu <- sapply( data, function (x) uco.rscu[ x, ] )
        # save(m3d.uco_rscu,file="data/Rdata/m3d.top10.codon_usage_freq.Rda")  
  # Generate matrices
  matrix.freq <- sapply( data.uco_freq, function(x) apply(x, 2, mean, na.rm=TRUE))
  matrix.rscu <- sapply( data.uco_rscu, function(x) apply(x, 2, mean, na.rm=TRUE))  
  
  return(list(matrix.freq,matrix.rscu))
}



boxcoxtransform <- function(x,lambda=0.5){
  min <- min(x,na.rm=T)
  if(lambda==0){
    tx <- logtransform(x, base=exp(1))
  } else {
    if(min<=0){ 
      constant <- abs(min) + 1
      tx <- ((x + constant )^(lambda) -1)/lambda 
    } else { tx <- (x^(lambda) -1)/lambda }
    
  }
  return(tx)
}

get.kurtosis_and_skewness <- function (x){
  require(fBasics, quietly=T, warn.conflicts=F)
  dagostino_test <- dagoTest(x)
  data.asymmetry <- data.frame( omnibus = as.numeric(dagostino_test@test$statistic[1]), 
                                skewness = as.numeric(dagostino_test@test$statistic[2]), 
                                kurtosis = as.numeric(dagostino_test@test$statistic[3]), 
                                p.omnibus = as.numeric(dagostino_test@test$p.value[1]),
                                p.skewness = as.numeric(dagostino_test@test$p.value[2]),
                                p.kurtosis = as.numeric(dagostino_test@test$p.value[3]) ,
                                row.names=NULL  )
  return( data.asymmetry ) 
}
# ---------------------------------------------- #
logtransform <- function(x, base=10, epsilon= 1e-20){
  min <- min(x,na.rm=T)
  constant <- abs(min) + 1
  if(min<=0){ tx <- log(x + constant, base= base ) } else { tx <- log(x + epsilon, base= base) } 
  #if(base==10){ if(min<=0){ tx <- log10(x + abs(min) + epsilon ) } else { tx <- log10(x + epsilon) } }
  return(tx)
}

# Rewrite variable's name in a more readable and natural format (replaces underscores by spaces)
reformat <- function (x,lab="\n", spt="_") { sapply(x, function (c) { paste(unlist(strsplit(as.character(c) , split=spt)),collapse=lab) }) }


boxcox_skewness_minimisation <- function (data=df, v="CV_ypd", range.boundaries=c(-4,4), by=0.1) {
  
  lambda_range <- seq(from=range.boundaries[1], to=range.boundaries[2], by=by)
  lambda_data <- ldply(lambda_range, function(l){ data.frame( variable = v, lambda= round(l,2), get.kurtosis_and_skewness( boxcoxtransform(x=data[,v], lambda=l))  ) })
  if (sum(is.na(lambda_data[,4:5])) > 0 ){ cat("WARNINGS [!] | NAs or NaN present indicating very high or very low values that cause rounding errors for some lambda.\n") }
  
  lambda_data <- subset(lambda_data, !is.na(skewness) )
  
  min.skewness <- subset(lambda_data, abs(skewness)==min(abs(skewness)), select=c(lambda,skewness) )
  
  
  d <- melt(lambda_data, id.vars=c("variable","lambda"))
  colnames(d)[3] <- "metric"
  d <- subset(d, metric %in% c("kurtosis","skewness"))
  S_m <-  c(min.skewness[[1]],min.skewness[[2]],subset(d,metric=="kurtosis" & round(lambda,2) == min.skewness[[1]])$value ) # Coordinates for minimum skewness
  S_0 = c(1, as.numeric(get.kurtosis_and_skewness(x=data[,v])[2:3])) # Coordinates of non-transformed data
  
  v2 <- paste0("transformed_",v)
  data[, v2]  <- boxcoxtransform( data[,v], lambda= S_m[[1]] )
  
    return(list(result = S_m, detail = lambda_data ))
}

refine.matrices <- function(selection,conditions_to_keep){
  require(som)
  M <-  generate.matrices(selection)
  if(missing(conditions_to_keep)){conditions_to_keep <- all()}
  M <- lapply( M, function(x) x[,conditions_to_keep] )
  M <- lapply( M, function(x) codon.subset(x, c("tgg", "atg")) )
  return(M)
}

rdc <- function(names){
# to print readable row/colnames in heat maps to be analysed
  require(seqinr)
  return(as.vector(readable.codons[names]))
}

# geometric men
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

##############################################################################################

# http://www.r-bloggers.com/extract-different-characters-between-two-strings-of-equal-length/
list.string.diff <- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE)
{
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case)
  {
    a <- toupper(a)
    b <- toupper(b)
  }
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
    ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
  if(!show.excluded) diff.info<-na.omit(diff.info)
  diff.info
}


type.pairing <- function( codon.nt = "a", anticodon.nt = "U" ){
  if( length(codon.nt) != length(anticodon.nt) ){ cat("Hey! provide equal-sized vectors for codon.nt and anticodon.nt!\n")} else {
   
    ref <- data.frame( on.codon = c("a","c","t","g","t","c","a","g"), 
                       on.anticodon = c("U","G","A","C","G","A","A","U"), 
                       complement = c("U","G","A","C","A","G","U","C"), 
                       type = c(rep("watson.crick",4), rep("wobble",4)) )
    
    result <- sapply( seq_along(codon.nt), function(u){
      is.pairing <- subset(ref, on.codon == tolower(codon.nt[[u]]) & on.anticodon == toupper(anticodon.nt[[u]]) )
      if(nrow(is.pairing)>0){ return( as.character(is.pairing$type )) } else {return("bad")}
    }        
    )
    
    return(result)
    
    
  }  
}

##############################################################################################
heatmap.3 <- function(matrix,a,b,main,KEY,ROWV,COLV,DENDRO,SCALE,SYM,SEPW,METHOD,LABROW,LABCOL,RowSideColors,ColSideColors, MIN,MAX,CENTRE, centered.break){
  require(marray)
  require(gplots)
  require(fBasics)
  require(grDevices)
  whisker <- boxplot.stats(matrix)$stats
  
  #op <- par(mar = par("mar")/2, oma=c(2,2,2,2))
  #pal0 <- maPalette(low="red", mid="white", high="#3C9DD0",50) #  
  #pal1 <- maPalette(low="black",mid="#6A94D4", high="#AFE7FF") #   pal1.1 <- maPalette(low="white", mid="gray", high="black",10) #   
  #pal2.0 <- maPalette(low="#6A94D4",mid="black", high="#FFCB73") #   pal2.1 <- c("#FFFFFF",maPalette(low="#6A94D4",mid="black", high="#FFCB73"))  
  palette <- maPalette(low="#18B2ED",mid="black", high="#F3F347") #   pal2.1 <- c("#FFFFFF",maPalette(low="#6A94D4",mid="black", high="#FFCB73"))  
  #pal2.0 <- rampPalette(n=60,name="cyan2magenta")
  
  if(missing(a)){a <- 1}; if(missing(b)){b <- 5}; if(missing(KEY)){KEY <- F}; if(missing(DENDRO)){DENDRO <- "both"}
  if(missing(SCALE)){SCALE <- "none"}; if(missing(SYM)){SYM <- F}; if(missing(SEPW)){SEPW <- NA} 
  if(missing(METHOD)){METHOD <- "ward"}; if(missing(main)){main <- ""}
  if(missing(LABROW)){LABROW <- colnames(matrix)}; if(missing(LABCOL)){LABCOL <- rownames(matrix)}  # not a mistake
  if(missing(ROWV)){ROWV <- T}; if(missing(COLV)){COLV <- T}
  if(missing(RowSideColors)){RowSideColors <- rep("white",nrow(matrix))}
  if(missing(ColSideColors)){ColSideColors <- rep("white",ncol(matrix))}  # not a mistake
  if(missing(MIN)){MIN <- whisker[1]}
  if(missing(MAX)){MAX <- whisker[5]}
  if(missing(CENTRE)){CENTRE <- median(matrix,na.rm=T) }
  if(missing(centered.break)){centered.break <- F}
  
  if(centered.break==T){
  A <- 2*abs(CENTRE - MIN)/length(palette); B <- 2*abs(MAX - CENTRE)/length(palette)
  BREAKS <- c( seq(from=MIN, to= (CENTRE - A) , by= A ) , CENTRE, seq( (CENTRE + B), MAX,  by= B ) )
  } else {
    BREAKS <- NULL
  }
  
  #if(KEY==T){LMAT <- rbind( c(0, 3), c(2,1), c(4,0) ); LHEI <- c(0.25, 3.25, 1.5 )}
  #if(KEY==F){LMAT <- rbind( c(0, 3), c(2,1), c(0,4) ); LHEI <- c(0.25, 7, 0.25 )}  
  
  dist2 <- function(x, ...) as.dist(1-cor(t(x), method="pearson"))
  
  h <- heatmap.2(
    matrix, Colv=COLV, Rowv=ROWV, dendrogram=DENDRO, key=KEY, keysize=1, 
    na.rm=TRUE, hclust=function(c){hclust(c, method=METHOD)}, distfun=dist2,
    main=main, labRow=LABCOL, labCol=LABROW,
    RowSideColors=RowSideColors,ColSideColors=ColSideColors,
    trace="none",notecex=0.25,scale=SCALE, symm=SYM,
    col=palette, breaks= BREAKS,
    density.info="density", denscol="white", 
    sepcolor='white', sepwidth=SEPW, colsep=1:ncol(matrix),
   # lmat= LMAT, lhei= LHEI,
    cexRow = 0.9*1/log10(nrow(matrix)), cexCol = 0.9*1/log10(ncol(matrix)), margins=c(a,b)
    ) 
 # par(op)
  
  h$hrow <- hclust(d = dist2(matrix), method = METHOD)
  h$hcol <- hclust(d = dist2(t(matrix)), method = METHOD)
  return(h)
}


codon.subset <- function(matrix,to.remove){
  matrix <- matrix[ -which(rownames(matrix) %in% to.remove), ]
  return( matrix )
}

experiment.subset <- function(matrix,to.remove,col){
  selection <- experiments[which(experiments[,2] %in% to.remove),col] 
  matrix <- matrix[ , -which(colnames(matrix) %in% selection) ]
  return( matrix )
}

color.clusters <- function(h){
  clusters.codons <- cutree(h$hcol, k=2)
  clusters.conditions <- cutree(h$hrow, k=2)
  # colors
  col <- list()
  col$codons <-     replace(clusters.codons, which(clusters.codons==1), "orange2"); col$codons <- replace(col$codons, which(col$codons==2), "skyblue2" )
  col$conditions <- replace(clusters.conditions, which(clusters.conditions==1), "red2"); col$conditions <- replace(col$conditions, which(col$conditions==2), "green2" )
  return(col)
}

highlights <- function(query,matrix){
  require(YET)
  query <- intersect(fORF(query), rownames(matrix))
  M.query <- as.matrix(matrix[query,]); head(M.query)
  rownames(M.query) <- fNAME(query); head(M.query)
  M.query[which(is.na(M.query))]  <- 0  
  return(M.query)
}

######## 

cai.sequence <- function(seq,w){
  require(seqinr)
  if(missing(w)){w <- caitab$sc; names(w) <- rownames(caitab)}
  out <- w[seq] 
  out <- as.data.frame(cbind(1:length(out), names(out), as.data.frame(out)))
  colnames(out) <- c("index","codons","CAI")
  return(out)
}

get.codons <- function(x){   
  require(seqinr)
  out <- c()
  for(i in seq(1,(3*floor(length(x)/3)),3) ){ out <- c(out, c2s(x[i:(i+2)])) } 
  return(out)
}

get.uco  <- function(biostring.seq,method="freq"){ uco(s2c(toString(biostring.seq)),index=method) }




get.ws <- function(tRNA, 
                   s = c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68), 
                   trim.above.ceiling = T,
                   ceiling = 1
                   )   
  # selective constraints (# optimised s-values)
  # super kingdom: 0-eukaryota, 1-prokaryota
  # I:U G:C U:A C:G 
  # G:U I:C I:A U:G
{
  # Anticodon (uppercase I,G,U,C) - codon (lower case, a,c,t,g)
  names(s) <- c("I:t","G:c","U:a","C:g","G:t","I:c","I:a","U:g")
  p = 1 - s

  # initialise w vector
  W = NULL  # don't confuse w (lowercase) and W (uppercase)
  S = NULL
  # obtain absolute adaptiveness values (Ws)
  ## [Guilhem: dos Reis generates the W scores 4 by 4, concatenating 4 new entries to W, 
  ## each of them being the sum of 1 conventional anticodon-codon pair and one non-conventional one
  ## since only 2 are possible per codon. If the anticodon doesn't exist at all, its contribution to the sum
  ## will be cancelled by having a tGCN of 0. It is thus normal to have weird inexisting anticodons in 
  ## Source/data_sources/tables/wobble_crick.txt ]
  for (i in seq(1, 61, by=4))
    W = c(W,
      p[1]*tRNA[i]   + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
      p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
      p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
      p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG

  
  # check methionine
  W[36] = p[4]*tRNA[36]
  
  # if bacteria, modify isoleucine ATA codon

  # get rid of stop codons (11, 12, 15) 
 # W = W[-c(11,12,15,36)]
  W[c(11,12,15)] <- NA
 
  # get ws #[G] new normalising factor:
  WN = gm_mean( sort(W,decreasing = T)[1:3], na.rm = T)
  w  = W/WN
 
  if( trim.above.ceiling == T){
    w[w> ceiling & !is.na(w) ] <- ceiling  
  }
  ##  w = W/max(W, na.rm=T)

  if(sum(w == 0, na.rm=T) > 0) {
    ws <- w[w != 0 & !is.na(w)] # zero-less ws
    gm <- exp(sum(log(ws))/length(ws)) # geometric mean
    w[w == 0] = gm # substitute 0-ws by gm
  }

  attr(w,"Wmax") <- WN
  return(w)
}  
  

gc.content <- function(string.vector, digits=3){
  out <- sapply( string.vector, function(i) {
    x <- unlist(strsplit(toupper(as.character(i)), split=""))
    as.numeric(round(sum(x %in% c("G","C"))/length(x),digits = digits))
  } )
  as.numeric(out)
}


############################################################################
# After calculating all w's, the tRNA adaptation index (tAI) can then be
# calculated. x is a matrix with all the codon frequencies per ORF in any
# analysed genome
#############################################################################

get.tai <- function(x,w.data, score="w") {
 # [G] does this ensure at all correspondance between x and w? No it doesn't!
 # [G] I will modify the code to ensure the ordering is correct  
  codons  <- colnames(x)
  w <- w.data[,score]
  names(w) <- w.data[,"codon"]
  w <- w[codons]
  w = log(w)      #calculate log of w
  n = apply(x,1,'*',w)		#multiply each row of by the weights
  n = t(n)                      #transpose
  n = apply(n,1,sum)		#sum rows
  L = apply(x,1,sum)		#get each ORF length
  tAI = exp(n/L)		#get tai
  return(tAI)
}

# compute metrics about the use of anticodons in CDS
compute.anticodon.enrichment <- function(  gene.names = as.character(unique(mRNA_molecule.copy.number$ORF)),
                                           codon.info = rosetta.codons,
                                           start = 1,
                                           l = "all",
                                           load.codon.eff = T
                                           ){
                                        
  codon.info <- subset(codon.info, select=c(aa, anticodon, codon))
  require(StingRay);data(SEQ)
  require(seqinr)
  require(data.table)
  
   if (l == "all"){
      codon.eff <- ldply( setNames(gene.names,gene.names), 
      # full ORF lenth minus STOP codon
                         function(x){ get.uco( SEQ$ORF_CDS[[x]][ (3*start-2):(length(SEQ$ORF_CDS[[x]]) - 3)], 
                                              method="eff") 
                                      } 
          )       
  
    } else { # assuming l is a number this time
      codon.eff <-  ldply( setNames(gene.names,gene.names), 
                            # segment of l codons
                            function(x){ get.uco(SEQ$ORF_CDS[[x]][ (3*start-2):((l*3) + 3*(start - 1) ) ], 
                                                 method="eff") 
                            } 
      )
    }
  
  codon.eff$N <- rowSums(codon.eff[,-1])
  codon.eff.m <- melt(codon.eff, id.vars=c(".id","N"))
   
  data <- data.table( merge( codon.info, codon.eff.m, by.x="codon" , by.y="variable" ), 
                      key = c(".id","N","anticodon","codon")
                      )
  

  n.synonyms.codons  <- ddply( codon.info, .(anticodon), summarise, n.codon = length(unique(codon)) ) 
  n.synonyms.anticodon <- merge( unique(codon.info[,c("aa","anticodon")]), ddply( codon.info, .(aa), summarise, n.syn.anticodon = length(unique(anticodon)) ) , by = "aa")[,-1]

  synonyms <-            data.table( merge( n.synonyms.anticodon, n.synonyms.codons, by="anticodon"), key = "anticodon" )
  n.anticodon         <- data[, list(n.anticodon = sum(value)), by=list(.id, N, anticodon)] 
  rel.freq.anticodon  <- data[, list(anticodon = anticodon, freq.anticodon = value/sum(value)), by=list(.id, N, aa)]
  # didn't manage to get an elegant one step solution, but anyway this is massively faster than with ddply
  rel.freq.anticodon  <- rel.freq.anticodon[, list(freq.anticodon = sum(freq.anticodon)), by=list(.id,N, aa, anticodon) ] 
  rel.freq.anticodon$freq.anticodon[which(is.nan(rel.freq.anticodon$freq.anticodon))] <- 0
  setkey( rel.freq.anticodon, anticodon )

  d <- rel.freq.anticodon[ synonyms, ]
  setkey(d, .id, N, anticodon)
  d <- d[ n.anticodon, ]
  setnames(d, ".id","ORF")

  return(d)
}


compute.codon.frequency.genes <- function( gene.names = unique(mRNA_molecule.copy.number$ORF), 
                                           seq.set = SEQ$ORF_CDS,
                                           l = "all", # n. of codons
                                           start = 1, # 2 to ignore the 1st atg codon
                                           method="freq"
                                           ){
  
   require(seqinr)
  
  if (l == "all"){
  codon.freq <-  ldply( setNames(gene.names,gene.names), 
                        # full ORF lenth minus STOP codon
                        function(x){ get.uco(seq.set[[x]][ (3*start-2):(length(seq.set[[x]]) - 3)  ], 
                                             method=method) 
                                     } 
                       )
  } else { # assuming l is a number this time
    codon.freq <-  ldply( setNames(gene.names,gene.names), 
                          # segment of l codons
                          function(x){ get.uco(seq.set[[x]][ (3*start-2):((l*3) + 3*(start - 1)) ], 
                                               method=method) 
                          } 
    )
  }
  # reproduce the ordering of codon used by dos Reis in dosReis2004
  nt <- c("t","c","a","g")
  codons <- apply(expand.grid(nt, nt, nt)[,3:1], 1, paste, collapse="") # == as.character(codon.table$codon)
  
  colnames(codon.freq)[1] <- "ORF"
  codon.freq.matrix <- as.matrix(codon.freq[,-1])
  rownames(codon.freq.matrix) <- codon.freq$ORF
  codon.freq.matrix <- codon.freq.matrix[,codons]
  return(codon.freq.matrix)
}


# negDistMat expSimMat linSimMat corSimMat linKernel
build.ap.clusters <- function(data = M.tRNAab$all.together ){
  require(apcluster)
  apclusters <- list( 
    negDistMat = apcluster( s =  negDistMat(data) ),
    expSimMat  = apcluster( s =  expSimMat(data) ),
    linKernel  = apcluster( s =  linKernel(data) )
  )
  return(apclusters)
}

identify.clusters <- function( apcluster.object, anticodon.info = rosetta.anticodons ){
  require(apcluster)
  require(plyr)
  clusters <- apcluster.object@clusters
  data <- ldply( 1:length(clusters), function(x){ data.frame( cluster.id = x, anticodon = names(clusters[[x]])  )  } )
  data <- merge(data, anticodon.info, by="anticodon")
  return(data)
}

pairwise.wilcox.test.gul <- function (x, g, p.adjust.method = p.adjust.methods, paired = FALSE, 
          ...) 
{
  p.adjust.method <- match.arg(p.adjust.method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  METHOD <- if (paired) 
    "Wilcoxon signed rank test"
  else "Wilcoxon rank sum test"
  
  compare.levels <- function(i, j) {
    xi <- x[as.integer(g) == i]
    xj <- x[as.integer(g) == j]
    w <- wilcox.test(xi, xj, paired = paired, ...)
    return( list( statistics = w$statistic, p.value = w$p.value))
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  ans
}

# need this to get the stupid W AND P-values in pairwise manner for the analysis
wilcox.test.batch <- function(x, grouping.vars = c("stress","faster.translation"), value = "value", adjust.method = "BY") { 
  x <- as.data.frame(x)
  require(coin)
  x$group  <-  if (length(grouping.vars)>1){ do.call(paste, x[,grouping.vars]) } else { as.character(x[,grouping.vars]) }
  x$value  <-  x[,value]
  c1 <- combn(unique(as.character(x$group)),2)
  test <- apply ( c1, 2, function(d){ sx <- subset( x, group %in% d );
                                      sx$group <- factor(sx$group)
                                      w <- wilcox_test( formula = value ~ group, data = sx, conf.int = T ) # was wilcox_text for some reason
                                      w.bis <- wilcox.test( formula = value ~ group, data = sx ) # was wilcox_text for some reason
                                      ns <- cast( data.frame( table (w@statistic@x) ), ~ Var1, value = "Freq")[,d]
                                      ref.n <- table(w@statistic@x)[1]
                                      medians <- round(cast( ddply(sx, .(group), summarise, median = median(value, na.rm=T)), ~ group, value = "median" )[,d],3)
                                      delta_medians <- as.numeric(medians[2] - medians[1])
                                      percent_median1 <- round(as.numeric(-delta_medians/medians[1]),2)
                                      percent_median2 <- round(as.numeric(delta_medians/medians[2]),2)
                                      colnames(medians) <- c("median.1","median.2")
                                      colnames(ns) <- c("n.1","n.2")
                                      U <- w@statistic@linearstatistic[[1]] - ref.n*(ref.n+1)/2
                                      # For the definition of U see my own answer on 
                                      # http://stats.stackexchange.com/questions/79843/is-the-w-statistic-outputted-by-wilcox-test-in-r-the-same-as-the-u-statistic
                                      data.frame( Z = prettyNum(w@statistic@teststatistic[[1]], digits=3), 
                                                  W = prettyNum(w@statistic@linearstatistic[[1]]), 
                                                  U = prettyNum(U), 
                                                #  expectation = prettyNum(w@statistic@expectation[[1]]),
                                                  effect_size = round(1 - 2*U / (ns[[1]] * ns[[2]]),2),
                                                 # r = round( abs(w@statistic@teststatistic[[1]]) / sqrt( sum(ns) ), 2),
                                                  medians,
                                                  delta = delta_medians,
                                                  percent_change.1 = percent_median1,
                                                  percent_change.2 = percent_median2,
                                                  n.total = nrow(sx),
                                                  ns, 
                                                  p.value = pvalue(w),
                                                  stringsAsFactors=F
                                                  )
                                     } )
  results <- data.frame( t(c1), do.call(rbind, test), stringsAsFactors=F )  
  colnames(results)[1:2] <- c("set.1","set.2")
  results$p.adjusted <-   noquote(format.pval(p.adjust(results$p.value, method = adjust.method), digits = 3)) 
  results$p.value <- noquote(format.pval(results$p.value, digits = 3)) 
  rownames(results) <- NULL
  return(results)
}

matrixify <- function(df, ref.n = 1){
  m <- as.matrix(df[, -ref.n])
  df[,ref.n] <- lapply(df[,ref.n, drop=F], as.character)
  rownames(m) <- apply( df[,ref.n, drop=F], 1 , paste, collapse="_")
  colnames(m) <- colnames(df)[-ref.n]
  return(m)
}

# pretty scientific notations, with threshold for use of x10^ form
scientific_10 <- function(x, threshold=10) {
  require(scales)
  parse(text=gsub("e[+?]", " %*% 10^", ifelse( abs(x)<=threshold, prettyNum(x), scientific_format()(x) ),perl = T ))
}


# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(p){
  tmp <- ggplotGrob(p)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


# http://stackoverflow.com/questions/10985224/r-heatmap-with-diverging-colour-palette
diverge.color <- function(data,pal_choice="RdGy",centeredOn=0,Min,Max){
  require(grDevices)
  require(classInt)
  require(RColorBrewer)
  nHalf=50
  N <- brewer.pal.info[pal_choice,"maxcolors"]
  if(missing(Min)){Min <- min(data,na.rm=TRUE)}
  if(missing(Max)){Max <- max(data,na.rm=TRUE)}
  Thresh <- centeredOn
  pal<-brewer.pal(n=N,pal_choice)
  rc1<-colorRampPalette(colors=c(pal[1],pal[2]),space="Lab")(10)
  for(i in 2:(10) ){
    tmp<-colorRampPalette(colors=c(pal[i],pal[i+1]),space="Lab")(10)
    rc1<-c(rc1,tmp)
  }
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  cuts <- classIntervals(data, style="fixed",fixedBreaks=rampbreaks)
  return(list(cuts,rc1))
}


gulmap <- function(data, centeredOn = 0, pal = "RdBu",Min, Max,...){
  require(pheatmap)
  require(RColorBrewer)
  if(missing(Min)){Min <- min(data,na.rm=TRUE)}
  if(missing(Max)){Max <- max(data,na.rm=TRUE)}
  pheatmap(data,
           color  = diverge.color(data  = data, centeredOn = centeredOn, pal_choice = pal, Min = Min, Max = Max)[[2]],
           breaks = diverge.color(data = data, centeredOn = centeredOn, Min = Min, Max = Max)[[1]]$brks,...
           )
}


recast3 <- function(data, formula, ..., id.var, measure.var) {
  if (any(c("id.vars", "measure.vars") %in% names(match.call()))) {
    stop("Use var, not vars\n")
  }
  
  molten <- melt(data, id.var, measure.var)
  dcast(molten, formula, ...)
}

# http://stackoverflow.com/questions/20372224/how-to-only-change-parameters-for-lower-plots-in-the-ggpairs-function-from-gga
add_p<-function(g,i,params){
  
  side=length(g$columns)                # get number of cells per side
  lapply(i,function(i){
    
    s<-as.character(g$plots[i])         # get existing call as a template
    l<-nchar(s)
    p<-paste0(substr(s,1,l-1),",",params,")")   # append params before last bracket
    r<-i%/%side+1                               # work out the position on the grid
    c<-i%%side
    
    array(c(p,r,c))                     # return the sub-plot and position data
    
  })
  
}

# http://codereview.stackexchange.com/questions/17905/compute-intersections-of-all-combinations-of-vectors-in-a-list-of-vectors-in-r
overlap <- function(l) {
  results <- lapply(l, unique)
  
  # combinations of m elements of list l
  for (m in seq(along=l)[-1]) {
    
    # generate and iterate through combinations of length m
    for (indices in combn(seq(length(l)), m, simplify=FALSE)) {
      
      # make name by concatenating the names of the elements
      # of l that we're intersecting
      name_1 <- paste(names(l)[indices[-m]], collapse="_")
      name_2 <- names(l)[indices[m]]
      name <- paste(name_1, name_2, sep="_")
      
      results[[name]] <- intersect(results[[name_1]], results[[name_2]])
      
    }
  }
  results
}

flip.list <- function(l){
  ll <- lapply(seq_along(l)[[1]], function(i) lapply( l, "[[", i));
  names(ll) <- names(l[[1]])
  return( ll )
}

vennit <- function(l){
  require(reshape2)
 d <- melt(l)
 colnames(d) <- c("name","ind")
 t <- arrange(ldply( reshape2::dcast( d, ind ~ name, value.var="ind" )[,-1], function(x){  data.frame( group = paste0(x[!is.na(x)],
                                                                                                    #ifelse(is.na(x), 0, 1), 
                                                                                                    collapse="&"));  }),
              group)   
 return(t)
}

split.venn <- function(t){
  t <-  dlply(t, .(group), function(x) x$.id)
  return(t)
}

findPalette <- function(x, N=8, palette="Dark2"){
  require(RColorBrewer)
  n <- length(unique(x))
  palette <- c(rep(brewer.pal(N, palette), n %/% N), brewer.pal(N, palette)[0:(n%%N)])
  return(palette)
}


#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
rbind_gtable_max <- function(...){
  require(gtable)
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}

# http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
cbind_gtable_max <- function(...){
  require(gtable)
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(nrow(x) == nrow(y))
    if (ncol(x) == 0) 
      return(y)
    if (ncol(y) == 0) 
      return(x)
    y$layout$l <- y$layout$l + ncol(x)
    y$layout$r <- y$layout$r + ncol(x)
    x$layout <- rbind(x$layout, y$layout)
    x$widths <- gtable:::insert.unit(x$widths, y$widths)
    x$colnames <- c(x$colnames, y$colnames)
    x$heights <- grid::unit.pmax(x$heights, y$heights)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  Reduce(bind2, gtl)
}



# automated linear regressions
batch.regressions = function(df, y = "total.tRNAab", x ="demand.mRNAab.up", display.mode = T ){
  f <- as.formula( paste(y, x, sep = " ~ ") )
  m = lm( formula = f, df);
  pears = cor.test( formula = as.formula( paste( "~", y, "+", x)), data = df, method = "pearson",exact = F);
  line1 <- substitute(italic(r)~"="~pears~" | "~alpha~"="~coeff, 
                      list( pears = format(pears$estimate, digits = 2),
                            coeff = format(summary(m)$coefficients[[2]], digits = 2)
                      )
  )
  line2 <- substitute(italic(R)^2~"="~R2,
                      list( R2 = format(summary(m)$r.squared, digits = 2))
  )
  if(display.mode == T ){
    as.character( as.expression( paste( c('atop(', line1, ',', line2,')'), collapse=" ") ) )  
  } else {
    data.frame( r = pears$estimate, coeff = summary(m)$coefficients[[2]], intercept = summary(m)$coefficients[[1]], R2 = summary(m)$r.squared )
    
  }
}
