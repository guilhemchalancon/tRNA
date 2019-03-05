require(XLConnect)
require(StingRay)
require(YET)
require(data.table)

# Data on tRNA export and synthesis regulators
wb <- loadWorkbook("~/Documents/MRC16/data/tRNA regulators/GENE_LIST.xlsx")
tRNA.regulators <- readWorksheet(wb, sheet = getSheets(wb),header=T)[["table"]]; head(tRNA.regulators)
tRNA.regulators$gene <- fORF(tRNA.regulators$name)

# Transcriptomic data (Gasch A, et al. 2000)
data(TRS)
Gasch.corrected.names <- read.table("~/Documents/StingRay/Source/data_sources/transcriptomes/gash.all.txt", check.names=F)
tRNA.Gasch <-  as.data.table(merge(tRNA.regulators, Gasch.corrected.names, by.x="gene",by.y=0, all.x=T)); tRNA.MD3


# Ad hoc commands
give <- function(x,pattern){
return(x[ grep(pattern,  x, ignore.case=T) ])
}
dist2 <- function(x, ...) as.dist(1 - cor(t(x), method = "pearson"))
 
selection <- list()
# Heat shock 
selection$heat.shock <- give(x=colnames(Gasch.corrected.names),pattern="heat_shock")[12:15]
# Hypo osmotic shock
selection$osmotic.shock <- give(x=colnames(Gasch.corrected.names),pattern="1M_sorbitol__")[4:10]
# Oxidative stress
selection$oxidative.stress <- give(x=colnames(Gasch.corrected.names),pattern="H2O2")[1:10]
# Stationary phase
selection$stationary.phase <- give(x=colnames(Gasch.corrected.names),pattern="stationary_phase")


# shorten labels to make them more readable and elegant by only catching the tag referring to time point
labels.selection <- reformat( sapply(do.call(c, selection), function(c){ #paste(unlist(
  string  <- as.character(c);
  pat <- '^.*__([0-9]+[hmind]+)\\.?.*$';
  sub(pat, '\\1', string[grepl(pat, string)]) 
}), lab=" ", spt="_")

labels.groups <- do.call(c,lapply( names(selection), function(x){
  rep(reformat(x, lab=" ", spt="\\."), times=length(selection[[x]]))
} )  )

# Build data set
  # combine tRNA lists with selected stress conditions
  data <- subset( tRNA.Gasch, select=c(colnames(tRNA.Gasch)[1:3], do.call(c, selection) ))
  data <- data[complete.cases(data),]
#     other.genes <- subset( TRS$Gasch2000_response.stress, select=c(heat.shocks, hypo.osmotic.shocks) )
#     mean.other.genes <- melt(other.genes)
  # molten data
  m.data <- data.table(melt(data, id=c("gene","name","functional_category")))
  m.data$stress <- factor(m.data$variable, levels= do.call(c, selection), labels = labels.groups )
  m.data$variable <- factor(m.data$variable, levels= do.call(c, selection), labels= labels.selection )
 # m.data$variable <- as.character(m.data$variable)  
  m.data$stress <- as.character(m.data$stress)
  #m.data$variable <- reformat(x=m.data$variable,lab=" ",spt="_")
  setnames(m.data, "value", "log2")
  # cluster value per functional_category etc.
  clustered.data <- ddply(m.data, .(functional_category), .fun=
               function(c){ y <- cast( c, formula= name ~ variable, value= "log2", mean ); 
                            z <- as.matrix(y[,-1]); 
                            rownames(z) <- y$name;
                            colnames(z) <- colnames(y)[-1];
                            d <- hclust(d= dist(x= as.matrix(z)),
                                        #dist2(x= as.matrix(z), na.rm=T),
                                        method="ward") 
                            return( data.frame(name=d$label, order=d$order) )
                            } )   


m.data$name <- factor(m.data$name, levels= as.character(arrange(clustered.data, functional_category, order )$name) )

# Construct heatmap # symetric color scale!
Ncol <- 50
coolColors <- colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(Ncol)
max_val <- max(abs(m.data$log2))
values <- seq(-max_val, max_val, length = Ncol)

save(m.data, file="results/tRNA.regulators/m.data_transcriptional.response.Rda")

ggplot(data=m.data, 
       aes(x=variable, y=name, fill=log2) ) + 
  ylab("gene") +  
  geom_tile(alpha=0.5) + geom_tile() + theme_bw() +
  scale_fill_gradientn(colours=coolColors, values = values, rescaler = function(x, ...) x, oob = identity ) + 
  facet_grid(functional_category ~ stress, scales="free", space="free", drop=T) +
  theme(axis.text.x=element_text(angle=-45, hjust=0, size=8), axis.text.y=element_text(size=7), legend.position="top") + 
  labs(fill="log2 ratio mRNA expression")

