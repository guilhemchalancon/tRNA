library("Biostrings") # notably to read FASTA files
library("reshape")
library(YET)
setwd("~/Documents/MRC16/")

# Protein properties
# http://www.yeastgenome.org/download-data/curation
protein.properties <- read.table(file="data/protein_properties.tab",sep="\t",row.names=1)
colnames(protein.properties) <- c("SGDID","MOLECULAR_WEIGHT","PI","CAI","PROTEIN_LENGTH","N_TERM_SEQ","C_TERM_SEQ","CODON_BIAS","ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","FOP_SCORE","GRAVY_SCORE","AROMATICITY_SCORE","Feature_type")
save(protein.properties,file="data/Rdata/protein.properties.Rda",safe=T)

# orf_coding_all.fasta.gz
# ORF Coding Sequences (CDS) only, without 5'-UTR, 3'-UTR, intron
# sequences, or bases not translated due to translational
# frameshifting.  For all ORFs.

orf_genes.coding <- readFASTA(file="data/orf_coding_all_sgn.fasta")
            desc <- sapply(orf_genes.coding, function(x) x$desc)
   orf_genes_seq <- sapply(orf_genes.coding, function(x) x$seq)
   orf_genes_seq <- as.list(orf_genes_seq)
            names(orf_genes_seq) <- desc
save(orf_genes_seq,file="data/Rdata/orf_genes.coding_all.Rda",safe=T)

# RNA sequences of yeast Open Reading Frames
# http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_dna.README
# rna.coding <- read.delim(gzfile("data/rna_coding.fasta.gz"), sep="\t") # careful, it is FASTA files, not a table
# RNA genes (rRNAs, tRNAs, snRNAs, snoRNAs, and ncRNAs)
rna_genes.coding <- readFASTA("data/rna_genes.coding_sgn.fasta")
            desc <- sapply(rna_genes.coding, function(x) x$desc)
   rna_genes_seq <- sapply(rna_genes.coding, function(x) x$seq)
   rna_genes_seq <- as.list(rna_genes_seq)
                    names(rna_genes_seq) <- desc
  save(rna_genes_seq,file="data/Rdata/rna_genes.coding.Rda",safe=T)


other <- read.table("~/Documents/MRC13/data/yg_s98_v3_Build_2_norm/expr.norm.tab",sep="\t")
boxplot(other)




# Gash et al., Mol Biol Cell (2000 and 2001)
# http://genome-www.stanford.edu/yeast_stress/
# http://genome-www.stanford.edu/yeast_stress/materials.pdf (Material and methods)
gash.stress <- read.table("data/Gash_and_Chen/gash.stress.tab",sep="\t",header=1,row.names=1,check.names=F);head(stress)
  save(gash.stress,file="data/Rdata/gash.stress.Rda")
# http://www-genome.stanford.edu/mec1
# The data represent normalized, background-corrected log2 values of the Red/Green ratio measured 
# on each microarray.
gash.dna <- read.table("data/Gash_and_Chen/gash.dna.tab",sep="\t",header=1,row.names=1,check.names=F); head(dna)
  save(gash.dna,file="data/Rdata/gash.dna.Rda")



# http://jura.wi.mit.edu/fink_public/h_chen/

# Many fungi undergo a developmental transition from a unicellular yeast form to an invasive 
# filamentous form in response to environmental cues. Here we describe a quorum-signaling pathway 
# that links environmental sensing to morphogenesis in Saccharomyces cerevisiae. Saccharomyces cells
# secrete aromatic alcohols that stimulate morphogenesis by inducing the expression of FLO11 through
# a Tpk2p-dependent mechanism. Mutants defective in synthesis of these alcohols show reduced 
# filamentous growth, which is partially suppressed by the addition of these aromatic alcohols. 
# The production of these auto-signaling alcohols is regulated by nitrogen: high ammonia restricts 
# it by repressing the expression of their biosynthetic pathway, whereas nitrogen-poor conditions 
# activate it. Moreover, the production of these aromatic alcohols is controlled by cell density and 
# subjected to positive feedback regulation, which requires the transcription factor Aro80p. 
# These interactions define a quorum-sensing circuit that allows Saccharomyces to respond both to cell
# density and the nutritional state of the environment. These same autoregulatory molecules do not 
# evoke the morphological switch in Candida albicans, suggesting these molecular signals are 
# species-specific.

chen.alcohol.original <- read.table("data/Gash_and_Chen/chen.alcohol.tab",sep="\t",header=1)

# Important: some genes are duplicated (obtained by different strains. I need to average them for now)
    # Get the list of identifier for duplicated ORFs
    duplicates.chen <- unique(chen.alcohol.original$ORF[duplicated(chen.alcohol.original$ORF)])
    # Subset unique genes from the original dataset
    unique.chen <- chen.alcohol.original[-which(chen.alcohol.original$ORF %in% duplicates.chen),]
    # For each duplicated, apply an average, for each data value
    averaged.duplicates <- lapply( duplicates.chen, function (x) { y <- as.character(x); names(y) <- "ORF"; c(y,apply(chen.alcohol.original[ which(chen.alcohol.original$ORF %in% x), 2:5],2,mean,na.rm=T) ) } )
    # Reformat the ex.duplicates into a data frame
    require(plyr); ex.duplicates.chen <- ldply(averaged.duplicates, fun = rbind)
    ex.duplicates.chen <- apply(ex.duplicates.chen, 1, function(x) {as.numeric(x[,2:5])})
    # Bind together the two subset into a new data frame
    chen.alcohol <- rbind(unique.chen,ex.duplicates.chen)
    chen.alcohol[,2:5] <- apply(chen.alcohol[,2:5], 1, as.numeric)
    # Construct the data frame in its final format
    rownames(chen.alcohol) <- chen.alcohol$ORF; chen.alcohol <- chen.alcohol[,-1] 
      
######## Combined datasets
gash.all <- merge(gash.stress,gash.dna,by.x=0,by.y=0); rownames(gash.all) <- gash.all[,1]; gash.all <- gash.all[,-1]
gash.all  <- merge(gash.all, chen.alcohol,by.x=0,by.y=0); rownames(gash.all) <- gash.all[,1]; gash.all <- gash.all[,-1]
    save(gash.all,file="data/Rdata/gash.all_sgd.outdated.Rda")
    # Remove the ORFs present in Gash et al. (2000, 2001) that have been deprecated in the latest releases of the SGD

gash.all <- gash.all[-which(rownames(gash.all) %in% setdiff(rownames(gash.all), names(orf_genes_seq))), ]
    save(gash.all,file="data/Rdata/gash.all.Rda")
  

## Checking the quantity of NAs per genes and per conditions
gash.colsumNAs <- apply( gash.all, 2, function (x) {length(which(is.na(x))) } )
gash.rowsumNAs <- apply( gash.all, 1, function (x) {length(which(is.na(x))) } )

pdf("results/gash_without_filtering.pdf")
plot(gash.rowsumNAs[order(gash.rowsumNAs)], ylab="number of conditions where missing", xlab="genes")
plot(gash.colsumNAs[order(gash.colsumNAs)], ylab="number of genes without data", xlab="conditions")
dev.off()

# Not sure if I should apply a filtering of the data yet
gash.genes.to.filter <- names(gash.rowsumNAs[gash.rowsumNAs>45])
gash.conditions.to.filter <- names(gash.colsumNAs[gash.colsumNAs>600]) # contains heat shocks


library(som)
# Matrix of M3D
m3d <- read.table("data/m3d.expr.original.tab",sep="\t",row.names=1)

# m3d.log <- normalize(m3d,byrow=F) # already normalised (affymetrix, hence single-channel)
m3d.log <- read.table("data/m3d.expr.norm.tab",sep="\t",row.names=1)
rownames(m3d.log) <- fORF(rownames(m3d.log))

    # Are there genes in m3d not found in orf_genes_seq?
    dismissed <- setdiff(rownames(m3d),names(orf_genes_seq))
    # Which genes from orf_genes_seq are absent in m3d?
    missing <- setdiff(names(orf_genes_seq),rownames(m3d))

    # Remove the two ORFs present in M3D that have been deprecated in the latest releases of the SGD
    m3d.log <- m3d.log[-which(rownames(m3d.log) %in% dismissed),]
    m3d <- m3d[-which(rownames(m3d) %in% dismissed),]

save(m3d,file="data/Rdata/m3d.Rda")
save(m3d.log,file="data/Rdata/m3d_log.transformed.Rda")

m3d.tRNAs <- read.table("data/Yeast_Genome_S98_v2_norm/cleaned.avg_Yeast_Genome_S98_v2_tRNAs.tab",sep="\t",row.names=1)
  save(m3d.tRNAs,file="data/Rdata/m3d.tRNAs.Rda")
# problems of duplicated names

pdf("results/normalisation.m3d.pdf")
boxplot(m3d,axes=F,main="normalised expression levels (248 conditions) - all ORFs")
box()
axis(2)
dev.off()

pdf("results/normalisation.m3d_log.tranform.pdf")
boxplot(m3d.log,axes=F,main="normalised expression levels, log transformed - all ORFs")
box()
axis(2)
dev.off()
