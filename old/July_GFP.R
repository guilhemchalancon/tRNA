#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#     OBJECTIVES    #--------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# control the execution of stochastic simulation of protein translation in batch, for 9 (or 13) distinct experimental conditions (1 normal, 4 stresses with 2 or 3 time points)


setwd("~/Documents/MRC16/")



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#      Commands     #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


converts.mRNAab.files <- function(
  path.perlscript = "data/SMoPT/utilities/convert.fasta.to.genom.pl",
  fasta.file = "data/SMoPT/batch/input_GFPs/S.cer.fasta",
  genetic.code.file = "data/SMoPT/utilities/genetic.code.tsv",
  experiment = "ox",
  time = "20",
  GFP,
  output.path = "data/SMoPT/batch/input_GFPs/genomes/"
){
 
    mRNAab.file <- paste0("data/SMoPT/batch/input_GFPs/mRNA.abundance/Scer.mRNAab.",experiment, if(!missing(GFP)){paste0("_",GFP)} ,"_",time,".txt") # note that the mRNA abundance files must preexist (another function to build)
    output.file <- paste0(output.path,"Scer.",experiment, if(!missing(GFP)){paste0("_",GFP)}, "_",time)
 
  paste("perl",path.perlscript, fasta.file, genetic.code.file, mRNAab.file, output.file )
}


#----------------------------------------------------------------------------------------------------------------------------------#

# get command to run the SMPoT binary
get.SMPoT.command <- function( 
  path.bin = "data/SMoPT/bin/SMoPT",
  duration = 1500,
  heating = 1000,
  n.ribosomes = 200000,
  n.tRNAs = 3.3e6,
  n.mRNAs = 3795,
  sequences.path = "data/SMoPT/batch/input_GFPs/genomes/",
  tRNA.info.path = "data/SMoPT/batch/input_GFPs/adjusted.tGCN/",
  seed = 1413,
  experiment = "normal",
  time = "0",
  GFP,
  output.path.root = "data/SMoPT/batch/output_GFPs/",
  options = c(1:5)
){
  sequences   = paste0(sequences.path, "Scer.", experiment, if(!missing(GFP)){paste0("_",GFP)}, "_", time, ".genom", collapse = "") # Scer.genome.normal_0.txt
  tRNA.info   = paste0(tRNA.info.path, "Scer.tRNA.",   experiment, "_", time, ".txt", collapse = "") # Scer.tRNA.normal_0.txt
  time.folder = sprintf("%03.0f", as.numeric(time)) # e.g. gets 0 written as 000, 20 as 020 and 120 as 120
  output.path = paste0(output.path.root,  experiment, "/", GFP, "/Scer" )
  paste( paste0("./", path.bin), 
         "-T", format(duration, scientific = F), 
         "-H", format(heating, scientific = F), 
         "-R", format(n.ribosomes, scientific = F),
         "-t", format(n.tRNAs, scientific = F),
         "-N", format(n.mRNAs, scientific = F),
         "-F", sequences,
         "-C", tRNA.info,
         "-s", seed,
         "-O", output.path,
         paste( paste0("-p", options), collapse=" ")
  )  
}



#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#        Data       #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# prepare all the data needed to process input files and exectute hte scripts
relative.GFP.level <- 0.5
GFPs <- c("eGFP_diauxic","eGFP_ox","eGFP_osm","eGFP_temp","eGFP_normal")
# 0 # define which conditions to compile
conditions <- arrange( data.frame( experiment = rep(c("diauxic","ox","osm","temp","normal"),each=5), 
                                          GFP = rep(GFPs,times=5),
                                          time = 20,
                                          stringsAsFactors=F ), experiment)
conditions[conditions$experiment == "normal","time"] <- 0

#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1 # gcn data

# stress-adjusted tGCN values (precomputed)
adjusted.tGCN <- arrange( read.table( "results/stress tAI/adjusted.tGCN.txt", header =1  ), anticodon)
adjusted.tGCN <- adjusted.tGCN[,c("anticodon","normal_0","diauxic_20","ox_20","osm_20","temp_20")]
adjusted.tGCN[,c("normal_0","diauxic_20","ox_20","osm_20","temp_20")] <- sapply(adjusted.tGCN[,c("normal_0","diauxic_20","ox_20","osm_20","temp_20")], 
                                                                                function(x){ x <- round(x); x[which(x==0)] <- 1; x } )
rosetta.codons <- read.table("data/info/rosetta.codons.txt",header=1)
# the output is as follows: codon number \t tRNA number \t adjusted tGCN \t wobbling parameter

# important: need to use the exact same codon and anticodon order than the SMoPT scripts
ordering.tRNA   <- data.frame( t.number=0:40, subset(read.table("data/SMoPT/utilities/S.cer.tRNA_copy_num", header=1), gcn >0 ))[,-3]
colnames(ordering.tRNA)[2] <- "anticodon"
ordering.codons <- data.frame( c.number=0:60, read.table("data/SMoPT/utilities/genetic.code.tsv", header=1) )[,-2]
ordering.wobble <- read.table("data/SMoPT/example/input/S.cer.tRNA",header=F)[,c(1,4)]
colnames(ordering.wobble) <- c("c.number","wobble")

# reformat anticodons to match with my own labelling system (the identifier being t.number)
ordering.tRNA$anticodon <- gsub( pattern ="T" , replacement ="U", x= as.character(ordering.tRNA$anticodon) )
# reformat codons to match with my own labelling system (the identifier being c.number)
ordering.codons$codon <- tolower( as.character(ordering.codons$codon) )

# "rosetta" stone for SMoPT. The ambiguity about codon and anticodon indexes is over!
template <- arrange( merge( merge( subset(rosetta.codons, select=c(codon,anticodon)), ordering.codons, by="codon" ), ordering.tRNA, by = "anticodon", all.x=T), c.number )
write.table( template, file="data/SMoPT/stochastic.simulation/template.txt",sep="\t", quote = F, row.names= F)


# Generate new Scer.tRNA_copy_num (without the need to run the perl script, in fact)
wide.table.tRNA <- merge( merge( template, adjusted.tGCN, by = "anticodon" ), ordering.wobble, by = "c.number")

# Write the tables in their respective input files
dlply( melt(wide.table.tRNA, id.vars = c("codon","anticodon", "c.number", "t.number", "wobble") ), 
       .(variable), 
       function(x) { y <- arrange( x[,c("c.number","t.number", "value","wobble")], c.number)
       write.table(y, file = paste0("data/SMoPT/batch/input_GFPs/adjusted.tGCN/Scer.tRNA.",unique(x$variable),".txt") , row.names=F, sep="\t", col.names=F, quote=F)
       } )


# 2 # total number of tRNAs

tRNA_molecule.copy.number <- compile.tRNAabundance.master.table()
#tRNA_molecule.copy.number <- read.table("data/tRNA abundance/tRNA_molecule.copy.number.txt", header=1, sep="\t")
# this is the opportunity to compute more accurate N.total.tRNA values than with the approache based on the ratio 'sum tGCNs / sum tGCNn'

total.tRNAs <- ddply(tRNA_molecule.copy.number, .(experiment, time), summarise, N.total.tRNA = sum(total.tRNAab) )




#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3 # mRNA abundance data

mRNA.abundance <- read.table("data/SMoPT/batch/input_GFPs/mRNA_abundance.txt",header=1)
IniProb        <- read.table("data/SMoPT/batch/input_GFPs/S.cer.mRNA.abndc.ini.tsv",header=1)
to.get <- c("diauxic_20","ox_20","osm_20","temp_20")

wide.table.mRNA  <- cbind(IniProb, mRNA.abundance[, to.get ])
colnames(wide.table.mRNA) <- c("ORF","IniProb","diauxic_20","normal_0","ox_20","osm_20","temp_20")

head(wide.table.mRNA,15)
tail(wide.table.mRNA,15)


# # melted table giving the abundance of all transcripts
# mRNA.abundance.m <- melt(wide.table.mRNA, id.vars = c("ORF","IniProb"))
# subset( mRNA.abundance.m, value == 0) # there are some transcripts that have an abundance of 0. Quick fix: set them to 1 # long term: set them to 0, adjust genom table and downstream routines to ensure the right genes are mapped correctly.
# mRNA.abundance.m$value[ mRNA.abundance.m$value == 0 ] <- 1

# NEW JULY 2015
wide.table.mRNA[,to.get] <- sapply(wide.table.mRNA[,to.get], function(x){x[x == 0] <- 1; return(x)}) 

#-------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#   Pre-processing  #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# write the tables in their respective input files
# dlply( mRNA.abundance.m, 
#        .(variable), 
#        function(x) { y <- arrange( x[,c("ORF","IniProb", "value")], ORF)
#        colnames(y)[3] <- "rand_mRNA"
#        write.table(y, file = paste0("data/SMoPT/batch/input/mRNA.abundance/Scer.mRNAab.",unique(x$variable),".txt") , quote=F,row.names=F, sep="\t")
#        } )


# NEW JULY 2015

apply(conditions, 1, function(x){
  condition <- paste0(x[["experiment"]],"_",as.numeric(x[["time"]]),collapse = "")
  tag <- paste0(x[["experiment"]],"_",x[["GFP"]],"_",as.numeric(x[["time"]]),collapse = "")
  mRNA.table<- wide.table.mRNA[,c("ORF","IniProb", condition)]
  # cancel the abundance of GFPs
  mRNA.table[ mRNA.table$ORF %in% GFPs, condition ] <- 0
  transcriptome.size <- sum(mRNA.table[,condition])
  expression.GFP <- round( relative.GFP.level * min(transcriptome.size,60000) )
  mRNA.table[ mRNA.table$ORF == x[["GFP"]], condition ] <- expression.GFP
  table <- mRNA.table
  # now each table has a specific expression level for only the given GFP expressed in the given condition (5 GFPs * 5 conditions)
  colnames(table)[3] <- "rand_mRNA"
  write.table(table, file = paste0("data/SMoPT/batch/input_GFPs/mRNA.abundance/Scer.mRNAab.", tag,".txt") , quote=F,row.names=F, sep="\t")
})


# execute perl scripts to format FASTA files into .genom files
apply( conditions, 1, function(x){ 
  cat(x[["experiment"]],x[["GFP"]],x[["time"]],"\n")
  system( converts.mRNAab.files(experiment = x[["experiment"]], time = x[["time"]], GFP = x[["GFP"]]) ) 
  }  
)


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#     Execution     #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# generates a list of commands to execute, one for each experimental condition... and run it.

write(
  c( "#!/bin/tcsh",
     apply( conditions, 1, function(x){
       #total.tRNAs
       get.SMPoT.command(n.tRNAs =        subset( total.tRNAs, experiment == x[["experiment"]] & time == x[["time"]] )$N.total.tRNA, 
                         experiment = x[["experiment"]], 
                         GFP = x[["GFP"]],
                         n.mRNAs = nrow(mRNA.abundance),
                         duration = 1750, # +200% increase
                         heating  = 1000,
                         time = as.numeric(x[["time"]]), options = 1:4 )
     })), file = "data/SMoPT/batch/script_GFPs.sh" 
  )

# ONE LINE COMMAND:
# cd ~/Documents/MRC16 && ./data/SMoPT/batch/script_GFPs.sh

eGFP.master.table <- compile.master.table.new(in.path = "data/SMoPT/batch/input_GFPs/", out.path = "data/SMoPT/batch/output_GFPs/", ANALYSIS = "GFP")
require(tidyr)
eGFP.master.table <- separate(eGFP.master.table, col = "stress", into = c("condition","GFP","construct"), sep="_" )
eGFP.master.table$construct <- toupper(eGFP.master.table$construct)
write.table(eGFP.master.table, file="~/Documents/MRC16/results/SMoPT/eGFP.master.table.txt",row.names=F, sep="\t",quote=F)

eGFP.data <- subset(eGFP.master.table, name %in% GFPs & n.events > 0, select = c(name, condition, construct, n.events, av.initiation_time, av.elongation_time) )
eGFP.data <- arrange( eGFP.data, condition, construct)
eGFP.data <- subset(eGFP.data, construct != "NORMAL")

write.table(eGFP.data, file="~/Documents/MRC16/results/SMoPT/eGFP.data_2.txt",row.names=F, sep="\t",quote=F)
eGFP.data <- read.table("~/Documents/MRC16/results/SMoPT/eGFP.data.txt", header=1)

elongation_GFPs <- reshape2::dcast(eGFP.data, condition ~ construct, value.var = "av.elongation_time")
initiation_GFPs <- reshape2::dcast(eGFP.data, condition ~ construct, value.var = "av.initiation_time")
n.events_GFPs   <- reshape2::dcast(eGFP.data, condition ~ construct, value.var = "n.events")

require(pheatmap)
pheatmap( mat = 1/matrixify(elongation_GFPs)*60, scale = "column" )
pheatmap( mat = 1/matrixify(initiation_GFPs)*60, scale = "column" )

