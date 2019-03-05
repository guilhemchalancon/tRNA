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
  fasta.file = "data/SMoPT/batch/input/S.cer.fasta",
  genetic.code.file = "data/SMoPT/utilities/genetic.code.tsv",
  experiment = "ox",
  time = "20",
  output.path = "data/SMoPT/batch/input/genomes/"
  ){
  mRNAab.file <- paste0("data/SMoPT/batch/input/mRNA.abundance/Scer.mRNAab.",experiment,"_",time,".txt") # note that the mRNA abundance files must preexist (another function to build)
  output.file <- paste0(output.path,"Scer.",experiment,"_",time)
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
  sequences.path = "data/SMoPT/batch/input/genomes/",
  tRNA.info.path = "data/SMoPT/batch/input/adjusted.tGCN/",
  seed = 1413,
  experiment = "normal",
  time = "0",
  output.path.root = "data/SMoPT/batch/output/",
  options = c(1:5)
){
  sequences   = paste0(sequences.path, "Scer.", experiment, "_", time, ".genom", collapse = "") # Scer.genome.normal_0.txt
  tRNA.info   = paste0(tRNA.info.path, "Scer.tRNA.",   experiment, "_", time, ".txt", collapse = "") # Scer.tRNA.normal_0.txt
  time.folder = sprintf("%03.0f", as.numeric(time)) # e.g. gets 0 written as 000, 20 as 020 and 120 as 120
  output.path = paste0(output.path.root,  experiment, "/", time.folder, "/Scer" )
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

# 0 # define which conditions to compile
conditions <- rbind( arrange( data.frame( experiment = c("diauxic","osm","ox","temp"), time = rep(c(20,120),each = 4), stringsAsFactors=F ), experiment), c(experiment = "normal", time = 0) )


#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1 # gcn data

  # stress-adjusted tGCN values (precomputed)
  adjusted.tGCN  <- arrange( read.table( "results/stress tAI/adjusted.tGCN.txt", header =1  ), anticodon)
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
                     write.table(y, file = paste0("data/SMoPT/batch/input/adjusted.tGCN/Scer.tRNA.",unique(x$variable),".txt") , row.names=F, sep="\t", col.names=F, quote=F)
        } )

# 2 # total number of tRNAs

tRNA_molecule.copy.number <- compile.tRNAabundance.master.table()
#tRNA_molecule.copy.number <- read.table("data/tRNA abundance/tRNA_molecule.copy.number.txt", header=1, sep="\t")
# this is the opportunity to compute more accurate N.total.tRNA values than with the approache based on the ratio 'sum tGCNs / sum tGCNn'

total.tRNAs <- ddply(tRNA_molecule.copy.number, .(experiment, time), summarise, N.total.tRNA = sum(total.tRNAab) )




#-------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3 # mRNA abundance data

  mRNA.abundance <- read.table("data/SMoPT/stochastic.simulation/mRNA_abundance.txt",header=1)
  IniProb        <- read.table("data/SMoPT/example/input/S.cer.mRNA.abndc.ini.tsv",header=1)
  to.get <- c("diauxic_20","diauxic_120","temp_20","temp_120","osm_20","osm_120","ox_20","ox_120")

  wide.table.mRNA  <- cbind(IniProb, mRNA.abundance[, to.get ])
  colnames(wide.table.mRNA) <- c("ORF","IniProb","normal_0","diauxic_20","diauxic_120","temp_20","temp_120","osm_20","osm_120","ox_20","ox_120")

head(wide.table.mRNA,15)




mRNA.abundance.m <- melt(wide.table.mRNA, id.vars = c("ORF","IniProb"))
subset( mRNA.abundance.m, value == 0) # there are some transcripts that have an abundance of 0. Quick fix: set them to 1 # long term: set them to 0, adjust genom table and downstream routines to ensure the right genes are mapped correctly.
mRNA.abundance.m$value[ mRNA.abundance.m$value == 0 ] <- 1


#-------------------------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#--------------------------------#   Pre-processing  #---------------------------------------------------------------------------------------------------
#--------------------------------#-#-#-#-#-#-#-#-#-#-#--------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# write the tables in their respective input files
dlply( mRNA.abundance.m, 
       .(variable), 
       function(x) { y <- arrange( x[,c("ORF","IniProb", "value")], ORF)
                     colnames(y)[3] <- "rand_mRNA"
                     write.table(y, file = paste0("data/SMoPT/batch/input/mRNA.abundance/Scer.mRNAab.",unique(x$variable),".txt") , quote=F,row.names=F, sep="\t")
       } )


# execute perl scripts to format FASTA files into .genom files
apply( conditions, 1, function(x){ 
   cat(x[[1]],x[[2]],"\n")
   system( converts.mRNAab.files(experiment = x[[1]], time = x[[2]]) ) }  )


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
apply( total.tRNAs, 1, function(x){
  get.SMPoT.command(n.tRNAs = x[["N.total.tRNA"]], 
                    experiment = x[["experiment"]], 
                    duration = 1500, # +200% increase
                    heating  = 1000,
                    time = as.numeric(x[["time"]]), options = 1:4 )
})), file = "data/SMoPT/batch/script.sh" )

