read_fasta <- function(filepath="data/Scer.GFP.fasta", type='DNA'){
  
  if(type=='DNA'){
    return(Biostrings::readDNAStringSet(filepath=filepath)) 
  }

  if(type == 'RNA'){
    return(Biostrings::readRNAStringSet(filepath=filepath))
  }
  
}