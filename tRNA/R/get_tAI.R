#'
#' @param seq a FASTA sequence
#' @param weights 
#' @param score
#' @return a data.table with stAI scores 
get_tAI <- function(seq, weights, score){
  codons  <- colnames(x)
  w <- weights[,score]
  names(w) <- weights[,"codon"]
  w <- w[codons]
  w = log(w)      #calculate log of w
  n = apply(x,1,'*',w)		#multiply each row of by the weights
  n = t(n)                      #transpose
  n = apply(n,1,sum)		#sum rows
  L = apply(x,1,sum)		#get each ORF length
  tAI = exp(n/L)		#get tai
  return(tAI)
}