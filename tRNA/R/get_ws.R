#'
#' @param tRNA
#' @param s
#' @param trim_above_ceiling
get_ws <- function(tRNA, 
                   s = c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68), 
                   trim_above_ceiling = T,
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
  
  if( trim_above_ceiling == T){
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