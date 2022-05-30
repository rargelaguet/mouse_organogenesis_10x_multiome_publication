#' Calculates similarity between two PFMs. 
#'
#' This function calculates the normalized motif correlation as a measure of motif
#' frequency matrix similarity. 
#'
#' This score is essentially a normalized version of the sum of column correlations
#' as proposed by Pietrokovski (1996). The sum is normalized by the average motif
#' length of m1 and m2, i.e. (ncol(m1)+ncol(m2))/2. Thus, for two idential motifs
#' this score is going to be 1. For unrelated motifs the score is going to be typically 
#' around 0. 
#'
#' Motifs need to aligned for this score to be calculated. The current implementation
#' tries all possible ungapped alignment with a minimal of two basepair matching, and 
#' the maximal score over all alignments is returned.
#'
#' Motif 1 is aligned both to Motif 2 and its reverse complement. Thus, the motif
#' similarities are the same if the reverse complement of any of the two motifs
#' is given.  
#'
#' @param m1 matrix with four rows representing the frequency matrix of first motif
#' @param m2 matrix with four rows representing the frequency matrix of second motif
#' @param trim bases with information content smaller than this value will be trimmed off both motif ends
#' @param self.sim if to calculate self similarity (i.e. without including offset=0 in alignment)
#' @references Pietrokovski S. Searching databases of conserved sequence regions by aligning protein multiple-alignments. Nucleic Acids Res 1996;24:3836-3845.
#' @export
#' @examples
#' 
#' if(requireNamespace("PWMEnrich.Dmelanogaster.background")){
#'    data(MotifDb.Dmel.PFM, package = "PWMEnrich.Dmelanogaster.background")
#'
#'    # calculate the similarity of tin and vnd motifs (which are almost identical)
#'    motifSimilarity(MotifDb.Dmel.PFM[["tin"]], MotifDb.Dmel.PFM[["vnd"]])
#'
#'    # similarity of two unrelated motifs
#'    motifSimilarity(MotifDb.Dmel.PFM[["tin"]], MotifDb.Dmel.PFM[["ttk"]])
#' }
motifSimilarity = function(m1, m2, trim=0.4, self.sim=FALSE){
  
  # remove any zero-variance columns
  colSd1 = apply(m1, 2, sd)
  colSd2 = apply(m2, 2, sd)
  if(any(colSd1 == 0)) m1 = m1[,colSd1!=0]
  if(any(colSd2 == 0)) m2 = m2[,colSd2!=0]
  
  # trim bases with small motif content
  bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
  
  ic1 <- colSums(m1*log2(m1/bg))
  ic2 <- colSums(m2*log2(m2/bg))
  
  r1 = range(which(ic1 >= trim))
  r2 = range(which(ic2 >= trim))
  
  if (sum(ic1>=trim)<=3)stop("Not enough information content in motif 1")
  if (sum(ic2>=trim)<=3)stop("Not enough information content in motif 1")
  
  m1 = m1[,r1[1]:r1[2]]
  m2 = m2[,r2[1]:r2[2]]
  
  # get max score
  if(self.sim){
    # we don't want offset=0 because then the score would always be 1
    score = tryAllMotifAlignments(m1, m2, exclude.zero=TRUE)
  } else {
    score = max(tryAllMotifAlignments(m1,m2), 
                tryAllMotifAlignments(m1,reverseComplement(m2)))
  }
  
  # normalize the max score by mean motif size
  return(score / mean(c(ncol(m1), ncol(m2))))
  
}



## functions to calculate motif similarity metrics

#' Returned the aligned motif parts
#'
#' This function takes the offset of first motif relative to second and 
#' chops off the end of both motifs that are not aligned. It returns a list
#' containing only the columns that align. 
#'
#' @param m1 frequency matrix of first motif
#' @param m2 frequency matrix of second motif
#' @param offset a number of nucleotides by which the first motif is offsetted compared to the second
#' @return a list of column-trimmed motifs m1, m2
maxAligned = function(m1, m2, offset){
  if(offset > 0){
    # need to cut off beginning of m2
    m2 = m2[,(offset+1):ncol(m2)]
  } else if(offset < 0){
    # net to cut off beginning of m1
    m1 = m1[,(abs(offset)+1):ncol(m1)]
  }
  
  # make sure we keep matrix 
  if(is.vector(m1))
    m1 = matrix(m1, ncol=1)
  if(is.vector(m2))
    m2 = matrix(m2, ncol=1)
  
  # maximal number of columns that align
  max.col = min(ncol(m1),ncol(m2))
  
  m1 = m1[,1:max.col]
  m2 = m2[,1:max.col]
  
  # make sure we keep it in matrix
  if(is.vector(m1))
    m1 = matrix(m1, ncol=1)
  if(is.vector(m2))
    m2 = matrix(m2, ncol=1)
  
  list(m1, m2)
}

#' Try all motif alignments and return max score
#'
#' This function tries all offsets of motif1 compared to motif2 and returns
#' the maximal (unnormalized) correlation score. 
#'
#' The correlation score is essentially the sum of correlations of individual aligned
#' columns as described in Pietrokovski (1996).
#'
#' @param m1 frequency matrix of motif 1
#' @param m2 frequency matrix of motif 2
#' @param min.align minimal number of basepairs that need to align
#' @param exclude.zero if to exclude offset=0, useful for calculating self-similarity
#' @references Pietrokovski S. Searching databases of conserved sequence regions by aligning protein multiple-alignments. Nucleic Acids Res 1996;24:3836-3845.
#' @return single maximal score
tryAllMotifAlignments = function(m1, m2, min.align=2, exclude.zero=FALSE){
  offset = (-ncol(m1)+min.align):(ncol(m2)-min.align)
  
  if(exclude.zero)
    offset = setdiff(offset, 0)
  
  scores = rep(0, length(offset))
  for(i in 1:length(offset)){
    # get aligned matrices
    aligned = maxAligned(m1, m2, offset[i])
    # convert to probabilities
    aligned.norm = lapply(aligned, function(x) divideRows(x, colSums(x)))
    # scores for first alignment (sum of column-wise correlations)
    scores[i] = sum(sapply(1:ncol(aligned.norm[[1]]), function(i) cor(aligned.norm[[1]][,i], aligned.norm[[2]][,i])))
    
  }
  
  max(scores)
}



#' Divide each row of a matrix with a vector
#'
#' @param m matrix to be divided
#' @param v the vector to use for division
divideRows = function(m, v){
  t(apply(m, 1, function(x) x/v))
}


