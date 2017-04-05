#' Coverage for single cells
#'
#'Provides Coverage for each cell in a library pool
#'@param bs bsseq object
#'@return vector of coverage for the cells in bs object
#'@examples coverage(bsseqObject)
#'@export

coverage <- function(bs) {

  covMatrix<-bsseq::getCoverage(bs)
  covVec<- colSums(covMatrix>0,na.rm=TRUE)
  return(covVec)
}
