#' Coverage statistics
#'
#'Provides Coverage metrics for the sample
#'@param takes bs object as input
#'@return vector of numbers as output
#'@examples coverage(bsseqObject)
#'@export

coverage <- function(bsObject) {
  load(bsObject)
  covMatrix<-getCoverage(bs)
  covVec<- colSums(covMatrix>0,na.rm=TRUE)
  return(covVec)
}
