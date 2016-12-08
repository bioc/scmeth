#' square a number
#'
#'Take in a numeric value and squares it
#'@param x A numeric input
#'@return The sqaure of the output
#'@example
#'@export

coverage <- function(bs) {
  covMatrix<-getCoverage(bs)
  covVec<- colSums(covMatrix>0,na.rm=TRUE)
  return(covVec)
}
