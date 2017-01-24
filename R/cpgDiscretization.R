#' Discretize the CpG methylation values
#' to align with single cell analysis
#'
#'Take in a numeric value and squares it
#'@param rda as an input
#'@return discretized methylation matrix and a plot showing how many of the CpGs are discarded
#'@example
#'@export
#

cpgDiscretization<-function(bs){
  covMatrix<-getCoverage(bs)
  methMatrix<-getCoverage(bs,type='M')
  methMatrix<-methMatrix/covMatrix
  tempMethylationMatrix<-methMatrix
  methMatrix[methMatrix<=0.2]<-0
  methMatrix[methMatrix>=0.8]<-1
  removedCpGs<-colSums(methMatrix>0.2 & methMatrix<0.8, na.rm=TRUE)
  removedCpGFrac<-(removedCpGs/(scmeth::coverage(bs)))*100
  returnList<-list('meth' = methMatrix, 'discard' = removedCpGs,'discard-perc' = removedCpGFrac)
  return(returnList)

}
