#' Mehtylation distribution function
#'
#'Plot the methylation distribution for a few of the cells
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@example
#'@export


methylationDist<-function(bs){
  covMatrix<-getCoverage(bs)
  methMatrix<-getCoverage(bs,type='M')/covMatrix
  d<-density(methMatrix[,sample(ncol(methMatrix),1)],na.rm=TRUE)
  methylationDensityPlot<-plot(d,main = "Methylation distribution for a single cell \n in the library")

  return(methylationDensityPlot)

}
