#' Mehtylation distribution function
#'
#'Plot the methylation distribution for a few of the cells
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@example
#'@export


methylationDist<-function(bs,all=FALSE){
  covMatrix<-getCoverage(bs)
  methMatrix<-getCoverage(bs,type='M')/covMatrix
  if (all== TRUE){

    methylationDensityPlot<-ggplot(methMatrix[,sample(ncol(methMatrix),1)])+geom_density()
    #d<-density(methMatrix[,sample(ncol(methMatrix),1)],na.rm=TRUE)
    #methylationDensityPlot<-plot(d,main = "Methylation distribution for a single cell \n in the library")
    return(methylationDensityPlot)

  }else{
    methylationDensityPlot<-ggplot(methMatrix[,sample(ncol(methMatrix),1)])+geom_density()
    #d<-density(methMatrix[,sample(ncol(methMatrix),1)],na.rm=TRUE)
    #methylationDensityPlot<-plot(d,main = "Methylation distribution for a single cell \n in the library")

    return(methylationDensityPlot)
  }

}

