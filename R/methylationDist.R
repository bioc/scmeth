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
  df <- as.data.frame(matrix(unlist(methMatrix), nrow = nrow(methMatrix)))

  if (all==TRUE){
    meltedDf<-melt(df)
    methylationDensityPlot<-ggplot(meltedDf)+geom_density(aes(x=value))+facet_wrap(~variable,ncol=3)+
                              ggtitle('Methylation distribution for all the cells')+xlab('Methylation')

    return(methylationDensityPlot)

  }else if (all==FALSE){
    indCell<-data.frame(x=df[,sample(ncol(df),1)])
    methylationDensityPlot<-ggplot()+geom_density(aes(x=x),data=indCell)+
                            ggtitle('Methylation Distribution for an arbitrary cell')+xlab('Methylation rate')


    return(methylationDensityPlot)
  }

}

