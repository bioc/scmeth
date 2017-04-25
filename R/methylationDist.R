#' Provide graphics for methylation distribution
#'
#'Plot the methylation distribution for the cells in bsseq object
#'@param bs bsseq object
#'@param all Indicator of whether to plot the distribution for all the cells or for
#'a single cell
#'@return plot of the methylation distribution
#'@examples
#'methylationDist(bsObject)
#'methylationDist(bsObject,all = TRUE)
#'@export


methylationDist<-function(bs,all=FALSE){
  covMatrix<-bsseq::getCoverage(bs)
  methMatrix<-bsseq::getCoverage(bs,type='M')/covMatrix
  df <- as.data.frame(matrix(unlist(methMatrix), nrow = nrow(methMatrix)))
  colnames(df)<-colnames(methMatrix)

  if (all==TRUE){
    meltedDf<-reshape2::melt(df)
    if (requireNamespace("ggplot2",quietly = TRUE)){
    methylationDensityPlot<-ggplot2::ggplot(meltedDf)+geom_density(aes(x=value))+facet_wrap(~variable,ncol=3)+
                              ggtitle('Methylation distribution for all the cells')+xlab('Methylation')

    return(methylationDensityPlot)
    }else{
      warning('ggplot2 required for plot rendering')
    }

  }else if (all==FALSE){
    indCell<-data.frame(x=df[,sample(ncol(df),1)])
    if (requireNamespace("ggplot2",quietly = TRUE)){
    methylationDensityPlot<-ggplot2::ggplot()+geom_density(aes(x=x),data=indCell)+
                            ggtitle('Methylation Distribution for an arbitrary cell')+xlab('Methylation rate')


    return(methylationDensityPlot)

    }else{
      warning('ggplot2 required for plot rendering')
    }

  }
}

