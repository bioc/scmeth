#' Provide graphics for methylation distribution
#'
#'Plot the methylation distribution for the cells in bsseq object
#'@param bs bsseq object
#'@param all Indicator of whether to plot the distribution for all the
#'cells or for
#'a single cell
#'@return plot of the methylation distribution
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'methylationDist(bs)
#'methylationDist(bs,all = TRUE)
#'@export


methylationDist<-function(bs,all=FALSE){
    covMatrix<-bsseq::getCoverage(bs)
    methMatrix<-bsseq::getCoverage(bs,type='M')/covMatrix
    df <- as.data.frame(matrix(unlist(methMatrix), nrow = nrow(methMatrix)))
    colnames(df)<-colnames(methMatrix)

    if (all==TRUE){
        meltedDf<-reshape2::melt(df)
    if (requireNamespace("ggplot2",quietly = TRUE)){
        g<-ggplot2::ggplot(meltedDf)
        g<-g+ggplot2::geom_density(ggplot2::aes_string(x='value'))
        g<-g+ggplot2::facet_wrap(~variable,ncol=3)
        g<-g+ggplot2::ggtitle('Methylation distribution for all the cells')
        g<-g+ggplot2::xlab('Methylation')

    return(g)
    }else{
        warning('ggplot2 required for plot rendering')
    }

  }else if (all==FALSE){
      indCell<-data.frame(x=df[,sample(ncol(df),1)])
      if (requireNamespace("ggplot2",quietly = TRUE)){
          g<-ggplot2::ggplot()+ggplot2::geom_density(ggplot2::aes_string(x='x'),data=indCell)
          g<-g+ggplot2::ggtitle('Methylation Distribution for an arbitrary cell')
          g<-g+ggplot2::xlab('Methylation rate')

      return(g)

    }else{
        warning('ggplot2 required for plot rendering')
    }

  }
}

