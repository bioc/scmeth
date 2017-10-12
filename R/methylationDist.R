#' Provide graphics for methylation distribution
#'
#'Plot the methylation distribution for the cells in bsseq object
#'@param bs bsseq object
#'@return plot of the methylation distribution
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'methylationDist(bs)
#'@importFrom bsseq getCoverage
#'@importFrom ggthemes theme_tufte
#'@export


methylationDist<-function(bs){
    covMatrix<-bsseq::getCoverage(bs)
    methMatrix<-bsseq::getCoverage(bs,type='M')/covMatrix
    nSamples<-ncol(methMatrix)

    methylationCutOff<-c(0,0.2,0.4,0.6,0.8,1.0)
    methylIntervals<-length(methylationCutOff)-1
    totCpGs<-DelayedArray::colSums(!is.na(methMatrix))


    methylationDistMatrix<-matrix(nrow=nSamples,ncol=methylIntervals)
    for (i in 1:methylIntervals){
      if (i==1){
        methylationDistMatrix[,i]<-DelayedArray::colSums(methMatrix>=methylationCutOff[i] &
                                                           methMatrix<methylationCutOff[i+1],na.rm=TRUE)
      } else {
      methylationDistMatrix[,i]<-DelayedArray::colSums(methMatrix>methylationCutOff[i] &
                                              methMatrix<=methylationCutOff[i+1],na.rm=TRUE)
      }
    }

    methylationDistMatrix<-apply(methylationDistMatrix,2,function(x) x/totCpGs)
    orderdMeth<-order(methylationDistMatrix[,1])
    methylationDistMatrix<-methylationDistMatrix[orderdMeth,]

    colnames(methylationDistMatrix)<-c('[0,0.2]','(0.2,0.4]','(0.4,0.6]','(0.6,0.8]','(0.8,1]')
    meltedMDistMatrix<-reshape2::melt(methylationDistMatrix)


    g <- ggplot2::ggplot(meltedMDistMatrix, ggplot2::aes_string(x='Var2', y='Var1',fill='value'))
    g <- g + ggplot2::geom_tile(color="white", size=0.1)
    g <- g + viridis::scale_fill_viridis(name="Fraction of \n CpGs observed", label=scales::comma)
    #g <- g + ggplot2::coord_equal() # Adds same aspect ratio
    g <- g + ggplot2::labs(x=NULL, y=NULL, title="Methylation Distribution")
    g <- g + ggthemes::theme_tufte(base_family="Helvetica")
    g <- g + ggplot2::theme(plot.title=ggplot2::element_text(hjust=0))
    g <- g + ggplot2::theme(axis.text.y=ggplot2::element_blank())
    g <- g + ggplot2::theme(axis.ticks=ggplot2::element_blank()) # get rid of tick marks
    g <- g + ggplot2::theme(axis.text=ggplot2::element_text(size=12)) # Change the font size
    g <- g + ggplot2::theme(legend.title=ggplot2::element_text(size=12))
    g <- g + ggplot2::theme(legend.text=ggplot2::element_text(size=10))

    return(g)

}

