#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param bs bsseq object
#'@return Plot showing the mapped and unmapped read information for each cell
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'readmetrics(bs)
#'@export
#'@importFrom utils read.delim


readmetrics<-function(bs){
    phenptypicData<-Biobase::pData(bs)
    pDcolNames<-colnames(phenptypicData)
    if (('totalReads' %in% pDcolNames & 'mappedReads' %in% pDcolNames)){
        dat<-data.frame(sample=rownames(phenptypicData),
                        total=phenptypicData$totalReads,
                        mapped=phenptypicData$mappedReads)
        dat$unmapped<-dat$total-dat$mapped
        o<-order(dat$total,dat$sample)
        sampleOrder<-dat$sample[o]
        m<-reshape2::melt(dat[,c("sample","mapped","unmapped")],
                        id.vars="sample",variable.name="Mapping_status")
        m$sample<-factor(m$sample,levels=sampleOrder)
        g<-ggplot2::ggplot(m,ggplot2::aes_string('sample','value',
                                                fill='Mapping_status'))
        g<-g+ggplot2::geom_bar(stat="identity")+ggplot2::coord_flip()
        g<-g+ggplot2::scale_y_continuous(name="Number of reads")
        g<-g+ggplot2::xlab("samples")
        g<-g+ggplot2::ggtitle("Number of reads in Samples")
        g<-g+ggplot2::theme_bw()
        return(g)

  }else{
      warning('Read information not provided in the phenotypic data')
  }
}


