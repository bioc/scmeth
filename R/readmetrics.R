#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param readData a .txt file providing mapped and unmapped reads for each sample
#'@return Plot showing the mapped and unmapped read information for each cell
#'@export
#'@importFrom utils read.delim


readmetrics<-function(readData){

  dat<-read.delim(readData,header=TRUE,sep='',row.names=1)
  dat$sample<-rownames(dat)
  dat$unmapped<-dat$total-dat$mapped
  o<-order(dat$total,dat$sample)
  sampleOrder<-dat$sample[o]
  m<-reshape2::melt(dat[,c("sample","mapped","unmapped")],id.vars="sample",variable.name="Mapping_status")
  m$sample<-factor(m$sample,levels=sampleOrder)
  if (requireNamespace("ggplot2",quietly = TRUE)){
    g<-ggplot2::ggplot(m,ggplot2::aes_string('sample','value',fill='Mapping_status'))+ggplot2::geom_bar(stat="identity")+ggplot2::coord_flip()+
      ggplot2::scale_y_continuous(name="Number of reads")+ggplot2::xlab("samples")+
      ggplot2::ggtitle("Number of reads in Samples")+ggplot2::theme_bw()
  }


    return(g)
  }



