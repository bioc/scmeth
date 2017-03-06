#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param readData a .txt file providing mapped and unmapped reads for each sample
#'@return Plot showing the mapped and unmapped read information for each cell
#'@examples
#'readmetrics(read.txt)
#'@export


readmetrics<-function(readData){

  dat<-read.delim(readData,header=TRUE,sep='',row.names=1)
  dat$sample<-rownames(dat)
  dat$unmapped<-dat$total-dat$mapped
  o<-order(dat$total,dat$sample)
  sampleOrder<-dat$sample[o]
  m<-melt(dat[,c("sample","mapped","unmapped")],id.vars="sample",variable.name="Mapping_status")
  m$sample<-factor(m$sample,levels=sampleOrder)
  g<-ggplot2::ggplot(m,aes(sample,value,fill=Mapping_status))+geom_bar(stat="identity")+coord_flip()+
    scale_y_continuous(name="Number of reads")+xlab("samples")+ggtitle("Number of reads in Samples")+
    theme_bw()


    return(g)
  }



