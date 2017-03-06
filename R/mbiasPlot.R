#' Mehtylation bias plot
#'
#'Plot the methylation at each position of the read to observe any biases in the methylation
#'based on the read position
#'
#'@param mbiasFile methylation bias file provided by methyldackel/pileometh
#'@return Plot showing the methylation across the read position in original top
#'and original bottom strand both in forward and reverse reads
#'@examples
#'mbiasplot(mbiasFile)
#'@export


mbiasplot<-function(mbiasFile){
  mbiasTable<-read.table(mbiasFile,header=TRUE)
  mbiasTable$methylation<-mbiasTable$nMethylated/(mbiasTable$nMethylated+mbiasTable$nUnmethylated)
  mbiasTable$Read<-as.factor(mbiasTable$Read)
  mbiasTable$Strand<-as.factor(mbiasTable$Strand)

  g<-ggplot2::ggplot(mbiasTable)+geom_line(aes(x=Position,y=methylation,colour=Read,linetype=Strand))+
    ylim(0,1)+ggtitle('Mbias Plot')

  return(g)


}

