#' Mehtylation distribution function
#'
#'Plot the methylation distribution for a few of the cells
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'mbiasplot(mbiasFile)
#'@export


mbiasplot<-function(mbiasFile){
  mbiasTable<-read.table(mbiasFile,header=TRUE)
  mbiasTable$methylation<-mbiasTable$nMethylated/(mbiasTable$nMethylated+mbiasTable$nUnmethylated)
  mbiasTable$Read<-as.factor(mbiasTable$Read)
  mbiasTable$Strand<-as.factor(mbiasTable$Strand)

  g<-ggplot(mbiasTable)+geom_line(aes(x=Position,y=methylation,colour=Read,linetype=Strand))+
    ylim(0,1)+ggtitle('Mbias Plot')

  return(g)


}

