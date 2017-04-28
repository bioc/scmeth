#' Methylation bias plot
#'
#'Plot the methylation at each position of the read to observe any biases in the methylation
#'based on the read position
#'
#'@param mbiasFile methylation bias file provided by methyldackel/pileometh
#'@return Plot showing the methylation across the read position in original top
#'and original bottom strand both in forward and reverse reads
#'@examples
#'mbiasplot(system.file("extdata",'16_trimmed_sorted.txt',package='scmeth'))
#'@export
#'@importFrom utils read.table


mbiasplot<-function(mbiasFile){
  mbiasTable<-read.table(mbiasFile,header=TRUE)
  mbiasTable$methylation<-mbiasTable$nMethylated/(mbiasTable$nMethylated+mbiasTable$nUnmethylated)
  mbiasTable$Read<-as.factor(mbiasTable$Read)
  mbiasTable$Strand<-as.factor(mbiasTable$Strand)
  g<-ggplot2::ggplot(mbiasTable)+ggplot2::geom_line(ggplot2::aes_string(x='Position',y='methylation',colour='Read',linetype='Strand'))+
    ggplot2::ylim(0,1)+ggplot2::ggtitle('Mbias Plot')

  return(g)


}

