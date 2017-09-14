#' Methylation bias plot
#'
#'Plot the methylation at each position of the read to observe any biases in
#'the methylation
#'based on the read position
#'
#'@param mbiasFile methylation bias file provided by methyldackel/pileometh
#'@return Plot showing the methylation across the read position in original top
#'and original bottom strand both in forward and reverse reads
#'@examples
#'mbiasFile<-'2017-04-21_HG23KBCXY_2_AGGCAGAA_TATCTC_pe.M-bias.txt'
#'mbiasplot(system.file("extdata",mbiasFile,package='scmeth'))
#'@importFrom utils read.table
#'@importFrom utils read.csv
#'@export



mbiasplot<-function(mbiasFile){
    MbiasFile<-readLines(mbiasFile)

    nvec<-length(MbiasFile)
    breaks<-which(!nzchar(MbiasFile))
    nbreaks<-length(breaks)
    if (breaks[nbreaks] < nvec) {
      breaks <- c(breaks, nvec + 1L)
      nbreaks <- nbreaks + 1L
    }

    if (nbreaks > 0L) {
      oracle <- mapply(function(a,b) paste(MbiasFile[a:b]),
                       c(1L, 1L + breaks[-nbreaks])+2,
                       breaks - 1L)
    }

    CpG_Mbias_Read1<-utils::read.csv(sep='\t',text=oracle[[1]])
    CpG_Mbias_Read1$read<-'read-1'
    CpG_Mbias_Read2<-utils::read.csv(sep='\t',text=oracle[[4]])
    CpG_Mbias_Read2$read<-'read-2'

    mbiasTable<-rbind(CpG_Mbias_Read1[,c('position','X..methylation','read')],
                      CpG_Mbias_Read2[,c('position','X..methylation','read')])


    g<-ggplot2::ggplot(mbiasTable)
    g<-g+ggplot2::geom_line(ggplot2::aes_string(x='position',
                            y='X..methylation',colour='read'))
    g<-g+ggplot2::ylim(0,100)+ggplot2::ggtitle('Mbias Plot')

    return(g)
}

