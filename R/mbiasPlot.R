#' Methylation bias plot
#'
#'Plot the methylation at each position of the read to observe any biases in
#'the methylation
#'based on the read position
#'
#'@param dir directory name with mbias files
#'@param mbiasFiles list of mbias files
#'@return Plot showing the methylation across the read position in original top
#'and original bottom strand both in forward and reverse reads
#'@examples
#'mbiasFile <- '2017-04-21_HG23KBCXY_2_AGGCAGAA_TATCTC_pe.M-bias.txt'
#'mbiasplot(mbiasFiles=system.file("extdata",mbiasFile,package='scmeth'))
#'@importFrom utils read.table
#'@importFrom utils read.csv
#'@importFrom stats sd
#'@importFrom stats aggregate
#'@export

mbiasplot <- function(dir=NULL,mbiasFiles=NULL){
    if (!is.null(dir)){
        mbiasFileList <- list.files(dir,pattern="*.M-bias.txt",full.names=TRUE)

    }else{
        mbiasFileList <- mbiasFiles
    }

    mbiasTableList <- list()
    nSamples <- length(mbiasFileList)
    for (i in seq_len(nSamples)){
        MbiasFile <- readLines(mbiasFileList[i])
        nvec <- length(MbiasFile)
        breaks <- which(!nzchar(MbiasFile))
        nbreaks <- length(breaks)
        if (breaks[nbreaks] < nvec) {
            breaks <- c(breaks, nvec + 1L)
            nbreaks <- nbreaks + 1L
        }

    if (nbreaks > 0L) {
        oracle <- mapply(function(a,b) paste(MbiasFile[a:b]),
                        c(1L, 1L + breaks[-nbreaks])+2,
                        breaks - 1L)
    }

    CpG_Mbias_Read1 <- utils::read.csv(sep='\t',text=oracle[[1]])
    CpG_Mbias_Read1$read <- 'read-1'
    CpG_Mbias_Read2 <- utils::read.csv(sep='\t',text=oracle[[4]])
    CpG_Mbias_Read2$read <- 'read-2'

    mbiasTable <- rbind(CpG_Mbias_Read1[,c('position','X..methylation','read')],
                        CpG_Mbias_Read2[,c('position','X..methylation','read')])
    colnames(mbiasTable)<-c('position','methylation','read')
    mbiasTableList[[i]] <- mbiasTable
    }

    mt <- reshape2::melt(mbiasTableList,
                        id.vars=c('position', 'methylation', 'read'))
    mt$read_rep <- paste(mt$read, mt$L1, sep="_")

    meanTable <- stats::aggregate( methylation ~ position+read, data=mt, FUN=mean)
    sdTable <- stats::aggregate( methylation ~ position+read, data=mt, FUN=sd)
    seTable <- stats::aggregate( methylation ~ position+read, data=mt, FUN=function(x){sd(x)/sqrt(length(x))})

    sum_mt<-data.frame('position'=meanTable$position,'read'=meanTable$read,
                       'meth'=meanTable$X..methylation, 'sdMeth'=sdTable$X..methylation,
                       'seMeth'=seTable$X..methylation)
    #sum_mt <- plyr::ddply(mt, .(position, read), plyr::summarise,
                    #meth = mean(X..methylation), sdMeth = sd(X..methylation),
                    #seMeth = sd(X..methylation)/sqrt(length(X..methylation)))
    #sum_mt <- mt %>% dplyr::group_by(position,read) %>%
                        #dplyr::summarise(meth = mean(X..methylation),
                           #     sdMeth=stats::sd(X..methylation))
    #sum_mt$seMeth <- sum_mt$sdMeth/sqrt(nSamples)
    sum_mt$upperCI <- sum_mt$meth+(1.96*sum_mt$seMeth)
    sum_mt$lowerCI <- sum_mt$meth-(1.96*sum_mt$seMeth)
    sum_mt$read_rep <- paste(sum_mt$read, sum_mt$position,sep="_")

    g <- ggplot2::ggplot(sum_mt)
    g <- g+ggplot2::geom_line(ggplot2::aes_string(x='position',y='meth',
                                                colour='read'))
    g <- g+ggplot2::geom_ribbon(ggplot2::aes_string(ymin = 'lowerCI',
                                ymax = 'upperCI', x='position',fill = 'read'),
                                alpha=0.4)
    g <- g+ggplot2::ylim(0,100)+ggplot2::ggtitle('Mbias Plot')
    g <- g+ggplot2::ylab('methylation')
    return(g)
}
