#' Methylation bias plot
#'
#'Plot the methylation at each position of the read to observe any biases in
#'the methylation
#'based on the read position
#'
#'@param dir directory name with mbias files
#'@param mbiasFiles list of mbias files
#'@return Returns a list containing the methylation across the read position in original top
#'and original bottom strand both in forward and reverse reads for multiple samples
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
    CpG_Mbias_Read1$methylation <- CpG_Mbias_Read1$count.methylated*100/
      (CpG_Mbias_Read1$count.methylated+CpG_Mbias_Read1$count.unmethylated)
    CpG_Mbias_Read2 <- utils::read.csv(sep='\t',text=oracle[[4]])
    CpG_Mbias_Read2$read <- 'read-2'
    CpG_Mbias_Read2$methylation <- CpG_Mbias_Read2$count.methylated*100/
      (CpG_Mbias_Read2$count.methylated+CpG_Mbias_Read2$count.unmethylated)


    mbiasTable <- rbind(CpG_Mbias_Read1[,c('position','methylation','read')],
                        CpG_Mbias_Read2[,c('position','methylation','read')])
    mbiasTableList[[i]] <- mbiasTable
    }

    return(mbiasTableList)
}
