#' Combine multiple bsseq rda object
#'
#'
#'Combine multiple bsseq rda object into one bs object
#'Preferably for a library pool of cells
#'
#'@param rdaList list of bsseq rda objects
#'@return Combined one rda object for all the cells
#'@examples
#'CpGBedGraphFile_1<-system.file("extdata",'sc-RRBS_zyg_01_chr1_CpG.bedGraph',package='scmeth')
#'readMetricsFile_1<-system.file("extdata",'sc-RRBS-zygote_01.read_metrics.txt',package='scmeth')
#'bsConversionFile_1<-system.file("extdata",'sc-RRBS-zygote_01.bsConv.txt',package='scmeth')
#'CpGBedGraphFile_2<-system.file("extdata",'sc-RRBS_zyg_02_chr1_CpG.bedGraph',package='scmeth')
#'readMetricsFile_2<-system.file("extdata",'sc-RRBS-zygote_02.read_metrics.txt',package='scmeth')
#'bsConversionFile_2<-system.file("extdata",'sc-RRBS-zygote_02.bsConv.txt',package='scmeth')
#'bsObject1<-createRDA(CpGBedGraphFile_1,readMetricsFile_1,bsConversionFile_1)
#'bsObject2<-createRDA(CpGBedGraphFile_2,readMetricsFile_2,bsConversionFile_2)
#'combineRDA(c(bsObject1, bsObject2))
#'@export


combineRDA<-function(rdaList){
    message("Combining bs object")
    bs<-Reduce(BiocGenerics::combine,rdaList)
    bs@parameters<-rdaList[[1]]@parameters

    return(bs)
}
