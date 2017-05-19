#' Coverage of CpGs in repeat masker region
#'
#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest, i.e. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Mmusculus.UCSC.mm10)
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'repMask(bs,Mmusculus,'mm10')
#'@export


repMask<-function(bs,organism,genome){
    hub <- AnnotationHub::AnnotationHub()
    repeatGr <- hub[[names(AnnotationHub::query(hub,
                      c("rmsk", GenomeInfoDb::organism(organism), genome)))]]
    rep <- GenomicRanges::countOverlaps(bs, repeatGr)>0
    cov<-bsseq::getCoverage(bs)
    covDf <- data.frame(coveredCpgs=colSums(cov[!rep,]>=1))
    return(covDf)
}
