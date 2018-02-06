#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest,
#'e.g. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Mmusculus.UCSC.mm10)
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'repMask(bs,Mmusculus,'mm10')
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


repMask <- function(bs,organism,genome){
    GenomeInfoDb::seqlevelsStyle(bs) <- "UCSC"
    hub <- AnnotationHub::AnnotationHub()
    repeatGr <- hub[[names(AnnotationHub::query(hub,
                        c("rmsk", GenomeInfoDb::organism(organism), genome)))]]
    rep <- GenomicRanges::countOverlaps(bs, repeatGr)>0
    cov <- bsseq::getCoverage(bs)
    covDf <- data.frame(coveredCpgs=DelayedArray::colSums(cov[!rep,]>=1))
    return(covDf)
}
