#' CpG covergae in each chromosome
#'
#'Provides Coverage metrics for each sample by the chromosome
#'@param bs bsseq object
#'@return matrix of chromsome covergae with
#'column and rows indicating the samples and the chromosome respectively
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'chromosomeCoverage(bs)
#'@importFrom bsseq getCoverage
#'@export


chromosomeCoverage <- function(bs) {
    bs<-GenomeInfoDb::keepStandardChromosomes(bs)
    covMatrix<-bsseq::getCoverage(bs)
    Granges<-GenomicRanges::granges(bs)
    standardChr <- GenomeInfoDb::seqnames(Granges)
    chrCov <- by(covMatrix>0, standardChr, colSums)
    chrCov <- do.call("rbind", chrCov)
    return(chrCov)

}
