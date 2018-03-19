#'Provides Coverage by the CpG density. CpG Density is defined as the number
#'of CpGs observed in certain base pair long region.
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest,
#'e.g. Mmusculus or Hsapiens
#'@param windowLength Length of the window to calculate the density
#'Default value for window length is 1000 basepairs.
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Hsapiens.NCBI.GRCh38)
#'directory <- system.file("extdata/bismark_data",package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'memory.limit(size=300)
#'cpgDensity(bs,Hsapiens,1000)
#'@import BSgenome
#'@importFrom bsseq getCoverage
#'@importFrom Biostrings DNAString
#'@importFrom Biostrings vmatchPattern
#'@export

cpgDensity <- function(bs,organism,windowLength=1000){

    #GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(organism)[1]
    cov <- bsseq::getCoverage(bs)
    gr <- GenomicRanges::granges(bs)
    cpg <- Biostrings::DNAString("CG")
    cpg_gr <- Biostrings::vmatchPattern(cpg, organism)
    cpg_gr <- GenomeInfoDb::keepStandardChromosomes(cpg_gr, pruning.mode= "coarse")

    r_cpg_gr <- GenomicRanges::resize(cpg_gr,width=(windowLength/2),fix='center')
    cpgd <- GenomicRanges::countOverlaps(gr,r_cpg_gr)

    #cpgd <- Repitools::cpgDensityCalc(gr, organism, window = windowLength)

    maxcpgd <- max(cpgd)
    cpgdCov <- sapply(seq_len(ncol(cov)), function(i) {
        cv = as.vector(cov[,i])
        cpgdCell <- cpgd[cv>0 ]
        tab <- table(cpgdCell)
        x <- rep(0, maxcpgd)
        x[as.numeric(names(tab))] <- tab
        x
    })

    rownames(cpgdCov) <- seq_len(maxcpgd)
    return(cpgdCov)
}
