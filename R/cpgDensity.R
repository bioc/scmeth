#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest,
#'i.e. Mus musculus or Homo sapiens
#'@param windowLength Length of the window to calculate the density
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Hsapiens.NCBI.GRCh38)
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'cpgDensity(bs,Hsapiens,1000)
#'@import BSgenome
#'@importFrom bsseq getCoverage
#'@export



cpgDensity<-function(bs,organism,windowLength=1000){

    cov<-bsseq::getCoverage(bs)
    gr <- granges(bs)
    GenomeInfoDb::seqlevelsStyle(gr) <-GenomeInfoDb::seqlevelsStyle(organism)[1]
    cpgd<-Repitools::cpgDensityCalc(gr, organism, window = windowLength)
    #cpgdBin<-cut(cpgd,seq(0,max(cpgd),5)

    # Need to find a better way to conduct this on-disk
    cpgdCov <- by(cov>0, cpgd, colSums)
    cpgdCov <- do.call("rbind", cpgdCov)
    return(cpgdCov)

}
