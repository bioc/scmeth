#' Generates an inclusive report on methylation analysis
#'
#' This function uses most of the functions in this package to
#' generate a report for the user
#'@param bsObj bsseq object
#'@param outdirectory name of the output directory where the
#'report will be saved
#'@param organism scientific name of the organism of interest,
#'e.g. Mmusculus or Hsapiens
#'@param genome reference alignment, e.g. mm10 or hg38
#'the report will have graphics on read information
#'@param mbiasDir Optional argument to provide directory name
#'that has the mbias files or the list of mbias files
#'@param subSample number of CpGs to subsample
#'Default value is 1000000.
#'@param offset how many CpGs to offset when subsampling
#'Default value is set to be 50000, i.e. first 50000 CpGs will
#'be ignored in subsampling.
#'@param small Indicator for a small dataset, cpg density is calculated more
#'@return Report will be an html file
#'@examples
#'library(BSgenome.Hsapiens.NCBI.GRCh38)
#'directory <- system.file("extdata/bismark_data", package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'mbiasDirectory=system.file("extdata", package='scmeth')
#'outDir <- system.file(package='scmeth')
#'report(bs, outDir, Hsapiens, 'hg38', mbiasDir=mbiasDirectory, small=TRUE)
#'@importFrom HDF5Array loadHDF5SummarizedExperiment
#'@import knitr
#'@import DT
#'@import SummarizedExperiment
#'@export
#

report <- function(bsObj, outdirectory, organism, genome, mbiasDir=NULL,
                   subSample=1e6, offset=50000, small=FALSE) {

    RmdFile <- system.file(".", 'qcReport.Rmd', package="scmeth")
    rmarkdown::render(RmdFile, params=list(outdir=outdirectory, samples=bsObj,
                                        organism=organism, genome=genome,
                                        mbias=mbiasDir, nCpGs=subSample,
                                        offset=offset, small=small),
                     knit_root_dir=getwd(),
                    output_file=paste0(outdirectory, "/qcReport.html"))

}
