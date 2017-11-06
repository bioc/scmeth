#' Generates an inclusive report on methylation analysis
#'
#' This function uses most of the functions in this package to
#' generate a report for the user
#'@param bsObj bsseq object
#'@param outdirectory name of the output directory where the
#'report will be saved
#'@param organism scientific name of the organism of interest,
#'i.e. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'the report will have graphics on read information
#'@param mbiasDir Optional argument to provide directory name
#'that has the mbias files or the list of mbias files
#'@return Report will be an html file
#'@examples
#'library(BSgenome.Hsapiens.NCBI.GRCh38)
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'mbiasDirectory=system.file("extdata",package='scmeth')
#'report(bs,'~',Hsapiens,'hg38',mbiasDir=mbiasDirectory)
#'@import knitr
#'@import DT
#'@import SummarizedExperiment
#'@export
#

report <- function(bsObj,outdirectory,organism,genome,mbiasDir=NULL) {
    RmdFile<-system.file(".",'qcReport.Rmd',package="scmeth")
    rmarkdown::render(RmdFile,params=list(outdir=outdirectory,samples=bsObj
                                        ,organism=organism,genome=genome,mbias=mbiasDir)
                    ,output_file=paste0(outdirectory,"/qcReport.html"))

}
