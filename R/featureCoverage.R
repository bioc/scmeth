#' Coverage based on the genomic feature
#'
#'Provides Coverage metrics for the sample by each genomic features provided
#'by the user
#'@param bs bsseq object
#'@param features list of genomic features, i.e. genes_exons, genes_introns,
#'cpg_islands, cpg_shelves
#'Names are based on the annotatr packages, so all the features provided by the
#'annotatr
#'package will be supported in this function
#'@param genomebuild reference alignment, i.e. mm10 or hg38
#'@return a data frame with genomic feature names and the number of
#'CpG covered in each feature
#'@examples
#'library(annotatr)
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'featureCoverage(bs,c('cpg_islands','genes_exons'),'hg38')
#'@importFrom DelayedArray rowSums
#'@importFrom GenomeInfoDb seqlevelsStyle
#'@importFrom annotatr builtin_genomes
#'@importFrom annotatr build_annotations
#'@importFrom annotatr annotate_regions
#'@importFrom annotatr summarize_annotations
#'@import GenomicRanges
#'@export


featureCoverage <-function(bs,features,genomebuild){


  annotationFeatures<-c()
  for (i in features){
    annotationFeatures<-c(paste0(genomebuild,'_',i),annotationFeatures)
  }


  annots_gr = annotatr::build_annotations(genome = genomebuild, annotations = annotationFeatures)
  GenomeInfoDb::seqlevelsStyle(bs)<-"UCSC"

  nSamples<-dim(bs)[2]

  sumAnnotMatrix<-matrix(nrow=length(features),ncol=nSamples)
  for (i in 1:nSamples){
    bsCell<-bs[,i]

    # CpGs that are observed
    coverageMatrix<-getCoverage(bsCell)
    ind<-DelayedArray::rowSums(coverageMatrix)>0
    # Intersect the regions with the reference annotations

    dm_annotated = annotatr::annotate_regions(
      regions = GenomicRanges::granges(bsCell)[ind,],
      annotations = annots_gr,
      ignore.strand = TRUE,
      quiet = TRUE)
    sumAnnot<-annotatr::summarize_annotations(dm_annotated,quiet=TRUE)
    sumAnnotMatrix[,i]<-sumAnnot$n/sum(ind)

  }
  colnames(sumAnnotMatrix)<-colnames(bs)
  rownames(sumAnnotMatrix)<-features

  return(sumAnnotMatrix)


}
