## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("scmeth")

## ---- eval=FALSE---------------------------------------------------------
#  directory<-system.file("extdata/bismark_data",package='scmeth')
#  bsObject<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)

## ---- eval=FALSE---------------------------------------------------------
#  library(scmeth)
#  scmeth::report(bsObject, '~/Documents',Hsapiens,"hg38")

## ----  warning=FALSE,message=FALSE,comment=FALSE-------------------------
library(scmeth)
directory<-system.file("extdata/bismark_data",package='scmeth')
bsObject<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)

## ------------------------------------------------------------------------
scmeth::coverage(bsObject)

## ----fig.width=6,fig.height=6--------------------------------------------
scmeth::readmetrics(bsObject)

## ---- warning=FALSE,message=FALSE,eval=FALSE-----------------------------
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  load(system.file("extdata",'bsObject.rda',package='scmeth'))
#  scmeth::repMask(bs,Mmusculus,"mm10")

## ---- warning=FALSE------------------------------------------------------
scmeth::chromosomeCoverage(bsObject)

## ---- warning=FALSE,message=FALSE----------------------------------------
library(annotatr)
featureList<-c('genes_exons','genes_introns')
DT::datatable(scmeth::featureCoverage(bsObject,features=featureList,"hg38"))



## ----warning=FALSE,message=FALSE-----------------------------------------
library(BSgenome.Hsapiens.NCBI.GRCh38)
DT::datatable(scmeth::cpgDensity(bsObject,Hsapiens,windowLength=1000))

## ----warning=FALSE-------------------------------------------------------
DT::datatable(scmeth::downsample(bsObject))

## ----warning=FALSE,message=FALSE,fig.width=6,fig.height=6----------------
methylationBiasFile<-'2017-04-21_HG23KBCXY_2_AGGCAGAA_TATCTC_pe.M-bias.txt'
scmeth::mbiasplot(mbiasFiles=system.file("extdata",methylationBiasFile,
                                        package='scmeth'))

## ----warning=FALSE,message=FALSE,fig.width=6,fig.height=6----------------
scmeth::methylationDist(bsObject)

## ----warning=FALSE,message=FALSE,fig.width=6,fig.height=6----------------
scmeth::bsConversionPlot(bsObject)

