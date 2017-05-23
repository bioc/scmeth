## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("scmeth")

## ---- eval=FALSE---------------------------------------------------------
#  load('methylationData.rda')

## ---- eval=FALSE---------------------------------------------------------
#  library(scmeth)
#  scmeth::report(bs, '~/Documents',Mmusculus,"mm10")

## ------------------------------------------------------------------------
library(scmeth)
load(system.file("extdata",'bsObject.rda',package='scmeth'))

## ------------------------------------------------------------------------
scmeth::coverage(bs)

## ---- warning=FALSE,message=FALSE----------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
scmeth::repMask(bs,Mmusculus,"mm10")

## ---- warning=FALSE------------------------------------------------------
scmeth::chromosomeCoverage(bs)

## ---- warning=FALSE,message=FALSE----------------------------------------
library(annotatr)
scmeth::featureCoverage(bs,features=c('genes_exons','genes_introns',
                                      'cpg_islands'),"mm10")



## ----warning=FALSE,message=FALSE-----------------------------------------
library(BSgenome.Mmusculus.UCSC.mm10)
scmeth::cpgDensity(bs,Mmusculus,windowLength=1000)

## ----warning=FALSE-------------------------------------------------------
scmeth::downsample(bs)

## ----warning=FALSE,message=FALSE,fig.width=6,fig.height=6----------------
library(ggplot2)
scmeth::mbiasplot(system.file("extdata",'16_trimmed_sorted.txt',package='scmeth'))

## ----warning=FALSE,message=FALSE,fig.width=6-----------------------------
#library(ggplot2)
scmeth::methylationDist(bs,all=TRUE)

## ------------------------------------------------------------------------
scmeth::bsConversionPlot(bs)

## ----warning=FALSE-------------------------------------------------------
CpGBedGraphFile_1<-system.file("extdata",
                      'sc-RRBS_zyg_01_chr1_CpG.bedGraph',package='scmeth')
readMetricsFile_1<-system.file("extdata",
                      'sc-RRBS-zygote_01.read_metrics.txt',package='scmeth')
bsConversionFile_1<-system.file("extdata",
                      'sc-RRBS-zygote_01.bsConv.txt',package='scmeth')
CpGBedGraphFile_2<-system.file("extdata",
                      'sc-RRBS_zyg_02_chr1_CpG.bedGraph',package='scmeth')
readMetricsFile_2<-system.file("extdata",
                      'sc-RRBS-zygote_02.read_metrics.txt',package='scmeth')
bsConversionFile_2<-system.file("extdata",
                      'sc-RRBS-zygote_02.bsConv.txt',package='scmeth')
CpGBedGraphFile_3<-system.file("extdata",
                      'sc-RRBS_zyg_03_chr1_CpG.bedGraph',package='scmeth')
readMetricsFile_3<-system.file("extdata",
                      'sc-RRBS-zygote_03.read_metrics.txt',package='scmeth')
bsConversionFile_3<-system.file("extdata",
                      'sc-RRBS-zygote_03.bsConv.txt',package='scmeth')


rda1<-createRDA(CpGBedGraphFile_1,readMetricsFile_1,bsConversionFile_1)
rda2<-createRDA(CpGBedGraphFile_2,readMetricsFile_2,bsConversionFile_2)
rda3<-createRDA(CpGBedGraphFile_3,readMetricsFile_3,bsConversionFile_3)

## ------------------------------------------------------------------------
combineRDA(c(rda1,rda2,rda3))

