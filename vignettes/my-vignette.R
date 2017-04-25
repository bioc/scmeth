## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("scmeth")

## ---- eval=FALSE---------------------------------------------------------
#  load('methylationData.rda')

## ---- eval=FALSE---------------------------------------------------------
#  library(scmeth)
#  scmeth::report(bs, '~/Documents',"Mus musculus","mm10")

## ------------------------------------------------------------------------
library(scmeth)
load(system.file("extdata",'bsObject.rda',package='scmeth'))

## ------------------------------------------------------------------------
scmeth::coverage(bs)

## ---- warning=FALSE,message=FALSE----------------------------------------
#library(GenomicRanges)
#scmeth::repMask(bs,"Mus musculus","mm10")

## ---- warning=FALSE------------------------------------------------------
scmeth::chromosomeCoverage(bs)

## ---- warning=FALSE,message=FALSE----------------------------------------
scmeth::featureCoverage(bs,features=c('genes_exons','genes_introns','genes_intergenic','cpg_islands'),"mm10")

## ----warning=FALSE,message=FALSE-----------------------------------------
scmeth::cpgDensity(bs,"Mus musculus",windowLength=1000)

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
rda1<-createRDA(system.file("extdata",'sc-RRBS_zyg_01_chr1_CpG.bedGraph',package='scmeth'))
rda2<-createRDA(system.file("extdata",'sc-RRBS_zyg_02_chr1_CpG.bedGraph',package='scmeth'))
rda3<-createRDA(system.file("extdata",'sc-RRBS_zyg_03_chr1_CpG.bedGraph',package='scmeth'))

## ------------------------------------------------------------------------
combineRDA(c(rda1,rda2,rda3))

