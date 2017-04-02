## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("scmeth")

## ---- eval=FALSE---------------------------------------------------------
#  load('methylationData.rda')

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::coverage(bs)

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::repmask(bs,"Mus musculus","mm10")

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::chromosomeCoverage(bs,"Mus musculus","mm10")

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::featureCoverage(bs,features=c('genes_exons','genes_introns','genes_intergenic','cpg_islands'),"mm10")

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::cpgDensity(bs,"Mus musculus",windowLength=1000)

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::downsample(bs)

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::mbiasplot(file)

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::methylationDist(file)

## ---- eval=FALSE---------------------------------------------------------
#  scmeth::report(bs, '~/Documents',"Mus musculus","mm10")

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

