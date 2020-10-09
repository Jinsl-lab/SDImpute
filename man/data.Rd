\name{data}
\alias{data}
\docType{data}
\title{Camp dataset}
\description{
  This data set is subset from Camp et al. (2017). Original data is

  RNA-seq data on 425 cells.Column names indicate five different cell types,
  row names indicate gene.

}
\examples{
library("SDImpute")
data(data)
imputed_data<-SDImpute(data,do.nor=TRUE,auto.k=FALSE,k=5,M=15,T=0.5)
}

