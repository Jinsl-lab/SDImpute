#' SDImpute: A statistical block imputation method based on cell-level and gene-level information for dropouts
#'           in single-cell RNA-seq data
#' @param data  A gene expression matrix,the rows correspond to genes and the columns correspond to cells.
#' @param do.nor  Logical. If TRUE, the data is Normalized by a global-scaling normalization method.
#' @param do.log  Logical. If TRUE, the data is logarithmic transformation (log2).
#' @param auto.k  Logical. If TRUE, k is estimated  by either the Calinski Harabasz index  or average silhouette width ;
#'            If FALSE, the parameter k need to be set manually.
#' @param criterion  One of "asw" or "ch". Determines whether average silhouette width or Calinski-Harabasz is applied.
#' @param krange  Integer vector. Numbers of clusters which are to be compared by the average
#'              silhouette width criterion. Note: average silhouette width and Calinski-Harabasz
#'              can't estimate number of clusters nc=1.
#' @param k   Integer. The number of cell clusters. This parameter can be determined based on prior knowledge or clustering result of raw data.
#' @param M   Integer. The number of nearest neighbors.When the number of nearest neighbors for each cell is small, the parameter M should not
#'            be too large to guarantee that it makes sense. In general, this parameter is set to an integer between 10 and 30.
#' @param T   Numeric between 0 and 1. The dropout probability candidate threshold which controls the degree of imputation to the gene expression
#'            matrix. The recommended value of parameter T is 0.5.
#' @export
#' @return  An imputation matrix
#' @import fpc
#' @import stats
#' @examples
#' library("SDImpute")
#' data(data)
#' imputed_data<-SDImpute(data,do.nor=TRUE,auto.k=FALSE,k=5,M=15,T=0.5)

SDImpute<-function(data,do.nor=TRUE,do.log=TRUE,auto.k=TRUE,criterion="asw",krange=c(5:15),k=5,M=15,T=0.5){

  print("Data preprocessing...")
  requireNamespace("stats")
  nor<-function(x){
    if(x==TRUE){
      ltpm <- t(t(data)/colSums(data))*1000000+1
    }
    else{ltpm=data
    }
  }
  ltpm<-nor(do.nor)
  logt<-function(x){
    if(x==TRUE){
      ltpm <- log2(ltpm+1)
    }
    else{ltpm=ltpm
    }
  }

  ltpm<-logt(do.log)


  print("Identification of dropouts...")

  # performing PCA and k-means.
  pca <- prcomp(t(ltpm))
  requireNamespace("fpc")
  #library()
  pkc<-function(x){
    if(x==TRUE){
      pkc0<-kmeansruns(pca$x,krange=krange,criterion=criterion,
                       iter.max=10,runs=10)
    }

    else{
      pkc0<-kmeans(pca$x,centers=k)
    }
  }
  pkc0<-pkc(auto.k)

  table(pkc0$cluster)

  result<-cbind(pca$x[,c(1,2)],pkc0$cluster)
  label<-t(result[,3])
  label<-c(label)
  ltpm<-as.matrix(data)
  labelltpm<-as.matrix(rbind(label,ltpm))


  #Calculating dropout rate.

  a<-matrix(0,ncol = length(pkc0$size) ,nrow = dim(ltpm)[1])
  p<-function(x){
    which(x==0)
    nzero<-length(which(x==0))
    return(nzero)
  }

  m<-matrix(0,ncol = length(pkc0$size) ,nrow = dim(ltpm)[1])
  for (c in 1: length(pkc0$size)) {
    cellid<-which(labelltpm[1,]==c)
    for (i in 1:dim(data)[1]) {
      m[i,c]=p(ltpm[i,cellid])
    }
  }
  N<-table(labelltpm[1,])
  N<-c(0,N)
  for (i in 1: length(pkc0$size)) {
    a[,i]=m[,i]/N[i+1]

  }

  aedata<-matrix(0,ncol = length(pkc0$size) ,nrow = dim(ltpm)[1])
  for (c in 1: length(pkc0$size)) {
    cellid<-which(labelltpm[1,]==c)
    aedata[,c]=rowMeans(ltpm[,cellid])
  }

  warnings('off')
  options(warn =-1)
  data0<-cbind(a,aedata)
  data0<-as.data.frame(data0)
  e<-matrix(0,ncol =2 ,nrow = length(pkc0$size) )
  for (c in 1: length(pkc0$size)) {
    data0$y[data0[,c]>T]<-1
    data0$y[data0[,c]<=T]<-0
    anes=glm(y~aedata[,c],family=binomial(link='logit'),data=data0,control=list(maxit=100))
    e[c,]<-anes$coefficients

  }
  options(warn=-1)

  pm<-matrix(0,ncol =dim(ltpm)[2],nrow =dim(ltpm)[1])
  clusterltpm<-labelltpm[2:dim(labelltpm)[1],]
  for (c in 1:length(pkc0$size)) {

    for (j in 1:dim(ltpm)[2]){

      for (i in 1:dim(ltpm)[1]){

        if(labelltpm[1,j]==c){

          pm[i,j]=exp(e[c,1]+e[c,2]*clusterltpm[i,j])/(1+exp(e[c,1]+e[c,2]*clusterltpm[i,j]))
        }
      }
    }
  }
  pm[is.na(pm)]<-0

  #Calculating dropout probability.
  cvltpm<-labelltpm[2:dim(labelltpm)[1],]
  cv <-matrix(0,ncol = length(pkc0$size) ,nrow = dim(cvltpm)[1])
  for (c in 1:length(pkc0$size)) {
    cellid<-which(labelltpm[1,]==c)
    for (i in 1:dim(cvltpm)[1]) {
      cv[i,c]=sqrt(var(cvltpm[i,cellid]))/((mean(cvltpm[i,cellid]))+0.001)

    }
  }

  cvnormal <-matrix(0,ncol = length(pkc0$size) ,nrow = dim(cvltpm)[1])
  for (c in 1:length(pkc0$size)) {
    for (i in 1:dim(cvltpm)[1]) {
      cvnormal[i,c]=atan(cv[i,c]*sqrt(var(cv[,c])))*2/pi

    }
  }



  dropoutp<-matrix(0,ncol =dim(ltpm)[2],nrow =dim(ltpm)[1])
  for (c in 1:length(pkc0$size)) {
    for (i in 1:dim(ltpm)[1]) {
      for (j in 1:dim(ltpm)[2]) {
        if(labelltpm[1,j]==c){
          dropoutp[i,j]=cvnormal[i,c]*pm[i,j]
        }
      }
    }
  }


  #Calculating Gaussian kernel coefficient matrix.
  print("Block imputation...")
  g<-pca$x[,1:2]
  dw<-as.matrix(dist(g))
  #dw<-as.matrix(dist(t(x)))#
  tknn<-c(rep(1,dim(ltpm)[2]))
  for (i in 1:dim(ltpm)[2]){
    rank=rank(dw[i,])
    J=(which(rank<="M"))
    tknn[i]=mean(dw[i,J])
  }


  gaussk=matrix(0,ncol = dim(ltpm)[2],nrow = dim(ltpm)[2])
  for (i in 1:dim(ltpm)[2]) {

    for (j in 1:dim(ltpm)[2]) {
      gaussk[i,j]=exp(-dw[i,j]/tknn[i])

    }
  }
  gaussk<-as.matrix(gaussk)

  #Performing imputation for dropouts.


  J= dim(ltpm)[2]
  I = dim(ltpm)[1]
  idwm<-matrix(0,ncol = J,nrow = I)
  imput<-function(x){
    for (c in 1:length(pkc0$size)) {
      clustercellid<-which(labelltpm[1,]==c)
      for (j in 1:J){
        if(labelltpm[1,j]==c){
          for (i in 1:I){
            if(dropoutp[i,j]>T){
              cellid<-which(dropoutp[i, clustercellid]<=T)
              a1<-c(ltpm[i,cellid])
              wt1<-c(gaussk[j,cellid])
              idwm[i,j] =weighted.mean(a1,wt1)
              if(idwm[i,j]== "NaN")
              {
                idwm[i,j]=0
              }
            }
            else{
              idwm[i,j]=ltpm[i,j]}
          }
        }
      }

    }
    return(idwm)
  }

  SDImputedata<-imput(ltpm)

  colnames(SDImputedata)<-colnames(data)
  rownames(SDImputedata)<-rownames(data)


  return(SDImputedata)


}
print("Done")
