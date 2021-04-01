library(ramify)
library(skmeans)
#spherical k-means based on cosine dissimilarity
clusterMeans=function(data,k=2,nruns=100){
  skmeans(data,k,method="pclust",control=list(nruns = nruns))$prototypes
} 
#######################################################################
################## spherical k-principal component clustering##########
#The main routine:
#clusterPC=function(data,k=2,nrep=100,tol=10^(-5),startFromMeans=FALSE)
#
#data is a matrix of observations, where each row is non-negative with unit Euclidean norm.
#k is the number of clusters
#nrep is a number of random restarts
#tol is used to stop searching for a local optimum for each restart
#startFromMeans=TRUE adds one restart using k-means centroids


#######################
#single iteration
#centroids is a k*d matrix with current proposals
clusterPC_iter=function(data, centroids){
  k=length(centroids[,1])
  n=length(data[,1])
  d=length(data[1,])
  M=data%*%t(centroids)
  #find current value
  v=mean(apply(M,1,max))
  gr=argmax(M,rows=T)
  for (i in 1:k){
    seldata=data[gr==i,]
    if (length(seldata)==d) #interpretation problem when just one vector
      seldata=t(seldata)
    Sig=t(seldata)%*%seldata/n
    res=eigen(Sig)
    centroids[i,]=abs(res$vectors[,1]) #use the first eigenvector, the entries fo which are necessarily positive
  }
  list(centroids,v)
}
#pick randomly the initial centers
clusterPCOnce=function(data,k,tol,startFromMeans=FALSE){
  val=0
  n=length(data[,1])
  if (startFromMeans)
    centroids=clusterMeans(data,k)  
  else{
    centroids=data[sample(1:n,k),]
    if (k==1)   #make sure it is a matrix
      centroids=t(as.matrix(centroids))
  }
  niter=0
  repeat{
    niter=niter+1
    res=clusterPC_iter(data, centroids)
    centroids=res[[1]]
    diff=res[[2]]-val
    val=res[[2]]
    if(diff<tol)
      break
  }
  #print(niter)
  list(centroids,val)
}
#iterate nrep times and pick the best
clusterPC=function(data,k=2,tol=10^(-5),nrep=100,startFromMeans=FALSE){
  maxval=0
  for( i in 1:nrep){
    res=clusterPCOnce(data,k,tol,startFromMeans && (i==1))
    if (res[[2]]>maxval){
      maxval=res[[2]]
      centroids=res[[1]]
    }
  }
  centroids
}



##################################################
####################Useful routines#################
#assign groups
getClusterIndex=function(data,centroids){
  M=data%*%t(centroids)
  gr=argmax(M,rows=T)
  gr
}
#compute dissimilarity cost
getCost=function(data,centroids,cosine=TRUE){
  M=data%*%t(centroids)
  products=apply(M,1,max)
  if(cosine)
    return(1-mean(products))
  else
    return(1-mean(products^2))
}

############################################
#================ some test ===============
d=3
n=100
data=matrix(runif(d*n),nrow=n)
nrms=sqrt(rowSums(data^2))
data=data/nrms

clusterMeans(data,k=1)
clusterPC(data,k=1,startFromMeans = TRUE)
clusterMeans(data,k=4)
clusterPC(data,k=4)
