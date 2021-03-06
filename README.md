################## spherical k-principal-component clustering ##########

############################### R language code ########################

described in
V. Fomichov and J. Ivanovs. Spherical clustering in detection of groups of concomitant extremes. Biometrika ???.
https://arxiv.org/abs/2103.11475

The main routine:
clusterPC(data,k=2,nrep=100,tol=10^(-5),startFromMeans=FALSE)

-data is a matrix of observations, where each row is non-negative with unit Euclidean norm.

-k is the number of clusters

-nrep is a number of random restarts

-tol is used to stop searching for a local optimum for each restart

-startFromMeans=TRUE adds one restart using k-means centroids
