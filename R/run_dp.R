#' cluster segments based on LAF_median and LRR_mean using dirichlet process   
#'
#' @param x x variable (LAF for each segment)
#' @param y y variable (LRR for each segment)
#' @param disp.param 
#' @param max.iter 
#' @param tolerance 
#'
#' @source http://statistical-research.com/dirichlet-process-infinite-mixture-models-and-clustering/
#' @export
runDP <- function(x, y, disp.param = 0.45, max.iter = 100, tolerance = .001){
  
  
  myData <- cbind( x, y )
  n <- dim(myData)[1]
  var.range<-c(0,0)
  var.range[1]<-max(myData[,1])-min(myData[,1])
  var.range[2]<-max(myData[,2])-min(myData[,2])
  small.scale<-which(var.range==min(var.range))
  
  myData[,small.scale]<-myData[,small.scale]*4
  pick.clusters <- rep(1, n)
  k <- 1
  
  mu <- matrix( apply(myData,2,mean), nrow=1, ncol=ncol(myData) )
  
  is.converged <- FALSE
  iteration <- 0
  
  ss.old <- Inf
  ss.curr <- Inf
  itr.cluster<-matrix(NA, nrow=n, ncol=0)
  while ( !is.converged & iteration < max.iter ) { 
    # Iterate until convergence
    iteration <- iteration + 1
    for( i in 1:n ) { 
      # Iterate over each observation and measure the distance each observation' from its mean center's squared distance from its mean
      distances <- rep(NA, k)
      
      for( j in 1:k ){
        distances[j] <- sum( (myData[i, ] - mu[j, ])^2 ) # Distance formula.
      }
      
      if( min(distances) > disp.param ) { # If the dispersion parameter is still less than the minimum distance then create a new cluster
        k <- k + 1
        pick.clusters[i] <- k
        mu <- rbind(mu, myData[i, ])
      } else {
        pick.clusters[i] <- which(distances == min(distances))
      } 
    }
    
    ####### Calculate new cluster means
    
    for( j in 1:k ) {
      if( length(pick.clusters == j) > 0 ) {
        mu[j, ] <- apply(subset(myData,pick.clusters == j), 2, mean)
      }
    }
    
    ####### Test for convergence
    
    ss.cur<-0
    for (i in 1:n){
      ss.cur<-ss.curr+sum((myData[i,]-mu[pick.clusters[i],])^2)
    }
    ss.diff<-ss.old-ss.curr
    ss.old<-ss.curr
    if (!is.nan(ss.diff) & ss.diff < tolerance){
      is.converged<-TRUE
    }
    
    
    
    itr.cluster<-cbind(itr.cluster,pick.clusters)
  }
  
  centers <- data.frame(mu)
  ret.val <- list("centers" = centers, "cluster" = factor(pick.clusters),
                  "k" = k, "iterations" = iteration, "itr_cluster" = itr.cluster)
  
  return(ret.val)
}
