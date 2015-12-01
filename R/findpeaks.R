# find peaks from density
peaks <- function(density,span=3){
  
  z <- embed(density, span)
  s <- span%/%2
  v <- max.col(z) == 1 + s
  result <- c(rep(FALSE,s),v)
  result <- result[1:(length(result)-s)]
  result
}



#' find the density peaks of minor allele copy number ratio (minCNR)
#'
#' @param y minCNR for all segments; the first element of the list returned by getASCN()
#' @return a list with two elements (fit_mat , peak_mat)  
#' # 1) fit_mat: Minor_Allelic_CN, density for each segments
#' 2) peak_mat: peak_value, peak_y for each identified peak
#' @export
#'
findPeaks <- function(y){
  
  
  
  # DPdensity parameters
  mcmc <- list(nburn = 1000, nsave = 1000, nskip = 1, ndisplay = 100)
  prior1 <- list(alpha = 1, m1 = rep(0, 1), psiinv1 = diag(0.5, 1), nu1 = 4, tau1 = 1,  tau2 = 100)
  prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.25,1)), nu1=1,nu2=4,tau1=1,tau2=1)
  prior3 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.25,1)), nu1=1,nu2=10,tau1=1,tau2=10)
  
  
  fit <- DPdensity(y = y, prior = prior2, mcmc = mcmc, state = state, status = TRUE)
  peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
  peak_value<-fit$x1[peak_pos]
  peak_y<-fit$dens[peak_pos]
  peak_mat<-as.data.frame(cbind(peak_value, peak_y))
  if ((max(peak_y)-min(peak_y)) > 0.8*(max(fit$dens)-min(fit$dens))){
    fit <- DPdensity(y = y, prior = prior3, mcmc = mcmc, state = state, status = TRUE)
    peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
    peak_value<-fit$x1[peak_pos]
    peak_y<-fit$dens[peak_pos]
    peak_mat<-as.data.frame(cbind(peak_value, peak_y))
    if ((max(peak_y)-min(peak_y)) > 0.8*(max(fit$dens)-min(fit$dens))){
      fit <- DPdensity(y = y, prior = prior1, mcmc = mcmc, state = state, status = TRUE)
      peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
      peak_value<-fit$x1[peak_pos]
      peak_y<-fit$dens[peak_pos]
      peak_mat<-as.data.frame(cbind(peak_value, peak_y))
    }
  }
  
  fit_mat<-as.data.frame(cbind(fit$x1, fit$dens))
  colnames(fit_mat)<-c("Minor_Allelic_CN","density")
  colnames(peak_mat)<-c("peak_value", "peak_y")
  return(list(fit_mat, peak_mat))
}

















# END