# all sub-functions in CLOSER

# convert BAF to LAF values
BAFtoLAF <- function(BAF.vec){
  
  # BAF.vec: a vector of BAF values (numeric)
  
  LAF.vec <- BAF.vec
  BAF.folding.index<-(1:length(BAF.vec))[BAF.vec>0.5]
  if (length(BAF.folding.index)>1){
    LAF.vec[BAF.folding.index]<-1-LAF.vec[BAF.folding.index]
  }
  return(LAF.vec)
}





#' calualte minCNR (minor copy number ratio), majCNR (major copy number ratio)
#'
#' @param LAF vector of LAF values
#' @param R ratio of read depth (tumor/normal) for each segments: usualy 2^LRR
#' @param scale scale factor for library size
#'
#' @return list with 2 elements (minCNR and majCNR for each segment)
#' @export
#'
getASCN <- function(LAF, R, scale=1){
  
  minCNR <- 2*LAF*R/scale
  majCNR <- 2*R*(1-LAF)/scale
  return(list(minCNR, majCNR))
}


#' find centroid of the CNR clsuters
#'
#' @param segs.cluster 
#' @param Ncluster number of clusters; third element of the list returned from runDP()
#'
#' @return centroid of each cluster
#' @export
#'
getCentroid<-function(segs.cluster, Ncluster){
  
  
  # 1. get the mean values of LAF and LRR for each cluster as the centroid
  LAF_centroid_vec<-c()
  logR_centroid_vec<-c()
  for (i in 1:Ncluster){
    points<-segs.cluster[segs.cluster[,1]==i,]
    LAF_centroid_vec<-c(LAF_centroid_vec,mean(points[,2]))
    logR_centroid_vec<-c(logR_centroid_vec, mean(points[,3]))
  }
  
  # 2. calculate minCNR and majCNR for centroid
  centroid_min_majCNR<-getASCN(LAF=LAF_centroid_vec, R=2^(logR_centroid_vec), scale=1)
  centroid<-as.data.frame(cbind(centroid_min_majCNR[[1]], centroid_min_majCNR[[2]], rep(0, length(Ncluster) ), rep(16, length(Ncluster)) ))
  colnames(centroid)<-c("minCNR","majCNR", "chr", "shape")
  return(centroid)
}


# get clonality from number of peaks
getClonality <- function(peak.mat){
  # paramters:
  # 1) peak.mat: peak_value, peak_y for each identified peak; the second element of the list returned from findPeaks()
  
  # return: a single value of clonality, numeric
  
  clonality<-dim(peak.mat)[1]-2
  if (clonality<=0){
    clonality<-0
  }
  return(clonality)
}


#' estimate the copy number status: normal, gainm loss, LOH, ambigous
#'
#' @param centroid 
#' @param CNR.mat 
#' @param segs 
#' @param segs.cluster 
#' @param Ncluster 
#'
#' @export
#'
getCNstatus <- function(centroid, CNR.mat, segs, segs.cluster, Ncluster){
  
  message("--> get the status for each segments based on cluster information ...")
  message("----> estimate the copy number status of each cluster centroid ...")
  normal_index<-(1:dim(centroid)[1])[centroid[,1]>=0.75 & centroid[,1]<=1.25 & centroid[,2]>=0.5 & centroid[,2]<=1.5]
  loss_loh_index<-(1:dim(centroid)[1])[centroid[,1]<=0.3 ]
  neutral_loh_index<-(1:dim(centroid)[1])[centroid[,1]<=0.3 & centroid[,2]>=1.8 & centroid[,2]<=2.2]
  gain_index<-(1:dim(centroid)[1])[centroid[,1]>=0.75 & centroid[,2]>=1.5]
  ambigu_index<-(1:dim(centroid)[1])[-c(normal_index, loss_loh_index, gain_index)]
  
  cluster_status<-rep("Normal", rep(dim(centroid)[1]))
  cluster_status[loss_loh_index]<-"Loss"
  cluster_status[neutral_loh_index]<-"CN Neutral LOH"
  cluster_status[gain_index]<-"Gain"
  cluster_status[ambigu_index]<-"Ambiguous"
  
  message("----> assign the CN status of cluster centroid to all segments in that cluster ...")
  seg_status_cluster<-rep("N", dim(segs)[1])
  for (i in 1:Ncluster){
    cluster_seg<-(1:dim(segs.cluster)[1])[segs.cluster[,1]==i]
    seg_status_cluster[cluster_seg]<-rep(cluster_status[i], length(cluster_seg))
  }
  
  
  message("--> assign status to each segments using segments infomation, not cluster ...")
  seg_normal_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]>=0.75 & CNR.mat[,1]<=1.25 & CNR.mat[,2]>=0.5 & CNR.mat[,2]<=1.5]
  seg_loss_loh_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]<=0.3]
  seg_neutral_loh_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]<=0.3 & CNR.mat[,2]>=1.8 & CNR.mat[,2]<=2.2]
  seg_gain_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]>=0.75 & CNR.mat[,2]>=1.5]
  seg_ambigu_index<-(1:dim(CNR.mat)[1])[-c(seg_normal_index, seg_loss_loh_index, seg_gain_index)]
  
  seg_status<-rep("Normal", rep(dim(CNR.mat)[1]))
  seg_status[seg_loss_loh_index]<-rep("Loss", length(seg_loss_loh_index))
  seg_status[seg_neutral_loh_index]<-rep("CN Neutral LOH", length(seg_neutral_loh_index))   
  seg_status[seg_gain_index]<-rep("Gain", length(seg_gain_index))
  seg_status[seg_ambigu_index]<-rep("Ambiguous", length(seg_ambigu_index))
  
  
  f<-cbind(paste("chr",segs[,5], sep=""), segs[,6:7],segs[,1:2], CNR.mat[,1:2], segs.cluster[,1], seg_status_cluster, seg_status)
  f2<-f[order(segs[,5], f[,2]),]
  colnames(f2)<-c("chr","start","end","LAF_median", "logR_mean", "minCNR", "majCNR", "cluster", "status_cluster", "status_seg")
  return(f2)
}


# get purity from peak values of minor allele CN ratio
getPurity<-function(peak.mat){
  
  peak.value <- peak.mat[,1]
  return(peak.value[length(peak.value)]-peak.value[1])
}








# prepare segmentations for CLOSE-R
segPrep<-function(Input){
  
  # paramters:
  # 1) Input:the Input file of CLOSER 
  
  # return: LAF_median and LRR_mean for each segment
  
  size_vec = region_length_vec = c()
  for (i in 1:dim(Input)[1]){
    region_length<-as.numeric(Input[i,3])-as.numeric(Input[i,2])+1
    region_length_vec<-c(region_length_vec, region_length)
    size<-as.numeric(region_length)/1000000
    if (size<1) {size<-1}
    if (size>8) {size<-8}
    size_vec<-c(size_vec,size)
  }
  
  mat<-cbind(Input[,4:5], region_length_vec,  size_vec, Input[,1:3] )
  return(mat)
  
}









