#' CLOSE-R main functions, estimate copy number based on CRP
#' calls all the sub functions inside of subFunc.R
#'
#' @param Input segments with LAF and LRR. See manual for more detail
#' @param codePath directory for subFunc.R 
#' @param sampleName output files prefix. "closer" is not specified
#' @import DPpackage
#' @export
CRP <- function(Input, sampleName="closer"){
  
  
  library("DPpackage")
  library(grid)
  library(ggplot2)
  
  set.seed(123)
  
  
  ############ segmentation preparation
  segs = segPrep(Input=Input)
  # output: LAF_median and LRR_mean for each segment
  
  
  ############ Dirichlet Process
  
  message("> running dirichlet process ...")
  DPcluster = runDP(x=segs[,1], y=segs[,2], disp.param = 0.45, max.iter = 30, tolerance = .001)
  # cluster segments based on LAF_median and LRR_mean
  segs.cluster = as.data.frame(cbind(DPcluster[[2]], segs))
  colnames(segs.cluster)<-c("cluster_id", "LAF_median","LRR_mean","length","size","chromosome","start","end")
  Ncluster = DPcluster[[3]]
  
  ############ calculate minCNR and majCNR
  
  message("> get ascn ...")
  CNR = getASCN(LAF=segs[,1], R=2^(segs[,2]), scale=1)
  CNR.mat = as.data.frame(cbind(CNR[[1]], CNR[[2]], factor(segs.cluster[,1]), segs.cluster[,5],segs.cluster[,6]))
  CNR.mat[,3] = factor(CNR.mat[,3])
  colnames(CNR.mat) = c("minCNR", "majCNR", "cluster2", "size", "chr")     
  
  ############ find density peaks
  
  message("> find peaks ...")
  CNR.fit = findPeaks(y=CNR[[1]])
  fit.mat = as.data.frame(CNR.fit[[1]])
  peak.mat = as.data.frame(CNR.fit[[2]])
  
  purity = getPurity(peak.mat = peak.mat)
  clonality = getClonality(peak.mat = peak.mat)
  
  
  ############ estimate copy number status

  message("> get centroid ...")
  centroid = getCentroid(segs.cluster = segs.cluster, Ncluster = Ncluster)
  
  message("> estimate cn status ...")
  CNstatus = getCNstatus(centroid, CNR.mat, segs, segs.cluster, Ncluster)
  
  message("> writing out ...")
  write.table(CNstatus, paste(sampleName,  ".CNstatus.txt", sep=""), col.names=colnames(CNstatus), row.names=F, quote=F)
  
  ############ plot CNR clusters
  
  message("> plotting ...")
  plotCNR(sampleName = sampleName, 
          segs.cluster = segs.cluster, 
          CNR.mat = CNR.mat, 
          fit.mat = fit.mat, 
          peak.mat = peak.mat, 
          purity = purity, 
          clonality = clonality)
  
  ret = list(CNstatus = CNstatus, 
             segs.cluster = segs.cluster, 
             CNR.mat = CNR.mat, 
             fit.mat = fit.mat, 
             peak.mat = peak.mat, 
             purity = purity, 
             clonality = clonality)
  
  invisible(ret)
}


