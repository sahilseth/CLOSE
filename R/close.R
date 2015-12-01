#Local CLOSE-R estimation
#xuefeng.wang@stonybrook.edu; xuefeng.wang@yale.edu


get.MCpoints <- function(p_grid=0.92,la=4){
  PP.matrix<-data.matrix(expand.grid(p_grid,la)); colnames(PP.matrix)<-c("purity","ploidy")
  BC.matrix<-BC.expand(x=6); colnames(BC.matrix)<-c("B","total")
  LRR_LAF.array<-array(dim=c(nrow(PP.matrix),nrow(BC.matrix),2))
  MM.array<-array(dim=c(nrow(PP.matrix),nrow(BC.matrix),2))
  for (i in 1:nrow(PP.matrix)) {
    for (j in 1:nrow(BC.matrix)) {
      #LRR_LAF.array[i,j,]<-ifelse(BC.matrix[j,2]>PP.matrix[i,2], NA, getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1]))  #if model loss only
      LRR_LAF.array[i,j,]<-getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1],bias=0)
      MM.array[i,j,]<-getBA(x=LRR_LAF.array[i,j,1],y=LRR_LAF.array[i,j,2],la=PP.matrix[i,2],p=1)
      
    }
  }
  name1=paste(PP.matrix[,1],PP.matrix[,2],sep="/")
  name2=paste(BC.matrix[,1],BC.matrix[,2],sep="/")
  #
  out=data.frame(Min=MM.array[,,1],Maj=MM.array[,,2],BC.matrix[,1],BC.matrix[,2])
  names(out)[3:4]<-c("B","C")
  return(out)
}

#get.MCpoints(p_grid=0.92,la=4)
#get.MCpoints(p_grid=0.95,la=2)

##main function to estimate local copy number given ploidy and purity (based on segmental ASCN coordinates)
##assume one clone
fitCLOSE1 <- function(data0,la,p){
  #data0: c1: chr; c2-c3: strat and end postion, c4, LAF, c5,LRR
  #slength: length of each segments
  #la=2;p=1
  MM.coor=get.MCpoints(p_grid=p,la=la)
  o=matrix(nrow=nrow(data0),ncol=4)
  for (i in 1:nrow(data0)) {
    #i=2
    OLL.BA=getBA(x=data0[i,5],y=data0[i,4],la=la,p=1)
    dist<-apply(MM.coor[,1:2],1,function(x)dist(rbind(x,OLL.BA)))
    o[i,]=c(OLL.BA, t(MM.coor[which.min(dist),3:4]))
  }
  close<-cbind(data0,o)
  names(close)[8:9]<-c("minCN","tCN")
  return(close)
}

#data<-read.table("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/yurifT.0.45.CNV.output.txt",head=T)
#data<-data[order(data[,1],data[,2]),]
#data0<-data[,1:5]
#foo1=fitCLOSE1(data0,la=2,p=1)


