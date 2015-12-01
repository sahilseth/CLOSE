
################# Testing Data ##################
#data<-read.table("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/yurolTM.0.45.CNV.output.txt",head=T)
data0 <- read.table("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/yurolTM.0.45.CNV.output.txt",head=T)
data <- data0[data0$logR_mean>-2&data0$logR_mean<2,]#remove outlier
OLL <- data.matrix(data[,c(5,4)])
slength<-(data$end-data$start)/1000000
##
p_grid<-seq(0.2,1,0.1); la_grid<-seq(1,5,1)  #la is ploidy
CanP<-Cpoints(p_grid,la_grid)
#fitPP(OLL,slength,LAF.wt=c(1,5))
fitPP(OLL,slength,LAF.wt=c(1,4))



##get ploidy estimates of 76 samples
sub <- read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sample_path.csv",head=T)[1:76,1:4]
tnames <- as.character(sub$Tumor)
results<-NULL
for (i in 1:length(tnames)) {
  fname<-paste("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/",tnames[i],".0.45.CNV.output.txt",sep="")
  data0<-read.table(fname,head=T)
  data<-data0[data0$logR_mean>-2&data0$logR_mean<2,]
  OLL<-data.matrix(data[,c(5,4)])
  slength<-(data$end-data$start)/1000000
  ##
  p_grid<-seq(0.1,1,0.01); la_grid<-seq(1,5,1)  #la is ploidy
  CanP<-Cpoints(p_grid,la_grid)
  results[i]<-as.character(fitPP(OLL,slength,LAF.wt=c(1,5))[1,1])
}
OUT<-cbind(tnames,results)
save(OUT,file="76ploidy1.Rdata")

sub<-read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sample_path.csv",head=T)[1:76,1:4]
tnames<-as.character(sub$Tumor)
results<-NULL
for (i in 1:length(tnames)) {
  fname<-paste("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/",tnames[i],".0.45.CNV.output.txt",sep="")
  data0<-read.table(fname,head=T)
  data<-data0[data0$logR_mean>-2&data0$logR_mean<2,]
  OLL<-data.matrix(data[,c(5,4)])
  slength<-(data$end-data$start)/1000000
  ##
  p_grid<-seq(0.1,1,0.01); la_grid<-seq(1,5,1)  #la is ploidy
  CanP<-Cpoints(p_grid,la_grid)
  results[i]<-as.character(fitPP(OLL,slength,LAF.wt=c(1,4))[1,1])
}
OUT<-cbind(tnames,results)
save(OUT,file="76ploidy2.Rdata")






###
par(mfrow=c(2,2))
BA<-matrix(nrow=nrow(data),ncol=2)
for (i in 1:nrow(data)) {
  #i=1
  BA[i,]=getBA(x=data[i,5],y=data[i,4],la=4,p=0.92)
}
data1<-cbind(data,BA)

names(data1)[11:12]<-c("min","maj")
plot(data1$min,data1$maj)
for (i in 1:nrow(data0)) {
  #i=1
  BA[i,]=getBA(x=data0[i,5],y=data0[i,4],la=2,p=0.9)
}
data1<-cbind(data0,BA)
names(data1)[11:12]<-c("min","maj")
plot(data1$min,data1$maj)



################# lab ##################

p_grid<-c(0.95, 0.95); la_grid<-c(2,4)  #la is relative average ploidy
foo=Cpoints(p_grid,la_grid)
par(mfrow=c(2,2))
plot(cbind(foo$cLL[,,2][1,], foo$cLL[,,1][1,]),xlim=c(0,0.5),ylim=c(-2,2))
plot(cbind(foo$cLL[,,2][4,], foo$cLL[,,1][4,]),xlim=c(0,0.5),ylim=c(-2,2))

