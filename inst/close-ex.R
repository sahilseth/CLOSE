
################# Testing Data ##################

#get ASCN estimates of 76 samples
sub <- read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sample_path.csv",head=T,as.is=T)[1:76,1:4]
tnames<-as.character(sub$Tumor)
load("~/Dropbox/exomeCNV_revision/CLOSE-results/76ploidy1_5p.Rdata"); 
PP<-OUT #load ploidy estimates
PP=cbind(PP,do.call(rbind,strsplit(PP[,2],"/")))

TOUT<-NULL
PP=PP[c(-26,-60),]; tnames=tnames[c(-26,-60)] #exclude yuquestT and then yugoeTM2

sublist<-which(PP[,4]=="2");
TCOR<-matrix(ncol=4,nrow=length(sublist))
TRate<-NULL
for (j in 1:length(sublist)) {
  #j=13
  i=sublist[j]
  # i=64 #yurif
  #i=13
  fname<-paste("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/",tnames[i],".0.45.CNV.output.txt",sep="")
  data<-read.table(fname,head=T,as.is=T)
  data<-data[order(data[,1],data[,2]),]
  data0<-data[,1:5]
  pp<-PP[PP[,1]==tnames[i],]
  pp1=as.numeric(strsplit(pp[2],"/")[[1]][1]); pp2=as.numeric(strsplit(pp[2],"/")[[1]][2])
  RES<-fitCLOSE1(data0,la=pp2,p=pp1)
  #read falcon results
  #tname_short<-strsplit(tnames[i],"T")[[1]][1]
  tname_short<-tnames[i]
  gname<-paste("~/Dropbox/xiaoqing_share/falcon_yusample/",tname_short,"/CLOSER_falcon.txt",sep="")
  gdata<-read.table(gname,head=T,as.is=T)[,c(1,2,11,12)]
  falcon_tCN<-gdata$falcon_minCN+gdata$falcon_majCN
  gdata<-cbind(gdata,falcon_tCN)
  #merge and calculate correlation
  OUT<-merge(RES,gdata,by=c("chr","start")); OUT<-OUT[order(OUT[,1],OUT[,2]),]
  OUT<-merge(OUT,data[,c(1,2,9)],by=c("chr","start"))
  close_status<-rep("NA",nrow(OUT))
  for (k in 1:nrow(OUT)){
    if (OUT$tCN[k]>2) {close_status[k]="G"}
    if (OUT$tCN[k]<2) {close_status[k]="L"}
    if (OUT$tCN[k]==2) {close_status[k]="N"}
    if (OUT$tCN[k]==2 & OUT$minCN[k]==0) {close_status[k]="CNL"}
  }
  OUT=data.frame(OUT,close_status)
  OUT1<-OUT[OUT$status_cluster!="A",]
  TRate[j]<-sum(OUT1$status_cluster==OUT1$close_status)/nrow(OUT1)
  OUT<-OUT[OUT$tCN<6,]    #remove outliers OUT$tCN>6
  #TRate[j]<-sum(OUT$status_cluster==OUT$close_status)/nrow(OUT)
  TCOR[j,]=c(cor(OUT$falcon_minCN,OUT$min,use="complete.obs"),cor(OUT$falcon_tCN,OUT$tCN,use="complete.obs"),
             cor(round(OUT$falcon_minCN),OUT$min,use="complete.obs"), cor(round(OUT$falcon_tCN),OUT$tCN,use="complete.obs")
  )
  TOUT[[j]]<-OUT; OUT<-NULL
}


################# Table 1 (diploid ~30) ##################
tab1<-data.frame(PP[sublist,],TRate, TCOR);tab1

tab1$V3<-as.numeric(as.character(tab1$V3)); tab1<-tab1[order(tab1$V3),]
par(mfcol=c(2,2))
boxplot(tab1[tab1$V3<0.6,"X1"],tab1[tab1$V3>0.6,"X1"],ylim=c(0,1))
boxplot(tab1[tab1$V3<0.6,"X2"],tab1[tab1$V3>0.6,"X2"],ylim=c(0,1))
boxplot(tab1[tab1$V3<0.6,"X3"],tab1[tab1$V3>0.6,"X3"],ylim=c(0,1))
boxplot(tab1[tab1$V3<0.6,"X4"],tab1[tab1$V3>0.6,"X4"],ylim=c(0,1))
boxplot(tab1[tab1$V3<0.6,"TRate"],tab1[tab1$V3>0.6,"TRate"],ylim=c(0,1))

##ggplot
#concorance Rate (CLOSE model based vs model free)
library(ggplot2)
pdf("fig1.pdf",width=5, height=5)
dat=rbind(data.frame(purity=">0.4",Rate=tab1[tab1$V3>0.4,"TRate"]),data.frame(purity=">0.7",Rate=tab1[tab1$V3>0.7,"TRate"]),
          data.frame(purity=">0.8",Rate=tab1[tab1$V3>0.8,"TRate"]),data.frame(purity=">0.9",Rate=tab1[tab1$V3>0.9,"TRate"])
)
names(dat)=c("purity","Con_rate")
ggplot(dat, aes(x=purity, y=Con_rate),fill=purity) + geom_boxplot() +ylim(0,1) +ylab("CN concordance rate")+xlab("")+
  scale_y_continuous(breaks = c(0,0.5,0.6,0.7,0.8,0.9,1))
dev.off()

#falcon (CLOSE_falcon)
pdf("fig2.pdf",width=6, height=5)
dat1=rbind(data.frame("minor CN","no","purity<0.6",tab1[tab1$V3<0.6,"X1"]),data.frame("minor CN","no","purity>0.6",tab1[tab1$V3>0.6,"X1"]))
names(dat1)=c("mm","rounding", "purity","cor")
dat2=rbind(data.frame("minor CN","yes","purity<0.6",tab1[tab1$V3<0.6,"X3"]),data.frame("minor CN","yes","purity>0.6",tab1[tab1$V3>0.6,"X3"]))
names(dat2)=c("mm","rounding", "purity","cor")
dat3=rbind(data.frame("total CN","no","purity<0.6",tab1[tab1$V3<0.6,"X2"]),data.frame("total CN","no","purity>0.6",tab1[tab1$V3>0.6,"X2"]))
names(dat3)=c("mm","rounding", "purity","cor")
dat4=rbind(data.frame("total CN","yes","purity<0.6",tab1[tab1$V3<0.6,"X4"]),data.frame("total CN","yes","purity>0.6",tab1[tab1$V3>0.6,"X4"]))
names(dat4)=c("mm","rounding", "purity","cor")
dat=rbind(dat1,dat2,dat3,dat4)
#ggplot(dat, aes(x=purity, y=cor)) + geom_boxplot() +ylim(0,1)+facet_wrap(~ rounding)
ggplot(dat, aes(x=purity, y=cor , fill=rounding)) + geom_boxplot() + scale_y_continuous(breaks = c(0,0.5,0.6,0.7,0.8,0.9,1))+ylab("correlation")+ xlab("") +facet_wrap(~ mm)
dev.off()

#pp correlation
seq76=read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sequenza_76.csv",as.is=T)
load("~/Google Drive/exomeCNV_revision/CLOSE-results/76ploidy1_5p.Rdata"); PP<-OUT #load ploidy estimates
PP=cbind(PP,do.call(rbind,strsplit(PP[,2],"/")))

seq76[c(37,64,66),]
PP[c(37,64,66),]


#sample diagnosis
OUT<-TOUT[[24]]
par(mfcol=c(2,2))
boxplot(OUT$falcon_minCN~OUT$min,xlab="minor allele",main=cor(OUT$falcon_minCN,OUT$min))
boxplot(OUT$falcon_tCN~OUT$tCN,xlab="total copy number",main=cor(OUT$falcon_tCN,OUT$tCN))


#S table 1
tab1<-data.frame(PP[sublist,],TRate, TCOR);tab1
seq76=read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sequenza_76.csv",as.is=T)
foo=merge(tab1,seq76,by.x="tnames",by.y="Tumor")
#dat=data.frame(close=as.numeric(as.character(foo$V3)),sequenza=foo$cellularity)
#require(ggplot2)
#ggplot(dat,aes(x=close,y=sequenza))+geom_point()+geom_line()
TABLE<-foo[,c(1,3,11,5)]