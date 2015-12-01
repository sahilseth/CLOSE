# plot clusters based on minCNR and majCNR 
plotCNR <- function(sampleName, segs.cluster, CNR.mat, fit.mat, peak.mat, purity, clonality){
  
  # plot clusters based on minCNR and majCNR 
  pdf(paste(sampleName, ".plotCNR.pdf", sep=""), width=14, height=8)
  
  print(ggplot(segs.cluster, aes_string(x = "LAF_median", y = "LRR_mean", color = "cluster_id", label="chromosome")) + 
          geom_point(aes_string(size="size"), alpha=0.3) +
          scale_size_continuous(range = c(5,15))+ 
          geom_text(aes_string(color="cluster_id")) + 
          xlim(0,0.5) + ylim(-2,2) + theme_bw() + 
          theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", 
                plot.title=element_text(face="bold",size=20)))
  
  print(ggplot(CNR.mat, aes_string(x = "minCNR", y = "majCNR", color = "cluster2", label="chr")) + 
          geom_point(aes_string(size="size"), alpha=0.3)+
          scale_size_continuous(range = c(5,15)) + 
          geom_text(aes_string(color="cluster2")) + 
          xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1]))+ 
          ylim(0,min(6, max(CNR.mat[,2])))+theme_bw() +
          theme(panel.grid.major = element_blank(), 
                panel.border = element_rect(colour = "black"),legend.position="none", 
                plot.title=element_text(face="bold", size=20))+
          geom_hline(aes_string(yintercept=1), linetype="dashed")+
          geom_vline(aes_string(xintercept=1), linetype="dashed"))
  
  print(ggplot(fit.mat, aes_string(x="Minor_Allelic_CN", y="density"))+
          geom_line( color="darkred")+xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1])) + 
          ylim(0,max(fit.mat[,2])+2) + 
          geom_text(data=peak.mat, x=peak.mat[,1], y=peak.mat[,2]+0.2, label=paste("peak ", signif(peak.mat[,1],3), sep="")) +
          geom_text(data=NULL, x=max(fit.mat[,1])*0.6, y=(max(fit.mat[,2])+2)*0.9, 
                    label=paste("Purity parameter = ",signif(purity,3), "\n", "Clonality parameter = ", clonality, sep="" ))+
          theme_bw() + ggtitle("Dirchlet Process Fit of Minor Allelic CN")+
          theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),
                legend.position="none", plot.title=element_text(face="bold", size=14)))
  
  #print(ggplot(CNR.mat, aes_string(x = "minCNR", y = "majCNR", color = "cluster3", label="chr")) +geom_point(aes_string(size="size"), alpha=0.3)+scale_size_continuous(range = c(5,15))+geom_text(aes_string(color="cluster3"))+xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1]))+ylim(0,min(max(CNR.mat[,2])))+theme_bw() +theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", plot.title=element_text(face="bold", size=20))+geom_hline(aes_string(yintercept=1), linetype="dashed")+geom_vline(aes_string(xintercept=1), linetype="dashed"))
  
  dev.off()
  
}


#' plot copy number status estimated from getCNstatus chr by chr
#'
#' @param CNstatus output from getCNstatus()
#' @param BAF three columns (chr, position, BAF values)
#' @param LRR three columns (chr, position, log2(read depth ratio) values)
#' @param sampleName string; output prefix
#'
#' @export
plotCNstatus.chr <- function(CNstatus, BAF, LRR, sampleName){
  
  chr.list<-as.data.frame(table(CNstatus[,1]))[,1]
  
  for (j in 1:length(chr.list)){
    LRR.chr<-LRR[LRR[,1]==as.character(chr.list[j]),]
    BAF.chr<-BAF[BAF[,1]==as.character(chr.list[j]),]
    seg.chr<-CNstatus[CNstatus[,1]==as.character(chr.list[j]),]
    seg.LRR.mat<-data.frame(start=seg.chr[,2], end=seg.chr[,3],LRR=seg.chr[,5],type=rep("LRR", dim(seg.chr)[1]))
    seg.BAF.mat<-data.frame(start=seg.chr[,2], end=seg.chr[,3],BAF=seg.chr[,4],type=rep("BAF", dim(seg.chr)[1]))
    
    
    if (dim(LRR.chr)[1]>0 && dim(BAF.chr)[1]>0){
      pdf(paste(sampleName,".",chr.list[j], ".pdf", sep=""), width = 8, height=5.5)
      
      status<-rep(0, dim(seg.chr)[1]+5)
      status[(1:dim(seg.chr)[1])[seg.chr[,9]=="Gain"]]<-1
      status[(1:dim(seg.chr)[1])[seg.chr[,9]=="Loss" | seg.chr[,9]=="CN Neutral LOH" ]]<--1
      status[length(status)-1]<-1.5
      status[length(status)]<- -1.5
      status.label<-c(as.character(seg.chr[,9]), "Normal","CN Neutral LOH","Ambiguous","Gain","Loss")
      
      combined<-data.frame(Position=c(as.numeric(BAF.chr[,2]), as.numeric(LRR.chr[,2]) ,  as.numeric(seg.chr[,2]),c(-100000001:-100000005)) , start.value=c(as.numeric(BAF.chr[,3]), as.numeric(LRR.chr[,3]) ,  as.numeric(status)), end=c(as.numeric(BAF.chr[,2]), as.numeric(LRR.chr[,2]) , as.numeric(seg.chr[,3]),c(-100000001:-100000005)), end.value=c(as.numeric(BAF.chr[,3]), as.numeric(LRR.chr[,3]) , as.numeric(status)), status.label=c(rep(NA,dim(BAF.chr)[1]), rep(NA, dim(LRR.chr)[1]), status.label), type=c(rep("BAF", dim(BAF.chr)[1]), rep("LRR", dim(LRR.chr)[1]), rep("CN Status", 5+dim(seg.chr)[1])))
      combined$type<-factor(combined$type, levels=c("CN Status", "LRR","BAF")) # fix the order of facet panels
      combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))
      
      g0<-ggplot(subset(combined,type=="CN Status"), aes(Position, start.value))+ coord_cartesian(xlim=c(0-1000000,  max(combined[,3])+1000000))+theme_bw() +theme(strip.text.y = element_text(size = 10, face="bold"), axis.title=element_text(face="bold", size="10"), axis.text.y=element_text(face="bold", size="12"))+geom_segment(aes(x = Position,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G")) 
      p0<-ggplotGrob(g0)
      
      g<-ggplot(combined, aes(Position, start.value))+ coord_cartesian(xlim=c(0-1000000,  max(combined[,3])+1000000))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top",  axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+labs(x = "Position", y = "")+facet_grid(type ~ . , scales = "free")
      
      g1<-g+geom_segment(data=subset(combined,type=="CN Status"),aes(x = Position,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) + scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))
      g2<-g1+ geom_point(data=subset(combined,type=="LRR"),color="darkgrey", size=1)+geom_segment(data=seg.LRR.mat,aes(x = start,y = LRR,xend = end,yend = LRR),size=1, color="blue")
      g3<-g2+ geom_point(data=subset(combined,type=="BAF"),color="darkgrey", size=1) +geom_segment(data=seg.BAF.mat,aes(x = start,y = BAF,xend = end,yend = BAF),size=1, color="blue")
      
      gg_table <- ggplot_gtable(ggplot_build(g3))
      gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
      gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
      grid.draw(gg_table)
      dev.off()                        
    }
    
  }
  
}


# plot LRR, BAF, and CNstatus for whole genome
plotCNstatus.WG<-function(CNstatus, BAF, LRR, sampleName){
  # paramters are the same as plotCNstatus.chr
  
  # prepare matrix for plot
  hg19.length<-c(249250621,243199373,198022430, 191154276, 180915260,171115067,
                 159138663,146364022, 141213431,135534747,135006516,133851895, 115169878,
                 107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895, 
                 51304566)
  
  cum.length<-cumsum(hg19.length)
  length(cum.length)
  cum.length<-c(0,cum.length)
  
  
  # remove chrX and chrY,chrM
  auto.chr<-paste("chr",1:22, sep="")
  LRR.auto<-subset(LRR, LRR[,1]%in%auto.chr)
  BAF.auto<-subset(BAF, BAF[,1]%in%auto.chr)
  
  
  chr<-apply(LRR.auto,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
  LRR.chr<-cbind(chr, LRR.auto[,-1])
  LRR.chr<-LRR.chr[order(LRR.chr[,1],LRR.chr[,2]),]
  row.names(LRR.chr)<-c(1:dim(LRR.chr)[1])
  new.pos<-unlist(apply(LRR.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))	
  LRR.plot.mat<-data.frame(pos=as.numeric(new.pos),log2R=LRR.auto[,3])
  
  
  chr<-apply(BAF.auto,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
  BAF.chr<-cbind(chr, BAF.auto[,2:3])
  BAF.chr<-BAF.chr[order(BAF.chr[,1],BAF.chr[,2]),]
  row.names(BAF.chr)<-c(1:dim(BAF.chr)[1])
  new.pos<-unlist(apply(BAF.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))	
  BAF.plot.mat<-data.frame(pos=as.numeric(new.pos),BAF=BAF.auto[,3])
  
  
  chr<-apply(CNstatus,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
  CNstatus.chr<-cbind(chr, CNstatus[,-1])
  CNstatus.chr<-CNstatus.chr[order(CNstatus.chr[,1],CNstatus.chr[,2]),]
  row.names(CNstatus.chr)<-c(1:dim(CNstatus.chr)[1])
  
  new.start<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))
  new.end<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[3])+cum.length[as.numeric(x[1])]}))
  
  status<-rep(0, dim(CNstatus)[1]+5)
  status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Gain"]]<-1
  status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Loss" | CNstatus.chr[,9]=="CN Neutral LOH" ]]<--1
  status[length(status)-1]<-1.5
  status[length(status)]<- -1.5
  status.label<-c(as.character(CNstatus.chr[,9]), "Normal","CN Neutral LOH","Ambiguous","Gain","Loss")
  
  CNS.plot.mat<-data.frame(start=c(as.numeric(new.start),c(-1:-5)), end=c(as.numeric(new.end),c(-1:-5)), status.start=as.numeric(status), status.end= as.numeric(status), status.label=status.label)	
  
  
  combined<-data.frame(start=c(as.numeric(BAF.plot.mat[,1]), as.numeric(LRR.plot.mat[,1]) , as.numeric(CNS.plot.mat[,1])) , start.value=c(as.numeric(BAF.auto[,3]), as.numeric(LRR.auto[,3]) , as.numeric(CNS.plot.mat[,3])), end=c(as.numeric(BAF.plot.mat[,1]), as.numeric(LRR.plot.mat[,1]) , as.numeric(CNS.plot.mat[,2])), end.value=c(as.numeric(BAF.auto[,3]), as.numeric(LRR.auto[,3]), as.numeric(CNS.plot.mat[,4])), status.label=c(rep(NA,dim(BAF.plot.mat)[1]), rep(NA, dim(LRR.plot.mat)[1]), status.label ), type=c(rep("BAF", dim(BAF.plot.mat)[1]), rep("LRR", dim(LRR.plot.mat)[1]), rep("CN Status", dim(CNS.plot.mat)[1])))
  combined$type<-factor(combined$type, levels=c("CN Status", "LRR","BAF")) # fix the order of facet panels
  combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))
  
  
  combined<-subset(combined, start.value>=-2)
  combined<-subset(combined, start.value<=2)
  
  ##### plot CN status, BAF, and LRR
  
  pdf(paste(sampleName,".","CNstatus.WG.pdf", sep=""), width = 8, height=5.5)
  g0<-ggplot(subset(combined,type=="CN Status"), aes(start, start.value))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), axis.title=element_text(face="bold", size="14"), axis.text.y=element_text(face="bold", size="12"))+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G"))
  p0<-ggplotGrob(g0)
  
  g<-ggplot(combined, aes(start, start.value))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+facet_grid(type ~ . , scales = "free")
  
  g1<-g+geom_segment(data=subset(combined,type=="CN Status"),aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) + scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+ labs(x = "", y = "")
  g2<-g1+ geom_point(data=subset(combined,type=="LRR"),color="grey", size=0.5) + labs(x = "", y = "")
  g3<-g2+ geom_point(data=subset(combined,type=="BAF"),color="grey", size=0.5) + labs(x = "", y = "")
  
  for (i in 1:22) {	
    g3<-g3+annotation_custom(grob = textGrob(i,gp=gpar(fontsize=12,fontface=2), hjust=1, rot=90),  xmin = cum.length[i], xmax = cum.length[i+1], ymin = -7, ymax = -7)
  }
  
  gg_table <- ggplot_gtable(ggplot_build(g3))
  gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
  grid.draw(gg_table)
  dev.off()
}

