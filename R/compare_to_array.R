



# compare copy number estimated by CLOSER to TCGA resutls
compareToArray<-function(CNstatus, SNParray,sampleName){
  # parameters:
  # 1) CNstatus: output from getCNstatus()
  # 2) SNParray: matrix with 4 columns 
  #(1) chr: integer
  #(2) start: numeric
  #(3) end: numeric
  #(4) copy number estimates
  # 3) sampleName: output prefix
  
  
  hg19.length<-c(249250621,243199373,198022430, 191154276, 180915260,171115067,159138663,146364022, 141213431,135534747,135006516,133851895, 115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895, 51304566)
  cum.length<-cumsum(hg19.length)
  length(cum.length)
  cum.length<-c(0,cum.length)
  
  new.start<-unlist(apply(TCGA, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))
  new.end<-unlist(apply(TCGA, 1, function(x) {as.numeric(x[3])+cum.length[as.numeric(x[1])]}))
  
  col<-rep("Normal", dim(TCGA)[1])
  col[TCGA[,4]< -0.05]<-"Loss"
  col[TCGA[,4]> 0.05]<-"Gain"
  seg_mean<-TCGA[,4]
  for (i in 1:length(seg_mean)){
    if (TCGA[i,4]>2){
      seg_mean[i]=2
    }else if (TCGA[i,4]< -2){
      seg_mean[i]=-2
    }
  }
  
  TCGA.mat<-data.frame(start=new.start, end=new.end, segment_mean.start=rep(0,length(seg_mean)), segment_mean.end=as.numeric(seg_mean),col=factor(col))
  
  
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
  
  
  combined<-data.frame(start=c(as.numeric(CNS.plot.mat[,1]), as.numeric(TCGA.mat[,1])) , start.value=c(as.numeric(CNS.plot.mat[,3]), as.numeric(TCGA.mat[,3])), end=c(as.numeric(CNS.plot.mat[,2]), as.numeric(TCGA.mat[,2])), end.value=c(as.numeric(CNS.plot.mat[,4]), as.numeric(TCGA.mat[,4])), status.label=c(status.label,col), type=c( rep("CLOSE-R", dim(CNS.plot.mat)[1]), rep("SNP Array", dim(TCGA.mat)[1])))
  combined$type<-factor(combined$type, levels=c("CLOSE-R", "SNP Array")) # fix the order of facet panels
  combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))
  
  
  
  pdf(paste(sampleName,".","compareToTCGA.pdf", sep=""), width = 8, height=5)
  g0<-ggplot(subset(combined,type=="CLOSE-R"), aes(start, start.value,color=status.label))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G"))
  p0<-ggplotGrob(g0)
  
  g_0<-ggplot(subset(combined,type=="SNP Array"), aes(start, start.value,color=status.label,fill=status.label))+ coord_cartesian(ylim=c(-2,2),xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+labs(x = "", y = "")+geom_rect( aes(xmin = start,xmax=end, ymin = 0, ymax=end.value, color=status.label,fill=status.label))
  p_0<-ggplotGrob(g_0)
  
  
  g<-ggplot(combined, aes(start, start.value, color=status.label,fill=status.label))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+scale_fill_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+ labs(x = "", y = "")+facet_grid(type ~ . , scales = "free")
  
  g1<-g+geom_segment(data=subset(combined,type=="CLOSE-R"),aes(x = start,y = start.value,xend = end,yend = end.value, color=status.label,fill=status.label),size=3) 
  g2<-g1+geom_rect(data=subset(combined,type=="SNP Array"), aes(xmin = start,xmax=end, ymin = 0, ymax=end.value, color=status.label,fill=status.label))
  for (i in 1:22) {	
    g2<-g2+annotation_custom(grob = textGrob(i,gp=gpar(fontsize=12,fontface=2), hjust=1, rot=90),  xmin = cum.length[i], xmax = cum.length[i+1], ymin = -2.2, ymax = -2.2)
  }
  
  gg_table <- ggplot_gtable(ggplot_build(g2))
  gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
  gg_table[["grobs"]][[3]]<-p_0[["grobs"]][[2]]
  
  grid.draw(gg_table)
  dev.off()
  
}



