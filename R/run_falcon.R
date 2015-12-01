

runFalcon<-function(input, sampleName, threshold=0.15){
  # function that runs Falcon
  # threshold: for a given segments,if for both alleles the abs(copy numbers  -1) < threshold, this segment will not be considered as CNV
  
  library(falcon)
  
  AN = as.numeric(input[,5])
  BN = as.numeric(input[,6])
  AT = as.numeric(input[,7])
  BT = as.numeric(input[,8])
  pos = as.numeric(input[,2])
  
  rdep = median(AT+BT)/median(AN+BN)
  
  result = c()
  
  for (i in 1:22){
    name = paste("chr",i,sep="")
    ids = which(input[,1]==name)
    cn = falcon::getASCN(AT[ids],BT[ids],AN[ids],BN[ids],rdep=rdep)
    #save(cn,file=paste(sampleName,"_chr",i,".Rdata",sep=""))
    png(paste(sampleName,"_chr",i,".png",sep=""),width=1200,height=750,pointsize=23)
    view(AT[ids],BT[ids],AN[ids],BN[ids],cn$tauhat,cn$ascn,pos[ids],rdep=rdep)
    dev.off()
    
    tauhat0 = c(1,cn$tauhat,length(ids)+1)
    m = length(tauhat0)
    start = tauhat0[1:(m-1)]
    end = tauhat0[2:m]-1
    temp = cbind(rep(i,m-1),ids[start], ids[end], pos[ids[start]], pos[ids[end]], round(t(cn$ascn),digits=3))
    result = rbind(result, temp)
  }
  
  removeid = c()
  
  for (j in 1:dim(result)[1]){
    a = result[j,6]
    b = result[j,7]
    if (abs(a-1)<threshold && abs(b-1)<threshold){
      removeid = c(removeid, j)
    }
  }
  
  colnames(result) = c("chr", "SNP.start", "SNP.end", "pos.start", "pos.end", "minor.CN", "major.CN")
  #return(result[-removeid,])
  return(result)
  # write.table(result[-removeid,], file=paste(samplename,"_ascn_table.txt",sep=""), quote=F, row.names=F, sep="\t")
  
}
