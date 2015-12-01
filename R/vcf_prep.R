VCFprep<-function(Normal.vcf, Tumor.vcf,  filter = FALSE){
  
  # calculate BAF and LRR from the .vcf files of Normal and Tumor sample
  
  # paramters:
  # 1) Normal.vcf and Tumor.vcf are the path to the two .vcf files
  # 2) filter: logical. Whether to filter the .vcf files based on FILTER (PASS)
  
  library(VariantAnnotation)
  
  N.vcf = readVcf(Normal.vcf, "hg19") 
  if (!("AO"%in%names(geno(N.vcf)) && "RO"%in%names(geno(N.vcf)))) {
    #N.valid=F
    stop("Normal vcf failed vcf foramt check: No RO/AO information!\n")
  }else{
    cat ("Normal.vcf passed vcf foramt check\n")
  }
  T.vcf = readVcf(Tumor.vcf, "hg19")
  if (!("AO"%in%names(geno(T.vcf)) && "RO"%in%names(geno(T.vcf)))) {
    #T.valid=F
    stop("Tumor vcf failed vcf foramt check: No RO/AO information!\n")
  }else{
    cat ("Tumor.vcf passed vcf foramt check\n")
  }
  
  
  T.FILTER = filt(T.vcf)
  N.FILTER = filt(N.vcf)
  
  
  if (filter) {
    T.pass = which(T.FILTER=="PASS")
    N.pass = which(N.FILTER=="PASS")
  } else {
    T.pass = 1:length(T.FILTER)
    N.pass = 1:length(N.FILTER)
  }
  
  T.FILTER = T.FILTER[T.pass]
  N.FILTER = N.FILTER[N.pass]
  
  T.POS = start(ranges(rowRanges(T.vcf)))[T.pass]
  N.POS = start(ranges(rowRanges(N.vcf)))[N.pass]
  
  T.CHR = as.character(seqnames(rowRanges(T.vcf)))[T.pass]
  N.CHR = as.character(seqnames(rowRanges(N.vcf)))[N.pass]
  
  
  #T.GT = t(sapply(geno(T.vcf)$GT[T.pass,], function(x) unlist(strsplit(as.character(x), split="/"))))
  #N.GT = t(sapply(geno(N.vcf)$GT[N.pass,], function(x) unlist(strsplit(as.character(x), split="/"))))
  
  T.GT= as.vector(geno(T.vcf)$GT[T.pass,])
  N.GT= as.vector(geno(N.vcf)$GT[N.pass,])
  
  T.REF = as.vector(ref(T.vcf))[T.pass]
  N.REF = as.vector(ref(N.vcf))[N.pass]
  
  T.ALT = unstrsplit(CharacterList(alt(T.vcf)), sep = ",")[T.pass]
  N.ALT = unstrsplit(CharacterList(alt(N.vcf)), sep = ",")[N.pass]
  
  T.AO = unstrsplit(CharacterList(geno(T.vcf)$AO), sep=",")[T.pass]
  N.AO = unstrsplit(CharacterList(geno(N.vcf)$AO), sep=",")[N.pass]
  
  T.RO = as.vector(geno(T.vcf)$RO)[T.pass]
  N.RO = as.vector(geno(N.vcf)$RO)[N.pass]
  
  T.mat = cbind(T.CHR, T.POS, T.REF, T.ALT, T.GT, T.AO, T.RO)
  N.mat = cbind(N.CHR, N.POS, N.REF, N.ALT, N.GT, N.AO, N.RO)
  
  T.DP = as.vector(geno(T.vcf)$DP)[T.pass]
  N.DP = as.vector(geno(N.vcf)$DP)[N.pass]
  
  normalize.factor = sum(as.numeric(T.DP))/(sum(as.numeric(N.DP)))
  ### match Normal and Tumor positions
  
  T.index = N.index = c()
  
  
  for (i in 1:22){
    name = paste("chr", i, sep="")
    T.id = which(T.CHR==name)
    N.id = which(N.CHR==name)
    T.pos = T.POS[T.id]
    N.pos = N.POS[N.id]
    comp = match(T.pos, N.pos)
    T.match = (1:length(T.id))[-which(is.na(comp))]
    N.match = comp[-which(is.na(comp))]
    T.index = c(T.index,T.id[T.match])
    N.index = c(N.index,N.id[N.match])
    
  }
  
  T.match.mat = T.mat[T.index,]
  N.match.mat = N.mat[N.index,]
  
  ### get the heterozygous SNPs in normal sample
  N.GT.mat = matrix(as.numeric(sapply(as.list(N.match.mat[,5]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)
  
  dif = N.GT.mat[,2] - N.GT.mat[,1]
  het = which(dif>0)
  
  N.matrix.het = N.match.mat[het,]
  T.matrix.het = T.match.mat[het,]
  mydata = cbind(N.matrix.het[,c(1,2,3)], T.matrix.het[,3], N.matrix.het[,4], T.matrix.het[,4], N.matrix.het[,5],T.matrix.het[,5],N.matrix.het[,6],T.matrix.het[,6], N.matrix.het[,7], T.matrix.het[,7])
  colnames(mydata) = c("chr", "pos", "N.ref", "T.ref", "N.alt", "T.alt", "N.GT", "T.GT", "N.AO", "T.AO", "N.RO","T.RO")
  
  ### remove the SNPs where Tumor and normal sample has different ref or alt alleles
  ref.dif = which(as.character(mydata[,3])!=as.character(mydata[,4]))
  alt.dif = which(as.character(mydata[,5])!=as.character(mydata[,6])) # ref.dif is uauslly a subset of alt.dif
  
  mydata1 = mydata[-unique(c(ref.dif, alt.dif)),]
  
  ### remove the SNPs that has more than 2 alternative alleles in normal sample
  mydata2 = mydata1[which(as.character(mydata1[,7])=="0/1" | as.character(mydata1[,7])=="0/2" | as.character(mydata1[,7])=="1/2"),]
  #N.GT.mat = matrix(as.numeric(sapply(as.list(mydata2[,7]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)
  #T.GT.mat = matrix(as.numeric(sapply(as.list(mydata2[,8]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)
  
  mydata2 = mydata2[!apply(mydata2,1,function(x) {max(unlist(strsplit(x[7], split="/"))) !=  max(unlist(strsplit(x[8], split="/")))}),]
  ### get number of reads for the two alleles A and B
  ### AN: number of reads for allele A
  ### BN: number of reads for allele B
  
  reads = matrix(0,dim(mydata2)[1],4)
  for (i in 1:dim(mydata2)[1]) {
    #cat("i is ", i, "\n")
    if (as.character(mydata2[i,7])=="0/1"){
      reads[i,]=c(as.numeric(mydata2[i,11]), as.numeric(mydata2[i,9]), as.numeric(mydata2[i,12]), as.numeric(mydata2[i,10]))
    }else {
      temp1 = as.numeric(unlist(strsplit(as.character(mydata2[i,9]), split=",")))
      temp2 = as.numeric(unlist(strsplit(as.character(mydata2[i,10]),split=",")))
      
      if (mydata2[i,7]=="1/2") {
        reads[i,] = c(temp1, temp2)
      } else {
        reads[i,] = c(as.numeric(mydata2[,11]),temp1[2], as.numeric(mydata2[,12]),temp2[2])
      }
    }
  }
  
  
  Afreq.N = round(reads[,1]/(reads[,1]+reads[,2]),3)
  Afreq.T = round(reads[,3]/(reads[,3]+reads[,4]),3)
  LRR = log(((reads[,3]+reads[,4])/normalize.factor)/(reads[,1]+reads[,2]),2)
  output = cbind(mydata2[,c(1,2,3,5)], reads, Afreq.N, Afreq.T, LRR)
  colnames(output) = c("chr", "pos", "ref", "alt", "AN", "BN", "AT", "BT", "Afreq.N", "Afreq.T", "LRR")
  
  return(output)
}
