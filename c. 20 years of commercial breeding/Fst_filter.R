#Rscript
#filename: Fst_filter.R
#data: 2024/4/7
#update: 2024/7/2
#author: Yibin Qiu
rm(list=ls());options(stringsAsFactors=FALSE)

args = commandArgs(T)
filepath = args[1]
pop1 = args[2]
pop2 = args[3]
value= args[4]
setwd(filepath)
#value=0.1
value=as.numeric(value)
library("data.table")
cat("*********************************************************","\n")
cat("reading file:",paste0(pop1,"_",pop2,".weir.fst"),"\n")
fst<-fread(paste0(pop1,"_",pop2,".fst"),header=T)
fst<-fst[!is.na(fst$Fst),]
cat("filter fst with " , value, " and write files:","\n")
fst_01<-fst[fst$Fst>=value,]
fst_01_ID<-fst_01[,c("SNP")]
fwrite(fst_01,paste0(pop1,"_",pop2,".fst_",value),append = FALSE,row.names = FALSE,col.names = TRUE,quote = FALSE,sep= "\t")
fwrite(fst_01_ID,paste0(pop1,"_",pop2,".fst_",value,".ID"),append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE,sep= "\t")
cat("all done","\n")
cat("*********************************************************","\n")