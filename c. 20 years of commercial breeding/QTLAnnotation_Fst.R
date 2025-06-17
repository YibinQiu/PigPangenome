#!/usr/bin/env Rscript
# Script: QTLAnnotation_Fst.R
#Date: 2024/9/2
#User: Yibin Qiu
#update: 2024/9/5

library(argparser, quietly = TRUE)
p <- arg_parser("QTL enrichment for S21 and W51 selectionPreference results", hide.opts = TRUE)

# Add command line arguments
p <- add_argument(p, "pop", help = "type DUR/LR", type = "character")
p <- add_argument(p, "pop1", help = "type ADUR/WLR", type = "character")
p <- add_argument(p, "distance", help = "10000/5000", type = "integer", default = "5000")
argv <- parse_args(p)

pop <- argv$pop
pop1 <- argv$pop1
distance <- argv$distance
distance<-as.integer(distance)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: path:", getwd(),"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: pop:", pop,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: distance:", distance,"\n")

#################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: loading Rpackages \n")
library(tidyverse)
library(reshape2)
library(GALLO)
library(ltc)
library(ggplot2)
library(data.table)
library(ggprism)
library(writexl, quietly = TRUE)
library(bedtoolsr)
library(foreach)
library(doParallel)
#library(bedtoolsr)
#library(furrr)
#plan(multisession, workers = availableCores())


############################################################################################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: loading QTLdatabase \n")
#QTL_database can obtained from Pigbiobank (https://pigbiobank.piggtex.bio/download)
QTL_database<-fread("QTL_database.txt",header = T,sep = "\t")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: QTLdatabase data.table to data.frame \n")
QTL_database<-as.data.frame(QTL_database)
QTL_database<-QTL_database[QTL_database$chr%in%paste0("chr",seq(1,18)),]
QTL_database<-QTL_database%>%select(chr,start_pos,end_pos,everything())
QTL_database$chr<-str_replace_all(QTL_database$chr,"chr","")
QTL_database$chr<-as.integer(QTL_database$chr)
QTL_database<-QTL_database%>%arrange(chr,start_pos,end_pos)
##########################################################

#############################################
no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading SV.SelectionPrefernece_sum file \n")

SV<-fread(paste0("SV.SelectionPrefernece_sum.txt"),header = T,sep = "\t")[,c(1,2,3,4,10,7,9)]
SV<-as.data.frame(SV)
SV<-reshape2::dcast(SV, chr+snp+a1+a2+Group_changes ~ Time)
colnames(SV)<-c("CHR","ID","ALT","REF","group","T1","T2","T3")
SV$CHR<-paste0("chr",SV$CHR)
SV$BP1<-str_split(SV$ID,":",simplify = T)[,2]
SV$BP1<-as.integer(SV$BP1)
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading 3.PigPanSV_modifiedSampleID",pop,".avinput file \n")
SV_info<-fread(paste0("3.PigPanSV_modifiedSampleID_",pop1,".avinput"),header = F,sep = "\t")[,1:5]
SV_info<-as.data.frame(SV_info)
colnames(SV_info)<-c("CHR","BP1","BP2","REF","ALT")
SV$MakerType<-"SV"
SV<-left_join(SV,SV_info)
SV<-SV%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType)
SV$BP1<-SV$BP1-0
SV$BP2<-SV$BP2+0
SV$BP1[SV$BP1<=0]<-1
SV$CHR<-str_replace_all(SV$CHR,"chr","")
SV$CHR<-as.integer(SV$CHR)
SV<-SV%>%arrange(CHR,BP1,BP2)
SV_QTL<-bedtoolsr::bt.intersect(SV,QTL_database,loj = "-loj")
colnames(SV_QTL)<-c("CHR","BP1","BP2","ID","REF","ALT","T1","T2","T3","group","MakerType",
                    "chr","start_pos","end_pos","database","QTL_type","QTL_ID",
                    "Name","Abbrev","PUBMED_ID",
                    "trait_ID","trait","traitType",
                    "FlankMarker","breed","VTO_name",
                    "PTO_name","Map_Type","Model",
                    "Test_Base", "significance",
                    "P-value","CMO_name" ,"Dominance_Effect",
                    "Additive_Effect","gene_ID","gene_IDsrc",
                    "Variance","peak_CM","F-Stat","LOD-score",
                    "Likelihood_Ratio","Bayes-value","LS-means","extra_info")
SV_QTL<-SV_QTL[SV_QTL$chr!=".",]
SV_QTL <- SV_QTL%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType,
                          chr,database,QTL_type,
                          start_pos,end_pos,QTL_ID,
                          Name,Abbrev,PUBMED_ID,
                          trait_ID,trait,traitType,
                          FlankMarker,breed,VTO_name,
                          PTO_name,Map_Type,Model,
                          Test_Base, significance,
                          `P-value`,CMO_name ,Dominance_Effect,
                          Additive_Effect,gene_ID,gene_IDsrc,
                          Variance,peak_CM,`F-Stat`,`LOD-score`,
                          Likelihood_Ratio,`Bayes-value`,`LS-means`)


cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading SNP.SelectionPrefernece_sum file \n")
SNP<-fread(paste0("SNP.SelectionPrefernece_sum.txt"),header = T,sep = "\t")[,c(1,2,3,4,10,7,9)]
SNP<-as.data.frame(SNP)
SNP<-reshape2::dcast(SNP, chr+snp+a1+a2+Group_changes ~ Time)
colnames(SNP)<-c("CHR","ID","ALT","REF","group","T1","T2","T3")
SNP$CHR<-paste0("chr",SNP$CHR)
SNP$BP1<-str_split(SNP$ID,":|_",simplify = T)[,2]
SNP$BP1<-as.integer(SNP$BP1)
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading 4.PigPanSNP_modifiedSampleID",pop,".avinput file \n")
SNP_info<-fread(paste0("4.PigPanSNP_modifiedSampleID_",pop1,".avinput"),header = F,sep = "\t")[,1:5]
SNP_info<-as.data.frame(SNP_info)
colnames(SNP_info)<-c("CHR","BP1","BP2","REF","ALT")
SNP_info$ALT[SNP_info$ALT=="0"]<-"*"
SNP$MakerType<-"SNP"
SNP<-left_join(SNP,SNP_info)
SNP<-SNP%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType)
SNP$BP1<-SNP$BP1-0
SNP$BP2<-SNP$BP2+0
SNP$BP1[SNP$BP1<=0]<-1
SNP$CHR<-str_replace_all(SNP$CHR,"chr","")
SNP$CHR<-as.integer(SNP$CHR)
SNP<-SNP%>%arrange(CHR,BP1,BP2)
SNP_QTL<-bedtoolsr::bt.intersect(SNP,QTL_database,loj = "-loj")
colnames(SNP_QTL)<-c("CHR","BP1","BP2","ID","REF","ALT","T1","T2","T3","group","MakerType",
                    "chr","start_pos","end_pos","database","QTL_type","QTL_ID",
                    "Name","Abbrev","PUBMED_ID",
                    "trait_ID","trait","traitType",
                    "FlankMarker","breed","VTO_name",
                    "PTO_name","Map_Type","Model",
                    "Test_Base", "significance",
                    "P-value","CMO_name" ,"Dominance_Effect",
                    "Additive_Effect","gene_ID","gene_IDsrc",
                    "Variance","peak_CM","F-Stat","LOD-score",
                    "Likelihood_Ratio","Bayes-value","LS-means","extra_info")
SNP_QTL<-SNP_QTL[SNP_QTL$chr!=".",]
SNP_QTL <- SNP_QTL%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType,
                          chr,database,QTL_type,
                          start_pos,end_pos,QTL_ID,
                          Name,Abbrev,PUBMED_ID,
                          trait_ID,trait,traitType,
                          FlankMarker,breed,VTO_name,
                          PTO_name,Map_Type,Model,
                          Test_Base, significance,
                          `P-value`,CMO_name ,Dominance_Effect,
                          Additive_Effect,gene_ID,gene_IDsrc,
                          Variance,peak_CM,`F-Stat`,`LOD-score`,
                          Likelihood_Ratio,`Bayes-value`,`LS-means`)



cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading INDEL.SelectionPrefernece_sum file \n")
INDEL<-fread(paste0("INDEL.SelectionPrefernece_sum.txt"),header = T,sep = "\t")[,c(1,2,3,4,10,7,9)]
INDEL<-as.data.frame(INDEL)
INDEL<-reshape2::dcast(INDEL, chr+snp+a1+a2+Group_changes ~ Time)
colnames(INDEL)<-c("CHR","ID","ALT","REF","group","T1","T2","T3")
INDEL$CHR<-paste0("chr",INDEL$CHR)
INDEL$BP1<-str_split(INDEL$ID,":|_",simplify = T)[,2]
INDEL$BP1<-as.integer(INDEL$BP1)
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading 4.PigPanINDEL_modifiedSampleID",pop,".avinput file \n")
INDEL_info<-fread(paste0("4.PigPanINDEL_modifiedSampleID_",pop1,".avinput"),header = F,sep = "\t")[,1:5]
INDEL_info<-as.data.frame(INDEL_info)
colnames(INDEL_info)<-c("CHR","BP1","BP2","REF","ALT")
INDEL_info$ALT[INDEL_info$ALT=="0"]<-"*"
INDEL$MakerType<-"INDEL"
INDEL<-left_join(INDEL,INDEL_info)
INDEL<-INDEL%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType)
INDEL$BP1<-INDEL$BP1-0
INDEL$BP2<-INDEL$BP2+0
INDEL$BP1[INDEL$BP1<=0]<-1
INDEL$CHR<-str_replace_all(INDEL$CHR,"chr","")
INDEL$CHR<-as.integer(INDEL$CHR)
INDEL<-INDEL%>%arrange(CHR,BP1,BP2)
INDEL_QTL<-bedtoolsr::bt.intersect(INDEL,QTL_database,loj = "-loj")
colnames(INDEL_QTL)<-c("CHR","BP1","BP2","ID","REF","ALT","T1","T2","T3","group","MakerType",
                     "chr","start_pos","end_pos","database","QTL_type","QTL_ID",
                     "Name","Abbrev","PUBMED_ID",
                     "trait_ID","trait","traitType",
                     "FlankMarker","breed","VTO_name",
                     "PTO_name","Map_Type","Model",
                     "Test_Base", "significance",
                     "P-value","CMO_name" ,"Dominance_Effect",
                     "Additive_Effect","gene_ID","gene_IDsrc",
                     "Variance","peak_CM","F-Stat","LOD-score",
                     "Likelihood_Ratio","Bayes-value","LS-means","extra_info")
INDEL_QTL<-INDEL_QTL[INDEL_QTL$chr!=".",]
INDEL_QTL <- INDEL_QTL%>%select(CHR,BP1,BP2,ID,REF,ALT,T1,T2,T3,group,MakerType,
                            chr,database,QTL_type,
                            start_pos,end_pos,QTL_ID,
                            Name,Abbrev,PUBMED_ID,
                            trait_ID,trait,traitType,
                            FlankMarker,breed,VTO_name,
                            PTO_name,Map_Type,Model,
                            Test_Base, significance,
                            `P-value`,CMO_name ,Dominance_Effect,
                            Additive_Effect,gene_ID,gene_IDsrc,
                            Variance,peak_CM,`F-Stat`,`LOD-score`,
                            Likelihood_Ratio,`Bayes-value`,`LS-means`)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: conbine SNP/INDEL/SV window and marker overlapped QTL file \n")
all_QTL<-rbind(SV_QTL,SNP_QTL)%>%rbind(INDEL_QTL)
write.table(all_QTL,paste0("variCombined.SelectionPrefernece_all_QTLoverlapped_",distance,"bp.txt"),col.names = T,row.names = F,quote = F,sep = "\t")

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: plot pie for variCombined.SelectionPrefernece_all_QTLoverlapped file \n")
pdf(paste0("variCombined.SelectionPrefernece_all_QTLoverlapped_",distance,"bp.pie.pdf"),height = 9,width = 14)
par(mar=c(1,30,1,1))
plot_qtl_info(all_QTL,qtl_plot = "qtl_type",cex=2,horiz =F)
dev.off()
########################################################################
##########

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: enrich for overlapped QTL file \n")
QTL_enrich<-qtl_enrich(qtl_db=QTL_database,qtl_file=all_QTL,qtl_type=c("Name"),enrich_type="genome",chr.subset=NULL,padj="fdr",nThreads=64)
QTL_enrich<-QTL_enrich%>%arrange(pvalue)
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: write SelectionPrefernece all enrich to xlsx \n")
stopImplicitCluster()
write_xlsx(
  list(variCombined=QTL_enrich),
  paste0("variCombined.SelectionPrefernece_all_QTLenrichment_",distance,"bp.xlsx"),
  col_names = TRUE,
  format_headers = TRUE
)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: all done \n")