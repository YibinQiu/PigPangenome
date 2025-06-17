#!/usr/bin/env Rscript
# Script: QTLAnnotation_XPCLR.R
#Date: 2024/9/1
#User: Yibin Qiu
#update: 2024/9/2
#update: 2024/9/4
#update: 2024/9/23
library(argparser, quietly = TRUE)
p <- arg_parser("QTL enrichment for xpclr and Fst intersect results", hide.opts = TRUE)

# Add command line arguments
p <- add_argument(p, "pop1", help = "type DUR/LR/LW", type = "character")
p <- add_argument(p, "pop2", help = "type 20MDUR/MLR/MLW", type = "character")
p <- add_argument(p, "distance", help = "10000/5000", type = "integer", default = "5000")
argv <- parse_args(p)

pop1 <- argv$pop1
pop2 <- argv$pop2
distance <- argv$distance
distance<-as.integer(distance)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: path:", getwd(),"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: pop1:", pop1,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: pop2:", pop2,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: distance:", distance,"\n")

#################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: loading Rpackages \n")
library(tidyverse)
library(GALLO)
library(ltc)
library(ggplot2)
library(ggprism)
library(data.table)
############################################################################################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: loading QTLdatabase \n")
#QTL_database can obtained from Pigbiobank (https://pigbiobank.piggtex.bio/download)
QTL_database<-fread("QTL_database.txt",header = T,sep = "\t")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: QTLdatabase data.table to data.frame \n")
QTL_database<-as.data.frame(QTL_database)
QTL_database<-QTL_database[QTL_database$chr%in%paste0("chr",seq(1,18)),]
##########################################################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reset col.pallet in plot_qtl_info() \n")
col.pallet<-as.character(ltc("minou"))
names(col.pallet)<-unique(QTL_database$QTL_type)
plot_qtl_info<-function(qtl_file,qtl_plot=c("qtl_type","qtl_name"), n="all",qtl_class=NULL,horiz=FALSE,...){
  #check method
  qtl_plot <- match.arg(qtl_plot)
  
  qtl_file=qtl_file
  if(qtl_plot=="qtl_type"){
    qtl_file<-qtl_file[order(qtl_file$QTL_type),]
    data.prop<-as.data.frame(table(qtl_file$QTL_type))
    data.prop$Var1<-gsub("_"," ",data.prop$Var1)
    labels.prop<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
    labels.comp<-paste(labels.prop, "%", sep="")
    #col.pallet<-RColorBrewer::brewer.pal(n=nrow(data.prop), name="Set3")
    col.pallet<-col.pallet
    pie(x=data.prop[,2], labels=labels.comp, col=col.pallet, ...)
    legend(x=-2.2,y=0.2, pch=15, col=col.pallet, legend=unique(data.prop$Var1),bty="n",horiz = horiz, xpd = TRUE,inset = c(0,0),y.intersp=1.2,xjust=0,yjust=0, ...)
    
  }
  
  if(qtl_plot=="qtl_name"){
    data.prop<-as.data.frame(table(qtl_file$Name))
    data.prop$percentage<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
    data.prop<-data.prop[order(data.prop$percentage, decreasing=TRUE),]
    if(length(which(qtl_class %in% "all"))!=0){
      
      if(n=="all"){
        n_final<-nrow(data.prop)
      }else{
        n_final<-n
      }
      barplot(data.prop$percentage[seq_along(1:n_final)], names.arg=data.prop$Var1[seq_along(1:n_final)],horiz=TRUE, xlab="QTL Names (%)", las=1, xlim=c(0,round((max(data.prop$percentage[seq_along(1:n_final)])+0.5),0)), ...)
      
    }
    
    if(length(which(qtl_class %in% "all"))==0){
      data.prop<-as.data.frame(table(qtl_file$Name))
      qtl_file<-qtl_file[-duplicated(qtl_file[,c("QTL_type","Name")]),]
      data.prop$qtl_type<-qtl_file[match(as.character(data.prop[,1]),as.character(qtl_file$Name)),"QTL_type"]
      data.prop$percentage<-c(round(((data.prop[,2])/sum(data.prop[,2]))*100,2))
      data.prop<-data.prop[order(data.prop$percentage, decreasing=TRUE),]
      data.prop<-data.prop[which(data.prop$qtl_type %in% qtl_class),]
      
      if(n=="all"){
        n_final<-nrow(data.prop)
      }else{
        n_final<-n
      }
      
      barplot(data.prop$percentage[seq_along(1:n_final)], names.arg=data.prop$Var1[seq_along(1:n_final)],horiz=TRUE, xlab="QTL Names (%)", las=1, xlim=c(0,round((max(data.prop$percentage[seq_along(1:n_final)])+0.5),0)),...)
      
    }
  }
}

#############################################
QTL_database<-QTL_database%>%select(chr,database,QTL_type,start_pos,end_pos,extra_info)
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading SV top5 file \n")
SV<-read.table(paste0(pop1,"_vs_",pop2,"_","SV","/",pop1,"_",pop2,"_top5_SV.avinput"),header = F,sep = "\t")[,1:5]
colnames(SV)<-c("CHR","BP1","BP2","REF","ALT")
SV$MakerType<-"SV"
SV_QTL<-find_genes_qtls_around_markers(db_file=QTL_database,marker_file=SV, method = "qtl", marker = "haplotype", interval = 5000, nThreads = 1)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading SNP intersect top5 window file \n")
SNP<-read.table(paste0(pop1,"_vs_",pop2,"_","SNP","/",pop1,"_",pop2,"_top5_intersect.avinput"),header = F,sep = "\t")[,1:5]
colnames(SNP)<-c("CHR","BP1","BP2","REF","ALT")
SNP$MakerType<-"SNP window"
SNP_QTL<-find_genes_qtls_around_markers(db_file=QTL_database,marker_file=SNP, method = "qtl", marker = "haplotype", interval = 5000, nThreads = 1)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: reading INDEL intersect top5 window file \n")
INDEL<-read.table(paste0(pop1,"_vs_",pop2,"_","INDEL","/",pop1,"_",pop2,"_top5_intersect.avinput"),header = F,sep = "\t")[,1:5]
colnames(INDEL)<-c("CHR","BP1","BP2","REF","ALT")
INDEL$MakerType<-"INDEL window"
INDEL_QTL<-find_genes_qtls_around_markers(db_file=QTL_database,marker_file=INDEL, method = "qtl", marker = "haplotype", interval = 5000, nThreads = 1)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: conbine SNP/INDEL/SV window and marker overlapped QTL file \n")
all_QTL<-rbind(SV_QTL,SNP_QTL)%>%rbind(INDEL_QTL)
write.table(all_QTL,paste0(pop1,"_",pop2,"_QTLoverlapped_",distance,"bp.txt"),col.names = T,row.names = F,quote = F,sep = "\t")

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: plot pie for overlapped QTL file \n")
pdf(paste0(pop1,"_",pop2,"_QTLoverlapped_",distance,"bp.pie.pdf"),height = 9,width = 14)
par(mar=c(1,30,1,1))
plot_qtl_info(all_QTL,qtl_plot = "qtl_type",cex=2,horiz =F)
dev.off()

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: enrich for overlapped QTL file \n")
QTL_enrich<-qtl_enrich(qtl_db=QTL_database,qtl_file=all_QTL,qtl_type=c("Name"),enrich_type="genome",chr.subset=NULL,padj="fdr",nThreads=1)
QTL_enrich<-QTL_enrich%>%arrange(pvalue)
write.table(QTL_enrich,paste0(pop1,"_",pop2,"_QTLenrichment_",distance,"bp.txt"),col.names = T,row.names = F,quote = F,sep = "\t")

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: all done \n")