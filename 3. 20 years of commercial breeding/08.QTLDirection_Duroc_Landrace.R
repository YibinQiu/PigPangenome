#!/usr/bin/env Rscript
# Script: 08.QTLDirection_Duroc_Landrace.R
#Date: 2024/9/1
#User: Yibin Qiu
#update: 2024/9/2
#update: 2024/9/4
#update: 2024/9/23
library(writexl)
library(data.table)
library(tidyverse)
setwd("./S21/")
S21_go<-readxl::read_xlsx("variCombined.SelectionPrefernece_all_QTLenrichment_QTLdirection_5000bp.xlsx",sheet=1)
S21_increase<-fread("variCombined.SelectionPrefernece_all_QTLoverlapped_Combined_5000bp.txt",header = T)
S21_increase$group[S21_increase$group%in%c("Increase","Up_Flat","Flat_Up")]<-"UP"
S21_increase$group[S21_increase$group%in%c("Decrease","Down_Flat","Flat_Down")]<-"DOWN"
S21_increase<-S21_increase[S21_increase$group%in%c("UP","DOWN"),]

setwd("./W51/")
W51_go<-readxl::read_xlsx("variCombined.SelectionPrefernece_all_QTLenrichment_QTLdirection_5000bp.xlsx",sheet=1)
W51_increase<-fread("variCombined.SelectionPrefernece_all_QTLoverlapped_Combined_5000bp.txt",header = T)
W51_increase$group[W51_increase$group%in%c("Increase","Up_Flat","Flat_Up")]<-"UP"
W51_increase$group[W51_increase$group%in%c("Decrease","Down_Flat","Flat_Down")]<-"DOWN"
W51_increase<-W51_increase[W51_increase$group%in%c("UP","DOWN"),]
########
common_QTLs<-intersect(S21_go$QTL,W51_go$QTL)
total<-NULL
for (QTL in common_QTLs) {
  #QTL="Total litter weight of piglets born alive"
  subcategory<-S21_go$subcategory[S21_go$QTL==QTL]
  maincategory<-S21_go$maincategory[S21_go$QTL==QTL]
  S21_increase_QTL<-S21_increase[S21_increase$Name==QTL,]
  W51_increase_QTL<-W51_increase[W51_increase$Name==QTL,]
  S21QTLnumber<-length(unique(S21_increase_QTL$QTL_ID))
  W51QTLnumber<-length(unique(W51_increase_QTL$QTL_ID))
  commonQTLnumber<-length(intersect(unique(S21_increase_QTL$QTL_ID),unique(W51_increase_QTL$QTL_ID)))
  S21_QTL_id<-S21_increase_QTL[S21_increase_QTL$ID%in%intersect(S21_increase_QTL$ID,W51_increase_QTL$ID),c("ID","group")]
  W51_QTL_id<-W51_increase_QTL[W51_increase_QTL$ID%in%intersect(S21_increase_QTL$ID,W51_increase_QTL$ID),c("ID","group")]
  S21_QTL_id<-S21_QTL_id[!duplicated(S21_QTL_id),]
  colnames(S21_QTL_id)<-c("ID","S21_group")
  W51_QTL_id<-W51_QTL_id[!duplicated(W51_QTL_id),]
  colnames(W51_QTL_id)<-c("ID","W51_group")
  
  data<-left_join(S21_QTL_id,W51_QTL_id)
  common_var<-nrow(data)
  consist_var<-sum(data$S21_group==data$W51_group)

  df<-tribble(~QTL,~S21QTLnumber,~W51QTLnumber,~commonQTLnumber,~common_var,~consist_var,~ratio,~subcategory,~maincategory,
              QTL,S21QTLnumber,W51QTLnumber,commonQTLnumber,common_var,consist_var,consist_var/common_var,subcategory,maincategory)
  total<-rbind(total,df)
}
total<-left_join(total,S21_go[,c("QTL","pvalue")])
colnames(total)<-c(colnames(total)[1:9],"S21pvalue")
total<-left_join(total,W51_go[,c("QTL","pvalue")])
colnames(total)<-c(colnames(total)[1:10],"W51pvalue")
#remian<-total[total$S21pvalue<=0.05 & total$W51pvalue<=0.05,]
remian<-total[!is.na(total$ratio),]
remian<-remian[remian$maincategory!="Adaptation",]
remian<-left_join(remian,S21_go[,c("QTL","N_QTLs_db")])
remian<-remian%>%select(QTL,N_QTLs_db,S21QTLnumber,W51QTLnumber,commonQTLnumber,common_var,consist_var,ratio,subcategory,maincategory,S21pvalue,W51pvalue)

con<-remian[remian$ratio>=quantile(remian$ratio,0.7),]
con%>%group_by(subcategory)%>%summarise(n=n())

dif<-remian[remian$ratio<=quantile(remian$ratio,0.3),]
dif%>%group_by(subcategory)%>%summarise(n=n())
write.table(remian,"S21&W51_QTL_diretion_compare.txt",row.names = F,,col.names = T,quote = F,sep = "\t")
write.table(con,"S21&W51_QTL_diretion_compare_convergent_Terms.txt",row.names = F,,col.names = T,quote = F,sep = "\t")
write.table(dif,"S21&W51_QTL_diretion_compare_divergent_Terms.txt",row.names = F,,col.names = T,quote = F,sep = "\t")
S21_spec<-setdiff(S21_go$QTL,W51_go$QTL)
W51_spec<-setdiff(W51_go$QTL,S21_go$QTL)
