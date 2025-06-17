#!/usr/bin/env Rscript
# Script: 09.Genetic_gain_per_year.R
# date: 2024/10/11
# update: 2024/10/21
# update: 2024/10/24
# user: Yibin Qiu
library(data.table)
library(tidyverse)
library(ggprism)
library(ggpubr)
library(ggbreak)
library(introdataviz)
library(scales)
library(patchwork)
setwd("./S21/")
snp_S21<-fread(input = "SNP.SelectionPrefernece_sum.txt")
indel_S21<-fread(input = "INDEL.SelectionPrefernece_sum.txt")
sv_S21<-fread(input = "SV.SelectionPrefernece_sum.txt")
all_S21<-rbind(snp_S21,indel_S21)%>%rbind(sv_S21)
all_dcast_S21<-all_S21%>%select(chr,snp,a1,a2,Group_changes,n,Time,maf)
all_dcast_S21<-reshape2::dcast(all_dcast_S21,chr+snp+a1+a2+Group_changes+n ~Time)
all_dcast_S21%>%group_by(Group_changes)%>%summarise(N=n())

up_flat_S21<-all_dcast_S21[all_dcast_S21$Group_changes%in%c("Up_Flat","Down_Flat"),]
up_flat_S21$gain_value<-abs(up_flat_S21$T2 - up_flat_S21$T1)/10
up_flat_S21$type<-"Phenotypic selection"

flat_up_S21<-all_dcast_S21[all_dcast_S21$Group_changes%in%c("Flat_Up","Flat_Down"),]
flat_up_S21$gain_value<-abs(flat_up_S21$T3 - flat_up_S21$T2)/7
flat_up_S21$type<-"Genomic selection"

#wilcox.test(up_flat_S21$gain_value,flat_up_S21$gain_value)


select_S21<-rbind(up_flat_S21,flat_up_S21)
#select2_S21<-select_S21%>%group_by(type)%>%summarise(mean_value=mean(gain_value))
select_S21$type<-factor(select_S21$type,levels = c("Phenotypic selection","Genomic selection"))
select_S21$pop<-"S21 Duroc"

setwd("./W51/")
snp_W51<-fread(input = "SNP.SelectionPrefernece_sum.txt")
indel_W51<-fread(input = "INDEL.SelectionPrefernece_sum.txt")
sv_W51<-fread(input = "SV.SelectionPrefernece_sum.txt")
all_W51<-rbind(snp_W51,indel_W51)%>%rbind(sv_W51)
all_dcast_W51<-all_W51%>%select(chr,snp,a1,a2,Group_changes,n,Time,maf)
all_dcast_W51<-reshape2::dcast(all_dcast_W51,chr+snp+a1+a2+Group_changes+n ~Time)


up_flat_W51<-all_dcast_W51[all_dcast_W51$Group_changes%in%c("Up_Flat","Down_Flat"),]
up_flat_W51$gain_value<-abs(up_flat_W51$T2 - up_flat_W51$T1)/9
up_flat_W51$type<-"Phenotypic selection"

flat_up_W51<-all_dcast_W51[all_dcast_W51$Group_changes%in%c("Flat_Up","Flat_Down"),]
flat_up_W51$gain_value<-abs(flat_up_W51$T3 - flat_up_W51$T2)/6
flat_up_W51$type<-"Genomic selection"

#wilcox.test(up_flat_W51$gain_value,flat_up_W51$gain_value)
mean(up_flat_S21$gain_value)
mean(flat_up_S21$gain_value)
mean(up_flat_W51$gain_value)
mean(flat_up_W51$gain_value)
select_W51<-rbind(up_flat_W51,flat_up_W51)
#select2_W51<-select_W51%>%group_by(type)%>%summarise(mean_value=mean(gain_value))
select_W51$type<-factor(select_W51$type,levels = c("Phenotypic selection","Genomic selection"))
select_W51$pop<-"W51 Landrace"

dataplot<-rbind(select_S21,select_W51)%>%select(snp,gain_value,type,pop)
dataplot$pop<-factor(dataplot$pop,levels = c("S21 Duroc","W51 Landrace"))
setwd("./genetic change per year/")
dataplot$vari[dataplot$snp%in%c(snp_S21$snp,snp_W51$snp)]<-"SNP"
dataplot$vari[dataplot$snp%in%c(indel_S21$snp,indel_W51$snp)]<-"INDEL"
dataplot$vari[dataplot$snp%in%c(sv_S21$snp,sv_W51$snp)]<-"SV"
p1<-ggplot(dataplot,aes(x = pop, y =gain_value ,fill= type))+
  geom_split_violin(trim =F,color = NA,adjust = 1.5)+
  guides(fill=guide_legend(title="group"))+
  stat_summary(fun.data = "mean_sd",position=position_dodge(0.15),geom = "errorbar",width = .1) +
  stat_summary(fun = "mean", geom = "point", position=position_dodge(0.15),show.legend = F)+
  stat_compare_means(aes(group = type), label = "p.signif",label.y=c(0.07,0.09), method="t.test",size=5)+
  scale_fill_manual(values=c("#788FCE","#E6956F"))+
  scale_y_continuous(breaks = c(0,0.02,0.04,0.06,0.08,0.1))+
  coord_cartesian(ylim = c(0,0.10))+
  labs(x=NULL,y="Genetic changes per year")+
  theme_classic(base_size = 16)+
  theme(
    axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
    panel.background = element_rect(fill = NA,color = NA),
    legend.key=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(color="black",size=10),
    legend.position = c(0,1), legend.justification = c(0,1),
    legend.background=element_blank())+guides(fill="none")

ggsave("genetic_gain_per_year.pdf",width = 2.5,height = 4)
ggsave("genetic_gain_per_year.png",width = 3,height = 6,dpi = 300)

dataplot2<-dataplot%>%group_by(pop,type)%>%summarise(N=n())

p2<-ggplot(dataplot2,aes(x = pop, y =N/100000 ,fill= type))+
  geom_bar(stat = "identity",position="dodge",width=0.5)+
  guides(fill=guide_legend(title="group"))+
  scale_fill_manual(values=c("#788FCE","#E6956F"))+
  scale_y_break(c(3,9),scales = "fixed",ticklabels = c(9))+
  labs(x=NULL,y=expression('Selection variant count (x 10'^5*')'))+
  theme_classic(base_size = 16)+
  theme(
    axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
    axis.title.y = element_text(hjust = 0.7),
    panel.background = element_rect(fill = NA,color = NA),
    legend.key=element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(color="black",size=10),
    legend.position = c(0,1), legend.justification = c(0,1),
    legend.background=element_blank())+guides(fill="none")
ggsave("selection_variant_count.pdf",width = 2.5,height = 4)
ggsave("selection_variant_count.png",width = 3,height = 6,dpi = 300)

p3<-p2|p1
ggsave(plot = p3,"selection_variant_total.pdf",width = 6,height = 6)
#dataplot3<-dataplot%>%group_by(pop,type,vari)%>%summarise(N=n())

#ggdraw()+ 
#  draw_plot(p1)+
#  draw_plot(p2,x =0.2, y = 0.6,width = 0.5, height = 0.5,
#            scale = 0.5) 

#ggsave("test.png",width = 8,height = 7,dpi = 300)
save.image(file = "genetic_gain_per_year.RData")
