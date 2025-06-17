#!/usr/bin/env Rscript
# Script: SelectionPrefernece_Fst.R
# date: 2024/4/11
# update: 2024/7/21
# update: 2024/7/24
# user: Yibin Qiu
library(argparser, quietly=TRUE)
p <- arg_parser("Selection preference from Fst freq table.",
                hide.opts = T )
# Add command line arguments
p <- add_argument(p, "fst", help="input: frq.strat", type="character" )
p <- add_argument(p, "deltaAF1", help="input: 0.2~0.7", type="numeric",default = 0.4 )
p <- add_argument(p, "deltaAF2", help="input: 0.2~0.7", type="numeric",default = 0.5 )
p <- add_argument(p, "T1_T2_fst01", help="for example: 03ADUR_13ADUR.fst_0.1.ID", type="character" )
p <- add_argument(p, "T2_T3_fst01", help="for example: 13ADUR_21ADUR.fst_0.1.ID", type="character" )
p <- add_argument(p, "outpre", help="output prefix", type="character", default = "Result" )
p <- add_argument(p, "exclude", help="exclude file, which SNP/INDEL/SV ID want to exclude?", type="character", default = "INDEL_large30bp.id" )
argv <- parse_args(p)
input <- argv$fst
deltaAF1 <- argv$deltaAF1
deltaAF2 <- argv$deltaAF2
T1_T2_fst01<-argv$T1_T2_fst01
T2_T3_fst01<-argv$T2_T3_fst01
pre <- argv$outpre
exclude<-argv$exclude
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: path:", getwd(),"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: input:", input,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: deltaAF1:", deltaAF1,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: deltaAF2:", deltaAF2,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: T1_T2_fst01:", T1_T2_fst01,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: T2_T3_fst01:", T2_T3_fst01,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: pre:", pre,"\n")
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: exclude:", exclude,"\n")

deltaAF1<-as.numeric(deltaAF1)
deltaAF2<-as.numeric(deltaAF2)
library(data.table)
library(reshape2)
library(tidyverse)
library(ggh4x)
library(ggprism)
## 数据格式
# CHR                 SNP     CLST   A1   A2      MAF    MAC  NCHROBS
# 1      chr1:13105_C_T   03ADUR    T    C    0.125      5       40
# 1      chr1:13105_C_T   13ADUR    T    C        0      0       60
# 1      chr1:13105_C_T   21ADUR    T    C  0.02542      3      118

fst<-fread(input,header = T)
ex_file<-fread(exclude,header = F)
fst<-fst[!fst$SNP%in%ex_file$V1,]
fst<-fst[,c(1:6)]
cp<-fst
fst<-fst%>%group_by(SNP)%>%mutate(MAF_ST=scale(MAF))
fst<-fst[,c(1,2,3,4,5,7)]
colnames(fst)<-c("CHR","SNP","CLST","A1","A2","MAF_ST")
fst<-as.data.frame(fst)
FST<-reshape2::dcast(fst, CHR+SNP+A1+A2 ~ CLST)

colnames(FST)<-c("chr","snp","a1","a2","T1","T2","T3")

DATA<-FST %>%
  mutate(Group="NULL") %>%
  mutate(Group=ifelse(T2-T1 > 0 & T3-T2 > 0,"Increase",Group)) %>%
  mutate(Group=ifelse(T2-T1 > 0 & T3-T2 > -0.5 & T3-T2 < 0.5,"Up_Flat",Group)) %>%
  mutate(Group=ifelse(T2-T1 > -0.5 & T2-T1 < 0.5 & T3-T2 > 0,"Flat_Up",Group)) %>%
  mutate(Group=ifelse(T2-T1 > 0.5 & T3-T2 < -0.5,"Up_Down",Group)) %>%
  mutate(Group=ifelse(T2-T1 < 0 & T3-T2 < 0,"Decrease",Group)) %>%
  mutate(Group=ifelse(T2-T1 < 0 & T3-T2 < 0.5 & T3-T2 > -0.5,"Down_Flat",Group)) %>%
  mutate(Group=ifelse(T2-T1 < 0.5 & T2-T1 > -0.5 & T3-T2 < 0,"Flat_Down",Group)) %>%
  mutate(Group=ifelse(T2-T1 < -0.5 & T3-T2 > 0.5,"Down_Up",Group)) %>%
  #filter(Group == "NULL")
  #filter(Group %in% c("Increase","Decrease")) %>%
  #filter(Group %in% c("Increase","Up_Down","Up_flat")) %>%
  filter(Group != "NULL")%>%
  mutate(Group=factor(Group, levels = c("Increase","Up_Flat","Flat_Up","Up_Down",
                                        "Decrease","Down_Flat","Flat_Down","Down_Up"))) %>%
  group_by(Group) %>%
  mutate(n=n()) %>%
  gather("Time","StdVarFreqChanges",5:7) %>%
  mutate(Time=factor(Time, levels = c("T1","T2","T3")))%>%ungroup()

DATA_plot<-as.data.table(DATA)
DATA_plot$Time <- as.character(DATA_plot$Time)
strip <- strip_themed(background_x = elem_list_rect(fill = rep("#F9D662", 8)))

if (all(c("03ADUR", "13ADUR", "21ADUR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 13,
    DATA_plot$Time == "T3", 21
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.4) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 13)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.01) +
    geom_line(data = DATA_plot[Time %in% c(13, 21)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.01) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    xlab("Years") +
    facet_wrap2(~ Group, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum.png", sep = "."), width = 7, height = 4, dpi = 400)
} else if (all(c("03WLR", "12WLR", "18WLR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 12,
    DATA_plot$Time == "T3", 18
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.5) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 12)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.01) +
    geom_line(data = DATA_plot[Time %in% c(12, 18)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.01) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    xlab("Years") +
    facet_wrap2(~ Group, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum.png", sep = "."), width = 7, height = 4, dpi = 400)
}


colnames(cp)<-c("chr","snp","Time","a1","a2","maf")
cp$Time[cp$Time=="03ADUR"]<-"T1"
cp$Time[cp$Time=="13ADUR"]<-"T2"
cp$Time[cp$Time=="21ADUR"]<-"T3"
cp$Time[cp$Time=="03WLR"]<-"T1"
cp$Time[cp$Time=="12WLR"]<-"T2"
cp$Time[cp$Time=="18WLR"]<-"T3"
DATA<-left_join(DATA,cp)
#write.table(DATA, file=paste(pre, "SelectionPrefernece_sum.txt", sep="."),sep = "\t",quote = F, row.names = F)

maf_select<-DATA[,c(1,2,3,4,5,6,7,9)]
maf_select<-reshape2::dcast(maf_select, chr+snp+a1+a2+Group+n ~ Time)
#maf_select$freq_changes<-maf_select$T2 - maf_select$T1
#sum(abs(maf_select$freq_changes)<=0.264226)
T1_T2_fst01<-fread(T1_T2_fst01,header = F)
T1_T2_fst01<-maf_select[maf_select$snp %in%T1_T2_fst01$V1,]
T1_T2_fst01$freq_changes<-abs(T1_T2_fst01$T2- T1_T2_fst01$T1)
T1_T2_fst01$group<-"T1_T2"
#hist(T1_T2_fst01$freq_changes)
T2_T3_fst01<-fread(T2_T3_fst01,header = F)
T2_T3_fst01<-maf_select[maf_select$snp%in%T2_T3_fst01$V1,]
T2_T3_fst01$freq_changes<-abs(T2_T3_fst01$T3- T2_T3_fst01$T2)
T2_T3_fst01$group<-"T2_T3"
#hist(T2_T3_fst01$freq_changes)
hist_data<-rbind(T1_T2_fst01,T2_T3_fst01)
hist_data$freq_changes<-abs(hist_data$freq_changes)
ggplot(data = hist_data,aes(x=freq_changes))+
  geom_density(aes(fill=group,color=group),alpha=0.3)+
  scale_fill_manual(values = c("#00afbb","#fc4e07"))+
  scale_color_manual(values = c("#00afbb","#fc4e07"))+
  xlab("Frequency changes")+ylab("Density")+
  theme_prism(base_size = 14)+
  guides(fill=guide_legend(override.aes = list(linewidth=1)),color=guide_legend(override.aes = list(linewidth=1)))+
  theme(legend.text  = element_text(face = "bold",size = 14),
        axis.title = element_text(face = "bold"),
        legend.title=element_blank(),
        legend.position = c(0.8, 0.8))
ggsave(paste(pre, "SelectionPrefernece_delta_AF_density.png", sep = "."), width =4, height = 4, dpi = 400)
T1_T2_top10_thread<-quantile(T1_T2_fst01$freq_changes,0.99)
T2_T3_top10_thread<-quantile(T2_T3_fst01$freq_changes,0.99)
T1_T2_top10<-T1_T2_fst01[T1_T2_fst01$freq_changes>=T1_T2_top10_thread,]
T2_T3_top10<-T2_T3_fst01[T2_T3_fst01$freq_changes>=T2_T3_top10_thread,]
T2_T1_top10_avg_delta_AF<-sum(abs(T1_T2_top10$T2 - T1_T2_top10$T1))/nrow(T1_T2_top10)
T3_T2_top10_avg_delta_AF<-sum(abs(T2_T3_top10$T3 - T2_T3_top10$T2))/nrow(T2_T3_top10)


T2_T1_avg_delta_AF<-sum(abs(T1_T2_fst01$T2 - T1_T2_fst01$T1))/nrow(T1_T2_fst01)
T3_T2_avg_delta_AF<-sum(abs(T2_T3_fst01$T3 - T2_T3_fst01$T2))/nrow(T2_T3_fst01)

df<-tribble(
  ~Gen, ~total_snp,~total_AF,~avg_delta_AF, ~top10_snp, ~top10_AF, ~top10_avg_delta_AF, ~max_AF,
  "T2_T1",nrow(T1_T2_fst01),sum(abs(T1_T2_fst01$T2 - T1_T2_fst01$T1)),T2_T1_avg_delta_AF,nrow(T1_T2_top10),sum(abs(T1_T2_top10$T2 - T1_T2_top10$T1)),T2_T1_top10_avg_delta_AF,max(T1_T2_fst01$freq_changes),
  "T3_T2",nrow(T2_T3_fst01),sum(abs(T2_T3_fst01$T3 - T2_T3_fst01$T2)),T3_T2_avg_delta_AF,nrow(T2_T3_top10),sum(abs(T2_T3_top10$T3 - T2_T3_top10$T2)),T3_T2_top10_avg_delta_AF,max(T2_T3_fst01$freq_changes),
)
write.table(df, file=paste(pre, "SelectionPrefernece_avg_delta_AF.txt", sep="."),sep = "\t",quote = F, row.names = F)
###########################
maf_select$Group<-as.character(maf_select$Group)
upflat<-maf_select[maf_select$Group%in%c("Up_Flat"),]
upflat_fix<-upflat[upflat$T2>=0.95 & upflat$T3>=0.95,]

downflat<-maf_select[maf_select$Group%in%c("Down_Flat"),]
downflat_fix<-downflat[downflat$T2<=0.05 & downflat$T3<=0.05,]

##flatup & flat down T1 and T2 frequency fixed may realted to introduction
#flatup<-maf_select[maf_select$Group%in%c("Flat_Up"),]
#flatup_fix<-flatup[flatup$T1<=0.05 & flatup$T2<=0.05,]

#flatdown<-maf_select[maf_select$Group%in%c("Flat_Down"),]
#flatdown_fix<-flatdown[flatdown$T1>=0.95 & flatdown$T2>=0.95,]

maf_select$deltaAF1<-abs(maf_select$T2 - maf_select$T1)
maf_select$deltaAF2<-abs(maf_select$T3 - maf_select$T2)
intro1<-maf_select$snp[maf_select$deltaAF1 > deltaAF1]
intro2<-maf_select$snp[maf_select$deltaAF2 > deltaAF2]
drift1<-maf_select$snp[maf_select$deltaAF1 < 0.1]
drift2<-maf_select$snp[maf_select$deltaAF2 < 0.1]
dirft_intersect<-intersect(drift1,drift2)
intro<-c(intro1,intro2,dirft_intersect)
##########################
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: start filter introduction \n")
DATA$Group_changes<-DATA$Group
DATA$Group_changes<-as.character(DATA$Group_changes)
DATA$Group_changes[DATA$snp%in%upflat_fix$snp]<-"Increase"
DATA$Group_changes[DATA$snp%in%downflat_fix$snp]<-"Decrease"
#DATA$Group_changes[DATA$snp%in%flatup_fix$snp]<-"Introduction"
#DATA$Group_changes[DATA$snp%in%flatdown_fix$snp]<-"Introduction"
DATA$Group_changes[DATA$snp%in%intro]<-"Introduction"
DATA$n[DATA$Group_changes=="Introduction"]<-sum(DATA$Group_changes=="Introduction")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Introduction loci:", sum(DATA$Group_changes=="Introduction")/3,"\n")
DATA$n[DATA$Group_changes=="Down_Up"]<-sum(DATA$Group_changes=="Down_Up")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Down_Up loci:", sum(DATA$Group_changes=="Down_Up")/3,"\n")
DATA$n[DATA$Group_changes=="Decrease"]<-sum(DATA$Group_changes=="Decrease")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Decrease loci:", sum(DATA$Group_changes=="Decrease")/3,"\n")
DATA$n[DATA$Group_changes=="Increase"]<-sum(DATA$Group_changes=="Increase")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Increase loci:", sum(DATA$Group_changes=="Increase")/3,"\n")
DATA$n[DATA$Group_changes=="Up_Flat"]<-sum(DATA$Group_changes=="Up_Flat")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Up_Flat loci:", sum(DATA$Group_changes=="Up_Flat")/3,"\n")
DATA$n[DATA$Group_changes=="Up_Down"]<-sum(DATA$Group_changes=="Up_Down")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Up_Down loci:", sum(DATA$Group_changes=="Up_Down")/3,"\n")
DATA$n[DATA$Group_changes=="Down_Flat"]<-sum(DATA$Group_changes=="Down_Flat")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Down_Flat loci:", sum(DATA$Group_changes=="Down_Flat")/3,"\n")
DATA$n[DATA$Group_changes=="Flat_Up"]<-sum(DATA$Group_changes=="Flat_Up")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Flat_Up loci:", sum(DATA$Group_changes=="Flat_Up")/3,"\n")
DATA$n[DATA$Group_changes=="Flat_Down"]<-sum(DATA$Group_changes=="Flat_Down")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Flat_Down loci:", sum(DATA$Group_changes=="Flat_Down")/3,"\n")
write.table(DATA%>%arrange(snp), file=paste(pre, "SelectionPrefernece_sum.txt", sep="."),sep = "\t",quote = F, row.names = F)
intro_per<-((sum(DATA$Group_changes=="Introduction")/3)/(nrow(DATA)/3))*100
names(intro_per)<-"%"
intro_per<-as.data.frame(intro_per)
write.table(intro_per, file=paste(pre, "SelectionPrefernece_intro_per.txt", sep="."),sep = "\t",quote = F, row.names = T)

############################################################################################################################################
DATA_plot<-as.data.table(DATA)
DATA_plot<-DATA_plot[DATA_plot$Group_changes!="Introduction",]

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: start delta_AF_density plot and delta_AF after filter introduction loci \n")
cp<-DATA_plot[,c("chr","snp","a1","a2","Time","Group_changes","n","maf")]
cp<-reshape2::dcast(cp, chr+snp+a1+a2+Group_changes+n ~ Time)
T1_T2_fst01<-cp[cp$snp %in% T1_T2_fst01$snp,]
T1_T2_fst01$freq_changes<-abs(T1_T2_fst01$T2- T1_T2_fst01$T1)
T1_T2_fst01$group<-"T1_T2"
#hist(T1_T2_fst01$freq_changes)
T2_T3_fst01<-cp[cp$snp%in%T2_T3_fst01$snp,]
T2_T3_fst01$freq_changes<-abs(T2_T3_fst01$T3- T2_T3_fst01$T2)
T2_T3_fst01$group<-"T2_T3"
#hist(T2_T3_fst01$freq_changes)
hist_data<-rbind(T1_T2_fst01,T2_T3_fst01)
hist_data$freq_changes<-abs(hist_data$freq_changes)
ggplot(data = hist_data,aes(x=freq_changes))+
  geom_density(aes(fill=group,color=group),alpha=0.3)+
  scale_fill_manual(values = c("#00afbb","#fc4e07"))+
  scale_color_manual(values = c("#00afbb","#fc4e07"))+
  xlab("Frequency changes")+ylab("Density")+
  theme_prism(base_size = 14)+
  guides(fill=guide_legend(override.aes = list(linewidth=1)),color=guide_legend(override.aes = list(linewidth=1)))+
  theme(legend.text  = element_text(face = "bold",size = 14),
        axis.title = element_text(face = "bold"),
        legend.title=element_blank(),
        legend.position = c(0.8, 0.8))
ggsave(paste(pre, "SelectionPrefernece_delta_AF_density_filter.png", sep = "."), width =4, height = 4, dpi = 400)
T1_T2_top10_thread<-quantile(T1_T2_fst01$freq_changes,0.99)
T2_T3_top10_thread<-quantile(T2_T3_fst01$freq_changes,0.99)
T1_T2_top10<-T1_T2_fst01[T1_T2_fst01$freq_changes>=T1_T2_top10_thread,]
T2_T3_top10<-T2_T3_fst01[T2_T3_fst01$freq_changes>=T2_T3_top10_thread,]
T2_T1_top10_avg_delta_AF<-sum(abs(T1_T2_top10$T2 - T1_T2_top10$T1))/nrow(T1_T2_top10)
T3_T2_top10_avg_delta_AF<-sum(abs(T2_T3_top10$T3 - T2_T3_top10$T2))/nrow(T2_T3_top10)


T2_T1_avg_delta_AF<-sum(abs(T1_T2_fst01$T2 - T1_T2_fst01$T1))/nrow(T1_T2_fst01)
T3_T2_avg_delta_AF<-sum(abs(T2_T3_fst01$T3 - T2_T3_fst01$T2))/nrow(T2_T3_fst01)

df<-tribble(
  ~Gen, ~total_snp,~total_AF,~avg_delta_AF, ~top10_snp, ~top10_AF, ~top10_avg_delta_AF, ~max_AF,
  "T2_T1",nrow(T1_T2_fst01),sum(abs(T1_T2_fst01$T2 - T1_T2_fst01$T1)),T2_T1_avg_delta_AF,nrow(T1_T2_top10),sum(abs(T1_T2_top10$T2 - T1_T2_top10$T1)),T2_T1_top10_avg_delta_AF,max(T1_T2_fst01$freq_changes),
  "T3_T2",nrow(T2_T3_fst01),sum(abs(T2_T3_fst01$T3 - T2_T3_fst01$T2)),T3_T2_avg_delta_AF,nrow(T2_T3_top10),sum(abs(T2_T3_top10$T3 - T2_T3_top10$T2)),T3_T2_top10_avg_delta_AF,max(T2_T3_fst01$freq_changes),
)
write.table(df, file=paste(pre, "SelectionPrefernece_avg_delta_AF_filter.txt", sep="."),sep = "\t",quote = F, row.names = F)

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: start plot SelectionPrefernece after filter introduction loci \n")
DATA_plot$Time <- as.character(DATA_plot$Time)
strip <- strip_themed(background_x = elem_list_rect(fill = rep("#F9D662", 8)))

DATA_plot$Group_changes<-factor(DATA_plot$Group_changes, levels = c("Increase","Up_Flat","Flat_Up","Up_Down",
                                                                    "Decrease","Down_Flat","Flat_Down","Down_Up"))
if (all(c("03ADUR", "13ADUR", "21ADUR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 13,
    DATA_plot$Time == "T3", 21
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.4) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 13)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.01) +
    geom_line(data = DATA_plot[Time %in% c(13, 21)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.01) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    xlab("Years") +
    facet_wrap2(~ Group_changes, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum_filterIntro.png", sep = "."), width = 7, height = 4, dpi = 400)
} else if (all(c("03WLR", "12WLR", "18WLR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 12,
    DATA_plot$Time == "T3", 18
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.5) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 12)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.01) +
    geom_line(data = DATA_plot[Time %in% c(12, 18)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.01) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    xlab("Years") +
    facet_wrap2(~ Group_changes, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum_filterIntro.png", sep = "."), width = 7, height = 4, dpi = 400)
}

######
DATA_plot<-as.data.table(DATA)
DATA_plot<-DATA_plot[DATA_plot$Group_changes=="Introduction",]
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: start plot SelectionPrefernece for introduction loci \n")
DATA_plot$Time <- as.character(DATA_plot$Time)
strip <- strip_themed(background_x = elem_list_rect(fill = rep("#F9D662", 8)))

DATA_plot$n[DATA_plot$Group=="Down_Up"]<-sum(DATA_plot$Group=="Down_Up")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Down_Up loci:", sum(DATA_plot$Group=="Down_Up")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Decrease"]<-sum(DATA_plot$Group=="Decrease")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Decrease loci:", sum(DATA_plot$Group=="Decrease")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Increase"]<-sum(DATA_plot$Group=="Increase")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Increase loci:", sum(DATA_plot$Group=="Increase")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Up_Flat"]<-sum(DATA_plot$Group=="Up_Flat")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Up_Flat loci:", sum(DATA_plot$Group=="Up_Flat")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Up_Down"]<-sum(DATA_plot$Group=="Up_Down")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Up_Down loci:", sum(DATA_plot$Group=="Up_Down")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Down_Flat"]<-sum(DATA_plot$Group=="Down_Flat")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Down_Flat loci:", sum(DATA_plot$Group=="Down_Flat")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Flat_Up"]<-sum(DATA_plot$Group=="Flat_Up")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Flat_Up loci:", sum(DATA_plot$Group=="Flat_Up")/3,"\n")
DATA_plot$n[DATA_plot$Group=="Flat_Down"]<-sum(DATA_plot$Group=="Flat_Down")/3
cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: total number of Flat_Down loci:", sum(DATA_plot$Group=="Flat_Down")/3,"\n")

na<-c("Increase","Up_Flat","Flat_Up","Up_Down",
      "Decrease","Down_Flat","Flat_Down","Down_Up")[!c("Increase","Up_Flat","Flat_Up","Up_Down",
                                                       "Decrease","Down_Flat","Flat_Down","Down_Up")%in%DATA_plot$Group]

for (i in na) {
  tmp<-tribble(
    ~chr, ~snp,~a1,~a2,~Group,~n,~Time,~StdVarFreqChanges,~maf,~Group_changes,
    NA,NA,NA,NA,i,0,"T1",NA,NA,"Introduction",
    NA,NA,NA,NA,i,0,"T2",NA,NA,"Introduction",
    NA,NA,NA,NA,i,0,"T3",NA,NA,"Introduction",
  )
  DATA_plot<-rbind(tmp,DATA_plot)
}

DATA_plot$Group<-factor(DATA_plot$Group, levels = c("Increase","Up_Flat","Flat_Up","Up_Down",
                                                    "Decrease","Down_Flat","Flat_Down","Down_Up"))
DATA_plot<-as.data.table(DATA_plot)
if (all(c("03ADUR", "13ADUR", "21ADUR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 13,
    DATA_plot$Time == "T3", 21
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.4) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 13)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(13, 21)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.05) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    xlab("Years") +
    facet_wrap2(~ Group, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum_forIntro.png", sep = "."), width = 7, height = 4, dpi = 400)
} else if (all(c("03WLR", "12WLR", "18WLR") %in% unique(fst$CLST))) {
  DATA_plot$Time <- fcase(
    DATA_plot$Time == "T1", 3,
    DATA_plot$Time == "T2", 12,
    DATA_plot$Time == "T3", 18
  )
  ggplot(DATA_plot) +
    #geom_violin(aes(x = Time, y = StdVarFreqChanges, fill = as.factor(Time)), trim = FALSE, width = 1.5) +
    geom_point(aes(x = Time, y = StdVarFreqChanges), size = 0.05, alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(3, 12)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#00afbb", alpha = 0.05) +
    geom_line(data = DATA_plot[Time %in% c(12, 18)], aes(x = Time, y = StdVarFreqChanges, group = snp), color = "#fc4e07", alpha = 0.05) +
    geom_text(aes(x = median(Time), y = 1, label = n)) +
    #scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", "-0.5", "0", "0.5", "1.0")) +
    #scale_fill_manual(values = c("#00afbb", "white", "#fc4e07")) +
    xlab("Years") +
    facet_wrap2(~ Group, ncol = 4, strip = strip) +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 13, colour = "black"), axis.text.y = element_text(colour = "black", size = 13)) +
    guides(fill = "none")
  ggsave(paste(pre, "SelectionPrefernece_sum_forIntro.png", sep = "."), width = 7, height = 4, dpi = 400)
}

cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"~~ info: all done \n")