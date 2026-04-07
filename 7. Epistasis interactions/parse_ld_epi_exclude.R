#!/usr/bin/env Rscript
# Script: parse_ld_epi_exclude.R
# Date: 2025/01/06
# Purpose: parse 1 MB r2 > 0.1 epi pairs and sort with pvalue and beta, and excluded loci where any of the three genotypes had fewer than 10 individuals.
# Author: Yibin Qiu
Sys.setenv(LANGUAGE = "en") # 设置提示信息为英文
# Load required libraries and suppress messages and warnings
suppressWarnings(suppressMessages({
  library(argparser)
  library(data.table)
  library(tidyverse)
}))
setDTthreads(threads=24)
# constants
.VERSION = "1.0"
################################################################################
writeLines(c("*********************************************************************",
             "* merge r2 and gwas P value -- a rscript: merge R2 and GWAS result",
             sprintf("* Version %s", .VERSION),
             "* Yibin Qiu (13422157044qyb@gmail.com)",
             "* South China Agricultural Univerisity",
             "*********************************************************************"
))
# 创建解析器并隐藏帮助信息
p <- arg_parser("This script was used to merge r2 > 0.1 epi pairs (1MB updownstream) and sort with epi1 pvalue and beta.", hide.opts = TRUE)

# 添加位置索引参数和 flag 参数
p <- add_argument(p, "--trait", help = "trait", type = "character", default = "carcass_length")
p <- add_argument(p, "--gwas", help = "gwas path", type = "character", default = "carcass_length")
p <- add_argument(p, "--epi", help = "epi path", type = "character", default = "carcass_length")
p <- add_argument(p, "--ld", help = "ld path", type = "character", default = "carcass_length")
p <- add_argument(p, "--include", help = "include ID file", type = "character", default = "./merge_SVtransMissing_retain_ID.txt")
p <- add_argument(p, "--output", help = "Output path", type = "character", default = "./")

# 解析命令行参数
argv <- parse_args(p)

trait<-argv$trait
gwas<-argv$gwas
epi<-argv$epi
ld<-argv$ld
include<-argv$include
output<-argv$output
######################### testing #############################################
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: read gwas for phenotype ---- ",trait,".")
gwas <- fread(paste0(gwas,"/",trait,"_GWAS.pval"),header=T)%>%as.data.frame%>%select(rs,P1df)
colnames(gwas) <- c("SNP1","P1df")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: read epi for phenotype ---- ",trait,".")
epi <- fread(paste0(epi,"/",trait,".epi.qt"),header=T)%>%select(SNP1,SNP2,BETA_INT,P)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: read ld for phenotype ---- ",trait,".")
ld <- fread(paste0(ld,"/pval_0.001_ld.ld"),header=T)%>%select(SNP_A,SNP_B,R2)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: read include ID for phenotype ---- ",trait,".")
include <- fread(include,header=F)
colnames(include)<-c("ID")
file_a_reversed <- copy(epi)
setnames(file_a_reversed, c("SNP1", "SNP2"), c("SNP2", "SNP1"))

# 合并两个文件，找到正向匹配的 R2 值
merged_data <- merge(epi, ld, by.x = c("SNP1", "SNP2"), by.y = c("SNP_A", "SNP_B"), all.x = TRUE)

# 合并文件，找到反向匹配的 R2 值
merged_data_reversed <- merge(file_a_reversed, ld, by.x = c("SNP1", "SNP2"), by.y = c("SNP_A", "SNP_B"), all.x = TRUE)

# 将反向匹配的 R2 值合并到原始数据中
merged_data[, R2_reversed := merged_data_reversed$R2]
merged_data[, R2 := fifelse(is.na(R2), R2_reversed, R2)]
merged_data[, R2_reversed := NULL]


merged_data <-as.data.frame(merged_data)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: left_join all and arrange for phenotype ---- ",trait,".")
merged_data <- left_join(merged_data,gwas)%>%arrange(P,desc(abs(BETA_INT)),P1df)
print(head(merged_data))
#                SNP1               SNP2 BETA_INT         P R2      P1df
#1 chr17:40232721_G_T chr17:41281017_G_A -2.98154 2.657e-20 NA 3.357e-04
#2 chr17:27959133_A_T chr17:28883215_C_T  3.55108 2.917e-20 NA 1.108e-06
#3 chr17:28883215_C_T chr17:27959133_A_T  3.55108 2.917e-20 NA 9.843e-05
#4 chr17:27959133_A_T chr17:28885052_A_G  3.58874 5.273e-20 NA 1.108e-06
#5 chr17:28885052_A_G chr17:27959133_A_T  3.58874 5.273e-20 NA 1.251e-04
#6 chr17:27959133_A_T chr17:28884422_C_T  3.54955 8.766e-20 NA 1.108e-06
merged_data<-as.data.frame(merged_data)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: exclude variants for phenotype ---- ",trait,".")
merged_data<-merged_data[merged_data$SNP1%in%include$ID,]
merged_data<-merged_data[merged_data$SNP2%in%include$ID,]

fwrite(merged_data,paste0(output,"/",trait,".epi_r2_p1df.txt"),row.names=F,col.names=T,sep="\t",quote=F,na="NA")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: done for phenotype ---- ",trait,".")