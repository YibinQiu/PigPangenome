#!/usr/bin/env Rscript
# Script: adjust_phe.R
# Date: 2024/10/24
# Purpose: Remove outliers (mean ± 3*SD) and adjust phenotype with fixed effects.
# Author: Yibin Qiu
Sys.setenv(LANGUAGE = "en") # 设置提示信息为英文
# Load required libraries and suppress messages and warnings
suppressWarnings(suppressMessages({
  library(argparser)
  library(data.table)
  library(tidyverse)
}))
# constants
.VERSION = "1.0"
################################################################################
writeLines(c("*********************************************************************",
             "* adjust phe -- a rscript: adjust raw phenotype for epstasis GWAS",
             sprintf("* Version %s", .VERSION),
             "* Yibin Qiu (13422157044qyb@gmail.com)",
             "* South China Agricultural Univerisity",
             "*********************************************************************"
))
# 创建解析器并隐藏帮助信息
p <- arg_parser("This script was used to remove outliers (mean ± 3*SD) and remove fixed effects.", hide.opts = TRUE)

# 添加位置索引参数和 flag 参数
p <- add_argument(p, "--onlysigfactor", help = "Adjust only for significant effects: TRUE/FALSE", type = "character", default = "TRUE")
p <- add_argument(p, "--input_export", help = "Input path for export_list.txt", type = "character", default = "./")
p <- add_argument(p, "--input_pca", help = "Input path for PCA", type = "character", default = "./")
p <- add_argument(p, "--input_qcov", help = "Input path for qcov", type = "character", default = "./")
p <- add_argument(p, "--input_cov", help = "Input path for cov", type = "character", default = "./")
p <- add_argument(p, "--input_trait", help = "Input path for trait", type = "character", default = "./")
p <- add_argument(p, "--export_list", help = "Path to export_list.txt", type = "character", default = "export_list.txt")
p <- add_argument(p, "--pcafile", help = "PCA file created by gcta", type = "character", default = "S21_S22_S23_merge_GRM_PCA.eigenvec")
#p <- add_argument(p, "qcovfile", help = "Qcov file", type = "character", default = "qcov_ceding_trait.txt")
#p <- add_argument(p, "covfile", help = "Cov file", type = "character", default = "cov_ceding_trait.txt")
p <- add_argument(p, "--traitfile", help = "Trait file", type = "character", default = "ceding_trait.txt")
p <- add_argument(p, "--output", help = "Output path for adjusted phenotype", type = "character", default = "./")

# 解析命令行参数
argv <- parse_args(p)

onlysigfactor<-argv$onlysigfactor
onlysigfactor<-as.logical(onlysigfactor)

input_export<-argv$input_export
input_pca<-argv$input_pca
input_qcov<-argv$input_qcov
input_cov<-argv$input_cov
input_trait<-argv$input_trait
export_list<-argv$export_list
pcafile<-argv$pcafile
#qcovfile<-argv$qcovfile
#covfile<-argv$covfile
traitfile<-argv$traitfile
output<-argv$output

######################### 1.Loading input files #############################################
message("######################### 1.Loading input files ###############################################")

message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Loading input files...")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: reading export_list from [ ",paste0(input_export,"/",export_list)," ].")
export_list<-fread(paste0(input_export,"/",export_list),header = F)%>%as.data.frame()
colnames(export_list)<-c("traitfile","index","trait","covfile","qcovfile","qcovfile2")
export_list<-select(export_list,"traitfile","index","trait","covfile","qcovfile")
export_list$index<-export_list$index+2
export_list<-export_list[export_list$traitfile==traitfile,]
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",nrow(export_list)," traits are included from export_list.")

message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: reading pcafile from [ ",paste0(input_pca,"/",pcafile)," ].")
pcafile<-fread(paste0(input_pca,"/",pcafile),header = F)%>%as.data.frame()
colnames(pcafile)<-c("fid","iid",paste0("pc",1:(ncol(pcafile)-2)))
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",nrow(pcafile)," individuals are included from pcafile.")

qcovfile<-unique(export_list$qcovfile)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: reading qcovfile from [ ",paste0(input_qcov,"/",qcovfile)," ].")
qcovfile<-fread(paste0(input_qcov,"/",qcovfile),header = F)%>%as.data.frame()
if (ncol(qcovfile)==2) {
  colnames(qcovfile)<-c("fid","iid")
} else {
  colnames(qcovfile)<-c("fid","iid",paste0("qcov",1:(ncol(qcovfile)-2)))
}
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",nrow(qcovfile)," individuals are included from qcovfile.")

covfile<-unique(export_list$covfile)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: reading covfile from [ ",paste0(input_cov,"/",covfile)," ].")
covfile<-fread(paste0(input_cov,"/",covfile),header = F)%>%as.data.frame()
if (ncol(covfile)==2) {
  colnames(covfile)<-c("fid","iid")
} else {
  colnames(covfile)<-c("fid","iid",paste0("cov",1:(ncol(covfile)-2)))
}
covfile<-covfile%>%mutate(across(-c(1,2),as.factor)) #actually, yaozekai do yc in hiblup, found high correaltion (r=0.98) between yc and adjusted phe in R (no as.factor, let cov as qcov)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",nrow(covfile)," individuals are included from covfile.")

message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: reading traitfile from [ ",paste0(input_trait,"/",traitfile)," ].")
traitfile<-fread(paste0(input_trait,"/",traitfile),header = F)%>%as.data.frame()
num_cols <- ncol(traitfile)
col_types <- c(rep("character", 2), rep("character", num_cols - 2))
for (col in 3:ncol(traitfile)) {
  traitfile[[col]] <- as.numeric(traitfile[[col]])
}
traitfile<-traitfile[,c(1,2,export_list$index)]
colnames(traitfile)<-c("fid","iid",export_list$trait)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",nrow(traitfile)," individuals are included from traitfile.")

commonSample<-intersect(qcovfile$iid,covfile$iid)%>%intersect(traitfile$iid)%>%intersect(pcafile$iid)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",length(commonSample)," individuals are in common in these files.")

message("######################### 1.Loading input files done ##########################################")
######################### 1.Loading input files #############################################

######################### 2.remove outliers #################################################
message("#########################  2.remove outliers ##################################################")

message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: remove outliers...")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: define outliers function (mean ± 3*SD) ...")
filter_outliers_df <- function(df) {
  # 确保输入是数据框
  if (!is.data.frame(df)) {
    stop("输入必须是数据框")
  }
  
  # 遍历数据框的每一列（从第三列开始）
  for (col in 3:ncol(df)) {
    # 获取该列数据
    column_data <- df[[col]]
    
    # 检查是否为数值型
    if (is.numeric(column_data)) {
      # 计算平均值和标准差
      mean_value <- mean(column_data, na.rm = TRUE)
      sd_value <- sd(column_data, na.rm = TRUE)
      
      # 计算上下界限
      lower_bound <- mean_value - 3 * sd_value
      upper_bound <- mean_value + 3 * sd_value
      
      # 将超出范围的值替换为 NA
      df[[col]][column_data < lower_bound | column_data > upper_bound] <- NA
    } else {
      warning(paste("列", colnames(df)[col], "不是数值型，已跳过"))
    }
  }
  
  return(df)
}
traitfile<-filter_outliers_df(traitfile)
traitfile<-traitfile[traitfile$iid%in%commonSample,]
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: ",length(commonSample)," individuals will used for furthur analysis.")

message("#########################  2.remove outliers done #############################################")
######################### 2.remove outliers #################################################

######################### 3.fit models #################################################
message("#########################  3.fit models #######################################################")

pheno<-suppressMessages(left_join(traitfile,pcafile)%>%left_join(qcovfile)%>%left_join(covfile))
fix_effext<-pheno[,!colnames(pheno)%in%c("fid","iid",export_list$trait)]
ids<-pheno[,c("fid","iid")]
for (trait in export_list$trait) {
  
  #trait <- "total_teat_number"
  message("###############################################################################")
  message("###############################################################################")
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: asjust phenotype ---- ",trait,".")
  
  df_trait <- as.data.frame(cbind(as.data.frame(pheno[,trait]),fix_effext))
  colnames(df_trait)<-c("rawphe",colnames(df_trait)[-1])
  df_original <- cbind(ids,df_trait[,c("rawphe")])
  fwrite(df_original,paste0(output,"/",trait,"_original.txt"),sep = "\t",na="NA",row.names = F,col.names = F,quote = F)
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: output original for phenotype ---- ",trait,".")
  # 拟合线性模型
  fit <- lm(rawphe ~ . , data = df_trait)
  
  # 提取系数
  coef <- as.data.frame(summary(fit)$coef)
  fwrite(coef,paste0(output,"/",trait,"_fitmodel_orignal.txt"),sep = ",",na="NA",quote = F,row.names = T,col.names = T)
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: coef for phenotype ---- ",trait,".")
  print(coef)
  fix_factor <- rownames(coef)
  # 是否只用显著的因子做矫正
  if( onlysigfactor ) {
    coef<-coef[coef$`Pr(>|t|)` <= 0.05,]
    message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: select coef p value <0.05.")
    fix_factor <- rownames(coef)
    # 如果除了截距项没有显著的因子，则直接输出原始表型
    if( nrow(coef[fix_factor!="(Intercept)",])==0 ){
      df_adjusted <- cbind(ids,as.data.frame(pheno[,trait]))
      fwrite(df_adjusted,paste0(output,"/",trait,"_adjusted.txt"),sep = "\t",na="NA",row.names = F,col.names = F,quote = F)
      warning(trait," has no significant fix effect, will return original phenotype value (removed outliers).")
      next
    }
  }
  remian_fix<-fix_factor[fix_factor!="(Intercept)"]
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: remian fix effects ~ ",remian_fix ," for phenotype ---- ",trait,".")
  df_trait<-df_trait[,c("rawphe",remian_fix)]
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: re fit for phenotype ---- ",trait,".")
  fit <- lm(rawphe ~ . , data = df_trait)
  # 提取系数
  coef <- as.data.frame(summary(fit)$coef)
  fwrite(coef,paste0(output,"/",trait,"_fitmodel_refit.txt"),sep = ",",na="NA",quote = F,row.names = T,col.names = T)
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: re coef for phenotype ---- ",trait,".")
  print(coef)
  fix_factor <- rownames(coef)
  # 提取非截距项因子
  est <- coef$Estimate[fix_factor!="(Intercept)"]
  #intercept <- 0
  
  # 提取截距项
  #if( "(Intercept)" %in% fix_factor ) {
  #  intercept<- coef$Estimate[fix_factor=="(Intercept)"]
  #}
  
  # 提取因子矩阵
  fix_matrix<-df_trait%>%select(any_of(fix_factor))
  
  # 计算剩余表型值(暂时用不减截距的版本)
  df_residual <- as.matrix(df_trait[, colnames(df_trait)=="rawphe"]) - as.matrix(fix_matrix) %*% as.matrix(c(est))
  # 计算剩余表型值
  #df_residual <- as.matrix(df_trait[, colnames(df_trait)=="rawphe"]) - as.matrix(fix_matrix) %*% c(est) - intercept
  
  colnames(df_residual) <- trait
  df_adjusted <- cbind(ids,df_residual)
  fwrite(df_adjusted,paste0(output,"/",trait,"_adjusted.txt"),sep = "\t",na="NA",row.names = F,col.names = F,quote = F)
  message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: adjust done for phenotype ---- ",trait,".")
}
message("#########################  3.fit models done ##################################################")
######################### 3.fit models #################################################

message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: all done")
