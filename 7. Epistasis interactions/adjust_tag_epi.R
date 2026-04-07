#!/usr/bin/env Rscript
# Script: adjust_tag_epi.R
# Date: 2024/12/31
# Purpose: adjust representative epistatic effect pairs 
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
             "* adjust p -- a rscript: adjust representative epistatic effect pairs",
             sprintf("* Version %s", .VERSION),
             "* Yibin Qiu (13422157044qyb@gmail.com)",
             "* South China Agricultural Univerisity",
             "*********************************************************************"
))

# Set up argument parser
p <- arg_parser("This script was used to pick representative epistatic effect pairs.", hide.opts = TRUE)
p <- add_argument(p, "--trait", help = "Input trait name", type = "character", default = "carcass_length")
p <- add_argument(p, "--epi_path", help = "Input epi_path", type = "character", default = "./")
p <- add_argument(p, "--pvalue", help = "Input filter pvalue", type = "character", default = "1e-6")
p <- add_argument(p, "--method", help = "fdr,bonferroni", type = "character", default = "bonferroni")
p <- add_argument(p, "--output", help = "Output path", type = "character", default = "./")

# Parse arguments
argv <- parse_args(p)
trait <- argv$trait
epi_path <- argv$epi_path
pvalue <- as.numeric(argv$pvalue)
method <- argv$method
output_path <- argv$output

# Load and filter epi data
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: loading represent epi table for ", trait,".....")
represent <- fread(file.path(epi_path, paste0(trait, "_representative_pair.txt")), header = T)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: before filter representative ", nrow(represent),".....")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: select method ", method,".....")
represent$lead_P_adjust<-p.adjust(represent$lead_P,method = method)
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: threadshold ", pvalue,".....")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: after filter remain", sum(represent$lead_P_adjust<=pvalue),".....")
remain<-represent[represent$lead_P_adjust<=pvalue,]
remain<-remain%>%arrange(lead_P_adjust,lead_P)
fwrite(remain,file.path(output_path, paste0(trait, "_representative_pair_adjusted.txt")),row.names = F,col.names = T,quote = F,sep = "\t")
remain<-remain%>%select(lead_SNP1,lead_SNP2)
fwrite(remain,file.path(output_path, paste0(trait, "_representative_pair_adjusted_pairname_forplot.txt")),row.names = F,col.names = F,quote = F,sep = "\t")
message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: all done.....")