#!/usr/bin/env Rscript
# Script: svmu2vcf.R
# Date: 2025/04/02
# Purpose: svmu2vcf.
# Author: Yibin Qiu
Sys.setenv(LANGUAGE = "en") # 设置提示信息为英文
# Load required libraries and suppress messages and warnings
suppressWarnings(suppressMessages({
  library(svmutools)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(Rsamtools)
  library(argparser)
  library(tidyverse)
}))
# constants
.VERSION = "1.0"
################################################################################
writeLines(c("*********************************************************************",
             "* svmu2vcf -- a rscript: svmu2vcf",
             sprintf("* Version %s", .VERSION),
             "* Yibin Qiu (13422157044qyb@gmail.com)",
             "* South China Agricultural Univerisity",
             "*********************************************************************"
))

# Set up argument parser
p <- arg_parser("This script was used to create svmu vcf.", hide.opts = TRUE)
p <- add_argument(p, "--sample_name", 
                  help = "The name of the sample. Will be used as the genotype column in the output VCF file and as a prefix for the VCF ID column.", 
                  type = "character", default = "Babraham")
p <- add_argument(p, "--svmu_file", 
                  help = "The path to the file containing the output of svmu (typically sv.txt)", 
                  type = "character", default = "sv.Babraham.txt")
p <- add_argument(p, "--output_file", 
                  help = "The VCF file to which output is directed.", 
                  type = "character", default = "Babraham.vcf")
p <- add_argument(p, "--min_size", 
                  help = "The minimum size of an SV for it to be kept for further processing", 
                  type = "integer", default = 30)
p <- add_argument(p, "--max_size", 
                  help = "The maximum size for an SV to be kept for further processing", 
                  type = "integer", default = 5000000)
#p <- add_argument(p, "--sv_types", help = "Output path", type = "character", default = "")
p <- add_argument(p, "--ref_fasta", 
                  help = "The path to the reference fasta file", 
                  type = "character", default = "W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa")
p <- add_argument(p, "--query_fasta", 
                  help = "Same as ref_fasta but for the query sequence", 
                  type = "character", default = "GCA_031225015.1_TPI_Babraham_pig_v1_genomic_simplifyID.fna")
p <- add_argument(p, "--cores", 
                  help = "The number of cores to use for extracting REF/ALT sequences in parallel. By default, cores = 1 which means that processing is not done in parallel.", 
                  type = "integer", default = 1)
p <- add_argument(p, "--logging", 
                  help = "Whether the user should be informed at various processing steps (TRUE by default).", 
                  type = "logical", default = TRUE)
# Parse arguments
argv <- parse_args(p)
sample_name <- argv$sample_name
svmu_file <- argv$svmu_file
output_file <- argv$output_file
min_size <- argv$min_size
max_size <- argv$max_size
#sv_types <- argv$sv_types
ref_fasta <- argv$ref_fasta
query_fasta <- argv$query_fasta
cores <- argv$cores
logging <- argv$logging


#setwd(dir = "/Users/yibinqiu/Documents/pangenome_20231110/进展/workflow/pangenome2/svmu/")

####redefine read_svmu function############
read_svmu <- function(filename, remove_duplicates = TRUE) {
  svmu_data <- read.table(filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE,na.strings = "")
  svmu_data <- svmu_data[!is.na(svmu_data$REF_CHROM),]
  # Removing the duplicates if remove_duplicates
  if(remove_duplicates) {
    svmu_data <- svmu_data[!duplicated(svmu_data[, c("REF_CHROM", "REF_START", "REF_END", "SV_TYPE", "ID", "LEN")]), ]
  }
  
  svmu_data
}

####replaced inverted_regions function############
inverted_regions <- function(filename, min_distance) {
  # Reading the data in and keeping only inversions
  svmu_data <- read_svmu(filename)
  inversions <- svmu_data[svmu_data$SV_TYPE == "INV", ]
  
  # Converting it to a GenomicRanges object
  inversions <- GenomicRanges::makeGRangesFromDataFrame(inversions,
                                                        ignore.strand = TRUE,
                                                        seqnames.field = "REF_CHROM",
                                                        start.field = "REF_START",
                                                        end.field = "REF_END")
  
  # Reducing the ranges by merging them if they are closer than min_distance
  inversions <- GenomicRanges::reduce(inversions, min.gapwidth = min_distance)
  
  inversions
}
###########################

# Read svmu file
svmu<-read_svmu(filename = svmu_file,remove_duplicates = T)
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Reading ", nrow(svmu), " svmu records.")

# filter svmu file
svmu_filter<-filter_svmu(svmu_data = svmu,
                         min_size = min_size,
                         max_size = max_size,
                         sv_types = c("DEL","INS"))
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: After filtering, ", nrow(svmu_filter), " svmu records remianed.")

# find INV ranges
inv_range<-inverted_regions(filename = svmu_file,min_distance = 1L)
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Grep INV records.")

# flag SVs
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Flag DELs or INSs in INVs.")
svmu_filter_inv<-flag_inversions(svmu_data = svmu_filter,inversions = inv_range,maxgap = 0)

if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: only keep autosome ref DELs or INSs.")
svmu_filter_inv<-svmu_filter_inv[svmu_filter_inv$REF_CHROM%in%paste0("chr",seq(1,18)),]

if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: remove DELs or INSs ranges larger than ref or query.")
# for (i in 1:nrow(svmu_filter_inv)) {
#   cat(i,"\n")
#   svmu_filter_addseq<-add_sequences(svmu_data = svmu_filter_inv[i,],
#                                     ref_fasta = ref_fasta,
#                                     query_fasta = query_fasta)
# }

ref_range<-as.data.frame(ranges(Rsamtools::scanFaIndex(ref_fasta)))
colnames(ref_range)<-c("ref_start","ref_end","ref_width")
query_range<-as.data.frame(ranges(Rsamtools::scanFaIndex(query_fasta)))
colnames(query_range)<-c("query_start","query_end","query_width")

ref_seq<-as.data.frame(seqnames(Rsamtools::scanFaIndex(ref_fasta)))
colnames(ref_seq)<-c("REF_CHROM")
ref_seq$REF_CHROM<-as.character(ref_seq$REF_CHROM)
query_seq<-as.data.frame(seqnames(Rsamtools::scanFaIndex(query_fasta)))
colnames(query_seq)<-c("Q_CHROM")
query_seq$Q_CHROM<-as.character(query_seq$Q_CHROM)

ref_seq<-cbind(ref_seq,ref_range)
query_seq<-cbind(query_seq,query_range)
test<-svmu_filter_inv
test<-dplyr::left_join(test,ref_seq)%>%dplyr::left_join(query_seq)

test$refture<-test$REF_START<=test$ref_end & test$REF_END <= test$ref_end
test$queryture<-test$Q_START<=test$query_end & test$Q_END <= test$query_end
test$conbimedture<- test$refture & test$queryture

svmu_filter_inv<-svmu_filter_inv[test$conbimedture,]

# add sequences for SVs (sequences for flag SVs will reverse)
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Add sequences for DELs or INSs.")
svmu_filter_addseq<-add_sequences(svmu_data = svmu_filter_inv,
                                  ref_fasta = ref_fasta,
                                  query_fasta = query_fasta)

# dataframe to vcf body
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Format to vcf body.")
vcf_records<-format_vcf(svmu_data = svmu_filter_addseq,sample_name = sample_name,logging = logging)

# reoder vcf body
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Reoder vcf body by ref_fasta.")
fai_seqnames <- as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(ref_fasta)))
vcf_records$CHROM <- factor(vcf_records$CHROM, levels = fai_seqnames)
vcf_records <- vcf_records[order(vcf_records$CHROM, vcf_records$POS), ]
vcf_records$CHROM <- as.character(vcf_records$CHROM)
colnames(vcf_records)<-c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name)

# create vcf header
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Format to vcf header.")
header<-vcf_header(ref_fasta =ref_fasta)
# Writing the VCF file to disk
if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: Writing ", nrow(vcf_records), " VCF records to ", output_file,".")

output_connection <- file(output_file, open = "w")
on.exit(close(output_connection), add = TRUE)

# Writing the meta-information lines
writeLines(header, con = output_connection)

# Writing the VCF records themselves, including the header line
write.table(vcf_records, file = output_connection, append = TRUE, quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE)

if(logging) message(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),". INFO: all done.")