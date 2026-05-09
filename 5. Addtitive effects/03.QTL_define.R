identify_QTLs <- function(input_file, output_file) {
  # Read GWAS results
  gwas_results <- read.table(input_file, header = TRUE)
  colnames(gwas_results) <- c("SNP", "allele1", "allele0", "Weight", "Zscore", "P", "Direction", "CHR", "POS","pve")
  
  # Set P-value significance threshold
  p_threshold <- 1/38136159
  
  # Filter significant SNPs
  significant_snps <- gwas_results[gwas_results$P < p_threshold, ]
  
  # Sort by chromosome and position
  significant_snps <- significant_snps[order(significant_snps$CHR, significant_snps$POS), ]
  
  # Initialize QTL list
  QTLs <- list()
  
  # Iterate to identify lead SNPs and QTLs
  for (chr in unique(significant_snps$CHR)) {
    chr_snps <- significant_snps[significant_snps$CHR == chr, ]
    
    while (nrow(chr_snps) > 0) {
      # Find the lead SNP with the smallest P-value on the current chromosome
      lead_snp <- chr_snps[which.min(chr_snps$P), ]
      
      # Define QTL boundaries
      left_boundary <- lead_snp$POS - 500000
      right_boundary <- lead_snp$POS + 500000
      
      # Find the farthest significant SNP within the left boundary
      left_snps <- chr_snps[chr_snps$POS <= lead_snp$POS & chr_snps$POS >= left_boundary, ]
      if (nrow(left_snps) > 0) {
        QTL_start <- min(left_snps$POS)
      } else {
        QTL_start <- left_boundary
      }
      
      # Find the farthest significant SNP within the right boundary
      right_snps <- chr_snps[chr_snps$POS >= lead_snp$POS & chr_snps$POS <= right_boundary, ]
      if (nrow(right_snps) > 0) {
        QTL_end <- max(right_snps$POS)
      } else {
        QTL_end <- right_boundary
      }
      
      qtl_snps <- chr_snps[chr_snps$POS >= QTL_start & chr_snps$POS <= QTL_end, ]
      
      # Add QTL to the list
      sv_count <- sum(grepl(">|UG", qtl_snps$SNP))  # Count SVs (containing ">" or "UG")
      snp_count <- sum(nchar(qtl_snps$allele1) == 1 & nchar(qtl_snps$allele0) == 1)  # Count SNPs
      indel_count <- nrow(qtl_snps) - sv_count - snp_count  # Count INDELs
      QTLs <- append(QTLs, list(data.frame(
        CHR = chr,
        QTL_start = max(QTL_start, 0),
        QTL_end = QTL_end,
        lead_SNP = lead_snp$SNP,
        lead_pve = lead_snp$pve,
        lead_P = lead_snp$P,
        sv_count = sv_count,
        indel_count = indel_count,
        snp_count = snp_count
      )))
      
      # Remove SNPs within the defined QTL region
      chr_snps <- chr_snps[chr_snps$POS > QTL_end | chr_snps$POS < QTL_start, ]
    }
  }
  
  # If no QTL is found, return an empty dataframe
  if (length(QTLs) == 0) {
    return(data.frame())
  }
  
  # Convert QTL list to dataframe
  QTLs_df <- do.call(rbind, QTLs)
  
  # Save QTL dataframe to output file
  write.table(QTLs_df, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Return QTL dataframe
  return(QTLs_df)
}

# Define a function to batch generate QTL files, skipping files without QTLs
batch_identify_QTLs <- function(input_dir, output_dir) {
  # Get the list of all P-value files in the input directory
  p_value_files <- list.files(input_dir, pattern = "*.txt", full.names = TRUE)
  
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Loop through each P-value file
  for (input_file in p_value_files) {
    # Extract file name (without path and extension)
    trait_name <- tools::file_path_sans_ext(basename(input_file))
    
    # Define output file path
    output_file <- file.path(output_dir, paste0(trait_name, "_QTL.txt"))
    
    # Generate QTL
    QTLs_df <- identify_QTLs(input_file, output_file)
    
    # Save file only if QTLs_df is not empty
    if (nrow(QTLs_df) > 0) {
      write.table(QTLs_df, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
      cat("Processed:", input_file, "->", output_file, "\n")
    } else {
      cat("No QTL found for:", input_file, "\n")
    }
  }
}

# Run the function to process all P-value files
input_dir <- "./additive/sig_site/"  # Input file directory
output_dir <- "./additive/QTL_resullt/"  # Output file directory

batch_identify_QTLs(input_dir, output_dir)

# Define a function to count the number of lines (excluding header) for each file and skip empty files
count_lines_in_files <- function(directory) {
  # Get the list of all files in the directory
  files <- list.files(directory, pattern = "*.txt", full.names = TRUE)
  
  # Initialize an empty dataframe to store results
  file_line_counts <- data.frame(File = character(), LineCount = integer(), stringsAsFactors = FALSE)
  
  # Loop through each file and count lines
  for (file in files) {
    total_lines <- length(readLines(file))
    line_count <- total_lines - 1  # Subtract header line
    # Ensure line count is not negative (avoid empty files)
    if (line_count > 0) {
      file_line_counts <- rbind(file_line_counts, data.frame(File = basename(file), LineCount = line_count))
    }
  }
  
  # Return dataframe with file names and line counts (excluding files with 0 lines)
  return(file_line_counts)
}

# Use the function to count lines in each file in the folder
directory <- "./additive/QTL_resullt/"  # Input file directory
line_counts <- count_lines_in_files(directory)
# Print results
print(line_counts)
write.table(line_counts, "./additive/QTL_number.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

directory <- "./additive/sig_site/"  # Input file directory
line_counts <- count_lines_in_files(directory)
print(line_counts)

# Define a function to combine all QTL files into one file
combine_QTL_files <- function(input_dir, output_file) {
  # Get the list of all QTL files in the input directory
  qtl_files <- list.files(input_dir, pattern = "*_QTL.txt", full.names = TRUE)
  
  # Initialize an empty dataframe to store all QTL data
  all_qtls <- data.frame()
  
  # Loop through each QTL file
  for (file in qtl_files) {
    # Read the current QTL file
    qtls <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Add file name column
    qtls$file_name <- basename(file)
    
    # Append to the total dataframe
    all_qtls <- rbind(all_qtls, qtls)
  }
  
  # Save the combined QTL dataframe to the output file
  write.table(all_qtls, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Return the combined QTL dataframe
  return(all_qtls)
}

# Use the function to combine all QTL files
combined_output_file <- "./additive/combined_QTLs2.txt"  # Output file path
combine_QTL_data <- combine_QTL_files(output_dir, combined_output_file)

combine_QTL_data <- combine_QTL_data %>%
  mutate(lead_variant_type = case_when(
    grepl(">|UG", lead_SNP) ~ "SV",
    nchar(gsub("[^A-Z]", "", lead_SNP)) >= 3 ~ "INDEL",
    TRUE ~ "SNP"
  ))
write.table(combine_QTL_data, file = "./additive/combined_QTLs.txt", row.names = FALSE, quote = FALSE, sep = "\t")
