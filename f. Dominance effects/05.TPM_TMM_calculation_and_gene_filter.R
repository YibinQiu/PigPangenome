# ============================================================
# TPM Calculation 
# ============================================================

# --------------------------
# Section 1: Data Loading
# --------------------------

# Read raw read count data from file
read_counts <- read.table("raw_read_count_data.txt",
                         header = FALSE)

# Transpose matrix to have samples as rows and features as columns
read_counts_t <- t(read_counts)

# --------------------------
# Section 2: Data Preparation
# --------------------------

# Set column names from first row of data
colnames <- read_counts_t[1, ]
colnames(read_counts_t) <- colnames

read_counts_t <- as.data.frame(read_counts_t[-1, ])

read_counts_t$gene_length <- as.numeric(read_counts_t$gene_length)

read_counts_t$gene_length_kb <- read_counts_t$gene_length / 1000

# Reorder columns to place gene_length_kb after gene_length
column_order <- c(colnames(read_counts_t)[1:2], 
                 "gene_length_kb", 
                 colnames(read_counts_t)[3:(ncol(read_counts_t) - 1)])
read_counts_t <- read_counts_t[, column_order]

# --------------------------
# Section 3: TPM Calculation
# --------------------------

for (i in 4:ncol(read_counts_t)) {
  read_counts_t[, i] <- as.numeric(read_counts_t[, i])
  read_counts_t[, i] <- read_counts_t[, i] / read_counts_t$gene_length_kb
}

column_sums <- colSums(read_counts_t[, 4:ncol(read_counts_t)])

for (i in 4:ncol(read_counts_t)) {
  sample_index <- i - 3  # Adjust index for column_sums
  read_counts_t[, i] <- read_counts_t[, i] / (column_sums[sample_index] / 1e6)
}

# --------------------------
# Section 4: Result Output
# --------------------------

# Write TPM values to output file
write.table(read_counts_t,
            "TPM.txt",
            sep = " ",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# ============================================================
# TMM Calculation 
# filters genes based on TPM and read count
# ============================================================

library(data.table)
library(edgeR)

# --------------------------
# Section 1: Data Loading and Preparation
# --------------------------

# Read raw count data
Counts <- fread("raw_read_count_data.txt",
               data.table = FALSE,
               header = FALSE)

Counts <- t(Counts)

# Set column names from first row
colnames <- Counts[1, ]
colnames(Counts) <- colnames

Counts <- as.data.frame(Counts[-1, ])

rownames(Counts) <- Counts[, 1]

Counts <- Counts[, -c(1, 2)]

# Convert all columns to numeric
Counts[] <- lapply(Counts, as.numeric)

# Read TPM data
TPM <- fread("TPM.txt",
            data.table = FALSE)

rownames(TPM) <- TPM[, 1]

# Remove first three columns (gene ID and gene_length and gene_length_kb)
TPM <- TPM[, -c(1, 2, 3)]

# Convert all columns to numeric
TPM[] <- lapply(TPM, as.numeric)

# --------------------------
# Section 2: TMM Normalization
# --------------------------

# Get sample IDs
samids <- colnames(Counts)

expr_counts <- Counts
expr <- DGEList(counts = expr_counts)

# Get dimensions
nsamples <- length(samids)  # Number of samples
ngenes <- nrow(expr_counts) # Number of genes

y <- calcNormFactors(expr, method = "TMM")

TMM <- cpm(y, normalized.lib.sizes = TRUE)

# --------------------------
# Section 3: Gene Filtering
# --------------------------

# Set expression thresholds
count_threshold <- 6          # Minimum read count
tpm_threshold <- 0.1          # Minimum TPM value
sample_frac_threshold <- 0.2  # Minimum fraction of samples that must pass thresholds

expr_tpm <- TPM[rownames(expr_counts), samids]

# Count number of samples passing TPM threshold per gene
tpm_th <- rowSums(expr_tpm >= tpm_threshold)

# Count number of samples passing count threshold per gene
count_th <- rowSums(expr_counts >= count_threshold)

# Create logical masks for filtering
ctrl1 <- tpm_th >= (sample_frac_threshold * nsamples)  # TPM threshold
ctrl2 <- count_th >= (sample_frac_threshold * nsamples) # Count threshold
mask <- ctrl1 & ctrl2                                  # Combined mask

# Apply filter to TMM matrix
TMM_pass <- TMM[mask, ]
TMM_pass <- as.data.frame(TMM_pass)

# --------------------------
# Section 4: Inverse Normal Transformation
# --------------------------

# Define inverse normal transformation function
inverse_normal_transform <- function(x) {
  qnorm(rank(x) / (length(x) + 1))
}

TMM_inv <- t(apply(TMM_pass, 
                  MARGIN = 1, 
                  FUN = inverse_normal_transform))
TMM_inv <- as.data.frame(TMM_inv)

TMM_inv_with_rownames <- cbind(gene_id = rownames(TMM_inv), TMM_inv)

# --------------------------
# Section 5: Output Results
# --------------------------

# Write transformed TMM data
write.table(TMM_inv_with_rownames,
           "TMM_gene_filter.txt",
           sep = " ",
           row.names = FALSE,
           col.names = TRUE,
           quote = FALSE)

# Write filtered TPM data
TPM_pass <- TPM[mask, ]
TPM_pass <- as.data.frame(TPM_pass)
TPM_pass_with_rownames <- cbind(gene_id = rownames(TPM_pass), TPM_pass)

write.table(TPM_pass_with_rownames,
           "TPM_gene_filter.txt",
           sep = " ",
           row.names = FALSE,
           col.names = TRUE,
           quote = FALSE)