# ============================================================
# PEER Analysis
# ============================================================

# --------------------------
# Section 1: Parameter Handling
# --------------------------

## Get arguments from command line
args <- commandArgs(trailingOnly = TRUE)

input_1 <- args[1]  # Path to bed files
input_2 <- args[2]  # Tissue type identifier
output_1 <- args[3] # Output directory

# --------------------------
# Section 2: Data Preparation
# --------------------------

# Load PEER library
library(peer)

# Read TMM-normalized expression data
bed <- read.table(paste0(input_1, "/", input_2, "_TSS_TMM.txt"),
                 sep = "\t",
                 header = TRUE,
                 comment.char = "",
                 check.names = FALSE)

# Set row names as gene IDs
rownames(bed) <- bed$gene_id

bed$`#chr` <- NULL    
bed$start <- NULL     
bed$end <- NULL       
bed$gene_id <- NULL   

# Prepare expression matrix for PEER analysis
expr_peer <- bed

# --------------------------
# Section 3: PEER Analysis
# --------------------------

k <- 10

# Initialize PEER model
model <- PEER()

# Set expression data in the model
PEER_setPhenoMean(model, as.matrix(t(expr_peer)))
print(dim(PEER_getPhenoMean(model)))
PEER_setNk(model, k)
PEER_setNmax_iterations(model, 1000)
print(PEER_getNk(model))
PEER_update(model)

# --------------------------
# Section 4: Result Extraction
# --------------------------

factors <- PEER_getX(model)

# Set appropriate row and column names
rownames(factors) <- colnames(expr_peer)  
colnames(factors) <- paste0("peer", c(1:k))  

Alpha <- PEER_getAlpha(model)

alpha0 <- data.frame(
    Peer = c(1:k),          
    Alpha = Alpha,          
    Relevance = 1.0 / Alpha 
)

residuals0 <- t(PEER_getResiduals(model))
residuals0 <- data.frame(
    GENE_ID = rownames(expr_peer),  
    residuals0                      
)
colnames(residuals0) <- colnames(expr_peer)
covariates0 <- data.frame(
    SampleID = colnames(expr_peer),  
    factors                         
)

# --------------------------
# Section 5: Result Output
# --------------------------

write.table(covariates0,
           paste0(output_1, "/", input_2, "_peer_covariates0.txt"),
           quote = FALSE,
           row.names = FALSE,
           sep = "\t")

write.table(alpha0,
           paste0(output_1, "/", input_2, "_peer_alpha0.txt"),
           quote = FALSE,
           row.names = FALSE,
           sep = "\t")

write.table(residuals0,
           paste0(output_1, "/", input_2, "_peer_residuals0.txt"),
           quote = FALSE,
           row.names = FALSE,
           sep = "\t")