#########################Construct the epistasis matrix######################
#!/usr/bin/Rscript
rm(list = ls())

# Read the G matrix file
G <- read.table(paste0(i,"_50K_G.txt"), header = FALSE)
G <- as.matrix(G)

# Calculate the matrix order
dims1 <- dim(G)
n <- dims1[1]
cat("Dimensions of matrix", dims1, "\n")

E <- G * G
# Compute the trace of the matrix E
trace_E <- sum(diag(E))
cat("Trace of matrix E:", trace_E, "\n")

E <- E/(trace_E/n)
output_file_name <- paste0(i, "_50K_E.txt")
write.table(E, output_file_name, row.names = FALSE, col.names = FALSE)
print(E[1:5, 1:5])
