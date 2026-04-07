#!/usr/bin/Rscript
#####################################
# Make phenotype predictions using BGLR
#
# Arguments: [1]   G_file
#            [2]   D_file
#            [3]   E_file
#            [4]   X_file
#            [5]   Y_file
#            [6]   cvs_file
#            [7]   pp
#
#####################################
rm(list=ls())
options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
#####################################

start.time <- Sys.time()
start.date <- format(Sys.Date(), format="%Y%m%d")
library(BGLR)

####################################

G_file <- args[1]   #G
D_file <- args[2]   #D
E_file <- args[3]   #E
X_file <- args[4]   #G_ID
Y_file <- args[5]   #Phenotype
cvs_file <- args[6] #CV group
pp <- args[7]       #save path
#####################################
# Read matrix file
G <- read.table(G_file,check.names=FALSE,header = F)
D <- read.table(D_file,check.names=FALSE,header = F)
E <- read.table(E_file,check.names=FALSE,header = F)
ID <- read.table(X_file,check.names=FALSE,header = F)
Y <- read.table(Y_file, row.names=1,check.names=FALSE,header = T) 
cvs <- read.csv(cvs_file, row.names=1,check.names=FALSE,header = T)

# Set the row name of the matrix
rownames(G)<-ID$V1
rownames(D)<-ID$V1
rownames(E)<-ID$V1
# Match the order of Y and G
Y <- Y[match(rownames(G), rownames(Y)), , drop = FALSE]

G <- as.matrix(G)
D <- as.matrix(D)
E <- as.matrix(E)

print("head(G)")
print(G[1:5,1:5],quote=FALSE)
print("head(Y)")
head(Y)

current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
cat("\n ===================== data is already :",current_time,"===================\n\n")

fold<-length(cvs)  # nfold
krep <-length(cvs) # ntime

setwd(pp)

 ETA=list(list(K=G,model='RKHS'),list(K=D,model='RKHS'),list(K=E, model='RKHS'))
 
# Cyclic trait
for (j in 1 : length(Y)){     #j <-1
  y=Y[, names(Y)[j]]          # The value of the jth trait was extracted
  yhat <- data.frame(cbind(y, yhat = 0))  # Initialize y and predict yhat
  yhat$yhat <- as.numeric(yhat$yhat)
  
# Storage accuracy and unbias
  accuracy_unbiased_list<-list()
  sum_accuracy <- 0
  mean_accuracy <- 0
  sum_unbiased <- 0
  mean_unbiased <- 0
  
# fold * krep
  for (z in 1 : (fold*krep)){
    if(z %% fold == 0){
      kk <- z/fold
      jj <- fold
    }else{
        kk <- as.integer(z/fold) + 1
        jj <- z %% fold}
    
    cat(paste("======================================= CV ",kk,"--",jj,"    ",as.character(Sys.time())),"\n")
    groupkk <- cvs[,kk]
    test <- which(groupkk==jj)
    train <- which(groupkk!=jj)
    
    yNA <- y  # yNA was used as the total set of training and test sets
    yNA[test] <- NA # Set the test set to NA
    
    #  Predicted values for the test set
    fm <- BGLR(y=yNA,ETA=ETA,verbose=FALSE,nIter=5000,burnIn=1000,thin=3)
    yhat$yhat[test] <- fm$yHat[test] 
    
    # Accuracy and unbias of calculation
    accuracy <- cor(yhat$y[test], yhat$yhat[test]) 
    unbiased <- lm(yhat$y[test] ~ yhat$yhat[test])$coefficients[2]
    cat(paste(accuracy,"---",unbiased,"\n"))
   
    # Calculate the average value
    sum_accuracy <- sum_accuracy + accuracy
    mean_accuracy <- sum_accuracy/z
    sum_unbiased <- sum_unbiased + unbiased
    mean_unbiased <- sum_unbiased/z
    
    # Store to accuracy_unbiased_list
    accuracy_unbiased_list[[z]] <- c(kk,jj,accuracy,mean_accuracy,unbiased,mean_unbiased)
  }
  #save predictive Ability
  write.csv(accuracy_unbiased_list,paste0("BLUP_ADI_",names(Y)[j],"_accuracy",".csv"),row.names = F)
} 
