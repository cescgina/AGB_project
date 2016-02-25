#!/usr/bin/env Rscript
# Makes a Naïve Bayes classifer for files of the TCGA

args <- commandArgs(trailingOnly=TRUE)
n_args <- length(args)
if (n_args <= 6){
  msg = "This program trains Naïve Bayes clasificator and work with the data.\n"
  msg2 = "-d Set the working directory.\n"
  msg3 = "-c Set the names of the cancers to analyze: brca kirc prad in lowercase
  and  space separated\n"
  msg4 = "-o Set the name of the output file with the testing result.\n"
  msg5 = "Example how to run it:$./first_script.R -d /home/user/ -c brca luad-o output \n"
  msg6 = "-p Flag to count with pseudocounts it should be on the last position"
  stop(msg, msg5, msg2, msg3, msg4, call.=FALSE)
} else if (args[n_args] == "-p"){
  wd <- args[2]
  tumors <- args[c(4:(length(args)-3))]
  output <- args[length(args)-1]
  pseudocount = TRUE
} else {
  wd <- args[2]
  tumors <- args[c(4:(length(args)-2))]
  output <- args[length(args)]
  pseudocount = FALSE
}


# TODO: Remove the follwing line:
# stop("Protect our computer to run the whole programm'by accident' ", call.=FALSE)

# Quitar valor

# Things to work with
# wd <- "../project_data"
# tumors <- c("brca", "kirc", "prad", "coad", "luad", "thca", "hnsc", "lusc")

# Setting work envrionment
#wd <- "../project_data"
setwd(wd)

# Initial options
final_txt <- "_gene_zscore_full-filtered.txt"

#tumors <- c("brca", "kirc", "prad", "coad", "luad", "thca", "hnsc", "lusc")

rownames_m <- c('up', 'no_change', 'down')

# Needed to classify the data.frame when we read it
factorize <- function(invec) {
  # Classify the values to up, down or no_change dependiing on their values
  # from: http://stackoverflow.com/a/35338274/
  factorized <- rep("no_change", length(invec))
  factorized[invec > 2]  <- "up"
  factorized[invec < -2]  <- "down"
  factor(factorized, c("down", "no_change", "up"))
}

pseudocounts <- function(x){
  # adds one to whatever is has to count frequencies with pseudocounts
  x+1
}

split_names <- function(dataset){
  # We stay with just the gene number just in case the letter  or ?
  # create some problems
  gene_names <- data.frame(do.call(rbind, strsplit(rownames(dataset), "|", fixed=TRUE)))[,2]
  rownames(dataset) <- gene_names
  dataset
}

probability <- function(x, n_patient){
  # Calculate the probability of x by dividing by the number of patients given
  x/n_patient
}


# Calculates the entropy of the labels (which type of cancer)
partial_entropy <- function(x){
  if (x != 0){
    - x * log2(x)
  }
}

## Functions to calculate the Mutual information
# Function to calculate the partial mutual information
funct3 <- function(stat, pr, ge, datase){
  if (datase[stat, ge] != 0){
    return(datase[stat, ge ] * log2(datase[stat, ge]/pr[stat, ge]))
  }
  return(0)
}

# Iterates over the states (up, down, no_change)
funct2 <- function(dataset, state, pro, gen){
  if (gen %in% colnames(dataset)){
    sum(sapply(state, funct3, pr = pro, ge = gen, datase = dataset)) 
  } else{
    0
  }
}

# First function iterate over the datasets
funct1 <- function(gene, dataset_list, states, prob){
  sum(sapply(dataset_list, funct2, state = states, pro = prob, gen = gene))
}

### Starting the analysis

for (cancer in tumors){
  
  # Reading the files if prepared files already exist then omit reading them
  # (it is very slow to prepare the data about 4 minuts and much memory RAM)
#   if (file.exists(paste0(cancer, "_training_RDS.bin"))){
#     next
#   }
  
  # Read file and prepare data
  file <- paste0(cancer, final_txt)
  cat("Reading file ", file, "\n")
  dataset <- read.csv(file, sep="\t", header=TRUE, nrows = 20531, row.names = 1)
  dataset <- na.omit(dataset)
  
  # Remove the name of the gene to work with just the numeric value
  dataset <- split_names(dataset)
  
  # Eliminate normal patients by keeping those ending with T
  cat("Removing healthy patients\n")
  tumor <- dataset[, grep("T$", names(dataset))]
  rm(dataset)
  
  # Classify all the values
  tumor[] <- lapply(tumor, factorize)
  
  # If we want to ensure randomness we should use:
  # minimum is 210 and then we would have 80% for training and 20% for testing 
  # for all genes
  selected<-sample(ncol(tumor), 210) 
  tumor_training <- tumor[,selected] 
  tumor_test <- tumor[,-selected]
  cat("Saving the testing data of ", cancer, "\n")
  saveRDS(tumor_test, file=paste0(cancer,"_test_RDS.bin"))
  rm(tumor)
  
  # Calculates the number of times each gene have each type of expression for 
  # each dataset
  cat("Calculating the frequencies of ", cancer, "\n")
  freq_training <- apply(tumor_training, 1, function(x)
    tabulate(factor(x, levels = rownames_m), length(rownames_m)))
  row.names(freq_training) <- rownames_m
  
  # Adds one as a pseudocount before doing any probability
  if (pseudocount){
    freq <- apply(freq_training, c(1,2), pseudocounts)
  }
  cat("Saving prepared training data of ", cancer, "\n")
  saveRDS(freq_training, file=paste0(cancer,"_training_RDS.bin"))
  
  #From time to time it is advaisable to restart the session to free memory of the
  #r session with ctrl + shift + F10, and continue where we left it
}

# Store the probability of being up, down or no_change regardless of which 
# type of cancer  it has been
data_g <- matrix(0, nrow = 3, ncol = 20531)

# Read the prepared data
for (cancer in tumors){
   assign(cancer, readRDS(paste0(cancer, "_training_RDS.bin")))

# Calculates the probability of being up, down or no_change regardless of which 
# type of cancer  it has been 
  cat("Calculating the probability of being up, down or no_change for each gene.\n")
  col_cancer <- colnames(get(cancer))
  mcol <- match(col_cancer, colnames(data_g))
  if (length(col_cancer) == 20531){
    data_g <- get(cancer) + data_g
  } else {
    data_g[, mcol] <- get(cancer) + data_g[, mcol]
  }
}

patient_gene <- apply(data_g, 2, sum)
P_g <- t(t(data_g) / patient_gene) # Calculates the probability of a gene to be up, 
# down or no_change taking into account all the tumor patients

# Total number of patients with tumor in the training set
 cat("Calculating the number of patient by cancer and the total number of patients\n")
patient_tumor <-  sapply(tumors, function(x)sum(get(x)[,1]))
total_n_patients <- sum(patient_tumor)

# Caclulates the probability of a gene to be up, down or no_change given it is 
# of a given tumor, ie the conditional probability of a gene
cat("Calculating the conditional probability of being up, down or no_change given a cancer\n",
    "Calculates the prior probability\n")
Pg_tumor <- c() # Probability of genes of a given tumor of being up, down or no_
P_cancer <- c() # Stores the prior probability 
for (cancer in tumors){
  n_patients <- sum(get(cancer)[,1])
  name_v <- paste0("Pg_", cancer)
  assign(name_v, apply(get(cancer), c(1,2), probability, 
                       n_patient = patient_tumor[cancer]))
  
  Pg_tumor <- c(Pg_tumor, c=name_v)
  P_cancer[cancer] <- n_patients/total_n_patients
}

# Calculates the probability of being up, down or no_change and of a tumor
# i.e it is not the conditional probability of a given cancer,  
# the joint probability

Pg_tumor_state <- c()

for (cancer in tumors){
  cat("Calculates the joint probability o f being up, down or no_change and being of tumor", cancer, ".\n")
  name_v <- paste0("Pg_", cancer, "_state")
  assign(name_v, apply(get(cancer), c(1,2), probability, 
                       n_patient = total_n_patients))
  Pg_tumor_state <- c(Pg_tumor_state, c=name_v)
}

cat("Calculating the entropy\n")
entropy <- sum(sapply(P_cancer,  partial_entropy))

# Calculate the partial entropy for all the genes 
results <- sapply(col_cancer, funct1, dataset_list = mget(Pg_tumor_state), 
                   states = rownames_m, prob = P_g)
cat("Calculating the mutual information of each gene\n")
MI <- results + entropy
saveRDS(MI, "mutual_information.bin")

# Select which genes should be taken to decide to build the classifier
quan <- quantile(MI)
threshold <- quan[4]+1.5*(quan[4]-quan[2])
relevant_genes <- names(MI[MI>threshold])

# Selecting with which genes are worth using
for (cancer in tumors) {
  temp_tumor = readRDS(paste0(cancer,"_test_RDS.bin"))
  relevant_genes<- intersect(relevant_genes,rownames(temp_tumor))
}
# For each patient it should return which label and which score has using the 
# information of the expressions of each gene.
#tumor_test has rownames=patients and colnames=genes
#Probabilites in Pg_cancername
predict <- function(patient, tumors, tumor_test, best_genes, cancer_probs){
  max <- -1e6
  max_cancer <- ''
  for (cancer in tumors){
    Probs <- get(paste0("Pg_",cancer)) # TODO: Everything should go as an argument!
    P <- log2(cancer_probs[cancer])
    for (genes in best_genes){
      p <- tumor_test[genes,patient]
      if (Probs[p,genes] == 0){
        logval <-  0
      }
      else {
        logval <-  log2(Probs[p,genes])
      }
      P <-  P + logval
    }
    if (P > max){
      max <-  P
      max_cancer <-  cancer
    }
  } 
  return(c(max, max_cancer,"target",patient))
}

# The dataset to store the predictions and score of the testing data 
output_dataset <-  data.frame(Score = numeric(),
                              Prediction = character(),
                              Label=character(),
                              Patient=character(),
                              stringsAsFactors = FALSE)

Pg_dataset <-  mget(Pg_tumor_state)
for (tumor_type in tumors){
  temp_tumor <-  readRDS(paste0(tumor_type,"_test_RDS.bin"))
  for (patient in colnames(temp_tumor)){
    prediction <-  predict(patient, tumors, temp_tumor, relevant_genes, P_cancer)
    prediction[3] <- as.character(tumor_type)
    output_dataset[nrow(output_dataset)+1,] <- prediction
  }
  cat("Finished with ",tumor_type)
}

saveRDS(output_dataset,file="output.bin")
write.table(ouptput_dataset,file=output,sep="\t",row.names = FALSE)


# To test the sensivity of the data:

#Calculate True Positives and False Positives
positives<-function(x,y){
  if (x == y){
    return(1)
  }
  else  {
    return(0)
  }
}

positv <- apply(output_dataset, 1, function(x) {positives(x[2], x[3])})

TP = sum(positv)
FP = length(positv) - TP
print(paste0("TP=",TP," : ",TP*100/(TP+FP),"%"))
print(paste0("FP=",FP," : ",FP*100/(TP+FP),"%"))

