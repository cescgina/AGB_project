# Makes a Naïve Bayes classifer for files of the TCGA

# Setting work envrionment
wd <- "/home/lluis/Documents/master/2TR/AGB/exercises"
setwd(wd)

# Initial options
final_txt <- "_gene_zscore_full-filtered.txt"
tumors <- c("brca", "kirc", "prad", "coad", "luad", "thca", "hnsc", "lusc")

# Needed to classify the data.frame when we read it
factorize <- function(invec) {
  # from: http://stackoverflow.com/a/35338274/
  factorized <- rep("no_change", length(invec))
  factorized[invec > 2]  <- "up"
  factorized[invec < -2]  <- "down"
  factor(factorized, c("down", "no_change", "up"))
}

pseudocounts <- function(x){
  x+1
}

for (cancer in tumors){
  
  if (file.exists(paste0(cancer, "_training_RDS.bin"))){
    next
  }
  # Read file and prepare data
  file <- paste0(cancer, final_txt)
  dataset <- read.csv(file, sep="\t", header=TRUE, nrows = 20531, row.names = 1)
  dataset <- na.omit(dataset)
  
  # We stay with just the gene number just in case the letter  or ?
  # create some problems
  gene_names <- data.frame(do.call(rbind, strsplit(rownames(dataset), "|", fixed=TRUE)))[,2]
  rownames(dataset) <- gene_names
  rm(gene_names)
  
  # Eliminate normal patients
  tumor <- dataset[, grep("T$", names(dataset))]
  rm(dataset)
  
  # Classify all the values
  tumor[] <- lapply(tumor, factorize)
  
  #Split for training and testing purposes
  mig <- dim(tumor)[2]/2
#   tumor_training <- tumor[,c(1:mig)]
#   tumor_test <- tumor[,c((mig+1):(mig*2))]
#   rm(tumor)
  # tumor_testing <- c(tumor_testing, tumor_test)
  
  # If we want to ensure randomness we should use:
  selected<-sample(ncol(tumor), mig) 
  tumor_training <- tumor[,selected] 
  tumor_test <- tumor[,-selected] 
  rm(tumor)
  
  # Calculates the number of times each gene have each type of expression
  rownames_m <- c('up', 'no_change', 'down')
  freq <- apply(tumor_training, 1, function(x)
    tabulate(factor(x, levels = rownames_m), length(rownames_m)))
  rownames(freq) <- rownames_m
  
  # Adds one as a pseudocount before doing any probability
  # freq <- apply(freq, c(1,2), function(x){x+1})
  saveRDS(freq, file=paste0(cancer,"_training_RDS.bin"))
  
  #From timet to time it is advaisable to restart the session to free memory of the
  #r session with ctrl + shift + F10, and continue where we left it

#   # Count how many of each type are and assign them to the name of the cancer
#   # This way we can calculate the probability of up, down or no_change in a tumor
#   assign(file, table(unlist(tumor_training)))
#   assign(file, table(unlist(tumor_test)))
#   
#   # Join the data with previous data
#   tumor_training <- rbind(tumor_training, get(file))
#   tumor_test <- rbind(tumor_test, get(file))
}

# Read the prepared data
for (cancer in tumors){
  assign(cancer, readRDS(paste0(cancer, "_training_RDS.bin")))
}

# Total number of patients with tumor in the training set
total_n_patients <- sum(sapply(tumors, function(x){get(x)[,1]}))

# Caclulates the probability of a gene to be up, down or no_change given it is 
# of a given tumor
Pg_tumor <- c() # Probability of genes of a given tumor of being up, down or no_
P_cancer <- c() # Stores the prior probability 
for (cancer in tumors){
  n_patients <- sum(get(cancer)[,1])
  name_v <- paste0("Pg_", cancer)
  assign(name_v, apply(get(cancer), c(1,2), function(x){x/ n_patients}))
  Pg_tumor <- c(Pg_tumor, name_v)
  
  
  cancer_n <- paste0("P_", cancer)
  P_cancer[cancer_n] <- n_patients/total_n_patients
}


# Calculates the probability of being up, down or no_change and of a tumor
# i.e it is not the conditional probability of a given cancer
Pg_tumor_state <- c()
for (cancer in tumors){
  name_v <- paste0("Pg_", cancer, "_state")
  assign(name_v, apply(get(cancer), c(1,2), function(x){x/ total_n_patients}))
  Pg_tumor_state <- c(Pg_tumor_state, get(name_v))
}

# Calculates the probability of being up, down or no_change regardless of which 
# type of cancer 
