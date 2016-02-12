# Makes a Naïve Bayes classifer for files of the TCGA

# Setting work envrionment
wd <- "/home/lluis/Documents/master/2TR/AGB/exercises"
setwd(wd)

# Initial options
final_txt <- "_gene_zscore_full-filtered.txt"
tumors <- c("brca", "kirc", "prad", "coad", "luad", "thca", "hnsc", "lusc")

# A list to store the dataframes of each cancer
tumor_testing <- data.frame()

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
  
  # We stay with just the gene number just in case the letter  or ?create some problems
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
  tumor_training <- tumor[,c(1:mig)]
  tumor_test <- tumor[,c((mig+1):(mig*2))]
  rm(tumor)
  # tumor_testing <- c(tumor_testing, tumor_test)
  
  # If we want to ensure randomness we should use:
  # selected<-sample(ncol(tumor), mig) 
  # tumor_training <- tumor[selected,] 
  # tumor_test <- tumor[-selected,] 
  
  # Calculates the number of times each gene have each type of expression
  occurences <- apply(tumor_training, 1, table)
  max_length <- max(sapply(occurences, length))
  freq <- sapply(occurences, function(x){
    c(x, rep(0, max_length - length(x)))
  })
  rm(occurences)
  rownames(freq) <- c("no_change", "up", "down")
  
  # Adds one as a pseudocount before doing any probability
  # freq <- apply(freq, c(1,2), function(x){x+1})
  saveRDS(freq, file=paste0(cancer,"_training_RDS.bin"))
  
  # Store the value
  assign(cancer, freq)
  
  
#   # Count how many of each type are and assign them to the name of the cancer
#   # This way we can calculate the probability of up, down or no_change in a tumor
#   assign(file, table(unlist(tumor_training)))
#   assign(file, table(unlist(tumor_test)))
#   
#   # Join the data with previous data
#   tumor_training <- rbind(tumor_training, get(file))
#   tumor_test <- rbind(tumor_test, get(file))
}
info <- mget(tumors)
assign(cancer, readRDS(paste0(cancer, "_training_RDS.bin")))
#From now and then it is advaisable to restart the session to free memory of the
#r session with ctrl + shift + F10
Pg_names <- c()
for (cancer in tumors){
  assign(cancer, readRDS(file=paste0(cancer,"_training_RDS.bin")))
  Pg_names <- c(Pg_names, paste0("Pg_", cancer))
}
# Total number of patients with tumor in the training set
total_n_patients <- sum(sapply(tumors, function(x){get(x)[,1]}))

# Caclulates the probability of a gene to be up, down or no_change given it is 
# of a given tumor
Pg_tumor <- c()
for (cancer in tumors){
  n_patients <- sum(get(cancer)[,1])
  assign(paste0("Pg_", cancer), apply(get(cancer), c(1,2), function(x){x/ n_patients}))
  Pg_tumor <- c(Pg_tumor, get(paste0("Pg_", cancer)))
}

# Calculates the probability of being up, down or no_change and of a tumor
# i.e it is not the conditional probability