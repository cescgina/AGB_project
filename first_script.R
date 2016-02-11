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

for (cancer in tumors){
  
  # Read file and prepare data
  file <- paste0(cancer, final_txt)
  dataset <- read.csv(file, sep="\t", header=TRUE, nrows = 205, row.names=1)
  dataset <- na.omit(dataset)
  
  # We stay with just the gene number just in case the letter  or ?create some problems
  gene_names <- data.frame(do.call(rbind, strsplit(rownames(dataset), "|", fixed=TRUE)))[,2]
  rownames(dataset) <- gene_names
  rm(genes_names)
  
  # Eliminate normal patients
  tumor <- dataset[, grep("T$", names(dataset))]
  rm(da)
  
  # Classify all the values
  tumor[] <- lapply(tumor, factorize)
  
  #Split for training and testing purposes
  mig <- dim(tumor)[2]/2
  tumor_training <- tumor[,c(1:mig)]
  tumor_test <- tumor[,c((mig+1):(mig*2))]
  tumor_testing <- c(tumor_testing, tumor_test)
  
  # Calculates the number of times each gene have each type of expression
  occurences <- apply(tumor_training, 1, table)
  max_length <- max(sapply(occurences, length))
  freq <- sapply(occurences, function(x){
    c(x, rep(0, max_length - length(x)))
  })
  rownames(freq) <- c("no_change", "up", "down")
  assign(cancer, freq)
  
#   # Count how many of each type are and assign them to the name of the cancer
#   # This way we can calculate the probability of up, down or no_change in a tumor
#   assign(file, table(unlist(tumor_training)))
#   assign(file, table(unlist(tumor_test)))
#   
#   # Join the data with previous data
#   tumor_training <- rbind(tumor_training, get(file))
#   tumor_test <- rbind(tumor_test, get(file))

  break
}
info <- mget(tumors)
first <- read.csv("brca_gene_zscore_full-filtered.txt", sep="\t", 
                  header=TRUE, nrows = 20, row.names=1)

# Remove all  normal samples
first_tumor <- first[, grep("T$", names(first))]

second <- read.csv("coad_gene_zscore_full-filtered.txt", sep="\t", 
                  header=TRUE, nrows = 20)
# Remove NA from the table
second <- na.omit(second)
# Remove all  normal samples
second_tumor <- second[, grep("T$", names(second))]
mig <- dim(second_tumor)[2]/2
second_tumor_training <- second_tumor[,c(1:mig)]
second_tumor_test <- second_tumor[,c((mig+1):(mig*2))]

# Change the values to it
# second_tumor_training[second_tumor_training$"X6323T" >= 2,] <- "up"
# second_tumor_training[second_tumor_training$"X6323T" <= -2,] <- "down"
# second_tumor_training[second_tumor_training$"X6323T" <= 2 & -2 >= second_tumor_training$"X6323T",] <- "no_change"




classifier <- function(col){
  colum <- c()
  for (x in col){
    if (x >= 2){
      colum <- c( colum, "up")
    }else if (x <= -2){
      colum <- c( colum, "down")
    }else {
      colum <- c( colum, "no_change")
    }
  }
  return(as.factor(colum))
}
table(cut(first_tumor, c(-Inf:-2, -2:2, 2:Inf), c("down", "no_change", "up")))

