wd <- "/home/lluis/Documents/master/2TR/AGB/exercises"
setwd(wd)
final_txt <- "_gene_zscore_full-filtered.txt"
tumors <- c("brca", "kirck", "prad", "coad", "luad", "thca", "hnsc", "lusc")
files <- paste0(tumors, final_txt)

# A list to store the dataframes of each cancer
tumor_training <- c()
tumor_test <- c()
for (file in files){
  
  # Read file and prepare data
  data <- read.csv(file, sep="\t", header=TRUE, nrows = 20531, row.names=1)
  data <- na.omit(data)
  
  # Eliminate normal patients
  tumor <- data[, grep("T$", names(data))]
  
  # Classify all the values
  tumor[] <- lapply(tumor, cut, c(-Inf, -2, 2, Inf), c("down", "no_change", "up"))
  
  #Split for training and testing purposes
  mig <- dim(tumor)[2]/2
  tumor_training <- c(tumor_training, tumor[,c(1:mig)])
  tumor_test <- c(tumor_test, tumor[,c((mig+1):(mig*2))])
  
}
training <- data.frame()
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

