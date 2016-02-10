wd <- "/home/lluis/Documents/master/2TR/AGB/exercises"
setwd(wd)

first <- read.csv("brca_gene_zscore_full-filtered.txt", sep="\t", 
                  header=TRUE, nrows = 20)
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
