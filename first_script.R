# Makes a Naïve Bayes classifer for files of the TCGA

# Setting work envrionment
wd <- "../project_data"
setwd(wd)

# Initial options
final_txt <- "_gene_zscore_full-filtered.txt"
tumors <- c("brca", "kirc", "prad", "coad", "luad", "thca", "hnsc", "lusc")
rownames_m <- c('up', 'no_change', 'down')

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
  # minimum is 210 and then we would have 80% for training and 20% for testing 
  # for all genes
  selected<-sample(ncol(tumor), 210) 
  tumor_training <- tumor[,selected] 
  tumor_test <- tumor[,-selected] 
  saveRDS(tumor_test, file=paste0(cancer,"_test_RDS.bin"))
  rm(tumor)
  
  # Calculates the number of times each gene have each type of expression for 
  # each dataset
  freq_training <- apply(tumor_training, 1, function(x)
    tabulate(factor(x, levels = rownames_m), length(rownames_m)))
  row.names(freq_training) <- rownames_m
  
  
  # Adds one as a pseudocount before doing any probability
  # freq <- apply(freq, c(1,2), function(x){x+1})
  saveRDS(freq_training, file=paste0(cancer,"_training_RDS.bin"))
  
  
  #From timet to time it is advaisable to restart the session to free memory of the
  #r session with ctrl + shift + F10, and continue where we left it

#   # Count how many of each type are and assign them to the name of the cancer
#   # This way we can calculate the probability of up, down or no_change in a tumor
#   assign(file, table(unlist(tumor_training)))
#   assign(file, table(unlist(tumor_test)))
}

# Read the prepared data
for (cancer in tumors){
  assign(cancer, readRDS(paste0(cancer, "_training_RDS.bin")))
}

# Calculates the probability of being up, down or no_change regardless of which 
# type of cancer  it has been 
data_g <- matrix(0, nrow = 3, ncol = 20531)
for (cancer in tumors){
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
patient_tumor <-  sapply(tumors, function(x)sum(get(x)[,1]))
total_n_patients <- sum(patient_tumor)

# Caclulates the probability of a gene to be up, down or no_change given it is 
# of a given tumor, ie the conditional probability of a gene
Pg_tumor <- c() # Probability of genes of a given tumor of being up, down or no_
P_cancer <- c() # Stores the prior probability 
probability <- function(x, n_patient){
  x/n_patient
}
for (cancer in tumors){
  n_patients <- sum(get(cancer)[,1])
  name_v <- paste0("Pg_", cancer)
  assign(name_v, apply(get(cancer), c(1,2), probability, 
                       n_patient = patient_tumor[cancer]))
  
  Pg_tumor <- c(Pg_tumor, name_v)
  P_cancer[cancer] <- n_patients/total_n_patients
}

# Calculates the probability of being up, down or no_change and of a tumor
# i.e it is not the conditional probability of a given cancer
Pg_tumor_state <- c()

for (cancer in tumors){
  name_v <- paste0("Pg_", cancer, "_state")
  assign(name_v, apply(get(cancer), c(1,2), probability, 
                       n_patient = total_n_patients))
  Pg_tumor_state <- c(Pg_tumor_state, name_v)
}

# Calculates the entropy of the labels (which type of cancer)
partial_entropy <- function(x){
  if (x != 0){
    - x * log2(x)
  }
}
entropy <- sum(sapply(P_cancer,  partial_entropy))

# Calculates the relative entropy of a given gene
gene <- "100130426" 
val <- 0
for (state in rownames_m){
  for (dataset in Pg_tumor_state){
    if (gene %in% colnames(get(dataset))){
      if (get(dataset)[state, gene] != 0){
        new_val <- get(dataset)[state, gene ] * 
          log2(get(dataset)[state, gene]/P_g[state, gene]) 
        val <- val + new_val
      }
    }
  }
}

# Calculates the mutual information of the given gene
result <- entropy + val


funct10 <- function(gene, dataset_list, states, prob){
  val <- 0
  for (dataset in dataset_list){
    if (gene %in% colnames(dataset)){
      for (state in states){
        if (dataset[state, gene] != 0){
          val <- val + dataset[state, gene ] * 
            log2(dataset[state, gene]/prob[state, gene])
        }
      }
    }
  }
  val
}


funct3 <- function(stat, pr, ge, datase){
  if (datase[stat, ge] != 0){
    return(datase[stat, ge ] * log2(datase[stat, ge]/pr[stat, ge]))
  }
  return(0)
}

funct2 <- function(dataset, state, pro, gen){
  if (gen %in% colnames(dataset)){
    sum(sapply(state, funct3, pr = pro, ge = gen, datase = dataset)) 
  } else{
    0
  }
}

funct1 <- function(gene, dataset_list, states, prob){
  sum(sapply(dataset_list, funct2, state = states, pro = prob, gen = gene))
}

results0 <- sapply(col_cancer, funct1, dataset_list = mget(Pg_tumor_state), 
                   states = rownames_m, prob = P_g)

results <- sapply(col_cancer, funct10, dataset_list = mget(Pg_tumor_state), 
                    states = rownames_m, prob = P_g)

MI <- results + entropy
saveRDS(MI, "mutual_information.bin")

# Select which genes should be taken to decide to build the classifier
quan <- quantile(MI)
threshold <- quan[4]+1.5*(quan[4]-quan[2])
relevant_genes <- names(MI[MI>threshold])

# For each patient it should return which label and which score has using the 
# information of the expressions of each gene.
classifier <- function(patient, Prior_prob, Cond_prob){

}

#tumor_test has rownames=patients and colnames=genes
#Probabilites in Pg_cancername
predict <- function(patient,tumors,tumor_test,best_genes,Pg_tumor_state){
  max=0
  max_cancer=''
  for (cancer in tumors){
    Probs = get(paste0("Pg_",cancer,"_state"))
    P = log2(cancer_probs[cancer])
    for (genes in best_genes){
      P = P + log2(Probs[tumor_test[genes,patient],genes])
    }
    if (P > max){
      max = P
      max_cancer = cancer
    }
  } 
  return(c(max,max_cancer,"target",patient))
}
output_dataset = data.frame(c("Score","Predicition","Lable","Patient"))
Pg_dataset = mget(Pg_tumor_state)
for (tumor_type in tumors){
  temp_tumor = readRDS(paste0(tumor_type,"_test_RDS.bin"))
  for (patient in colnames(temp_tumor)){
    prediction = predict(patient,tumors,temp_tumor,best_genes,Pg_dataset)
    prediction[3]=tumor_type
    rbind(output_dataset,prediction)
  }
  