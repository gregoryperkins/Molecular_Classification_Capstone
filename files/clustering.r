# This file will demomstrate clustering in R on the
# Golub (1999) data set

source("custom_k_means.r")

# load in the data

clean_train_data = read.csv("clean_gene_train_data.csv", header=F)

clean_test_data = read.csv("clean_gene_test_data.csv", header=F)

# extracting the accession number that identifies
# the patient that the gene came from

train_accession_nums = clean_train_data[1,2:ncol(clean_train_data)]

test_accession_nums = clean_test_data[1,2:ncol(clean_test_data)]

clean_train_data = clean_train_data[2:nrow(clean_train_data),]

clean_test_data = clean_test_data[2:nrow(clean_test_data),]

# Hold onto the feature names

train_feature_names = clean_train_data[2:nrow(clean_train_data),2]

test_feature_names = clean_test_data[2:nrow(clean_test_data),2]

# Splice the data feature names out of the data

clean_train_data = clean_train_data[,2:ncol(clean_train_data)]

clean_test_data = clean_test_data[,2:ncol(clean_test_data)]


# read in the actual class labels

actual_class = read.table("gene-expression-data/actual.csv", header=TRUE, sep=",")

actual_class = t(actual_class[2])

# Form the data into a data frame

clean_train_data = data.frame(t(clean_train_data))

clean_test_data = data.frame(t(clean_test_data))

# Name the rows and columns appropriately

rownames(clean_train_data) = as.character(train_accession_nums) 
rownames(clean_test_data) = as.character(test_accession_nums)

colnames(clean_train_data) = train_feature_names
colnames(clean_test_data) = test_feature_names
# perfrom clustering


train_cluster = kmeans(clean_train_data, centers=2, nstart=40)
test_cluster = kmeans(clean_test_data, centers=2, nstart=40)

# find the actual class for each training and testing
# sample for verbose output

actual_train_classes = c()

for(index in seq(length(clean_train_data[,1]))) {
  class = actual_class[as.numeric(train_accession_nums[index])]  
  actual_train_classes = append(actual_train_classes, class)
}

actual_test_classes = c()

for(index in seq(length(clean_train_data[,1]))) {
  class = actual_class[as.numeric(test_accession_nums[index])]  
  actual_test_classes = append(actual_test_classes, class)
}

# output the results of the clustering

train_accession_nums = as.numeric(train_accession_nums)
test_accession_nums = as.numeric(test_accession_nums)

train_cluster_results = data.frame(accession_num = train_accession_nums, class = actual_train_classes, cluster = train_cluster$cluster)

test_cluster_results = data.frame(accession_num = test_accession_nums, class = actual_test_classes, cluster = test_cluster$cluster)

write.table(train_cluster_results, file="train_cluster_results.csv", sep=",", row.names=F, col.names=c("Accession Num", "Type", "Cluster"))

write.table(test_cluster_results, file="test_cluster_results.csv", sep=",", row.names=F, col.names=c("Accession Num", "Type", "Cluster"))

# Print the percentage of each type in each training cluster

  ALL_1 = dim(subset(train_cluster_results, train_cluster_results$class == "ALL" & train_cluster_results$cluster == 1))[1]

  ALL_2 = dim(subset(train_cluster_results, train_cluster_results$class == "ALL" & train_cluster_results$cluster == 2))[1]

  AML_1 = dim(subset(train_cluster_results, train_cluster_results$class == "AML" & train_cluster_results$cluster == 1))[1]

  AML_2 = dim(subset(train_cluster_results, train_cluster_results$class == "AML" & train_cluster_results$cluster == 2))[1]

  train_cluster_makeup = data.frame(ALL=c(ALL_1/(ALL_1+ALL_2), ALL_2/(ALL_1+ALL_2)), AML=c(AML_1/(AML_1+AML_2), AML_2/(AML_1+AML_2)))

  train_cluster_makeup

# Print the percentage of each type in each testing cluster

ALL_1 = dim(subset(test_cluster_results, test_cluster_results$class == "ALL" & test_cluster_results$cluster == 1))[1]

ALL_2 = dim(subset(test_cluster_results, test_cluster_results$class == "ALL" & test_cluster_results$cluster == 2))[1]

AML_1 = dim(subset(test_cluster_results, test_cluster_results$class == "AML" & test_cluster_results$cluster == 1))[1]

AML_2 = dim(subset(test_cluster_results, test_cluster_results$class == "AML" & test_cluster_results$cluster == 2))[1]

test_cluster_makeup = data.frame(ALL=c(ALL_1/(ALL_1+ALL_2), ALL_2/(ALL_1+ALL_2)), AML=c(AML_1/(AML_1+AML_2), AML_2/(AML_1+AML_2)))

test_cluster_makeup


# Now lets try clustering with principle components

train_pcs = read.csv("train_principle_components.csv", header=T)

test_pcs = read.csv("test_principle_components.csv", header=T)

# Splice row names out of data

train_pcs = train_pcs[,2:ncol(train_pcs)]
test_pcs = test_pcs[,2:ncol(test_pcs)]

train_cluster = kmeans(train_pcs, centers=2, nstart=40)
test_cluster = kmeans(test_pcs, centers=2, nstart=40)

train_cluster_results = data.frame(accession_num = train_accession_nums, class = actual_train_classes, cluster = train_cluster$cluster)

test_cluster_results = data.frame(accession_num = test_accession_nums, class = actual_test_classes, cluster=test_cluster$cluster)

# Caluclate the cluster makeup for training data

ALL_1 = dim(subset(train_cluster_results, train_cluster_results$class == "ALL" & train_cluster_results$cluster == 1))[1]

ALL_2 = dim(subset(train_cluster_results, train_cluster_results$class == "ALL" & train_cluster_results$cluster == 2))[1]

AML_1 = dim(subset(train_cluster_results, train_cluster_results$class == "AML" & train_cluster_results$cluster == 1))[1]

AML_2 = dim(subset(train_cluster_results, train_cluster_results$class == "AML" & train_cluster_results$cluster == 2))[1]

train_cluster_makeup = data.frame(ALL=c(ALL_1/(ALL_1+ALL_2), ALL_2/(ALL_1+ALL_2)), AML=c(AML_1/(AML_1+AML_2), AML_2/(AML_1+AML_2)))

train_cluster_makeup

# Calculate the cluster makeup for testing data

ALL_1 = dim(subset(test_cluster_results, test_cluster_results$class == "ALL" & test_cluster_results$cluster == 1))[1]

ALL_2 = dim(subset(test_cluster_results, test_cluster_results$class == "ALL" & test_cluster_results$cluster == 2))[1]

AML_1 = dim(subset(test_cluster_results, test_cluster_results$class == "AML" & test_cluster_results$cluster == 1))[1]

AML_2 = dim(subset(test_cluster_results, test_cluster_results$class == "AML" & test_cluster_results$cluster == 2))[1]

test_cluster_makeup = data.frame(ALL=c(ALL_1/(ALL_1+ALL_2), ALL_2/(ALL_1+ALL_2)), AML=c(AML_1/(AML_1+AML_2), AML_2/(AML_1+AML_2)))

test_cluster_makeup


