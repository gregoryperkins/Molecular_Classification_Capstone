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

custom_k_means(clean_train_data, k=2)

