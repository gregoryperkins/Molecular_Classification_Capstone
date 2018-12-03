# load in the data
gene_train_data = read.csv(file="gene-expression-data/data_set_ALL_AML_train.csv", header=FALSE, sep=",",stringsAsFactors = FALSE, quote = "")

gene_test_data = read.csv(file="gene-expression-data/data_set_ALL_AML_independent.csv", header=FALSE, sep=",",stringsAsFactors = FALSE, quote = "")


# Hold onto the feature names

train_feature_names = gene_train_data[2:nrow(gene_train_data),2]

test_feature_names = gene_test_data[2:nrow(gene_test_data),2]


# clean the data by removing all cols with header "call"

clean_train_data = gene_train_data[,-which(gene_train_data 
== "call", arr.ind=T)[,2]]

clean_test_data = gene_test_data[,-which(gene_test_data == "call", arr.ind=T)[,2]]

# Hold onto the accession numbers

train_accession_nums = clean_train_data[1,3:ncol(clean_train_data)]

test_accession_nums = clean_test_data[1,3:ncol(clean_test_data)]

# clean the data into just features

clean_train_data = clean_train_data[2:nrow(clean_train_data),3:ncol(clean_train_data)]

clean_test_data = clean_test_data[2:nrow(clean_test_data),3:ncol(clean_test_data)]

# set row and column names to true names

rownames(clean_train_data) = c(train_feature_names)
colnames(clean_train_data) = as.character(train_accession_nums)

rownames(clean_test_data) = c(test_feature_names)
colnames(clean_test_data) = as.character(test_accession_nums)

# write this new cleaned data to file

write.csv(clean_train_data, "clean_gene_train_data.csv")

write.csv(clean_test_data, "clean_gene_test_data.csv")

