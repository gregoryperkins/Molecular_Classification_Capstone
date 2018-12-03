biclustering_test = function(){ 

  library("biclust")
  library(multtest)

  clean_train_data = read.csv("clean_gene_train_data.csv", header=F)

  clean_test_data = read.csv("clean_gene_test_data.csv", header=F)

  # extracting the accession number that identifies
  # the patient that the gene came from

  train_accession_nums = as.numeric(clean_train_data[1,2:ncol(clean_train_data)])

  test_accession_nums = as.numeric(clean_test_data[1,2:ncol(clean_test_data)])

  clean_train_data = clean_train_data[2:nrow(clean_train_data),]

  clean_test_data = clean_test_data[2:nrow(clean_test_data),]

  # Hold onto the feature names

  train_feature_names = clean_train_data[,1]

  test_feature_names = clean_test_data[,1]

  # Splice the data feature names out of the data

  clean_train_data = clean_train_data[,2:ncol(clean_train_data)]

  clean_test_data = clean_test_data[,2:ncol(clean_test_data)]


  # read in the actual class labels

  actual_class = read.table("gene-expression-data/actual.csv", header=TRUE, sep=",")

  actual_class = t(as.character(actual_class[,2]))

  # Form the data into a data frame

  clean_train_data = t(clean_train_data)

  clean_test_data = t(clean_test_data)

  # take a subset of the features

  
   
  clean_train_data = clean_train_data[,1:30]

  # perform clustering

  res = biclust(clean_train_data, method = BCCC(), delta=.3, 
                alpha=.5, number=3)

  drawHeatmap(clean_train_data, res, 1, local=T)

}

