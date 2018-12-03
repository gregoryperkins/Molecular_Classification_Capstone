load_data = function(verbose=F) {

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

  # transpose the data

  clean_train_data = t(clean_train_data)

  clean_test_data = t(clean_test_data)

  # Name the rows and columns appropriately

  rownames(clean_train_data) = as.character(train_accession_nums)
  rownames(clean_test_data) = as.character(test_accession_nums)

  colnames(clean_train_data) = train_feature_names
  colnames(clean_test_data) = test_feature_names

  # Associating each sample with its class

  actual_train_classes = c()

  for(index in seq(length(clean_train_data[,1]))) {
    class = actual_class[as.numeric(train_accession_nums[index])]  
    actual_train_classes = append(actual_train_classes, class)
  }

  actual_test_classes = c()

  for(index in seq(length(clean_test_data[,1]))) {
    class = actual_class[as.numeric(test_accession_nums[index])]  
    actual_test_classes = append(actual_test_classes, class)
  }

  # partition the data off into a set of each cancer
  # type

  ALL_train_set = clean_train_data[actual_train_classes == "ALL", ]
  ALL_test_set = clean_test_data[actual_test_classes == "ALL", ]
  ALL_set = rbind(ALL_train_set, ALL_test_set)

  AML_train_set = clean_train_data[actual_train_classes == "AML", ]
  AML_test_set = clean_test_data[actual_test_classes == "AML", ]
  AML_set = rbind(AML_train_set, AML_test_set)


  # make all required variables global
  
  assign("clean_train_data", clean_train_data, 
                              envir = .GlobalEnv)

  assign("actual_train_classes", actual_train_classes,
                                    envir = .GlobalEnv)


  assign("clean_test_data", clean_test_data, 
                          envir = .GlobalEnv)

  assign("actual_class", actual_class, envir=.GlobalEnv)


  assign("actual_test_classes", actual_test_classes,
                                    envir = .GlobalEnv)

  assign("AML_set", AML_set, envir=.GlobalEnv)

  assign("ALL_set", ALL_set, envir=.GlobalEnv)

  assign("ALL_train_set", ALL_train_set, envir=.GlobalEnv)

  assign("ALL_test_set", ALL_test_set, envir=.GlobalEnv)

  assign("AML_train_set", AML_train_set, envir=.GlobalEnv)

  assign("AML_test_set", AML_test_set, envir=.GlobalEnv)


  if(verbose) {

    cat("\n\nVariables clean_train_data, clean_test_data, actual_train_classes, actual_test_classes,AML_set, ALL_set have all been loaded\n\n\n")
  
  }
}
