kmeans_pca_clustering = function() {
  
  # load in the custom kmeans script and output 
  # generation script

  source("custom_kmeans.r")
  source("generate_cluster_output.r")

  # load in the principle components

  train_pcs = read.csv("train_principle_components.csv",
        header=F)

  test_pcs = read.csv("test_principle_components.csv", 
        header=F)
    
  # grab the accession numbers for each sample

  train_accession_nums = t(train_pcs[2:nrow(train_pcs),1])
  test_accession_nums = t(test_pcs[2:nrow(test_pcs),1])

  # remove the accession number column and 
  # principle component row 
  
  train_pcs = train_pcs[2:nrow(train_pcs), 2:ncol(train_pcs)]

  test_pcs = test_pcs[2:nrow(test_pcs), 2:ncol(test_pcs)]

  # convert into the appriate type

  train_pcs = apply(train_pcs, 2, 
    function(x) return(as.numeric(as.character(x))))

  test_pcs = apply(test_pcs, 2, 
    function(x) return(as.numeric(as.character(x))))

  # read in the actual class labels

  actual_class = read.table("gene-expression-data/actual.csv", header=TRUE, sep=",")

  actual_class = t(as.character(actual_class[,2]))

  # Associating each sample with its class

  actual_train_classes = c()

  for(index in seq(length(train_pcs[,1]))) {
    class = actual_class[as.numeric(train_accession_nums[index])]  
    actual_train_classes = append(actual_train_classes, class)
  }

  actual_test_classes = c()

  for(index in seq(length(test_pcs[,1]))) {
    class = actual_class[as.numeric(test_accession_nums[index])]  
    actual_test_classes = append(actual_test_classes, class)
  }


  # cluster the data by principle components iteratively, 
  # removing a PC each run and keeping track of the ratio

  num_pcs = ncol(train_pcs)
  ratios = rep(0, num_pcs)

  for(i in 1:(num_pcs-1)){

    cat("\n using PCS 1 to ", ncol(train_pcs)-(i-1), " \n")
    
    #cluster = custom_kmeans(train_pcs[,1:(ncol(train_pcs)-(i-1))], metric="manhattan", 2, nstart=1000)

    #ratios[i] = cluster[1,2]

    cluster = kmeans(train_pcs[,1:(ncol(train_pcs)-(i-1))], 
                centers=2, nstart=1000)

    ratios[i] = cluster$betweenss / 
            (cluster$betweenss + cluster$tot.withinss)

    cluster = data.frame(class = actual_train_classes,
              cluster=cluster$cluster)

    generate_cluster_output(cluster)

  }

  # now for only the first principle component

  cat("\n using PC 1 \n")

  pcs = cbind(train_pcs[,1], rep(0, length(train_pcs[,1])))

  #cluster = custom_kmeans(pcs, 2, nstart=1000)

  #ratios[i+1] = cluster[1,2]
  
  cluster = kmeans(pcs, centers=2, nstart=1000)

  ratios[i+1] = cluster$betweenss / 
            (cluster$betweenss + cluster$tot.withinss)

  cluster = data.frame(class = actual_train_classes,
              cluster=cluster$cluster)

  plot(1:length(ratios), rev(ratios), xlab = "# of PCs used", ylab = "between_ss / total_ss")

  generate_cluster_output(cluster)  

}


