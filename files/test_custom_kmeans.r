# This file will demomstrate clustering in R on the
# Golub (1999) data set


test_custom_kmeans = function(k=2){
  source("custom_kmeans.r")
  source("generate_cluster_output.r")

  # load in the data

  source("load_data.r")
  load_data()

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


  # perfrom clustering


  custom_cluster = custom_kmeans(clean_train_data, centers=2, metric="euclidean", nstart=100)$cluster

  std_cluster = kmeans(clean_train_data, centers=2,
                    nstart=100)

  ratio = std_cluster$betweenss / 
      (std_cluster$tot.withinss + std_cluster$betweenss)

  #cat("\nbetween_ss / total_ss = ", ratio*100, "%\n\n")

  std_cluster = std_cluster$cluster


  # interpreting reults

  custom_cluster_results = data.frame(class =
   actual_train_classes, cluster=custom_cluster)

  std_cluster_results = data.frame(class = 
    actual_train_classes, cluster=std_cluster)

  generate_output(custom_cluster_results)

}

generate_output = function(cluster_results){

  unique_classes = unique(cluster_results$class)

  num_classes = length(unique_classes)

  num_clusters = length(unique(cluster_results$cluster))

  output_matrix = matrix(nrow=num_clusters, ncol=num_classes)

  colnames(output_matrix) = unique_classes

  for(i in 1:num_clusters){
    for(j in 1:num_classes){
      set = subset(cluster_results, cluster_results$class ==
           unique_classes[j] & cluster_results$cluster == i)
      num_cases = dim(set)[1]
      output_matrix[i,j] = num_cases
    }
  } 

  output_matrix[,1] = output_matrix[,1] / 
                        sum(output_matrix[,1])

  output_matrix[,2] = output_matrix[,2] / 
                        sum(output_matrix[,2])

  cat("Cluster Results: \n")
  print(output_matrix)

}


