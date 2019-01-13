generate_cluster_output = function(cluster_results){

  unique_classes = unique(cluster_results$class)

  num_classes = length(unique_classes)

  num_clusters = length(unique(cluster_results$cluster))

  output_matrix = matrix(nrow=num_classes, ncol=num_clusters)

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
