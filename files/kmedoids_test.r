cat("\n\n function kmedoids_test() loaded \n\n\n")

kmedoids_test = function() {

  library(MLInterfaces)
  source("generate_cluster_output.r")

  # load in the data
    
  source("load_data.r")
  load_data()  

  # perform the clustering
  
  cluster = pam(clean_train_data, k=2, metric="euclidean")

  cluster = data.frame(class = actual_train_classes,
              cluster=cluster$cluster)

  cat("\nResults for Euclidean \n")

  generate_cluster_output(cluster)

  cluster = pam(clean_train_data, k=2, metric="manhattan")

  cluster = data.frame(class = actual_train_classes,
              cluster=cluster$cluster)

  cat("\nResults for Manhattan \n")

  generate_cluster_output(cluster)

}
