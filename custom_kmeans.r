# Test data

# This file will implement k means in r, allowing for 
# a few different distance metrics

# input data should be of the standard form [samples x features]
custom_kmeans = function(data, centers, metric="euclidean", nstart=20, one_point_cluster=T){

  # optimize for operations, use a matrix

  data = as.matrix(data) 

  # we must set the approporiate distance metric

  if(metric == "euclidean"){ distance = euclidean_distance
  } else if(metric == "manhattan") { distance =   manhattan_distance
  } else if(metric == "cosine") {distance = cosine_distance
  } else {
    return("distance measure ", metric, " not defined \n")
  }

  # intializing variables

  num_data_points = length(data[,1])
  num_features = length(data[1,])
  cluster_assignments_list = c()
  cluster_ratio_list = rep(0,nstart)
  within_cluster_distances = matrix(nrow=nstart, ncol=centers)	
  cluster_assignments = rep(-1,num_data_points)
  distances =  matrix(nrow=num_data_points, 
                        ncol=centers)
  cluster_centers = matrix(nrow=centers, ncol=num_features)
  chosen_indices = rep(0,2)

  # This algorithm is most stable when random starting points
  # are chosen many times, and the best clustering is selected

  for(iteration in 1:nstart){

    prev_cluster_assignments = rep(0, num_data_points)
    
    # first we will need to randomly choose data points in the  
    # to be the starting cluster centers

    for (i in 1:centers){
      index = ceiling(runif(1)*num_data_points)
      while(index %in% chosen_indices){
          index = ceiling(runif(1)*length(data[,1]))
      }
      chosen_indices[i] = index
      cluster_centers[i,] = data[index,]
    }
  
    # the k-means algorithm iterates until the data points do
    # not switch clusters

    while(!all(prev_cluster_assignments == cluster_assignments)){  

      # set the previous cluster assignments to current cluster
      # assignments before finding new ones

      prev_cluster_assignments = cluster_assignments
  
      # calculate the distance of every data point
      # to the center of the clusters

      for(i in 1:num_data_points){
        for(j in 1:centers){
          distances[i,j] = distance(data[i,], 
                                    cluster_centers[j,])
        }
      }

      # assign each data point to its closest cluster
      
      for(i in 1:num_data_points){
        cluster_assignments[i] = which.min(distances[i,])
      }

      # calculate new cluster centers
     
      for(i in 1:centers){
        cluster = data[cluster_assignments == i,]

        # if there is more than one data point in the cluster, 
        # averge them, if not, the cluster center will not change

        if(length(cluster) > num_features) {
          cluster_centers[i,] = colMeans(cluster)
        } else if(length(cluster) == num_features) {
          cluster_centers[i,] = cluster
        } 
      }
    }

    # square the distances for calculation of 
    # between_ss and within_ss

    distances = distances^2
    
   # sum the distances squared within the clusters 
  
    for(i in 1:centers){
  
      # grab all data points in the cluster

      cluster = distances[cluster_assignments == i, ]
     
      # if there is more than one data point in the cluster, 
      # sum all distances from the points to the cluster center,
      # otherwise the only distance is from that one point to
      # the cluster center

      if(length(cluster) > centers){  
        within_cluster_distances[iteration, i] = 
                sum(apply(cluster, 1, function(x) return(min(x))))
      } else if(length(cluster) == centers){
        within_cluster_distances[iteration, i] = min(cluster)
      } 
    }

    # if option is set correctly, do not include results
    # that contain only one data point in a cluster

    if(one_point_cluster == F){
      within_cluster_distances[within_cluster_distances==0] =
          max(within_cluster_distances)
    }

    # add the clustering to a list to choose the best one 
    # after n iterations

    cluster_assignments_list = 
      rbind(cluster_assignments_list, cluster_assignments)


    # sum the distances between the clusters 

    between_ss =
      sum(apply(distances, 1,
              function(x) return(sum(x)-min(x))))

    # sum all distances within cluters
  
    within_ss =
      sum(within_cluster_distances[iteration, ])
  
    # add the ratio of between dist / total dist 
    # for informative output 

    cluster_ratio_list[iteration] = 
      between_ss /  (between_ss + within_ss)
  }

  # choose the best clustering and return it

  total_within_cluster_distances = 
    apply(within_cluster_distances, 1, 
    function(x) return(sum(x)))                      

  best_cluster_index = which.min(total_within_cluster_distances)

  cluster_assignments =  
    cluster_assignments_list[best_cluster_index,]
  
  ratio = cluster_ratio_list[best_cluster_index]

  # return the optimal cluster assignments

  cat("\n between_ss / total_ss = ", ratio*100, "%\n\n")

  cluster = data.frame(cluster_assignments, ratio)
   
  return(cluster)
}


euclidean_distance = function(p1, p2){
  if(length(p1) != length(p2)){
    return("Arguments not of the same length")
  }  

  distance = sqrt(sum((p1-p2)^2))
    
  return(distance)
}

manhattan_distance = function(p1,p2){
   if(length(p1) != length(p2)){
    return("Arguments not of the same length")
  }  

  distance = sum(abs(p1-p2))
    
  return(distance)

}
