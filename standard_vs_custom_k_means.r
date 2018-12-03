# This file will demomstrate clustering in R on the
# Golub (1999) data set

source("custom_k_means.r")

# load in the data

source("load_data.r")
load_data()

# Compare clustering

custom_results = c()
standard_results = c()

for(i in 1:100){
  cluster = kmeans(clean_train_data, centers=2, nstart=20)
  standard_results = rbind(standard_results, cluster)

  cluster = custrom_k_means(clean_train_data, k=2, 

}




custom_k_means(clean_train_data, k=2)

