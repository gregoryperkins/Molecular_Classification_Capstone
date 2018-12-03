# Variance distribution of genes in the data

source("load_data.r")

load_data()

all_data = rbind(ALL_set, AML_set)

num_genes = dim(all_data)[2]

variance_list = rep(1, num_genes)

for (i in 1:num_genes) {

  variance_list[i] = var(all_data[,i])

}


boxplot(variance_list)
