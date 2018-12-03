# This script contains functions for performing gene-wise t-tests
# with multiple correction types, as well as a function that includes
# permutation testing with the genes


select_de_by_gene_wise_t_test = function(correction="bonferroni", alpha=0.05) {


  cat("T-test correction: ", correction, "\n\n")

  source("load_data.r")
  load_data()

  # Now apply the t-test 

  num_genes = dim(AML_set)[2]

  chosen_indices = c()
  
  p_vals = rep(1, num_genes)

  names(p_vals) = colnames(ALL_set)

  # For each gene, perform a t-test with its specific
  # standard deviation
  for (i in 1:num_genes) {

    p_vals[i] = t.test(AML_set[,i], ALL_set[,i])$p.value

  }

  if(correction == "bonferroni" || correction == "bonf") {

    de_genes = names(p_vals[p_vals < alpha/num_genes])

  } else if(correction == "holms") { 

    # sort the p-values

    sorted_p_vals = sort(p_vals, decreasing=F)
  
    m = length(p_vals)

    cur_divisor = m

    cur_ind = 1

    while(sorted_p_vals[cur_ind] < alpha/cur_divisor) {

      # incriment the current index and divisor

      cur_ind = cur_ind + 1
      cur_divisor = cur_divisor - 1 
  
    } 
    
    de_genes = names(sorted_p_vals[1:cur_ind-1])

  } else if(correction == "fdr") {

    sorted_p_vals = sort(p_vals, decreasing=F)

    m = length(p_vals)

    cur_ind = 1

    while(sorted_p_vals[cur_ind] < ((cur_ind / m) * alpha)){

      # incriment the current index and divisor

      cur_ind = cur_ind + 1

    }

    de_genes = names(sorted_p_vals[1:cur_ind-1])

  }
  
  assign("t_test_selected_genes", de_genes, envir=.GlobalEnv)

}



select_de_by_gene_wise_t_test_with_permutations = function(alpha=0.05, correction="fdr") {

  # For permutation testing, we need to calculate p-values for nperm amount
  # of random permuatations of the data sets, then calculate our actual
  # p-value based on the distribution of p-values over the permuatations

  # We will use the sample() function in R to draw a random permutation of 
  # the list of indices / values of the gene expression data 

  source("load_data.r")
  load_data()

  # Now apply the t-test 

  num_genes = dim(AML_set)[2]

  chosen_indices = c()
  
  actual_stats = rep(1, num_genes)

  # For each gene, perform a t-test with its specific
  # standard deviation

  for (i in 1:num_genes) {

    actual_stats[i] = t.test(AML_set[,i], ALL_set[,i])$statistic

  }

  
  # load in the stats

  # require(data.table)

  cat("\nreading in permutation statistics ... \n\n")

  # stats = fread("permutation_stats.csv") 

  # now we can calculate a p-value for each gene

  p_vals = rep(0, num_genes)

  for (i in 1:num_genes) {

    # gene_stats = stats[,i]
  
    gene_stats = t(stats[,..i])

    # How many of the randomly permuted statistics
    # are equal to or above the actual statistic?

    actual_stat = actual_stats[i]
  
    p_vals[i] = ((length(which(gene_stats >= actual_stat))) 
                    /   length(gene_stats))

  }

  # now we have to adjust the p-vals like in a usual t-test
  # we will use the fdr correction as bonferroni and holms
  # are not suitable with 7000 tests

  names(p_vals) = colnames(ALL_set)

  sorted_p_vals = sort(p_vals, decreasing=F)

  m = length(p_vals)

  cur_ind = 1

  while(sorted_p_vals[cur_ind] < ((cur_ind / m) * alpha)){

    # incriment the current index and divisor

    cur_ind = cur_ind + 1

  }

  # de_genes = names(sorted_p_vals[1:cur_ind-1])

  de_genes = names(sorted_p_vals[1:500])


  assign("permutation_test_selected_genes", de_genes, envir=.GlobalEnv)
 
}


generate_permutation_statistics = function(nperms=100000) {

  # Define the run size, the amount of runs to store in ram while
  # running, the larger the run size, the more ram is used, though
  # the smaller the ram size the more writes are performed

  run_size = 100


  # We must permute patients between groups 

  # We can use data.table for faster reading and writing
  # with fread and fwrite, though accessing 

  require(data.table)

  # Load the data

  source("load_data.r")

  load_data()

  num_genes = dim(ALL_set)[2]

  stats = matrix(1, nrow=run_size, ncol=num_genes)

  all_data = rbind(ALL_set, AML_set)

  labels = rep(1, dim(all_data)[1])

  start = dim(ALL_set)[1]+1
  end = length(labels)

  labels[start:end] = 2

  for (i in 1:nperms) {

    cat("Permutation: ", i ,"\n\n")
  
    labels = sample(labels)
  
    for (j in 1:num_genes) {


      stats[i%%run_size,j] = t.test(all_data[labels==1,j], all_data[labels==2,j])$statistic

    }  
    
    if(i %% run_size == 0){
  
      fwrite(data.table(stats), file="permutation_stats.csv", append=T)

    } 
  }

  if(nperms %% run_size != 0){

    fwrite(data.table(stats), file="permutation_stats.csv", append=T)

  }
}


## Code graveyard

#permute_cols = function(data) {
#
#  num_cols = dim(data)[2]
#
#  nperms = num_cols/2
#
#  config = 1:num_cols
#
#  for (i in 1:nperms) {
#
#    # switch two cols
#
#    col_1 = ceiling(runif(n=1, min=0, max=num_cols))
#
#    col_2 = ceiling(runif(n=1, min=0, max=num_cols))
#
#    while(col_1 == col_2){
#       
#      col_2 = ceiling(runif(n=1, min=0, max=num_cols))
#
#    }
#
#    tmp = data[,col_1]
#
#    data[,col_1] = data[,col_2]
#
#    data[,col_2] = tmp
#  }
#
#  return(data)  
#}

