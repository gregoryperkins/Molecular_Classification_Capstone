# This script will include functions for selecting
# differentially expressed genes including by fold change, 
# unusual ratio and hypothesis testing

cat("\n\n functions main_test(), select_de_by_fold_change(), select_de_by_unusual_ratio() loaded \n\n\n")

main_test = function() {

  # load in limma to do venn diagrams

  library(limma)

  # load in the moderated t statistic source

  source("moderated_t_statistic.r")

  # load in the gene wise t-test source

  source("gene_wise_t_test.r")

  # load in the data to test the different methods

  source("load_data.r")
  load_data()

  # run all methods on the data

  select_de_by_fold_change()

  select_de_by_unusual_ratio()

  select_de_by_mod_t_stat()

  select_de_by_gene_wise_t_test(correction="fdr")

  select_de_by_sam()

  select_de_by_gene_wise_t_test_with_permutations()

  # print a venn digram of the selected genes
  
  # all_de_id = sam_selected_genes
    
  # all_de_id = append(all_de_id, fold_change_selected_genes)

  all_de_id = permutation_test_selected_genes

  all_de_id = append(all_de_id, moderated_t_stat_selected_genes)

  all_de_id = append(all_de_id, sam_selected_genes)

  all_de_id = unique(all_de_id)

  # tmp1 = all_de_id %in% fold_change_selected_genes

  # tmp2 = all_de_id %in% sam_selected_genes

  tmp2 = all_de_id %in% sam_selected_genes

  tmp3 = all_de_id %in% moderated_t_stat_selected_genes

  tmp4 = all_de_id %in% permutation_test_selected_genes

  venn_data = cbind(tmp2, tmp3, tmp4)

  venn_count = vennCounts(venn_data)
  
  vennDiagram(venn_count, names = c("sam", 
            "moderated t stat", "permutation t-test"), cex=1, counts.col="red")


}

select_de_by_fold_change = function(verbose=F) {

  # load the data
  
  source("load_data.r")
  load_data()

  # create the ratios

  create_ratios()

  # select the DE genes

  gene_color = rep("black", length(mean_expression))

  gene_color[ratios > 2] = "red"

  gene_color[ratios < -2] = "blue"

  # set the names to identify which genes have been chosen

  names(gene_color) = names(ratios)


  # plot the fold change vs the intensity of the gene
  # expression level 

  if(verbose){

    plot(mean_expression, ratios, xlab="mean intensity", 
                     ylab="log(|fold change|)", col=gene_color)
    
    abline(h=c(2, -2), col="purple", lty="dashed")
  
  }

  de_genes = gene_color[gene_color != "black"]

  num_DE_genes = length(de_genes)

  if(verbose){

    cat("\n Number of DE genes: ", num_DE_genes, "\n\n")

    cat ("Press [enter] to see the distribution of ratios")
    line = readline()

    hist_dens(ratios, scale=30)

  }

  assign("fold_change_selected_genes", names(de_genes),
                                       envir=.GlobalEnv)

}

select_de_by_unusual_ratio = function(verbose=F){

  # load the data
  
  source("load_data.r")
  load_data()

  # generate the ratios

  create_ratios()

  # z score the ratios

  ratios = (ratios-mean(ratios))/sd(ratios)

  # select the genes that are 2 std devs away from the mean

  lower_bound = mean(ratios) - (2*sd(ratios))

  upper_bound = mean(ratios) + (2*sd(ratios))

  de_genes = ratios[ratios > upper_bound |
                    ratios < lower_bound]

  if(verbose){

    cat("Number of DE genes: ", length(de_genes), "\n\n")

    cat ("press [enter] to see the distribution of ratios")
    line = readline()

    hist_dens(ratios, scale=30)

    abline(v=lower_bound)
    abline(v=upper_bound)

  }

  assign("unusual_ratio_selected_genes", names(de_genes), 
                                      envir=.GlobalEnv)
  
}


select_de_by_sam = function(delta=1.2, nperms=1000) {

  library(samr)

  training_data = rbind(ALL_train_set[1:11,], AML_train_set)

  outcomes = append(rep(1, 11), rep(2, 11))
  

  data = list(x = t(training_data), y = outcomes, 
                genenames=colnames(training_data), logged2=T)

  samr.obj = samr(data, resp.type="Two class unpaired", nperms=nperms, random.seed=123)

  delta.table = samr.compute.delta.table(samr.obj, 
                                    min.foldchange=1.5)

  de_genes = samr.compute.siggenes.table(
      samr.obj, del=delta, data=data, delta.table=delta.table,
      min.foldchange=1.5)

  up_genes = de_genes$genes.up[, "Gene ID"]

  down_genes = de_genes$genes.lo[, "Gene ID"]

  de_genes = append(up_genes, down_genes)


  assign("sam_selected_genes", de_genes, 
                          envir=.GlobalEnv)

}



# helper function to pre-process data

create_ratios = function(){

  # caluclate the the mean expression levels for all genes
  # in the AML and ALL set

  mean_AML = colMeans(AML_set)

  mean_ALL = colMeans(ALL_set)

  median_AML = apply(AML_set, 2, function(x) median(x))

  median_ALL = apply(ALL_set, 2, function(x) median(x))

  # now calculate the ratio of each gene
  # between the two conditions

  ratios = log2(abs(mean_ALL / mean_AML))

  # calculate the mean expression level of each gene
   
  all_data = rbind(AML_set, ALL_set)

  mean_expression = colMeans(all_data)

  median_expression = apply(all_data, 2, 
                        function(x) median(x))

  assign("ratios", ratios, envir=.GlobalEnv)

  assign("mean_expression", mean_expression, envir=.GlobalEnv)

}

# helper function to plot density of ratios

hist_dens = function(data, num_bin=200, scale=37, xlab="fold change", ylab="density") {

  dens = density(data, n=num_bin)
  hist(data, xlim=range(dens$x), breaks=dens$x, xlab=xlab, 
                    ylab=ylab, freq=F, main=NULL)

  lines(dens$x, dens$y*dens$bw*scale, col="red", lty=2)

}

# helper function to filter genes before hypothesis
# testing 

filter_genes = function(){

  # load the data

  source("load_data.r")
  load_data()
  
  # since we are unable to reference the original values, 
  # we do not know which genes have significant values

  # we will now filter for genes that are variable between 
  # AML and ALL. this will be done by making sure that the
  # max intesnity / min intesnity for each gene is atleast 
  # 1.5

  data = rbind(ALL_set[1:25, ], AML_set)

  min_val = apply(data, 2, function(x) min(x))

  max_val = apply(data, 2, function(x) max(x))

  data = data[,max_val/min_val >= 1.5]  

  cat("\n\nVariable filtered_data now available \n\n\n")

  assign("filtered_data", data, envir=.GlobalEnv)

}
