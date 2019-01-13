# This script will implement the moderated t statistic
# approach in R

# First we need the limma (Linear Models for MicroArrays)

select_de_by_mod_t_stat = function() {

  library(limma)

  # Load in data

  source("load_data.r")
  load_data()
  data = t(rbind(ALL_set[1:25,], AML_set))

  # Form the Linear Model of 2 groups, with 3 elements in 
  # each group, with the appropriate labels

  lin_mod = gl(2,25,labels=c("control", "treatment"))

  # Form the design matrix, with a treatment-contrast
  # parameterization

  design_mat = model.matrix(~0+lin_mod)

  colnames(design_mat) = levels(lin_mod)

  cat("Design matrix: \n")

  print(design_mat)

  # Use the lmFit function to estimate the parameters of 
  # the linear model

  fit = lmFit(data, design_mat)

  # Make the contrast 

  cont.mat = makeContrasts(contrast=treatment-control, 
                       levels=design_mat)

  cat("contrast Matrix: \n")

  print(cont.mat)

  fit = contrasts.fit(fit, cont.mat)

  # Build the hierarchical Bayes model and estimate the 
  # moderated statistics and associated p-values as evidence 
  # of differential expression

  fit = eBayes(fit)

  results = decideTests(fit)

  # Choose all genes that have an adjusted pvalue less than .05

  de_genes = topTable(fit, coef="contrast", number=nrow(data), 
                    adjust="fdr", p.value=.05, lfc=log2(2), 
                    sort.by = "logFC")

  de_genes = rownames(de_genes)


  assign("moderated_t_stat_selected_genes", de_genes, 
                                      envir=.GlobalEnv)



}
