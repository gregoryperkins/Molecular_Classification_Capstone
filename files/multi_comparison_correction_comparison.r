# This script will generate a plot of the behavior
# of different multiple comparison correction techniques

source("load_data.r")

load_data()

num_genes = dim(AML_set)[2]

for (i in 1:num_genes) {

  p_vals = t.test(AML_set[,1], ALL_set[,1])

}


alpha = 0.05

bonf_vals = rep(alpha/num_genes, num_genes)

holms_vals = rep(1, num_genes)

fdr_vals  = rep(1, num_genes)


for (i in 1:num_genes) {

  holms_vals[i] = (1/(num_genes - i + 1)) * alpha

  fdr_vals[i] = (i/num_genes) * alpha

}

plot(bonf_vals, type="o")

lines(holms_vals)

lines(fdr_vals)

