#!/usr/bin/env Rscript

# This script will be used to produce a t-test and order 
# the genes by significance of p-value

# First we will randomly generate a set of 5 tumor samples
# and 5 healthy tissue samples, with 20 genes each

set.seed(100)
ngenes = 20
ncases = 5
rownames = paste("g",1:ngenes,sep="")
colnames = c(paste("T",1:ncases,sep=""), paste("C",1:ncases,sep=""))
l = list(rownames,colnames)
data = matrix(rnorm(ngenes*2*ncases), nrow=ngenes, ncol=2*ncases, dimnames=l)

data 

# Now we will perform a t-test

labels = c(rep("tumor", ncases), rep("control", ncases))

mytest <- function(x){
  tumors = x[labels=="tumor"]
  controls = x[labels=="control"]
  p = t.test(tumors, controls, var.equal=TRUE)$p.value
  foldchange = mean(tumors) - mean(controls)
  return(c(foldchange,p))
}

results = t(apply(data,1,mytest))
rownames(results)<-rownames(data)
colnames(results)<-c("logFC", "p")
results<-results[order(results[,2]),]

results
