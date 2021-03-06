# In order to prepare the data for ML techniques we must
# select differentially expressed gene then perform pca
# on the this set of genes. The first few PCs will be 
# used to represent each patient in the classification 
# of ALL vs AML


library(class)
library(plotly)

Sys.setenv("plotly_username"="jake_sauter")
Sys.setenv("plotly_api_key"="not_public_key")

# first we must load the data 

source("load_data.r")

load_data()

# now we have AML_train_set, AML_test_set
# ALL_train_set, ALL_test_set

# now we will select differentiall expressed
# genes using SAM

#source("selecting_DE_genes.r")

#select_de_by_sam(1.2)

# now we have sam_selected_genes, being the name
# of all de genes selected by sam

de_indices = which(colnames(AML_train_set) %in% sam_selected_genes)

# now we can do pca on the training and testing data

#gene_pca = prcomp(clean_train_data, scale=T, center=T)

gene_pca = prcomp(clean_train_data[, de_indices], scale=T, center=T)

# We will use PCs 2-5 to classify the data, as the first
# PC is just expression level based, which we would 
# expect to vary and carriest no information

train_data = gene_pca$x

train_data = data.frame(train_data[,1:7])

train_labels = as.factor(actual_train_classes)

test_data = predict(gene_pca, clean_test_data[,de_indices])

#test_data = predict(gene_pca, clean_test_data)

test_data = data.frame(test_data[,1:7])

test_labels = as.factor(actual_test_classes)



#------------------KNN------------------#

knn_output = knn(train = train_data[1:3], test = test_data[1:3], cl = as.factor(train_labels), k=5)

table(test_labels, knn_output)


#------------Logistic-Regression-----------#

#train_table = cbind(as.factor(train_labels), train_data)

#colnames(train_table)[1] = "labels"

#model = glm(labels~.,family=binomial(link='logit'),data=train_table)

# By using the type=response setting in predict, we will
# be given probabilites of class membership of each 
# sample

#log_output = predict(model, test_data)

# now we will convert our outputs to classes
# in the usual way for logisitc regession 


#log_output[log_output > 0] = "AML"

#log_output[log_output <= 0] = "ALL"

#table(log_output, test_labels)

#----------Support-Vector-Machine----------#

#library(e1071)


#model = svm(labels ~ ., data=train_table)

#test_table = cbind(as.factor(test_labels), test_data)

#colnames(test_table)[1] = "labels"

#svm_output = predict(model, test_data)

#table(test_labels, svm_output)

# can be used to tune svm params when kernels are
# more complex 

#svm_tune = tune(svm, train.x = train_data, train.y = as.factor(train_labels), kernel="linear", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

#-------------Decision tree-----------------#

#library(tree)

#model = tree(labels ~ ., data.frame(train_table))

#tree_output = predict(model, test_data, type="class")

#table(tree_output, test_labels)

#-------------Random Forest-----------------#

#library(randomForest)

#model = randomForest(x = train_data, y = train_labels, ntree = 1000)

#forest_output = predict(model, data.frame(test_table), type="class")

#table(forest_output, test_labels)



# if we want to plot the PCs, we can do so with 
# plotly

plot = T

if(plot==T){

# factors are limited class type variables that allow
# us to have categorical colors

train_group = as.factor(train_labels)

test_group = as.factor(test_labels)

# add PCs to data frame

test_PC1 = test_data[,1]
test_PC2 = test_data[,2]
test_PC3 = test_data[,3]
test_PC4 = test_data[,4]

df_test = data.frame(test_PC1, test_PC2, test_PC3, test_PC4, test_group)


#Plot the testing data for pc1 pc2 pc3

p <- plot_ly(df_test, x = ~test_PC1, y = ~test_PC2, z = ~test_PC3, color = ~test_group, colors = c('#BF382A', '#0C4B8E')) %>%#
  add_markers()

api_create(p, type="scatter3d", filename="pca_all_genes_1_2_3")

p <- plot_ly(df_test, x = ~test_PC2, y = ~test_PC3, z = ~test_PC4, color = ~test_group, colors = c('#BF382A', '#0C4B8E')) %>%#
  add_markers()

api_create(p, type="scatter3d", filename="pca_all_genes_2_3_4")

}
