library(ggplot2)
library(gg3D)
library(plotly)

Sys.setenv("plotly_username"="jake_sauter")
Sys.setenv("plotly_api_key"="vHqi1mDjclWV7xwrpZGa")

# load in the data
gene_train_data = read.csv(file="gene-expression-data/data_set_ALL_AML_train.csv", header=FALSE, sep=",",stringsAsFactors = FALSE, quote = "")

gene_test_data = read.csv(file="gene-expression-data/data_set_ALL_AML_independent.csv", header=FALSE, sep=",",stringsAsFactors = FALSE, quote = "")


# Hold onto the feature names

train_feature_names = gene_train_data[2:nrow(gene_train_data),2]

test_feature_names = gene_test_data[2:nrow(gene_test_data),2]


# clean the data by removing all cols with header "call"

clean_train_data = gene_train_data[,-which(gene_train_data 
== "call", arr.ind=T)[,2]]

clean_test_data = gene_test_data[,-which(gene_test_data == "call", arr.ind=T)[,2]]

# Hold onto the accession numbers

train_accession_nums = clean_train_data[1,3:ncol(clean_train_data)]

test_accession_nums = clean_test_data[1,3:ncol(clean_test_data)]

# clean the data into just features

clean_train_data = clean_train_data[2:nrow(clean_train_data),3:ncol(clean_train_data)]

clean_test_data = clean_test_data[2:nrow(clean_test_data),3:ncol(clean_test_data)]


# read in the actual class labels

actual_class = read.table("gene-expression-data/actual.csv", header=TRUE, sep=",")

actual_class = t(actual_class[2])

# Form the data into a data frame

clean_train_data = data.frame(t(clean_train_data))

clean_test_data = data.frame(t(clean_test_data))

# Name the columns appropriately so pca pred will know
# what features to use

colnames(clean_train_data) = train_feature_names

colnames(clean_test_data) = test_feature_names

# perform pca

gene_pca = prcomp(clean_train_data, scale=T, center=T)

# Use pca$x[,principle compenent num] to access principle componenents for each 

# plot the test data on the pca plot generated from training data

pred = predict(gene_pca, clean_test_data)


# Now we need to get the data into the form that the principle
# components can be plotted for each point with different
# colors for the class points

train_group = c()

for(index in seq(length(clean_train_data[,1]))) {
  if(actual_class[as.numeric(train_accession_nums[index])] == "ALL") {
    train_group = append("blue", train_group) } else {
    train_group = append("red", train_group)}
}

test_group = c()

for(index in seq(length(clean_test_data[,1]))) {
  if(actual_class[as.numeric(test_accession_nums[index])] == "ALL") {
    test_group = append("blue", test_group) } else {
    test_group = append("red", test_group)}
}

# factors are limited class type variables that allow
# us to have categorical colors

train_group = as.factor(train_group)

test_group = as.factor(test_group)

# add PCs to data frame

train_PC1 = gene_pca$x[,1]
train_PC2 = gene_pca$x[,2]
train_PC3 = gene_pca$x[,3]
train_PC4 = gene_pca$x[,4]

test_PC1 = pred[,1]
test_PC2 = pred[,2]
test_PC3 = pred[,3]
test_PC4 = pred[,4]

df_train = data.frame(train_PC1, train_PC2, train_PC3, train_PC4, train_group)

df_test = data.frame(test_PC1, test_PC2, test_PC3, test_PC4, test_group)

# Plot the training data for pc1 pc2 pc3
p <- plot_ly(df_train, x = ~train_PC1, y = ~train_PC2, z = ~train_PC3, color = ~train_group, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers()

api_create(p, type="scatter3d", filename="pca_train_pc1_pc2_pc3")


# Plot the training data for pc2 pc3 pc4

p <- plot_ly(df_train, x = ~train_PC2, y = ~train_PC3, z = ~train_PC4, color = ~train_group, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers()

api_create(p, type="scatter3d", filename="pca_train_pc2_pc3_pc4")

# plot the predicted for pc1 pc2 PC3

p <- plot_ly(df_test, x = ~test_PC1, y = ~test_PC2, z = ~test_PC3, color = ~test_group, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers()

api_create(p, type="scatter3d", filename="pca_test_pc1_pc2_pc3")

# plot the predicted for pc2 pc3 PC4

p <- plot_ly(df_test, x = ~test_PC2, y = ~test_PC3, z = ~test_PC4, color = ~test_group, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers()

api_create(p, type="scatter3d", filename="pca_test_pc2_pc3_pc4")





#-------------Useful Bits of Code----------------#

# write this new cleaned data to file

#write.table(clean_train_data, "clean_gene_train_data.csv", sep=",", row.names=F, col.names=train_accession_nums)

#write.table(clean_test_data, "clean_gene_test_data.csv", sep=",", row.names=F, col.names=test_accession_nums)

# read in the cleaned data file

#clean_train_data = read.csv("clean_gene_train_data.csv", header=F)

#clean_test_data = read.csv("clean_gene_test_data.csv", header=F)

# extracting the accession number that identifies
# the patient that the gene came from

# train_accession_nums = clean_train_irdata[1,]

# test_accession_nums = clean_test_data[1,]

# clean_train_data = clean_train_data[2:nrow(clean_train_data),]

# clean_test_data = clean_test_data[2:nrow(clean_test_data),]

# We want to identify which gene is differentialbe
# between the cases, so we need col vectors of
# gene expression levels


# 2d plot
#ggplot(data = d, aes_string(x = "PC2", y = "PC3", col="group")) + geom_point(size = 3)

# 3d plot

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col="group")) +
#  theme_void() +
#  axes_3D() +
# stat_3D()

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  theme_void() +
#  axes_3D(theta=30) +
# stat_3D(theta=30)library(ggplot2)


# 2d plot
#ggplot(data = d, aes_string(x = "PC2", y = "PC3", col="group")) + geom_point(size = 3)

# 3d plot

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col="group")) +
#  theme_void() +
#  axes_3D() +
# stat_3D()

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  theme_void() +
#  axes_3D(theta=30) +
# stat_3D(theta=30)

# We can do many angles

#g_angle1 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group")) +
#  axes_3D(theta=100) +
#  stat_3D(theta=100)
#
#g_angle2 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D() +
#  stat_3D()
#
#g_angle3 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D(theta=170) +
#  stat_3D(theta=170)
#
#leg = get_legend(g_angle1)
#no_leg=theme(legend.position = "none")
#plot_grid(
#  g_angle1+theme_void()+no_leg,
#  g_angle2+theme_void()+no_leg,
#  g_angle3+theme_void()+no_leg,
#  leg, ncol=4)
#
# 3D interactive plotly chart


#plot_ly(x=gene_pca$x[, 1], y=gene_pca$x[, 2], z=gene_pca$x[, 3], type="scatter3d", mode="markers", color=group)



# We can do many angles

#g_angle1 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group")) +
#  axes_3D(theta=100) +
#  stat_3D(theta=100)
#
#g_angle2 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D() +
#  stat_3D()
#
#g_angle3 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D(theta=170) +
#  stat_3D(theta=170)
#
#leg = get_legend(g_angle1)
#no_leg=theme(legend.position = "none")
#plot_grid(
#  g_angle1+theme_void()+no_leg,
#  g_angle2+theme_void()+no_leg,
#  g_angle3+theme_void()+no_leg,
#  leg, ncol=4)
#

# 2d plot
#ggplot(data = d, aes_string(x = "PC2", y = "PC3", col="group")) + geom_point(size = 3)

# 3d plot

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col="group")) +
#  theme_void() +
#  axes_3D() +
# stat_3D()

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  theme_void() +
#  axes_3D(theta=30) +
# stat_3D(theta=30)library(ggplot2)


# 2d plot
#ggplot(data = d, aes_string(x = "PC2", y = "PC3", col="group")) + geom_point(size = 3)

# 3d plot

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col="group")) +
#  theme_void() +
#  axes_3D() +
# stat_3D()

#ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  theme_void() +
#  axes_3D(theta=30) +
# stat_3D(theta=30)

# We can do many angles

#g_angle1 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group")) +
#  axes_3D(theta=100) +
#  stat_3D(theta=100)
#
#g_angle2 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D() +
#  stat_3D()
#
#g_angle3 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D(theta=170) +
#  stat_3D(theta=170)
#
#leg = get_legend(g_angle1)
#no_leg=theme(legend.position = "none")
#plot_grid(
#  g_angle1+theme_void()+no_leg,
#  g_angle2+theme_void()+no_leg,
#  g_angle3+theme_void()+no_leg,
#  leg, ncol=4)
#
# 3D interactive plotly chart


#plot_ly(x=gene_pca$x[, 1], y=gene_pca$x[, 2], z=gene_pca$x[, 3], type="scatter3d", mode="markers", color=group)



# We can do many angles

#g_angle1 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group")) +
#  axes_3D(theta=100) +
#  stat_3D(theta=100)
#
#g_angle2 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D() +
#  stat_3D()
#
#g_angle3 = ggplot(d, aes(x=PC1, y=PC2, z=PC3, col=group)) +
#  axes_3D(theta=170) +
#  stat_3D(theta=170)
#
#leg = get_legend(g_angle1)
#no_leg=theme(legend.position = "none")
#plot_grid(
#  g_angle1+theme_void()+no_leg,
#  g_angle2+theme_void()+no_leg,
#  g_angle3+theme_void()+no_leg,
#  leg, ncol=4)
#
# 3D interactive plotly chart

# Plot the training data for pc1 pc2 pc3

#plot_ly(x=gene_pca$x[, 1], y=gene_pca$x[, 2], z=gene_pca$x[, 3], type="scatter3d", mode="markers", color=group)
