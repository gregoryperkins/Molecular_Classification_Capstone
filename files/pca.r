library(ggplot2)
library(gg3D)
library(plotly)

Sys.setenv("plotly_username"="jake_sauter")
Sys.setenv("plotly_api_key"="I0DCLRoNDsgs8MtxqUeZ")

# load in the data

source("load_data.r")
load_data()

# perform pca

gene_pca = prcomp(clean_train_data, scale=T, center=T)

# Use pca$x[,principle compenent num] to access principle componenents for each 

# plot the test data on the pca plot generated from training data

pred = predict(gene_pca, clean_test_data)

# Write the principle components to file to be used later

train_pcs = data.frame(gene_pca$x)
#rownames(train_pcs) = train_accession_nums
write.csv(train_pcs, file="train_principle_components.csv")

test_pcs = data.frame(pred)
#rownames(test_pcs) = test_accession_nums
write.csv(test_pcs, file="test_principle_components.csv")

# Now we need to get the data into the form that the principle
# components can be plotted for each point with different
# colors for the class points

train_group = as.factor(actual_train_classes)

test_group = as.factor(actual_test_classes)

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
