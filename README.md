# Molecular_Classification_Capstone

Project for capstone in Applied Mathematics at SUNY Oswego. Classification of cancer type (ALL/AML) by genomic mircoarray samples (expression levels of genes). Most techniques were applied to data from [Golub (1999)](https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999))

## Week 1 (Cell Biology)

The focus on week 1 was on the introduction to cell biology to have adequate background on the data. 

[Microarray Introduction Presentation](files/Microarrays_Introduction.pdf)


## Week 2 (Reliability and Reproducibility of Microarray Data)

The focus of week 2 was on the reliability and reproducibility of Microarray data.

[Microarray Reliability Presentation](files/Microarrays_Reliability.pdf)


## Week 3 (Multiple Comparisons)

The focus of week 3 was on the problems that arise when performing hypothesis tests on large data sets, with many comparisions happening.

[Muliple Comparisons Presentation](files/Microarrays_Mulitple_Comparisons.pdf)


## Week 4 (Principle Componenent Analysis)

Week 4 is the first we that we are hands-on data. The first step of data analysis on datasets with a large amount of features and a small amount of examples is dimensionality reduction. Principal Compenenet Analysis (PCA) is a form of dimenstionality reduction that aims to keep the covariance structure of the data sound, while reducting as many dimensions as possible.

[PCA Presentation](files/PCA.pdf)

[clean_data.r](clean_data.r) --> [clean_gene_train_data.csv](files/clean_gene_train_data.csv) and [clean_gene_test_data.csv](files/clean_gene_test_data.csv)

[load_data.r](load_data.r)

[pca.r](pca.r) --> [train_principle_components.csv](files/train_principle_components.csv) and [test_principle_components.csv](files/test_principle_components.csv)

Please view the PCA visualizations at the following links: 

[The plotting of training data on PC1, PC2 and PC3](https://plot.ly/~jake_sauter/113)

<img src="/images/cluster_1.png" alt="drawing" width="300"/>


[The plotting of training data on PC2, PC3, and PC4](https://plot.ly/~jake_sauter/115)

<img src="/images/cluster_2.png" alt="drawing" width="300"/>

[The plotting of testing data on PC1, PC2 and PC3 generated from training data](https://plot.ly/~jake_sauter/117)

<img src="/images/cluster_3.png" alt="drawing" width="300"/>

[The plotting of testing data on PC2, PC3 and PC4 generated from training data](https://plot.ly/~jake_sauter/119)

<img src="files/images/cluster_4.png" alt="drawing" width="300"/>


## Week 5 (Cluster Analysis)

Week 5 was all about clustering, many disatance metrics were reviewed as well as a wide range of clustering algorithms. Kmeans, Kmedoids, Hierarchical clustering, Biclustering and other algorithms were covered as well as metrics for determining the goodness of clusters. Simple Clustering algorithms were implemented in R and results were shown. 

[clustering.r](files/clustering.r)

[k_means.r](files/k_means.r)


[generate_cluster_results.r](files/generate_cluster_results.r)

[Cluster Analysis Presentation](files/Clustering.pdf)


## Week 6 (Clustering Continued)

The concern of week 6 was to really get a hold of all of the clustering algorithms. KMeans, being so popular and widely used, was implemented from scratch in R. This also provided the capability of using different distance metrics for KMeans, as this is not an option in the R kmeans() function. KMedoids and Hierarchical clustering was also explored more and implemented though library functions in R on the [Golub (1999) data](https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)).

[custom_kmeans.r](files/custom_kmeans.r)

[k_means_pca_clustering.r](files/kmeans_pca_clustering.r)

[hierarchical_clustring.r](files/hierarchical_clustering.r)

[biclustering_test.r](files/biclustering_test.r)

[test_cluster_results.csv](files/test_cluster_results.csv)

[train_cluster_results.csv](files/train_cluster_results.csv)

[Clustering Continued Presentation](files/Clustering_Continued.pdf)


## Week 7 (Selecting DE Genes)

Selecting differentially expressed (DE) genes is a key step in the analysis of microarray data. Ensuring that only DE genes need to be evaluated later for classification ensures that negligable genes do not increase the problem space dimensionality and greatly reduce the complexity of the problem. Many methods were explored, as well as fold change and unusal ratio techniques being implemented from scratch in R.

[Selecting Differentially Expressed Genes Presentation](files/Selecting_DE_Genes.pdf)

[selecting_DE_genes.r](files/selecting_DE_genes.r)

## Week 8 (Selecting DE Genes Continued)

During the first week of selecting differentially expressed genes between the ALL and AML clinical cancer groups I was not able to implemement many technqiues, so this week surved that purpose. During the course of this week I was able to implement the  Holms and FDR Family Wise Error Rate corrections for t-test, use SAM in R. I was also drawn to the power of permutation testing and wanted to gain more insight and intuition about the process, so I also implemented permutation testing for t-test in R from scratch. 

[Selecting Differenentially Expressed Genes Continued Presentation](files/Selecting_DE_Genes_Continued.pdf)

[selecting_DE_genes.r](files/selecting_DE_genes.r)

[gene_wise_t_test.r](files/gene_wise_t_test.r)

## Week 9 (Selecting DE Genes: Modterated t-Statistic)

The topic of selecting differentially expressed genes seemed very important to the project goal of developing the most robust classifier for these clinical groups, and the book used for this project noted a seemigly powerful technique that was not outlined very well. In this situation my advisor and I chose to take another week to go over this method. To do so, I reviewed the paper where the method originated from and gathered other resources to make sense of the techniques this method implemented. 

In brief the moderated t-statistic makes the assumption that all genes have equal variance, and implements a hierachical model and pooled variance to produce a more reliable estimate for the variance of each gene based on this assumption. I later found this assumption to be questionable, and decided to go with SAM for the final selection of DE genes before PCA for feature construction. 

[Selecting Differentially Expressed Genes: Moderated t-Statistic Presentation](files/Moderated_t_statistic.pdf)

[selecting_DE_genes.r](files/selecting_DE_genes.r)

## Week 10 (Machine Learning)

This week I reviewed many different machine learning algorithms, metrics and concpets, as well as implemented a few. Due to selecting differentially expressed genes before PCA for feature construction, amazing results were obtained using very simple KNN classifier.

[Machine Learning Presentation](files/Machine_Learning.pdf)

[ml.r](files/ml.r)

## Week 11 (Machine Learning Continued)

For the final week of the project I wrapped up loose ends by implementing a few more machine learning implimentations as well as diving deeper into the mathematical intuition of SVMs

[Machine Learning Continued Presentation](files/Machine_Learning_Continued.pdf)

[ml.r](files/ml.r)

## Complete Project 

[Concatenation of All Presentations](files/Capstone_Complete_Presentation.pdf)

[Capstone Talk](files/Capstone_Talk.pdf)

[Capstone Paper](files/Capstone_Final_Paper.pdf)
