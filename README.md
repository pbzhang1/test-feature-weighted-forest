# test-feature-weighted-forest

Biological data often exhibit natural clustering and batch effects. In this study, we investigate the impact of ensembling Random Forest learners trained on data clusters derived from a single dataset, accounting for feature distribution heterogeneity and feature mapping heterogeneity.

Our findings reveal that constructing ensembles of forests, trained on machine-learning-determined clusters, with each point weighted according to test dataset feature information, leads to substantial improvements in accuracy and generalizability compared to traditional Random Forest algorithms.

This work is built upon Maya Ramchandran's 'Cross-cluster weighted forest.' My modification involves incorporating information from the features of the test dataset into the weighting of individual data points during the stacked regression step. This adaptation enhances Maya's algorithm, ensuring robust performance even in the presence of feature mapping heterogeneity.





sample run:

load every function and package in the F1simsep.r file. and run the var_levels functions using the following parameters


F1imbalance investigates how imbalance influences the performance of the data

F1highdimension investigates how far away the true clusters affects the performance of algorithm. It uses the average silhouette score implemented in checksilhouette.r for cluster distance. 

F1outcome.r investigates heterogeneity in quadratic outcome 

F1separation.r  investigates the separation of underlying clusters

F1signalstrength.r investigates the signal strength( the norm of coefficient vector)

F1underlyingclusters.r investigates the number of underlying true clusters 

MayaCrossClusterWeightedForest.r 
includes Maya's cross cluster weighted forest 

create_combined.r 
is how we implement our test data weighted forest 

automatedplot.r
generates plots according to the directory of the CSV files 


simsep = 1 
is the separation of underlying true clusters, 1 means very well separated, 0 means not separated at all

ncoef = 20
is the dimensionality of features

nchoose = 10
the number of active features

inbalance = 1
how imbalanced the test data is distributed in the underlying cluster, 1 means a giant centroid, 0 means the test data are evenly distributed among all clucsters

hetero = 1
the heterogeneity of the feature mapping. 0 indicates no feature mapping heterogeneity


k_list
is a list of clusters we choose when we use GMM to cluster the data


ndat = 5
at the simulation step, we have 5 underlying true clusters 
ntest = 1
we only consider one test data set. 

we use random forest to fit and model and predict the outcome
modfit = randomforestfit
modpred =  randomforestpredict




k_list  <- seq(10, 50, 5)

k_list = c(2:5,k_list)


k_cluster = 5
#ndat and ntest is simulation level
ndat = 5
ntest = 1
out_str = "crossclusterk"

k = kcluster
modfit = randomforestfit
modpred =  randomforestpredict


cs <- vary_levels(reps = 50, var_list = k_list, modfit = randomforestfit, modpred
                  = randomforestpredict, ndat = ndat, ntest = ntest,ncoef = ncoef, nchoose = nchoose, out_str = "cs", k = kcluster, hetero = hetero, inbalance = inbalance, simsep = simsep)


Remark:

When comparing how the number of underlying clusters affects the performance of the algorithm we use the oneCentroid simplification.


