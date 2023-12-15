library(tibble)
#we do not have study specific regressions here
library(locfit)
library(fungible)
library(mvtnorm)
library(ClusterR)
#install.packages("ClusterR")
library(cluster)
library(MLmetrics)
#install.packages("MLmetrics")
library(randomForest)
library(glmnet)
#install.packages("glmnet")
library(nnet)
library(nnls)
#install.packages("nnls")

library(doParallel)
#install.packages("doParallel")

#install.packages("BiocGenerics")

library(foreach)
library(plyr)
#library(BiocGenerics)
library(clusterGeneration)
#install.packages("clusterGeneration")
#install.packages("ranger")
library(ranger)
library(wavethresh)
#install.packages("Bioconductor")
#install.packages("BiocGenerics")

library(BiocGenerics)
#install.packages('party')
library(party)
library(Rborist)


setwd("C:/Users/17347/Desktop/newsampling")





rcomb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}
norm_vec <- function(x) sqrt(sum(x^2))

absnorm <- function(vec, max.norm = FALSE){
  sgn <- sign(vec)
  vec <- abs(vec)
  if(max.norm){
    mvec <- max(vec)	
  } else {
    mvec <- min(vec)
  }
  
  vec <- abs(vec - mvec)
  sgn*(vec/sum(vec))
}


#using Rborist package for fast randomforest 

randomforestfit <- function(data, ...){
  md1 <- rfArb( x =  data[,-1], y = data$y , ntree = 200, impPermute = 1)
  

}
randomforestpredict <- function(data, newdata){
  AA = predict(data, newdata=newdata)
  as.vector(AA$yPred)
}

#split the training and test data, get rid of the labels
trainTestDecoupling <- function(studies_list,nstudies,ntest){
  
  
  
  liststudy = studies_list
  
  merged <- do.call(rbind, studies_list[1:length(studies_list)])
  
  testdata <- merged[merged[,2]<0,][-2]
  
  for (i in 1:nstudies){
    data_i = studies_list[[i]] 
    liststudy[[i]] = data_i[data_i[,2]>0,][,-2]
  }
  
  liststudy[[nstudies+1]] = testdata
  return(liststudy)
  
}


sampleID  <-  read.csv("sampleID.csv",skip = 1, header = F)





normal_or_neoplastic  <-  read.csv("normal_or_neoplastic.csv",skip = 1, header = F)


patientlist = unique(sampleID$V2)
typelist = unique(normal_or_neoplastic$V2)
typelist = typelist[1:2]



create_list <- function(patientName,cellType ){
  
  
  
  #patienName is one of the factors in patientlist
  #cellType is one of typelist
  
  topgenes <-  read.csv( "topgenes.csv",skip = 1, header = F)
  clusters  <-  read.csv("cluster.csv",skip = 1, header = F)
  
  
  
  X_pcaNoCCND1  <-  read.csv("X_pcaNoCCND1.csv",skip = 1, header = F)
  
  sampleID  <-  read.csv("sampleID.csv",skip = 1, header = F)
  
  
  
  
  
  normal_or_neoplastic  <-  read.csv("normal_or_neoplastic.csv",skip = 1, header = F)
  
  
  patientlist = unique(sampleID$V2)
  typelist = unique(normal_or_neoplastic$V2)
  
  typelist = typelist[1:2]
  
  
  dims = dim(X_pcaNoCCND1)
  
  N = dims[1]
  
  testlabel = rep(1,N)
  
  #we choose half of specified cells from specified patient as test data 
  
  N2 = length(testlabel[sampleID$V2 == patientName & normal_or_neoplastic$V2 == cellType])
  N3 = floor(N2/2)
  
  labs = which(sampleID$V2 == patientName & normal_or_neoplastic$V2 == cellType)
  testlabels = sample(labs,N3)
  
  
  
  
  #label -1 as test data, label 1 as training data
  testlabel[testlabels ] = -1
  
  studies_list = add_column(topgenes, testlabel, .after = "V2")
  
  drop = c('V1')
  
  studies_list = studies_list[,!(names(studies_list) %in% drop)]
  
  
  dims = dim(studies_list)
  nvals = dims[2]-2
  colnames(studies_list) <- c("y", paste0("V", 1:(nvals+1)))
  
  X_pcaNoCCND1 = X_pcaNoCCND1[,!(names(X_pcaNoCCND1) %in% drop)]
  
  #X_pcaNoCCND1
  
  pca_list = X_pcaNoCCND1
  
  return(list(studies_list = studies_list,pca_list =pca_list ))
  
} 





#temp = X_pcaNoCCND1[,c(1:4)]

create_combined_realdata <- function(studies_list,pca_list, ntest, k){
  #minimum number of training data in a cluster, if smaller than this number, we merge this cluster to nearby clusters
  minTraining = 5
  
  
  merged <- studies_list
  
  
  #cluster without using y, or the labels
  
  
  k3 = GMM(pca_list[,c(1:21)], k, dist_mode  = "eucl_dist", seed_mode = "random_subset", 10, 10)
  
  #training data log-likelihood
  kk3 = k3$Log_likelihood[merged[,2]>0,]
  #numvars = length(kk3)/k
  
  # the clustering step
  #k2 = max.col(k3$Log_likelihood, 'first')
  
  #k2 is label
  testdata = merged[merged[,2]<0,]
  
  h3 = k3$Log_likelihood
  
  dims = dim(merged)
  #number of variables of the data (i.e. dimension of the observation)
  nvals = dims[2]-2
  colnames(merged) <- c("y", paste0("V", 1:(nvals+1)))
  
  #each column of l3 is the order
  #l3[i,ind], ind = ind of variables
  l3 = apply(-h3,1,order)
  L3 = t(l3)
  l4 = data.frame(L3)
  colnames(l4) <- c(paste0("C", 1:k))
  
  
  #aa = l4[[5]]
  #h3 = k3$Log_likelihood
  h4 = data.frame(h3)
  colnames(h4) <- c(paste0("L", 1:k))
  
  #we put the data together with its clusters
  m2 = cbind(merged,l4)
  #m2 = cbind(merged,h4)
  
  nrows =nrow(m2)
  
  
  #this is very important, it changes the row names of m2 to 1:length, so that it is trackable
  row.names(m2) =  as.character( c(1:nrows))
  
  #training data belongs to cluster1 of GMM algorithm
  
  deletedclusters = c(0)
  
  finishMerging = 0
  
  while (finishMerging == 0){
    finishMerging = 1
    
    for (i in 1:k){
      
      if (i %in% deletedclusters){
        next
      }else{
        #number of training data in 
        a3 = nrow(m2[which(m2$V1 >0 & m2$C1== i),])
        #print(a3)
        if (a3 < minTraining){
          #then we merge this cluster
          deletedclusters = append(deletedclusters,i)
          
          #temp = m2[which(m2$V1 >0 & m2$C1== i),]
          
          
          # We first find each data that belongs to cluster i
          
          
          #i changed this m2[which(m2$V1 >0 & m2$C1== i),]
          trainMove = rownames(m2[which(m2$C1== i),])
          temp4 = as.integer(trainMove)
          len1 = length(temp4)
          for (j in 1:len1){
            #the index of the data
            ind = temp4[j]
            
            #the first index of cluster 
            clustInd = nvals+2+1
            
            
            #then we put it in the nearest cluster that is not deleted
            while (m2[ind,clustInd] %in% deletedclusters){
              clustInd = clustInd +1
            }
            
            m2[ind,nvals+2+1] = m2[ind,clustInd]
            
            
            #check
            if( m2[ind,]$C1 !=  m2[ind,nvals+2+1]){
              print("wrong!")
            }
          }
          finishMerging = 0
        }
      }
    }
  }
  
  
  
  #minClusterSize = 5
  
  k2 = m2$C1
  
  
  
  clusters_list <- lapply(split(seq_along(k2), k2), 
                          function(m, ind) m[ind,], m = merged)
  
  #clusters_list <- lapply(split(seq_along(k2), k2), 
  #                        function(m, ind) m[ind,], m = m2)
  
  
  #we do not need to permutate the vectors 
  #clusters_list <- lapply(split(seq_along(k2), k2), 
  #                        function(m, ind) m[ind,], m = merged)[order(unique(k2))]
  
  len_clust <- sapply(clusters_list, function(i) nrow(i))
  
  #exclude clusters less than 2 pts
  wlc <- which(len_clust <= 4)
  if (any(wlc)){
    clusters_list <- clusters_list[-wlc]
    #kk3 <- kk3[,-wlc]
  }
  
  
  nclusters <- length(clusters_list)
  
  
  # we only leave the log-likelihood of those that matters
  deletedclusters = as.integer(deletedclusters)
  deletedclusters = deletedclusters[deletedclusters != 0]
  nntest = rep(1,k) 
  nntrain = rep(1,k)
  for (i in 1:  k){
    nntest[i] = nrow(m2[which(m2$V1 <0 & m2$C1== i),])
    #nntest[i] = nrow(clsti[clsti[,2]<0,])
    #need modification
    nntrain[i] = nrow(m2[which(m2$V1 >0 & m2$C1== i),])
    
  }
  
  if (any(deletedclusters)){
    h5 = h4[,-deletedclusters]
    nntrain = nntrain[-deletedclusters]
    nntest = nntest[-deletedclusters]
  } else{
    h5 = h4
  }
  
  
  #we record the likelihood
  
  m1 = cbind(merged,h5)
  
  Dims = dim(m1)
  ncols = Dims[2]
  
  kk5 = lapply(split(seq_along(k2), k2), 
               function(m, ind) m[ind,], m = m1)
  
  kk6 = kk5
  for (i in 1:nclusters){
    temp = kk5[[i]]
    temp2 = temp[which(temp$V1>0),]
    kk6[[i]] = temp2[,(nvals+2+1):ncols]
    
    if (nrow(kk6[[i]])!=nntrain[i] ){
      print("wrong loglikelihood")
    }
    
    if (ncol(kk6[[i]])!= nclusters ){
      print("wrong loglikelihood, again")
    }
    
    
  }
  
  kk7 = do.call(rbind, kk6[1:nclusters])
  # #whether the algorithm give the correct number of clusters
  # ktest = length(deletedclusters)+nclusters
  # if (ktest != k){
  #   print("wrong")
  # }else{
  #   print("right")
  # }
  
  
  #the test data set
  #  for (t in 1:ntest){
  #    clusters_list[[(nclusters + t)]] <- testdata
  #}
  
  #throws away entire clusters, including the test data 
  s1 = trainTestDecoupling(clusters_list,nclusters,ntest)
  #len_clust2 <- sapply(s1, function(i) nrow(i))
  
  #kk3 is training data loglikelihood
  #we want the test data loglikelihood
  return(list(clusters_list = s1,nntest = nntest,nntrain = nntrain,gmmfit = kk7))
}





#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating sub-studies, 3 if just running on the original set of simulated studies
clusters_fit_realdata <- function(modfit, modpred, ntest, studies_list, pca_list,cluster_ind, nsep){
  
  #The studies list have already been decoupled here for ind = 1,2,3
  if (cluster_ind == 1){ #K-means clustering
    cc <- create_clusters(studies_list,pca_list, ntest, k)
    edat <- cc$clusters_list
    len_clust <- sapply(edat, function(i) nrow(i))
    
  }else if (cluster_ind == 4){ #combining the clustering with test data sets
    cc <- create_combined_realdata(studies_list, pca_list,ntest, nsep)
    edat <- cc$clusters_list
    #edat <- studies_list
    numtest <- cc$nntest
    gmmfit <- cc$gmmfit
    len_clust <- sapply(edat, function(i) nrow(i))
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  # wlc2 <- which(len_clust < 3)
  # ntrain3  = ntrain-length(wlc2) 
  # trainingclusters2 = c(1:ntrain)
  # 
  # if (any(wlc2)){
  #   trainingData2 <- edat[-wlc2]
  #   trainingclusters2 = trainingclusters2[-wlc2]http://127.0.0.1:44767/graphics/plot_zoom_png?width=1200&height=900
  #   
  # }else{
  #   trainingData2 <- edat
  # }
  # 
  # 
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  #mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  #mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  data1 = as.data.frame(matstack[sample(nrow(matstack)), ])
  mod0 <- rfArb( x =  data1[,-1], y = data1$y , ntree = 500, impPermute = 1)
  
  
  
  #this is the percentage increase in MSE
  #aa = mod0$importance[,1]
  
  #m1 <- cforest(y~.,  data =as.data.frame(matstack[sample(nrow(matstack)), ]),control=   cforest_unbiased(mtry=5,ntree=500))
  
  mods = vector("list", ntrain)
  if (cluster_ind == 4){
    #used to train models
    #trainingData = edat
    # wlc <- which(len_clust <= 5)
    # ntrain2  = ntrain-length(wlc) 
    #trainingclusters = c(1:ntrain)
    
    
    # if (any(wlc)){
    #   trainingData <- trainingData[-wlc]
    #   trainingclusters = trainingclusters[-wlc]
    #   
    # }
    
    
    for (j in 1:ntrain){
      mods[[j]] <- modfit(edat[[j]])
      
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[j]], newdata = x[, -1])})
      
      
      curpreds <- lapply(preds, as.numeric)
      allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
      
    }
    
  }else{
    
    for (j in 1:ntrain){
      mods[[j]] <- modfit(edat[[j]])
      
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[j]], newdata = x[, -1]) 
      })
      #mses[j,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
      #  newdata = x[, -1]
      #  preds <- modpred(mods[[j]], newdata = newdata) 
      #  mean((preds - x[,"y"])^2)}
      #))
      
      curpreds <- lapply(preds, as.numeric)
      allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
    }
    
    
  }
  
  
  #diag(mses) <- NA
  
  # CS Weights
  
  #tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  #weights <- absnorm(sqrt(tt), max.norm = TRUE)
  #nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  #nwts <- absnorm(nk[1:ntrain])
  
  predstack <- do.call(rbind, allpreds)
  
  len_total = length(matstack$y)
  
  w = rep(1,len_total)
  
  #specific treatment
  if (cluster_ind == 4){
    loglike = gmmfit
    ll = apply(loglike, 1, max)
    loglike = loglike-ll
    likelihood = exp(loglike)
    normll = apply(likelihood,1,norm_vec)
    likelihood = as.matrix(likelihood)
    likelihood = likelihood/normll
    
    numtest = numtest/sum(numtest)
    #weights 
    
    w = likelihood %*% numtest
    
    #number of data points that are weighted significantly in weighted regression
    #sum(w[w>0.1])
  }else{
    w = rep(1,len_total)
  }
  
  
  # Regression: stacked (intercept and no intercept)
  
  X = as.numeric(sqrt(w)) * predstack
  y = as.numeric(sqrt(w)) * as.numeric(as.character(matstack$y))
  
  
  X2 = as.numeric(sqrt(w)) * cbind(rep(1,nrow(predstack)),predstack)
  #coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  #coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  
  coefs_stack_noint <- nnls::nnls(X, y)$x
  coefs_stack_int <- nnls::nnls(X2, y)$x
  
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), weights = w, alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), weights = w, alpha = 0, lower.limits = 0, intercept = T)))
  
  
  
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0
  
  coefs_stack_noint_norm <- absnorm(coefs_stack_noint)
  
  
  
  outmat <- matrix(NA, ntest, 7)
  colnames(outmat) <- c("Merged", 
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "Stack_lasso", "Stack_ridge","just_t_clus")
  
  #edat[[25]][,-1]
  
  
  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1])
    
    merged <- as.vector(sapply(merged, as.numeric))
    allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
    allmod <- apply(allmod, 2, as.numeric)
    # Unweighted average
    #unweighted <- colMeans(allmod)
    
    # sample size weighted
    #sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})
    # cross-study weighted 
    #cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})
    
    # regression: stacked (noint, int, each normed) + lasso
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    
    #sample_wtd<-stack_noint
    #cs_wtd <-sample_wtd
    # # regression: study_specific (noint, int, noint normed) + lasso
    # ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
    # ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
    # ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})
    # ss_lasso <- apply(allmod, 2, function(x){coefs_ss_lasso[1] + sum(coefs_ss_lasso[-1]*x)})
    # ss_ridge <- apply(allmod, 2, function(x){coefs_ss_ridge[1] + sum(coefs_ss_ridge[-1]*x)})
    # 
    # 
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    
    testclus = which.max(numtest)
    
    just_t_clus =  modpred(mods[[testclus]], newdata = edat[[i]][,-1])
    
    
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2),
                                  mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                  mean((cury - stack_lasso)^2), 
                                  mean((cury - stack_ridge)^2),mean((cury - just_t_clus)^2)))
    
    mergederror = mean((cury - merged)^2)
  }
  
  
  #calculate the importance, weighted by the stacked regression weights 
  importanceStack = rep(0,ncol(edat[[1]])-1)
  
  for (j in 1:ntrain){
    temp1 = mods[[j]]$importance$mse*abs(coefs_stack_ridge[j])
    importanceStack = importanceStack+temp1
    
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  
  
  
  return(list(outmat = colMeans(outmat), mses = mses,allmod = allmod,mergederror = mergederror,importanceStack=importanceStack))
}

#which.max(nntest)
patientName = "SMM-9"

cellType = typelist[1]
cluster_ind = 4
nsep = 25
ntest = 1
modfit = randomforestfit
modpred =  randomforestpredict

aa = create_list(patientName,cellType)

studies_list= aa$studies_list
pca_list = aa$pca_list
k = 25

errors_combined_fast <-clusters_fit_realdata(modfit, modpred, ntest, studies_list, pca_list, cluster_ind = 4, nsep = k)


m1 = 0



