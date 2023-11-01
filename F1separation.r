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
#install.packages('wavethresh')

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
randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}
randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=newdata))
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


#####
# Setup 1: Simulating data using the monte function - non-gaussian clusters
#more inbalance means test distribution is more skewed
#simsep how much the cluster are seperated at generation step



findTestAllocation <- function(nstudies, inbalance,iters = 40000){
  #find the test data proportion in each simulation cluster according to given inbalance
  n = nstudies
  k = sqrt(n)
  mm = 1
  a = rep(0,iters)
  for (i in 1:iters){
    x = runif(n)
    x = x/sum(x)
    #this ensures x sum up to 1
    #a[i] has range 0 to 1
    
    a[i] = sd(x)*k
    if (abs(a[i]-inbalance)<mm){
      y =x
      mm = abs(a[i]-inbalance)
    }
  }
  return (y)
  
}

sim_data <- function(nstudies, ncoef,nchoose = 10, ntest,hetero,inbalance,simsep){
  studies_list  <- vector("list", nstudies)
  #uncomment if want no sparsity:
  #nchoose <- ncoef
  
  oneCentroid = 1
  
  
  #linear term heterogeniety
  
  
  #general predictor-outcome rule:
  coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  vars <- sample(1:ncoef, nchoose)
  
  #the function g, we make the coefficients the same magnitude as f
  vars2 <-sample(1:ncoef, nchoose)
  coefs2 <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  coefs2 <- coefs2*(mean(abs(coefs)))/mean(abs(coefs2))
  #vars = vars+1
  icoefs <- c(3.5, 0.7)
  #If want to normalize coefficients, uncomment the following
  #norm.betas <- norm_vec(c(coefs,icoefs))
  #coefs <- coefs/norm.betas
  #cicoefs <- icoefs/norm.betas
  
  #cormat <- matrix(.80, ncoef, ncoef)
  #diag(cormat) <- rep(1, ncoef)
  #in.cor.list <- replicate(nstudies, cormat, simplify=FALSE)
  
  #separation 
  
  #how much the cluster are seperated at generation step
  #sep = sep1*simsep
  sep = runif(ncoef, 0.9, 1.1)*simsep
  minV = 0
  maxV = 0.95
  sep = sapply(sep, function(y) min(max(y,minV),maxV))
  
  ## one big centroid
  
  
  
  
  
  #  nlength = length(m$data[,1])
  
  #  ntestdata = floor(nlength/5)
  
  #  testcoeff = sample(1:nlength,ntestdata)
  #  testdata = m$data[testcoeff,]
  cluSize1 = floor(runif(nstudies, 600, 700))
  totalsize = sum(cluSize1)
  
  
  #randomly select data from each cluster for the test data 
  
  
  
  
  
  #more inbalance means test distribution is more skewed
  
  numtests = findTestAllocation(nstudies,inbalance,iters = 40000)
  #numtests = numtests/sum(numtests)
  n2 = runif(1,0.95,1.05)
  
  #roughly 20% of data is test data 
  studyportion = 0.2
  
  testsize = totalsize*studyportion*n2
  
  
  # this decides whether we only have test data in 1 cluster
  if (oneCentroid == 1 ){
    cent = which.max(numtests)
  }
  
  
  
  numtests = floor(numtests*testsize)
  
  
  
  porportion = runif(nstudies,1,10)
  
  temp1 = as.integer(1/studyportion)
  
  
  #here we decide whats the porportion of testdata in each cluster
  #we let 
  porportion = porportion/sum(porportion)*(temp1-2)
  clus_size = floor((porportion+2)*numtests+runif(nstudies,100,500))
  
  #each cluster minimum 200 points
  clus_size = sapply(clus_size,function(y) max(y,400))
  
  #print(sum(clus_size)/sum(numtests))
  #print(numtests)
  m <- monte(seed = sample(1:1000, 1), nvar = ncoef, nclus = nstudies,
             clus.size = clus_size, eta2 = sep,
             cor.list = NULL, random.cor = FALSE, skew.list = NULL,
             kurt.list = NULL, secor = NULL, compactness = NULL,
             sortMeans = FALSE)
  
  
  gamma = runif(nstudies,0,3)
  gamma = scale(gamma)
  
  
  if (oneCentroid == 1 ){
    cent = which.max(numtests)
  }
  
  for (i in 1:nstudies){
    
    #linear term heterogeneity
    #curcoefs <- sapply(coefs, function(x){runif(1, x -hetero, x + hetero)})
    
    curcoefs <- coefs
    #quadratic term heterogeneity
    #icoefs <- sapply(icoefs, function(x){runif(1, x - .7, x + .7)})
    
    #scale means center the column vectors to their mean(shifting the coordinate system
    
    
    #cut off the labels
    #data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ][,-1]))
    
    
    #data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ]))
    
    #scaling it
    #data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ][,-1]))
    
    data_i <- as.data.frame(m$data[m$data[,1] == i, ][,-1])
    
    #not scaling
    data_j<- as.data.frame(m$data[m$data[,1] == i, ])
    #data_j <- scale(as.data.frame(m$data[m$data[,1] == i, ]))
    
    #To binarize the data for part b, uncomment the following line:
    #data_i <- sapply(1:ncol(data_i), function(x) {ifelse(data_i[,x] <= quantile(data_i[,x], runif(1, .2, .8)), 0, 1)})
    
    #generate outcome
    #scaled here, but the original variables are left unscaled for the clustering step
    #baseline: linear model
    
    #we didnt cutoff the labels here
    y2 <- gamma[i]*hetero*(as.matrix((data_i[,vars2])) %*% coefs2)
    
    
    y <- as.matrix((data_i[,vars])) %*% curcoefs +
      #quadratic:
      icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      #interactions:
      # + icoefs[1]*data_i[,vars[1]]*data_i[,vars[2]]
      #- icoefs[2]*data_i[,vars[1]]*data_i[,vars[3]] +
      y2+
      #add noise and hetero term
      cbind(rnorm(nrow(data_i))/3)
    
    
    #testdata 
    
    
    scalednumtests = min(clus_size[i],numtests[i])
    
    
    if (oneCentroid == 1){
      if (i != cent){
        scalednumtests = 0
      }
    }
    
    if (scalednumtests >0){
      testlabels = sample(1:clus_size[i],scalednumtests)
      data_j[testlabels,1] = -i
      data_j[-testlabels,1] = i
    }else{
      data_j[,1] = i
    }
    
    
    
    #??????? why you scale the data?
    
    #studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    
    studies_list[[i]] <- as.data.frame(cbind(y, data_j))
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:(ncoef+1)))
    
  }
  
  
  return(list(studies_list = studies_list))
}


#nstudies is the number of true clusters
#ncoef is the sparsity of the coefficients
#ntest is the unused number of test data points?


#####




create_combined <- function(studies_list, ntest, k){
  #minimum number of training data in a cluster, if smaller than this number, we merge this cluster to nearby clusters
  minTraining = 5
  
  
  merged <- do.call(rbind, studies_list[1:length(studies_list)])
  merged <- merged[sample(nrow(merged)), ]
  
  #cluster without using y, or the labels
  
  #k2 <- kmeans(merged[,-c(1,2)], centers = k, nstart = 25)$cluster
  
  #clustering without using y and the true cluster membership 
  k3 = GMM(merged[,-c(1,2)], k, dist_mode  = "eucl_dist", seed_mode = "random_subset", 10, 10)
  
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
  
  #each column of l3 is the order
  #l3[i,ind], ind = ind of variables
  l3 = apply(-h3,1,order)
  L3 = t(l3)
  l4 = data.frame(L3)
  colnames(l4) <- c(paste0("C", 1:k))
  
  
  #aa = l4[[5]]
  h3 = k3$Log_likelihood
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
        print(a3)
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

#old way of clustering



create_clusters <- function(studies_list, ntest, k){
  
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  k2 <- kmeans(merged[,-1], centers = k, nstart = 25)$cluster
  clusters_list <- lapply(split(seq_along(k2), k2), #split indices by a
                          function(m, ind) m[ind,], m = merged)[order(unique(k2))]
  len_clust <- sapply(clusters_list, function(i) nrow(i))
  
  
  
  wlc <- which(len_clust <= 2)
  if (any(wlc)){
    clusters_list <- clusters_list[-wlc]
  }
  nclusters <- length(clusters_list)
  for (t in 1:ntest){
    clusters_list[[(nclusters + t)]] <- studies_list[[(length(studies_list) - ntest + t)]]
  }
  return(list(clusters_list = clusters_list))
}






#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating sub-studies, 3 if just running on the original set of simulated studies
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep){
  
  #The studies list have already been decoupled here for ind = 1,2,3
  if (cluster_ind == 1){ #K-means clustering
    cc <- create_clusters(studies_list, ntest, nsep)
    edat <- cc$clusters_list
    len_clust <- sapply(edat, function(i) nrow(i))
    
  }else if (cluster_ind == 4){ #combining the clustering with test data sets
    cc <- create_combined(studies_list, ntest, nsep)
    edat <- cc$clusters_list
    #edat <- studies_list
    numtest <- cc$nntest
    gmmfit <- cc$gmmfit
    len_clust <- sapply(edat, function(i) nrow(i))
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  

  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  #mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  
  
  mods = vector("list", ntrain)
  if (cluster_ind == 4){

    
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

      
      
      curpreds <- lapply(preds, as.numeric)
      allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
    }
    
    
  }

  
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
    

  }
  
  
  # Regression: stacked (intercept and no intercept)
  
  X = as.numeric(sqrt(w)) * predstack
  y = as.numeric(sqrt(w)) * as.numeric(as.character(matstack$y))
  
  
  X2 = as.numeric(sqrt(w)) * cbind(rep(1,nrow(predstack)),predstack)

  
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
  
  
  
  outmat <- matrix(NA, ntest, 6)
  colnames(outmat) <- c("Merged", 
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "Stack_lasso", "Stack_ridge")
  
  
  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1])
    
    merged <- as.vector(sapply(merged, as.numeric))
    allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
    allmod <- apply(allmod, 2, as.numeric)

    
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    

    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2),
                                  mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                  mean((cury - stack_lasso)^2), 
                                  mean((cury - stack_ridge)^2)))
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  
  
  
  return(list(outmat = colMeans(outmat), mses = mses,allmod = allmod))
}


rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef,nchoose, k,hetero,inbalance,simsep){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  num.threads <- as.integer(10)
  threads <- makeCluster(num.threads, outfile=logfile, setup_timeout = 0.5)
  registerDoParallel(threads)
  
  getDoParWorkers()
  #######
  results <- foreach(i = 1:reps, .combine = 'rcomb', .multicombine = TRUE, .export = ls(globalenv())) %dopar% {
    library(locfit)
    library(fungible)
    library(mvtnorm)
    library(ClusterR)
    library(cluster)
    library(MLmetrics)
    library(randomForest)
    library(glmnet)
    library(nnet)
    library(nnls)
    library(plyr)
    #library(BiocGenerics)
    library(clusterGeneration)
    
    sd <- sim_data(ndat, ncoef = ncoef,nchoose = nchoose, ntest = ntest,hetero = hetero,inbalance = inbalance,simsep = simsep)$studies_list
    
    sdDecoupled <- trainTestDecoupling(sd,ndat,ntest =1)
    
    
    #for cluster: 
    cf <- clusters_fit2(modfit, modpred, ndat, ncoef, ntest, studies_list = sdDecoupled, cluster_ind = 1, nsep = k)
    errors_cluster <- cf$outmat
    
    
    
    
    #our method 
    errors_combined <-clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 4, nsep = k)$outmat
    
    return(list(errors_cluster = errors_cluster,  
                errors_combined = errors_combined ))
  }
  closeAllConnections()
  
  errors_cluster = results[[1]]
  #errors_random = results[[2]]
  #errors_multi = results[[3]]
  errors_combined =results[[2]]
  
  
  
  
  colnames(errors_cluster) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                                "Stack_noint", "Stack_noint_norm", "Stack_int",
                                "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  errors_cluster = errors_cluster[,-c(2,3,4,8,9,10,12,14)]
  colnames(errors_combined) <- c("Merged",
                                 "Stack_noint", "Stack_noint_norm", "Stack_int",
                                 "Stack_lasso",  "Stack_ridge")
  
  #means of error
  # means_multi <- colMeans(errors_multi)
  means_cluster <- colMeans(errors_cluster)
  #means_random <- colMeans(errors_random)
  means_combined <- colMeans(errors_combined)
  
  #sd means standard deviation
  # sds_multi <- apply(errors_multi, 2, sd)
  sds_cluster <- apply(errors_cluster, 2, sd)
  # sds_random <- apply(errors_random, 2, sd)
  sds_combined <- apply(errors_combined, 2, sd)
  
  return(list( means_cluster = means_cluster,means_combined = means_combined,
               sds_cluster = sds_cluster,sds_combined = sds_combined,
               errors_cluster = errors_cluster,errors_combined = errors_combined))   
}


vary_levels <- function(reps, var_list, modfit, modpred, ndat, ntest,ncoef,nchoose,out_str,k,hetero,inbalance,simsep){
  ptm = proc.time()
  colnames_total <- c("Merged",
                      "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "Stack_lasso", "Stack_ridge")
  total_means_cluster  <-total_means_combined <- array(0, c(length(var_list), 6))
  total_sds_cluster  <- total_sds_combined <- array(0, c(length(var_list), 6))
  
  colnames(total_means_cluster)  <-colnames(total_means_combined) <- colnames_total 
  colnames(total_sds_cluster)<-colnames(total_sds_combined)  <- colnames_total
  
  
  # every i is one specific parameter of interest
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    
    level_rep <- rep.clusters_fit(reps, modfit, modpred, ndat, ntest, ncoef = ncoef, nchoose = nchoose,k = level, hetero = hetero,inbalance=inbalance,simsep = simsep)
    #means
    #total_means_multi[i,] <- level_rep$means_multi
    total_means_cluster[i,] <- level_rep$means_cluster
    #total_means_random[i,] <- level_rep$means_random
    total_means_combined[i,] <- level_rep$means_combined
    #sds
    #total_sds_multi[i,] <- level_rep$sds_multi
    total_sds_cluster[i,] <- level_rep$sds_cluster
    #total_sds_random[i,] <- level_rep$sds_random
    total_sds_combined[i,] <- level_rep$sds_combined
    
    
    # write.table(level_rep$errors_multi, paste0(out_str,"_errors_multi",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_cluster, paste0(out_str,"_errors_cluster",level,".csv"), sep = ",", col.names = colnames_total)
    #write.table(level_rep$errors_random, paste0(out_str,"_errors_random",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_combined, paste0(out_str,"_errors_combined",level,".csv"), sep = ",", col.names = colnames_total)
  }
  
  
  colnames(total_means_cluster) <- colnames(total_means_combined) <- colnames_total 
  colnames(total_sds_cluster)  <- colnames(total_sds_combined)  <- colnames_total
  
  #write.table(total_means_multi, paste0(out_str,"_means_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_cluster, paste0(out_str,"_means_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  #write.table(total_means_random, paste0(out_str,"_means_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_means_combined, paste0(out_str,"_means_combined.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  #write.table(total_sds_multi, paste0(out_str,"_sds_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_cluster, paste0(out_str,"_sds_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  #write.table(total_sds_random, paste0(out_str,"_sds_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_sds_combined, paste0(out_str,"_sds_combined.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  #write.table(indices_mat, paste0(out_str,"_numclusters.csv"), sep = ",", row.names = var_list, col.names = c("mean", "sd"))
  
  return(list( total_means_cluster = total_means_cluster, total_means_combined =total_means_combined,
               total_sds_cluster = total_sds_cluster,total_sds_combined=total_sds_combined))
}

#For gaussian simulations, make sure the correct sim_data is uncommented, and comment the non-gaussian version
#Uncomment the specified lines to create a linear, quadratic, or binary outome.



#test runs 
ncoef = 20
nchoose = 10
inbalance = 1
hetero = 1
ndat = 5
ntest = 1

k_list  <- seq(10, 60, 5)
k_list = c(1:5,k_list)
k_list = k_list[-1]

kcluster = 5
k = kcluster
modfit = randomforestfit
modpred =  randomforestpredict
simsep = 0
cs <- vary_levels(reps = 10, var_list = k_list, modfit = randomforestfit, modpred
                  = randomforestpredict, ndat = ndat, ntest = ntest,ncoef = ncoef, nchoose = nchoose, out_str = "cs", k = kcluster, hetero = hetero, inbalance = inbalance, simsep = simsep)





