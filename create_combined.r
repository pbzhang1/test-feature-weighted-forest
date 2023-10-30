

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
