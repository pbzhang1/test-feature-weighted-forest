
create_clusters2 <- function(studies_list, ntest, ncoef, nsep){
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  k2 <- kmeans(merged[,-1], centers = k, nstart = 25)$cluster
  clusters_list <- lapply(split(seq_along(k2), k2), 
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
clusters_fit2 <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep){
  if (cluster_ind == 1){ #K-means clustering
    cc <- create_clusters2(studies_list, ntest, ncoef, nsep)
    edat <- cc$clusters_list
  }
  else if (cluster_ind == 2){ #random 'pseudo-studies'
    edat <- create_random(studies_list, ntest, ncoef, nsep)
  }
  else if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    edat <- studies_list
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  #mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  
  
  for (j in 1:ntrain){
    mods[[j]] <- modfit(edat[[j]])
    
    preds <- lapply(edat[1:ntrain], function(x){
      modpred(mods[[j]], newdata = x[, -1]) 
    })
    mses[j,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
      newdata = x[, -1]
      preds <- modpred(mods[[j]], newdata = newdata) 
      mean((preds - x[,"y"])^2)}
    ))
    
    curpreds <- lapply(preds, as.numeric)
    allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
  }
  diag(mses) <- NA
  
  # CS Weights
  tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  weights <- absnorm(sqrt(tt), max.norm = TRUE)
  nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  nwts <- absnorm(nk[1:ntrain])
  
  # Regression: stacked (intercept and no intercept)
  predstack <- do.call(rbind, allpreds)
  coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 0, lower.limits = 0, intercept = T)))
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0
  
  coefs_stack_noint_norm <- absnorm(coefs_stack_noint)
  
  # Regression: study-specific (intercept and no intercept)
  coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)
  coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)
  coefs_ss_lasso <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 1, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_ridge <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 0, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_lasso <- colMeans(do.call(rbind, coefs_ss_lasso), na.rm = T)
  coefs_ss_ridge <- colMeans(do.call(rbind, coefs_ss_ridge), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_lasso[which(is.na(coefs_ss_lasso))] <- 0
  
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  outmat <- matrix(NA, ntest, 14)
  colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1])
    
    merged <- as.vector(sapply(merged, as.numeric))
    allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
    allmod <- apply(allmod, 2, as.numeric)
    # Unweighted average
    unweighted <- colMeans(allmod)
    
    # sample size weighted
    sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})
    
    # cross-study weighted 
    cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})
    
    # regression: stacked (noint, int, each normed) + lasso
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    # regression: study_specific (noint, int, noint normed) + lasso
    ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
    ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
    ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})
    ss_lasso <- apply(allmod, 2, function(x){coefs_ss_lasso[1] + sum(coefs_ss_lasso[-1]*x)})
    ss_ridge <- apply(allmod, 2, function(x){coefs_ss_ridge[1] + sum(coefs_ss_ridge[-1]*x)})
    
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
                                  mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                  mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), 
                                  mean((cury - ss_int)^2), mean((cury - stack_lasso)^2), mean((cury - ss_lasso)^2), 
                                  mean((cury - stack_ridge)^2), mean((cury - ss_ridge)^2)))
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  return(list(outmat = colMeans(outmat)))
}
