






#calculate the average silhouette score given a data frame that has the cluster membership labeled by variable V1

checksilhouette <- function(studies_list){
  

  merged <- do.call(rbind, studies_list[1:(length(studies_list))])
  
  clus = abs(merged$V1)
  
  k4 = silhouette(clus, dist(merged[,c(-1,-2)]))[,3]
  ss <- mean(k4)
  return(ss)
}






#calculate the average silhouette score given a specific simulation scheme, This shows that even for the same parameters, the silhouette score can vary a lot
CSil <- function(simsep){
  
  sd <- sim_data(ndat, ncoef = ncoef,nchoose = nchoose, ntest = ntest,hetero = hetero,inbalance = inbalance,simsep = simsep)$studies_list
  
  merged <- do.call(rbind, sd[1:(length(sd))])
  
  clus = abs(merged$V1)
  
  k4 = silhouette(clus, dist(merged[,c(-1,-2)]))[,3]
  ss <- mean(k4)
  return(ss)
}


