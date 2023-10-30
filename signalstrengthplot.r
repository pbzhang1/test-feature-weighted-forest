



colnames_total <- c("varlist", "Merged", 
                    "Stack_noint", "Stack_noint_norm", "Stack_int",
                    "Stack_lasso",  "Stack_ridge")


compare3 <-function(filename,n=30){
  
  #means_cluster <- read.csv( "cs_means_cluster.csv",skip = 1, header = F)
  #sds_cluster <- read.csv("cs_sds_cluster.csv", skip = 1, header = F)
  
  #means_combined <- read.csv( "cs_means_combined.csv",skip = 1, header = F)
  #sds_combined <- read.csv( "cs_sds_combined.csv", skip = 1, header = F)
  
  
  
  
  
  means_cluster <-  read.csv( paste0("justiceForMaya/", filename, "/cs_means_cluster.csv"),skip = 1, header = F)
  means_combined <-  read.csv(paste0("justiceForMaya/", filename, "/cs_means_combined.csv"),skip = 1, header = F)
  
  
  x = means_cluster$V1
  sds_combined <-  read.csv(paste0("justiceForMaya/", filename, "/cs_sds_combined.csv"),skip = 1, header = F)
  
  sds_cluster <- read.csv(paste0("justiceForMaya/", filename, "/cs_sds_cluster.csv"),skip = 1, header = F)
  
  
  cluster_avg = means_cluster$V7
  cluster_std = sds_cluster$V7*1.96/sqrt(n)
  
  combined_avg = means_combined$V7
  combined_std = sds_combined$V7*1.96/sqrt(n)
  plot(x, cluster_avg,
       ylim=range(c(combined_avg-combined_std, cluster_avg+cluster_std)),
       pch=19, xlab="norm of coefficient a & b ", ylab="% change RMSE",
       col = 'red'
  )
  arrows(x, cluster_avg-cluster_std, x, cluster_avg+cluster_std, length=0.05, angle=90, code=3,col = 'black')
  
  
  
  
  
  
  combined_avg = means_combined$V7
  combined_std = sds_combined$V7*1.96/sqrt(n)
  
  points(x, combined_avg,
         ylim=range(c(combined_avg-combined_std, cluster_avg+cluster_std)),
         pch=19, xlab="norm of coefficient a & b ", ylab="% change RMSE",
         col = 'blue'
  )
  arrows(x, combined_avg-combined_std, x, combined_avg+combined_std, length=0.05, angle=90, code=3,col = 'black')
  
  
  legend(x = "topright", legend=c("unweighted", "weighted"), 
         fill = c("red","blue")
  )
  
  
  
}
means_cluster$Stack_ridge - 1.96*sds_cluster$Stack_ridge



compare1(1)

