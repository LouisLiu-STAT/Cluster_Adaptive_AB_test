################################################################################
#################   Load Function, Facebook pages data   #######################
################################################################################

source("/CCAS/home/louisliu/FB_Exp/NETWORK_FUNCTIONS.R")

# source("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/NETWORK_FUNCTIONS_Loc.R")

###################################################################################################################################
###################################################################################################################################  
#########################         1. Function to  Generate Clusters           #####################################################
###################################################################################################################################
###################################################################################################################################
CL.Generation <- function(K, k.min, alpha, d.min, d.max, od ){
  K.nodes <- rpldis(K, k.min, alpha)                               ### generate cluster sizes from power law distribution
  n.nodes <- sum(K.nodes)                                          ### number of nodes
  Clu.lab <- sapply(1:K, function(x) rep(x,K.nodes[x] )) %>%unlist ### clusters
  Adj.mat <- matrix(0, nrow = n.nodes, ncol = n.nodes)             ### Adj Matrix
  ### maximum probability of connectiveness 
  d.k     <- runif(K,d.min,d.max)                                  ### generate the probability of the connectiveness for each cluster
  Be.max  <- 1/max(K.nodes)^(od)/K                                 ### maximum probability between clusters 
  CL.ID   <- 1:K                                                   ### cluster seq
  for( k in c(1:K)){
    In.clp  <- d.k[k]+ runif(K.nodes[k], -0.075, 0.075 )                              ### Make the degrees different     
    In.adj  <- matrix(0, nrow =K.nodes[k] , ncol =K.nodes[k]  )                   ### Generate the within cluster adjacency matrix
    In.adj[1, 2:K.nodes[k] ] <-   rbinom(K.nodes[k]-1,1, In.clp[1] )
    for ( i in 2:(K.nodes[k]-1) ) {
      p <- (In.clp[i]*K.nodes[k]- sum(In.adj[c(1:(i-1)),i] ))/ (K.nodes[k]-(i-1) )
      p <- ifelse( p>1,1,p)
      p <- ifelse( p<0,0,p)
      In.adj[i, c( (i+1):K.nodes[k] )   ] <-  rbinom(K.nodes[k]-i,1,p   )
    }
    clk.id <-   which(Clu.lab==k)
    Adj.mat[clk.id ,  clk.id] <- In.adj
    for( x in CL.ID[-k]  ){
      kx.p    <- runif(1,0,Be.max)
      clkx.id <-which( Clu.lab== x)
      Adj.mat[ clk.id,clkx.id] <- sapply(1:length( clk.id),function(x) rbinom(clkx.id,1,    kx.p )   ) %>%t()
    }
  }
  Adj.mat       <- Adj.mat + t(Adj.mat)
  result        <- list(   K.nodes , Clu.lab , Adj.mat  )
  names(result) <- c( "Cluster_sizes" , "memberships", "Adj.matrix"   )
  return(result)
}

 
###################################################################################################################################
###################################################################################################################################
##########################            parameter settings           ################################################################
###################################################################################################################################
###################################################################################################################################
K         <- 20                                                    ### number of clusters
k.min     <- 10                                                    ### minimal number of nodes in each cluster
alpha     <- 2.8                                                   ### parameter for the power-law distribution, i.e., x^-alpha
d.min     <- 0.3                                                   ### minimal probability of connectiveness 
d.max     <- 0.5
od        <- 1
B         <- 1000
M         <- 100

                                                 ### number of clusters


Clusters2  <- lapply(1:M, function(x) CL.Generation(K,k.min, alpha,d.min,d.max,od) )
CL.stat2   <- lapply(1:M, function(x)  Used.STAT(Clusters2[[x]]$Adj.matrix, Clusters2[[x]]$memberships))


K         <- 50                                                    ### number of clusters


Clusters5  <- lapply(1:M, function(x) CL.Generation(K,k.min, alpha,d.min,d.max,od) )
CL.stat5   <- lapply(1:M, function(x)  Used.STAT(Clusters5[[x]]$Adj.matrix, Clusters5[[x]]$memberships))



#########################################################################################################################################################################################
#########################################################################################################################################################################################

B         <- 1000
sigma     <- 1
mu        <- c(1.6,1)
zeta      <- c(1.4,1)


#########################################################################################################################################################################################
#########################################################################################################################################################################################


noncores  <- detectCores()
registerDoParallel(cores =noncores)
cl              <- makeCluster(noncores, type="FORK")
 CR.weight.Eval2      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,   "CR",CL.stat2[[k]]))
CAR.weight.Eval2      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,  "CAR",CL.stat2[[k]]))
NEW.weight.Eval2      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,  "NEW",CL.stat2[[k]]))
stopCluster(cl)




 CR.weight2 <- lapply(1:M, function(m)   apply(    CR.weight.Eval2[[m]] , 1, function(x) mean( I(x==1))  ) )
CAR.weight2 <- lapply(1:M, function(m)   apply(   CAR.weight.Eval2[[m]] , 1, function(x) mean( I(x==1))  ) )                             
NEW.weight2 <- lapply(1:M, function(m)   apply(   NEW.weight.Eval2[[m]] , 1, function(x) mean( I(x==1))  ) )                          




registerDoParallel(cores =noncores)
cl             <- makeCluster(noncores, type="FORK")
 CR.RE2      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters2[[k]]$memberships,"CR"  ,CL.stat2[[k]],  CR.weight2[[k]]) )
CAR.RE2      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters2[[k]]$memberships,"CAR" ,CL.stat2[[k]], CAR.weight2[[k]]))
NEW.RE2      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters2[[k]]$memberships,"NEW" ,CL.stat2[[k]], NEW.weight2[[k]]))
stopCluster(cl)
20

#########################################################################################################################################################################################
#########################################################################################################################################################################################


noncores  <- detectCores()
registerDoParallel(cores =noncores)
cl              <- makeCluster(noncores, type="FORK")
 CR.weight.Eval5      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,   "CR",CL.stat5[[k]]))
CAR.weight.Eval5      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,  "CAR",CL.stat5[[k]]))
NEW.weight.Eval5      <- sapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr ,  "NEW",CL.stat5[[k]]))
stopCluster(cl)


 CR.weight5 <- lapply(1:M, function(m)   apply(    CR.weight.Eval5[[m]] , 1, function(x) mean( I(x==1))  ) )
CAR.weight5 <- lapply(1:M, function(m)   apply(   CAR.weight.Eval5[[m]] , 1, function(x) mean( I(x==1))  ) )                             
NEW.weight5 <- lapply(1:M, function(m)   apply(   NEW.weight.Eval5[[m]] , 1, function(x) mean( I(x==1))  ) )                          


registerDoParallel(cores =noncores)
cl             <- makeCluster(noncores, type="FORK")
 CR.RE5      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters5[[k]]$memberships,"CR"  ,CL.stat5[[k]],  CR.weight5[[k]]) )
CAR.RE5      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters5[[k]]$memberships,"CAR" ,CL.stat5[[k]], CAR.weight5[[k]]))
NEW.RE5      <- lapply(1:M, function(k)  foreach(itr = 1:1000,.combine = cbind) %dopar% TreatImb.Eval(itr,Clusters5[[k]]$memberships,"NEW" ,CL.stat5[[k]], NEW.weight5[[k]]))
stopCluster(cl)
save.image(file = "SIMU_Multi_CLuster.Rdata")
