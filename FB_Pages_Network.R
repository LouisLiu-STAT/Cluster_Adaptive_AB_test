################################################################################
#################   Load Function, Facebook pages data   #######################
################################################################################

load("/CCAS/home/louisliu/FB_Exp/FB_New.RData")
source("/CCAS/home/louisliu/FB_Exp/NETWORK_FUNCTIONS.R")

# load("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/FB_New.RData")
#source("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/NETWORK_FUNCTIONS_Loc.R")





mu        <- c(1.6,1)
zeta      <- c(1.4,1)
sigma     <- 1
B         <- 1000
N         <- 10000

noncores  <- detectCores()
registerDoParallel(cores =noncores)
cl             <- makeCluster(noncores, type="FORK")
 CR.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"CR"  , FB.STAT    )
CAR.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"CAR" , FB.STAT    )
NEW.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"NEW" , FB.STAT    )
stopCluster(cl)

###################################################################################################################################
###################################################################################################################################    
######################                  Evaluate Covariates Imbalance for the network                           ###################
###################################################################################################################################
###################################################################################################################################      



 CR.Weight   <-  apply(  CR.weight.Eval, 1, function(x) mean( I(x==1)) )
CAR.Weight   <-  apply( CAR.weight.Eval, 1, function(x) mean( I(x==1)) )
NEW.Weight   <-  apply( NEW.weight.Eval, 1, function(x) mean( I(x==1)) )

N


registerDoParallel(cores =noncores)
cl         <- makeCluster(noncores, type="FORK")
 CR.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , FB.mbs , "CR"  , FB.STAT ,   CR.Weight  ) 
CAR.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , FB.mbs , "CAR" , FB.STAT ,  CAR.Weight  ) 
NEW.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , FB.mbs , "NEW" , FB.STAT ,  NEW.Weight  ) 
stopCluster(cl)

save.image(file = "FB_NEW_Simulation.RData")