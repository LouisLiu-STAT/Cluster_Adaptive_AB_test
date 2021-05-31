################################################################################
#################   Load Function, MIT data   #######################
################################################################################

source("/CCAS/home/louisliu/FB_Exp/NETWORK_FUNCTIONS.R")

# load("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/MIT_CL.Rdata")
# source("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/NETWORK_FUNCTIONS_Loc.R")

load("/CCAS/home/louisliu/MIT_EXP/MIT_CL.Rdata")

mu        <- c(1.6,1)
zeta      <- c(1.4,1)
sigma     <- 1
B         <- 1000
N         <- 10000



noncores  <- detectCores()
registerDoParallel(cores =noncores)
cl             <- makeCluster(noncores, type="FORK")
 CR.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"CR"  , MIT.STAT    )
CAR.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"CAR" , MIT.STAT    )
NEW.weight.Eval      <- foreach(itr = 1:1000,.combine = cbind) %dopar% Eval.Inclusion(itr,"NEW" , MIT.STAT    )
stopCluster(cl)


#  CR.weight.Eval<-sapply(1:1000, function(itr) Eval.Inclusion(itr,"CR"  , MIT.STAT    ))
# CAR.weight.Eval<-sapply(1:1000, function(itr) Eval.Inclusion(itr,"CAR" , MIT.STAT    ))
# NEW.weight.Eval<-sapply(1:1000, function(itr) Eval.Inclusion(itr,"NEW" , MIT.STAT    ))


 CR.Weight   <-  apply(  CR.weight.Eval, 1, function(x) mean( I(x==1)) )
CAR.Weight   <-  apply( CAR.weight.Eval, 1, function(x) mean( I(x==1)) )
NEW.Weight   <-  apply( NEW.weight.Eval, 1, function(x) mean( I(x==1)) )



registerDoParallel(cores =noncores)
cl         <- makeCluster(noncores, type="FORK")
    CR.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , MIT.mbs , "CR"  , MIT.STAT ,   CR.Weight  ) 
   CAR.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , MIT.mbs , "CAR" , MIT.STAT ,  CAR.Weight  ) 
   NEW.Treat  <- foreach(itr = 1:N,.combine = cbind) %dopar% TreatImb.Eval( itr , MIT.mbs , "NEW" , MIT.STAT ,  NEW.Weight  ) 
stopCluster(cl)

save.image(file = "MIT_Results.RData")

