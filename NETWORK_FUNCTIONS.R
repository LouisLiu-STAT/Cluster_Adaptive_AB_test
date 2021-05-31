###################################################################################################################################
###################################################################################################################################  
#########################                     0. Load packages                #####################################################
###################################################################################################################################
###################################################################################################################################
address.ser    <- "/CCAS/home/louisliu/"

SER.packages   <- c( "dplyr"       , "Rcpp"       , "parallel"  )
Rib_packages   <- c( "iterators"   , "foreach"    , "registry"   , "pkgmaker" ,
                     "rngtools"    , "doRNG"      , "doParallel" )

library(igraph)
library(poweRlaw)
library(Matrix)

LD.ser         <- function(x) require(x,lib.loc="/CCAS/home/louisliu/R_lib", character.only = TRUE)
# LD.ser         <- function(x) require(x, character.only = TRUE)


lapply(Rib_packages, LD.ser  )
lapply(SER.packages, LD.ser )
library(RcppArmadillo)
Rcpp::sourceCpp(paste(  address.ser,"/NWcpp/NW_CL_Rcpp.cpp", sep = "" ))

#Rcpp::sourceCpp("C:/Users/louis/Dropbox/2021 SPRING/Network_Cluster_Covariates_Balance/Large_Network_Proj/Network_Rcpp.cpp")
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
##########################          2. Functions for design the experiments           #############################################
###################################################################################################################################
###################################################################################################################################    
###################################################################################################################################
###################################################################################################################################    
######################          Functions to   calculate the cluster covariates related to the inner nodes    #####################
###################################################################################################################################
################################################################################################################################### 
Cal.NW.STAT <- function( Adj.mat = "Adjacancy matrix " ,  mbs = "clusters memberships"){
  n.cluster  <- length(unique(mbs))
  CL.Edges   <- lapply(1:n.cluster, function(k) sapply(1:n.cluster ,  function(x) rowSums( Adj.mat[which(mbs==k), which(mbs==x)] )))
  Inner.ID   <- lapply(1:n.cluster, function(k) 1*(rowSums(CL.Edges[[k]][,-k])==0))
  #### the number of the inner nodes or outer nodes can be <= 1, to avoid the error of  apply. 
  Cal.CL.STAT <-function(k){
    In.num        <- sum(Inner.ID[[k]])
    In.deg.ave    <- ifelse( In.num == 0 , 0 ,  sum(    Inner.ID[[k]]*rowSums(CL.Edges[[k]]))/In.num   )
    De.num        <- sum(Inner.ID[[k]]==0)
    De.deg.ave    <- ifelse( De.num == 0 , 0 ,  sum((1-Inner.ID[[k]])*rowSums(CL.Edges[[k]]))/De.num   )
    re<- c(   In.num    ,      In.deg.ave  ,  De.num ,     De.deg.ave )
    return(re)
  }
  ### 1. inner nodes; 2. inner degrees; 3. outer nodes;  4. outer degrees
  CL.covariates  <- t( sapply( 1 : n.cluster , Cal.CL.STAT ) )
  result         <- list( CL.Edges,  Inner.ID , CL.covariates )
  names(result)  <- c(    "Edges" , "Inner.ID", "CL.Cov"      )
  return( result )
}
###################################################################################################################################
###################################################################################################################################    
Used.STAT <- function( Adj.mat = "Adjacancy matrix " ,  mbs = "clusters memberships"){
  n.nodes          <- nrow(Adj.mat)
  CL.STAT           <- Cal.NW.STAT(Adj.mat, mbs)
  n.cluster        <- nrow(CL.STAT[[3]])
  result         <- list( n.nodes, n.cluster , CL.STAT[[1]] ,  CL.STAT[[2]]  ,  CL.STAT[[3]])
  names( result )       <- c("n.nodes", "n.cluster","CL.Edges","Inner.ID", "CL.Cov")
  return( result )   
}
###################################################################################################################################
###################################################################################################################################    
######################    Pairwise sequential design for balancing covariates of the network characteristics    ###################
###################################################################################################################################
###################################################################################################################################    
PSR.NT  <- function( Adj.mat = "Adjacancy matrix " ,  mbs = "clusters memberships" ){
  CL.STAT    <- Cal.NW.STAT( Adj.mat , mbs )
  n.cluster <- nrow(CL.STAT[[3]])
  RE.Covar   <- matrix(0, nrow = n.cluster, ncol = 4)  
  CL.ID      <- sample( 1:n.cluster , n.cluster  , replace = F )
  RE.Covar   <- CL.STAT[[3]][CL.ID,]
  RE.Treat   <- PSR(  RE.Covar ,0.85 )
  CL.Treat   <- Match(RE.Treat,CL.ID)
  return(CL.Treat)
}
###################################################################################################################################
###################################################################################################################################    
######################    Proposed design for balancing both of the two kinds of the network characteristics    ###################
###################################################################################################################################
###################################################################################################################################    
###################################################################################################################################
###################################################################################################################################    
New.Design    <- function( Adj.mat = "Adjacancy matrix "     ,  mbs = "clusters memberships"   , 
                           B = "Times of Simulation"   ,  Ass = "Assumption: Full or Frac" ){
  CL.STAT     <-  Cal.NW.STAT( Adj.mat , mbs )
  n.cluster  <- nrow(CL.STAT[[3]])
  FS.Rand     <- function(  x  , CL.STAT ){
    CL.ID     <- RE.Covar <- RE.Treat <- rep(0,n.cluster)
    RE.Covar  <- matrix(0, nrow = n.cluster, ncol = 4)  
    CL.ID     <- sample( 1:n.cluster , n.cluster  , replace = F )
    RE.Covar  <- CL.STAT[[3]][CL.ID,]
    RE.Treat  <- PSR(   RE.Covar ,  0.85 )
    CL.Treat  <- Match( RE.Treat , CL.ID )
    Out.use.cov    <- CAL_TREAT_Cov( CL.Treat , CL.STAT[[2]] , CL.STAT[[1]] , Ass , frac )
    all.cov    <- cbind( Out.use.cov , CL.STAT[[3]] )
    diff.mean   <- CL.Treat %*%  all.cov   / sum( CL.Treat ) - (1-CL.Treat) %*% all.cov  / sum(1-CL.Treat)
    result    <- list( CL.Treat , diff.mean      )
    names(result)   <- c(    "Treat"  , "diff-in-mean" )
    return( result ) 
    gc()
  }
  #### Calculate the first kind of covariates 
  Treat.ALL  <- lapply( 1:B ,  function(x)    FS.Rand( x ,  CL.STAT ))  #### Evaluate the treatment assignments
  ######################################################################################################
  ########  Imbalance evaluation: calculate the Mah dist for each of the treatment assignments  ########
  ######################################################################################################
  DIM  <- sapply( 1:B , function(b) Treat.ALL[[b]][[2]]) 
  DIM  <- t(DIM)
  Mah.dist  <- Cal_diff( DIM ) 
  Treat      <- Treat.ALL[[sample(order(Mah.dist)[1:(B*0.05)],1)]][[1]]  #### assign the treatment
  return(Treat)
}
###################################################################################################################################
###################################################################################################################################    
######################                   Function for the three types of the designs                            ###################
###################################################################################################################################
###################################################################################################################################    
Design.func <- function( Adj.mat = "Adjacancy matrix "    ,
                         mbs = "clusters memberships" ,
                         Types = "Design types"               ){
  if(Types == "CR"){
    n.cluster        <- length(unique(mbs))
    CL.Treat         <- rep(0,n.cluster)
    CL.ID            <- sample(1:n.cluster, floor( n.cluster/2)  , replace =  F  )
    trt              <- rbinom(1,1,0.5)
    CL.Treat[CL.ID]  <- trt
    CL.Treat[-CL.ID] <- 1-trt
  }else if( Types == "CAR"  ){
    CL.Treat         <- PSR.NT( Adj.mat , mbs )
  }else if(Types == "NEW"){
    CL.Treat         <- New.Design(Adj.mat,mbs, 1000, Ass) 
  }  
  return(CL.Treat)
  gc()
}

CAL.imb.stat <- function(x){
  Mean        <-mean(x)
  Bias       <-  mean(abs(x))
  Variance   <-  var(x)
  return(c(Mean,Bias, Variance))
}
###################################################################################################################################
###################################################################################################################################    
######################                   Function for the three types of the designs                            ###################
###################################################################################################################################
###################################################################################################################################
########################################################################################################
########################################################################################################
################ Function to calculate weights for each of the randomization designs  ##################
########################################################################################################
########################################################################################################
Eval.Inclusion <- function( itr,  Types = "Design types"             , 
                            CL.stat = "Cluster Statistics"  ){
  FS.Rand       <- function(  x  ,n.cluster , CL.Edges , Inner.ID , CL.cov){
    CL.ID       <- RE.Covar <- RE.Treat <- rep(0,n.cluster)
    RE.Covar    <- matrix(0, nrow = n.cluster, ncol = 4)  
    CL.ID       <- sample( 1:n.cluster , n.cluster  , replace = F )
    RE.Covar    <- CL.cov[CL.ID,]
    RE.Treat    <- PSR(   RE.Covar ,  0.85 )
    CL.Treat    <- Match( RE.Treat , CL.ID )
    Out.use.cov  <- CAL_TREAT_Cov( CL.Treat , Inner.ID , CL.Edges )
    all.cov      <- cbind( Out.use.cov , CL.cov )
    diff.mean   <- CL.Treat %*%  all.cov   / sum( CL.Treat ) - (1-CL.Treat) %*% all.cov  / sum(1-CL.Treat)
    result      <- list( CL.Treat , diff.mean      )
    names(result)   <- c(    "Treat"  , "diff-in-mean" )
    return( result ) 
  }
  
  n.nodes          <- CL.stat[[1]]
  n.cluster        <- CL.stat[[2]]
  CL.Edges         <- CL.stat[[3]]
  Inner.ID         <- CL.stat[[4]]
  First.cov        <- CL.stat[[5]]
  
  if(Types == "CR"){
    CL.Treat         <- rep(0,n.cluster)
    CL.ID            <- sample(1:n.cluster, floor( n.cluster/2)  , replace =  F  )
    trt              <- rbinom(1,1,0.5)
    CL.Treat[ CL.ID] <- trt
    CL.Treat[-CL.ID] <- 1-trt
  }else if( Types == "CAR"  ){
    CL.Treat       <- rep(0,n.cluster)
    CL.ID          <- sample( 1:n.cluster , n.cluster  , replace = F )
    RE.Covar       <- First.cov[CL.ID,]
    RE.Treat       <- PSR(   RE.Covar ,  0.85 )
    CL.Treat       <- Match( RE.Treat , CL.ID )
  }else if(Types == "NEW"){
    Treat.ALL      <- lapply( 1:B ,  function(x)    FS.Rand( x , n.cluster,  CL.Edges , Inner.ID ,First.cov  ))  #### Evaluate the treatment assignments
    DIM      <- sapply( 1:B , function(b) Treat.ALL[[b]][[2]]) 
    DIM      <- t(DIM)
    Mah.dist     <- Cal_diff( DIM ) 
    CL.Treat    <- Treat.ALL[[sample(order(Mah.dist)[1:(B*0.05)],1)]][[1]]  #### assign the treatment
  }   
  ngb.frac   <-   CAL_NHTreat(n.nodes, CL.Treat ,  Inner.ID , CL.Edges)
  return( ngb.frac )
  gc()
}
###################################################################################################################################
###################################################################################################################################    
######################                     Evaluation Treatment effect                                 ############################
###################################################################################################################################
###################################################################################################################################      



TreatImb.Eval <- function( itr ,  mbs     = "clusters memberships"    ,
                           Types     = "Design types"            ,  
                           CL.stat     = "Cluster Statistics"      ,  
                           HT.weights  = "HT weights"              ){
  FS.Rand       <- function(  x  , n.cluster , CL.Edges , Inner.ID , CL.cov ){
    CL.ID       <- RE.Covar <- RE.Treat <- rep(0,n.cluster)
    RE.Covar    <- matrix(0, nrow = n.cluster, ncol = 4)  
    CL.ID       <- sample( 1:n.cluster , n.cluster  , replace = F )
    RE.Covar    <- CL.cov[CL.ID,]
    RE.Treat    <- PSR(   RE.Covar ,  0.85 )
    CL.Treat    <- Match( RE.Treat , CL.ID )
    Out.use.cov  <- CAL_TREAT_Cov( CL.Treat , Inner.ID , CL.Edges )
    all.cov      <- cbind( Out.use.cov , CL.cov )
    diff.mean   <- CL.Treat %*%  all.cov   / sum( CL.Treat ) - (1-CL.Treat) %*% all.cov  / sum(1-CL.Treat)
    result      <- list( CL.Treat , diff.mean      )
    names(result)   <- c(    "Treat"  , "diff-in-mean" )
    return( result ) 
  }
  
  n.nodes          <- CL.stat[[1]]
  n.cluster        <- CL.stat[[2]]
  CL.Edges         <- CL.stat[[3]]
  Inner.ID         <- CL.stat[[4]]
  First.cov        <- CL.stat[[5]]
  
  n.deg            <- lapply( CL.Edges , rowSums )
  n.deg            <- unlist( n.deg )
  Ave.edges        <- sum(  n.deg       )/n.nodes
  epsilon          <- rnorm( n.nodes, 0 , sigma  )
  Epsilon          <- lapply( 1 : n.cluster , function(k) epsilon[mbs==k])
  
  sd.deg           <-  n.deg/Ave.edges
  PO.Y0            <-  mu[2] +  zeta[2] *  sd.deg  + epsilon
  PO.Y1            <-  mu[1] +  zeta[1] *  sd.deg  + epsilon
  
  CL.PO.barY1  <-  sapply( 1:n.cluster  , function(k)  mean(   PO.Y1[which(mbs==k)] ) )
  CL.PO.barY0  <-  sapply( 1:n.cluster  , function(k)  mean(   PO.Y0[which(mbs==k)] ) )
  CL.PO.sumY1  <-  sapply( 1:n.cluster  , function(k)   sum(   PO.Y1[which(mbs==k)] ) )
  CL.PO.sumY0  <-  sapply( 1:n.cluster  , function(k)   sum(   PO.Y0[which(mbs==k)] ) )
  
  IN.TAU        <-  mean(PO.Y1-PO.Y0)
  CL.TAU        <-  mean(CL.PO.barY1- CL.PO.barY0)
  
  HT.wgts     <-  lapply(1:n.cluster, function(k)   HT.weights[mbs==k] )          
  
  if(Types == "CR"){
    CL.Treat         <- rep(0,n.cluster)
    CL.ID            <- sample(1:n.cluster, floor( n.cluster/2)  , replace =  F  )
    trt              <- rbinom(1,1,0.5)
    CL.Treat[ CL.ID] <- trt
    CL.Treat[-CL.ID] <- 1-trt
    Out.cov          <- CAL_TREAT_Cov( CL.Treat , Inner.ID , CL.Edges )
    all.cov          <- cbind( Out.cov , First.cov )
    diff.mean        <- CL.Treat %*%  all.cov / sum( CL.Treat ) - (1-CL.Treat) %*% all.cov  / sum(1-CL.Treat)
  }else if( Types == "CAR"  ){
    CL.Treat        <- RE.Treat <- rep(0,n.cluster)
    CL.ID           <- sample( 1:n.cluster , n.cluster  , replace = F )
    RE.Covar        <- First.cov[CL.ID,]
    RE.Treat        <- PSR(   RE.Covar ,  0.85 )
    CL.Treat        <- Match( RE.Treat , CL.ID )
    Out.cov          <- CAL_TREAT_Cov( CL.Treat , Inner.ID , CL.Edges )
    all.cov          <- cbind( Out.cov , First.cov )
    diff.mean         <- CL.Treat %*%  all.cov / sum( CL.Treat ) - (1-CL.Treat) %*% all.cov  / sum(1-CL.Treat)
  }else if(Types == "NEW"){
    CL.Treat             <- RE.Treat <- rep(0,n.cluster)
    Treat.ALL          <- lapply( 1:B ,  function(x)    FS.Rand( x , n.cluster,CL.Edges , Inner.ID , First.cov  ) )  #### Evaluate the treatment assignments
    DIM              <- sapply( 1:B , function(b) Treat.ALL[[b]][[2]]) 
    DIM              <- t(DIM)
    Mah              <- Cal_diff( DIM )
    Full.ID           <- sample(order(Mah)[1:floor((B*0.05))], 1 )
    CL.Treat        <- Treat.ALL[[Full.ID]][[1]]   #### assign the treatment
    diff.mean         <- Treat.ALL[[Full.ID]][[2]]   
  }   
  
  In.Indic          <- Inc_Full( CL.Treat , Inner.ID ,  CL.Edges        )
  CL.TREAT.PY       <- Gen_response(  mu , zeta , CL.Treat , CL.Edges , Epsilon ,  Ave.edges )
  Adj.wgts           <- Adj_Wgts( In.Indic$Out_Full , Inner.ID )
  
  
  Adj.Y.sum  <-   sapply(1:n.cluster, function(k) sum( Inner.ID[[k]]*CL.TREAT.PY[[k]]) + ifelse( Adj.wgts[k] ==0,0,sum(In.Indic$Out_Full[[k]]*CL.TREAT.PY[[k]])/Adj.wgts[k] )       )  
  Adj.Y.bar  <-   Adj.Y.sum / sapply(1:n.cluster, function(k) length(Inner.ID[[k]]) )
  
  #### Calculate the estimates using the correct outcomes and
  HT.Y.sum  <-   sapply(1:n.cluster, function(k) sum( ifelse(  HT.wgts[[k]]  ==0,0, In.Indic$In_Full[[k]]*CL.TREAT.PY[[k]]/ HT.wgts[[k]]   ) ) ) 
  HT.Y.bar  <-    HT.Y.sum /  sapply(1:n.cluster, function(k) length(Inner.ID[[k]]) )
  
  HT.In <-  sum( CL.Treat* HT.Y.sum  -(1-CL.Treat )*HT.Y.sum )/ (0.5*n.nodes)
  HT.CL <-  sum( CL.Treat* HT.Y.bar - (1-CL.Treat )*HT.Y.bar) / (0.5*n.cluster) 
  
  #### Calculate the new proposed estimator
  
  Adj.In <-   sum( CL.Treat* Adj.Y.sum - (1-CL.Treat)*Adj.Y.sum  )/(0.5*n.nodes  )
  Adj.CL <-   sum( CL.Treat* Adj.Y.bar - (1-CL.Treat)*Adj.Y.bar  )/(0.5*n.cluster)
  
  
  Est.re  <-  c( IN.TAU, HT.In , Adj.In , CL.TAU, HT.CL , Adj.CL )
  names(Est.re) <- c( "Tau.IN","HT.IN",  "Adj.IN","Tau.CL","HT.CL",  "Adj.CL")
  names(diff.mean) <- c(paste("X",1:6,sep = ""))  
  re <- c( Est.re , diff.mean )
  return( re )
}  