# Cluster_Adaptive_AB_test

This repository contains  the R codes used for the numerical studies in  "Adaptive  A/B Test on  Networks with Cluster Structures".

## Preparations

 Please install the following packages: "dplyr"       , "Rcpp"       , "parallel" , "iterators"   , "foreach"    , "registry"   , "pkgmaker" ,
                     "rngtools"    , "doRNG"      , "doParallel" ,"igraph", "randnet","dplyr".
            
## Instructions
1. The Rcpp file "NW_CL_Rcpp.cpp" contains the Rcpp fuctions used for the design and esimation.
2. The R file "NETWORK.FUNCTIONS.R" contains the functions for different designs, the functions used to evaluate the weights for the Horivitz Thompson estimators, and the evaluation of the ATEs. 
3. The R file "2TypeCLsResults.R" contains the codes for the four clusters example.
4. The R file "Large_Numerical_Clusters.R" contains the code for the numerical study for the hypothestical network.
5. The R file "MIT_NETWORK.R" contains the code for the numerical study for the MIT Phone call network.   
6. The R file "FB_Pages_Network.R" contains the code for the numerical study for the Facebook page to page network.                     
7. The Rdata files "FB_NEW.Rdata" and "MIT_CL.Rdata" consist of the information of the adjacancy matrices, and the labels of the clusters of the Facebook page to page network and the MIT phone call network used in the numerical studies.  
