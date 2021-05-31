#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List  CAL_Inner( List   CAL_Edge   ){
  int               n_CL   =  CAL_Edge.size()             ;
  List                        Inner_ID( n_CL  )           ;
  NumericMatrix                  CL_cov(n_CL,4)           ;
  for( int k = 0 ;  k < n_CL ; k++ ){
    NumericMatrix  CLE  =  CAL_Edge[k]                    ;
    int             nk  = CLE.nrow()                      ;
    NumericVector    Inner_K( nk )                        ;
    NumericVector     degree( nk )                                  ;
    NumericVector   Inner_deg( nk )                                  ;
    for(int i = 0; i<nk;i++){
      NumericVector Odegree(n_CL)                                    ;
      degree[i] = sum(CLE(i,_))                                ;
      for(int j = 0; j<n_CL; j++){
        if(j==k){
          Odegree[j] =  0                                            ;
        }else{
          Odegree[j] = CLE(i,j)                                      ;
        }
      }
      int Odeg   = sum(Odegree)                                        ;
      if(Odeg==0){
        Inner_K[i] = 1                                        ;
      }else{
        Inner_K[i] = 0                                        ;
      }
      Inner_deg[i]= Inner_K[i]*degree[i]                     ;
    }   
    Inner_ID[k] = Inner_K                                  ; 
    CL_cov(k,0) = sum(Inner_K)                             ;
    CL_cov(k,1) = sum(Inner_deg)/ CL_cov(k,0)              ;
    CL_cov(k,2) = nk-sum(Inner_K)                          ;
    
    if( CL_cov(k,2)==0 ){
      CL_cov(k,3) = 0 ;
    }else{
      CL_cov(k,3) = (sum(degree)-sum(Inner_deg))/CL_cov(k,2) ;
    }
    
  }
  
  List  L = List::create( Named("Inner_ID") = Inner_ID  ,  Named("CL_Cov") = CL_cov );
  return L;
}

// [[Rcpp::export]]
NumericMatrix   CAL_TREAT_Cov( NumericVector  CL_Treat ,
                               List        Inner_ID    ,
                               List           CL_Edge   ){
  int                  n_CL   =  CL_Treat.size()              ;
  NumericMatrix               De_Cov( n_CL , 2 )              ;
  for( int k = 0 ; k < n_CL ;  k++ ){
    
    NumericVector           in_id   = Inner_ID[k]               ;
    int                        nk   = in_id.size()              ;
    NumericVector                      N_trt(nk)                ;
    NumericVector                  Outer_num(nk)                ;
    NumericVector                  Outer_deg(nk)                ;
    NumericMatrix            k_Edge =  CL_Edge[k]               ; 
    
    for(int i = 0 ; i<nk; i++){
      int degree = sum(k_Edge(i,_) )                        ;
      if( in_id[i] == 1 ) {
        N_trt[i] = sum(k_Edge(i,_) )                    ;
        Outer_num[i] =  0                                   ;
        Outer_deg[i] =  0                                   ;
      }else{
        int     num_trt = 0                                 ;
        int         trt = CL_Treat[k]                       ;
        if(trt==1){
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + CL_Treat[j]*k_Edge(i,j)    ; 
          } 
        }else{
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + (1-CL_Treat[j])*k_Edge(i,j)   ; 
          } 
        }
        N_trt[i] = num_trt  ;
        
        if( num_trt == degree ){ 
          Outer_num[i] = 1                                      ;
          Outer_deg[i] = degree                                 ;      
        }else{
          Outer_num[i] = 0                                      ;
          Outer_deg[i] = 0                                      ;  
        }
        
      }
    }
    int       n_use  = sum( Outer_num )                            ;
    De_Cov( k , 0 )  = n_use                                       ;
    if( n_use == 0 ){
      De_Cov( k , 1 )  = 0                                         ;
    }else{
      De_Cov( k , 1 )  = sum( Outer_deg ) / n_use                  ;        
    }
  }
  return De_Cov;
}

// [[Rcpp::export]]
NumericVector   CAL_NHTreat(  int              n_node      ,
                              NumericVector   CL_Treat     ,
                              List            Inner_ID     ,
                              List            CL_Edge   ){
  int                  n_CL   =  CL_Treat.size()                 ;
  NumericVector                 Num_Treat(n_node)                         ;
  int                    id   = 0                                ;
  for( int k = 0 ; k < n_CL ;  k++ ){
    
    NumericVector           in_id   = Inner_ID[k]                ;
    int                        nk   =    in_id.size()            ;
    NumericMatrix            k_Edge =    CL_Edge[k]              ; 
    int                         trt =    CL_Treat[k]             ;
    for(int i = 0 ; i<nk; i++){
      int degree = sum(k_Edge(i,_) )                             ;
      
      if( in_id[i] == 1 ) {
        Num_Treat(id) = 1                                        ;
      }else{
        int     num_trt = 0                                      ;
        if(trt==1){
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + CL_Treat[j]*k_Edge(i,j)         ; 
          } 
        }else{
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + (1-CL_Treat[j])*k_Edge(i,j)     ; 
          } 
        }
        Num_Treat(id)  =  num_trt/ degree                        ;
      }
      id  = id + 1          ;
    }
  }
  return Num_Treat ;
}

// [[Rcpp::export]]
NumericVector Match(  NumericVector  CLtreat ,
                      NumericVector  CLid      ) {
  int  n  = CLtreat.length()      ;
  NumericVector TREAT(n)          ; 
  for(int i =0; i<n; i++){
    int id   =  CLid[i]-1;
    TREAT[id] = CLtreat[i]  ; 
  }
  return TREAT;
}

// [[Rcpp::export]]
List Gen_response(   arma::vec    mu            ,
                     arma::vec  alpha           ,
                     arma::vec  CL_Treat        ,
                     List       CL_Edge         , 
                     List       epsilon         ,
                     double    Ave_edges
){
  int                n_CL   =  CL_Edge.size()                         ;
  List                     CL_res( n_CL )                             ;
  
  for( int k = 0 ; k < n_CL ; k ++ ){
    arma::mat  CL_Adj  = CL_Edge[k]                  ;
    int     nk   = CL_Adj.n_rows               ;
    arma::vec            resp(nk)                    ;
    arma::vec    eps_k = epsilon[k]                  ;
    int    treat = CL_Treat[k]                 ; 
    arma::vec    Adj_Treat = CL_Adj * CL_Treat       ;
    arma::vec    Adj_Contr = CL_Adj * (1-CL_Treat)   ;
    
    resp  =  treat*mu[0] + (1-treat)*mu[1] +  (Adj_Treat*alpha[0]+Adj_Contr*alpha[1])/Ave_edges + eps_k ; 
    CL_res[k] = resp;
  }
  
  return CL_res ;
}



// [[Rcpp::export]]
arma::vec  Cal_diff( arma::mat Diff   ){
  int                     n    =  Diff.n_rows                 ;
  arma::mat              DCov  =  cov(Diff)                   ;
  arma::mat              Dinv  =  pinv(DCov)                   ;  
  arma::vec              Imb(n)                               ;  
  
  for( int i  = 0 ; i < n  ; i ++ ){
    arma::rowvec             diff_i = Diff.row(i)             ;
    arma::vec                 Q     = diff_i*Dinv*diff_i.t()  ; 
    Imb(i)   = Q(0)                                           ;
  }
  return Imb     ; 
}

// [[Rcpp::export]]
arma::vec  PSR(      arma::mat X    , double  P   ){
  int                     n    =  X.n_rows        ;
  int                   halfn  =  floor( n/2)     ;
  arma::mat              XCov  = cov(X)           ;
  arma::mat              Xinv  = pinv(XCov)        ;
  arma::vec             Treat(n)                  ;
  double                 ran   = R::rbinom(1,1/2) ; 
  if (ran ==1){ 
    Treat(0) = 0;
    Treat(1) = 1;
  }else  { 
    Treat(0) = 1;
    Treat(1) = 0;
  }
  for( int i  = 1 ; i <= halfn  ; i++){
    arma::colvec          TREAT0(i*2)            ;
    arma::colvec          TREAT1(i*2)            ;
    for( int t=0; t<i;t++)  {     
      TREAT0(2*t)      = Treat(2*t)              ;
      TREAT0(2*t+1)    = Treat(2*t+1)            ;
      TREAT1(2*t)      = Treat(2*t)              ;
      TREAT1(2*t+1)    = Treat(2*t+1)            ;
    }
    TREAT0(2*i-2) = 0;
    TREAT0(2*i-1) = 1;
    TREAT1(2*i-2) = 1;
    TREAT1(2*i-1) = 0;
    arma::mat     Xsele = X.rows(0,2*i-1)       ;
    arma::rowvec  Xdiff0  = (TREAT0.t()*Xsele - (1-TREAT0.t())*Xsele )/i ;
    arma::rowvec  Xdiff1  = (TREAT1.t()*Xsele - (1-TREAT1.t())*Xsele )/i ;   
    arma::rowvec   Imb0   = Xdiff0* Xinv*Xdiff0.t()      ;
    arma::rowvec   Imb1   = Xdiff1* Xinv*Xdiff1.t()      ;
    if(Imb1(0)>Imb0(0)){
      double           ran   = R::rbinom(1,P);
      
      if (ran ==1){ 
        Treat(2*i-2) = 0;
        Treat(2*i-1) = 1;
      }else  { 
        Treat(2*i-2) = 1;
        Treat(2*i-1) = 0;
      }
    }else if(Imb0(0)>Imb1(0)){
      double           ran   = R::rbinom(1,P);
      if (ran ==1){ 
        Treat(2*i-2) = 1;
        Treat(2*i-1) = 0;
      }else  { 
        Treat(2*i-2) = 0;
        Treat(2*i-1) = 1;
      }
    }else if(Imb0(0)==Imb1(0)){
      double           ran   = R::rbinom(1,1/2);
      if (ran ==1){ 
        Treat(2*i-2) = 1;
        Treat(2*i-1) = 0;
      }else  { 
        Treat(2*i-2) = 0;
        Treat(2*i-1) = 1;
      }
    }
    
    
    
  }
  if(n> 2*halfn){
    arma::colvec          TREAT0(n)            ;
    arma::colvec          TREAT1(n)            ;
    for( int t=0; t<(n-1);t++)  {     
      TREAT0(t)      = Treat(t)              ;
      TREAT1(t)      = Treat(t)              ;
    }
    TREAT0(n-1) = 0;
    TREAT1(n-1) = 1;
    arma::rowvec  Xdiff0  = (TREAT0.t()*X - (1-TREAT0.t())*X )/n ;
    arma::rowvec  Xdiff1  = (TREAT1.t()*X - (1-TREAT1.t())*X )/n ;   
    arma::rowvec   Imb0   = Xdiff0* Xinv*Xdiff0.t()      ;
    arma::rowvec   Imb1   = Xdiff1* Xinv*Xdiff1.t()      ;
    if(Imb1(0)>Imb0(0)){
      double           ran   = R::rbinom(1,P);
      if (ran ==1){ 
        Treat(n-1) = 0;
      }else  { 
        Treat(n-1) = 1;
      }
    }else if(Imb0(0)>Imb1(0)){
      double           ran   = R::rbinom(1,P);
      if (ran ==1){ 
        Treat(n-1) = 1;
      }else  { 
        Treat(n-1) = 0;
      }
    }else if(Imb0(0)==Imb1(0)){
      double           ran   = R::rbinom(1,1/2);
      if (ran ==1){ 
        Treat(n-1) = 1;
      }else{ 
        Treat(n-1) = 0;
      }
      
    }
  }
  return Treat       ; 
}

// [[Rcpp::export]]
arma::vec  Adj_Wgts( List Out_In  ,  List Inner_ID   ){
  int           n_CL = Inner_ID.size()      ;
  arma::vec     Weights(n_CL)               ;
  
  for(int k = 0 ; k< n_CL; k ++ ){
    arma::vec k_Inner = Inner_ID[k]      ; 
    double    n_outer = sum(1- k_Inner)  ;
    arma::vec    k_In = Out_In[k]        ;
    
    if(n_outer == 0 ){
      Weights[k] = 0 ;
    }else{
      Weights[k] = sum(k_In)/n_outer ;
    }
  } 
  return Weights ;
}


// [[Rcpp::export]]
List   Inc_Full(    NumericVector     CL_Treat    ,
                    List           Inner_ID       ,
                    List              CL_Edge    ){
  
  int                  n_CL   =  CL_Treat.size()                 ;
  List                           CL_In_Full(n_CL)                ;
  List                           CL_In_Full_Out(n_CL)            ;
  
  
  for( int k = 0 ; k < n_CL ;  k++ ){
    
    NumericVector           in_id   = Inner_ID[k]                ;
    int                        nk   = in_id.size()               ;
    NumericMatrix            k_Edge = CL_Edge[k]                 ; 
    int                         trt = CL_Treat[k]                ;
    NumericVector                   k_In_Full(nk)                ;
    NumericVector                   k_In_Full_Out(nk)            ;
    
    
    for(int i = 0 ; i<nk; i++){
      double degree = sum(k_Edge(i,_) )                          ;
      double Frac_Treat                                         ;  
      if( in_id[i] == 1 ) {
        k_In_Full[i]  = 1                              ;
        k_In_Full_Out[i] = 0                              ;
      }else{
        int     num_trt      = 0                                 ;
        if(trt==1){
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + CL_Treat[j]*k_Edge(i,j)         ; 
          } 
        }else{
          for( int j = 0 ; j< n_CL ; j++ ){
            num_trt  = num_trt + (1-CL_Treat[j])*k_Edge(i,j)     ; 
          } 
        }
        Frac_Treat =  num_trt / degree            ;
        if( Frac_Treat == 1 ) {
          k_In_Full[i]     = 1         ;
          k_In_Full_Out[i] = 1         ;
        }else{
          k_In_Full[i]     = 0         ;
          k_In_Full_Out[i] = 0         ;
        }
        
      }
    }
    CL_In_Full[k]     = k_In_Full                         ;   
    CL_In_Full_Out[k] = k_In_Full_Out                     ;
    
  }
  
  List  L = List::create( Named("In_Full")  = CL_In_Full , Named("Out_Full") = CL_In_Full_Out) ;
  return  L                                      ;
}
