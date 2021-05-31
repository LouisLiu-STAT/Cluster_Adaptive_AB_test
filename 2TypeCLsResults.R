########################Chapter 2###############################################
library(igraph)
################################################################################
Adj.matrix<-function( n_node_cluster, In_cluster_p,Be_cluster_p){
  In.cluster.edges<-  function(n_cluster, In_cluster_p){
    Adj_Incluster <-  rbinom(n_cluster^2, 1, In_cluster_p)   %>%  matrix(nrow = n_cluster, ncol = n_cluster)
    Adj_Incluster[lower.tri(Adj_Incluster,diag = T) ]<- 0 
    Adj_Incluster <-   Adj_Incluster + t(Adj_Incluster)
    return(Adj_Incluster)
  }
  n_node     <- sum(n_node_cluster)
  adj.matrix <- rbinom(n_node^2,1,Be_cluster_p) %>%matrix(nrow= n_node, ncol = n_node )
  adj.matrix[lower.tri(adj.matrix,diag = T) ]<- 0 
  adj.matrix <-   adj.matrix + t(adj.matrix)
  n_node_cluster<-c(0,n_node_cluster )
  n_node_cluster <- n_node_cluster%>%cumsum()
  for (k in c(1:( length(n_node_cluster)-1))) {
    ID <- c((n_node_cluster[k]+1) : n_node_cluster[k+1])
    adj.matrix[ c(ID), c(ID)  ] <- In.cluster.edges(n_node_cluster[k+1]- n_node_cluster[k], In_cluster_p[k])
  }  
  return(adj.matrix)
  
}


set.seed(19960503)
## Two groups with not only few connection between groups
n     <-  50    # the largest number of nodes per cluster
k1    <-  0.8
k2    <-  0.6 
p1    <-  0.5
p2    <-  0.3
n_node_cluster <-  c( k1 , k1 , k2 , k2)*n
n.node         <-  sum(n_node_cluster)    
n_nodes.c      <-  cumsum(c(0,n_node_cluster) )
n.cluster      <-  length(n_node_cluster)

In_cluster_p   <- c( p1 , p1 , p2 , p2 )
Be_cluster_p   <- 1/180
cluster   <- c( rep(1,n_node_cluster[1]) ,  rep(2,n_node_cluster[2]) ,   rep(3,n_node_cluster[3]) ,    rep(4,n_node_cluster[4]))
color     <- c("#cbd5e8","#f4cae4","#e6f5c9","#fff2ae")
color.cluster  <- color[cluster]
adj.mat  <- Adj.matrix(n_node_cluster, In_cluster_p, Be_cluster_p)
G    <- graph_from_adjacency_matrix( adj.mat, mode =  "undirected") 
l    <- layout.kamada.kawai(G)

pdf("2types_network.pdf",  width=10, height=10)
### Plot the clusters           
plot(G  , vertex.color = color.cluster     , vertex.size = 8 , layout= l,  vertex.label = "")  
####### Generate the correct potential outcomes
dev.off()


###########################################################################################
################### Evaluate  Estimation  #################################################
###########################################################################################
mu        <-  c(1.6,1)
alpha     <-  c(1.4,1)
sig       <-  1

CL.adj    <- lapply( 1 : n.cluster ,  function(k) sapply(1: n.cluster , function(x) rowSums( adj.mat[ which(cluster==k),which(cluster==x)] ) ) )
n.deg     <- unlist(lapply(CL.adj, rowSums) )
Ave.edges <- sum(  n.deg  )/n.node
epsilon   <- rnorm(n.node,0, sig)
sd.deg    <- n.deg/Ave.edges
PO.Y0     <- mu[2] + alpha[2]* sd.deg + epsilon
PO.Y1     <- mu[1] + alpha[1]* sd.deg + epsilon

CL.barY1  <-  sapply( 1:n.cluster  , function(k)  mean(   PO.Y1[which(cluster==k)] ) )
CL.barY0  <-  sapply( 1:n.cluster  , function(k)  mean(   PO.Y0[which(cluster==k)] ) )


IN.TAU        <-   mean(PO.Y1-PO.Y0)
CL.TAU        <-   mean(CL.barY1-CL.barY0)
Tau           <-  c( IN.TAU, CL.TAU ) 
names(Tau)    <-  c("In","Cl" )
################################################

CL.adj        <- lapply( 1 : n.cluster ,  function(k) sapply(1: n.cluster         , function(x) apply( adj.mat[ which(cluster==k),which(cluster==x)], 1, sum ) ) ) 
node.In         <- lapply( 1 : n.cluster ,  function(k) sapply(1: nrow(CL.adj[[k]]) , function(x)     1*( sum(CL.adj[[k]][x,-k]  )==0       )))
node.de         <- lapply( 1 : n.cluster ,  function(k) sapply(1: nrow(CL.adj[[k]]) , function(x) ifelse( sum(CL.adj[[k]][x,-k]>0)==1, 1 ,0 )) )
Inclusion       <- lapply( 1 : n.cluster ,  function(k) matrix(0, nrow =sum(cluster==k) , ncol= n.cluster))

for (k  in c(1:n.cluster)  ) {
  Inclusion[[k]][,k]<-node.In[[k]]
  Inclusion[[k]][which(node.de[[k]]==1),-k] <- 1*(CL.adj[[k]][which(node.de[[k]]==1),-k]>0)
}



##assign 6 possible treatment for the clusters
Treat       <- cbind( c(1,1,0,0) , c(0,0,1,1) ,  c(1,0,1,0) ,
                      c(0,1,0,1) , c(1,0,0,1) ,  c(0,1,1,0) )

    In.nodes <- lapply(node.In, sum) %>%unlist()
    Ot.nodes <- lapply(node.de, sum) %>%unlist()
In.ave.edges<-  sapply(1:n.cluster , function(x)   mean( apply(CL.adj [[x]][which(node.In[[x]]==1),],1,sum  )))
de.ave.egess<-  sapply(1:n.cluster , function(x)   mean( apply(CL.adj [[x]][which(node.de[[x]]==1),],1,sum )))
  
  
  
Treat.eval <- function(x){
  CL.Treat        <- Treat[,x]
  #### Generate Individual level treatment
  Res.gen <- function(k){
    Adj.Treat <-  CL.adj[[k]]%*%  CL.Treat
    Adj.Contr <-  CL.adj[[k]]%*%(1-CL.Treat) 
      P.Y     <- mu[1]*CL.Treat[k] + mu[2]*(1-CL.Treat[k]) +  ( Adj.Treat*alpha[1]+Adj.Contr*alpha[2])/Ave.edges + epsilon[cluster==k]  
      return(P.Y)
  }
  #### Generate response according to the observed treatment
  P.Y       <- lapply(1:n.cluster, Res.gen)
  T.cluster <-   CL.Treat *c(1:n.cluster)
  T.cluster <- T.cluster[ which(T.cluster!= 0) ]
  C.cluster <- (1-CL.Treat)*c(1:n.cluster)
  C.cluster <- C.cluster[ which(C.cluster!= 0) ]
  HT.ind    <- rep(0,n.cluster)
  HT.cl     <- rep(0,n.cluster)
  
  for (k in c(1: n.cluster)) {
    P.Yx    <-   P.Y[[k]]
    P.Yx1   <-   P.Yx[ which( node.In[[k]]==1 )]   
    if(  CL.Treat[k]==1){
      if (sum( Inclusion[[k]][,T.cluster[which(T.cluster!=k)]]==1 ) ==0 ) {
        P.Yx2   <- 0  
      }else{
        P.Yx2   <-   P.Yx[ which( Inclusion[[k]][,T.cluster[which(T.cluster!=k)]]==1)  ]
      }
    }else{ 
      if(sum( Inclusion[[k]][,C.cluster[which(C.cluster!=k)] ]==1)==0){
        P.Yx2  <- 0
            }else{
        P.Yx2  <-    P.Yx[ which( Inclusion[[k]][,C.cluster[which(C.cluster!=k)] ]==1)]
            }
    }
    HT.ind[k]  <-   sum(P.Yx1)*2 + sum(P.Yx2)*6
    HT.cl[k]   <-   HT.ind[k]  /  length(P.Yx)
  }
  HT.est.ind <-  ( sum( CL.Treat*HT.ind ) - sum( (1-CL.Treat)*HT.ind )  ) / n.node
  HT.est.cl  <-  ( sum( CL.Treat*HT.cl  ) - sum( (1-CL.Treat)*HT.cl  )  ) / n.cluster
    re  <- c(     HT.est.ind  , HT.est.cl   )
  names(re) <- c("HT.In","HT.CL")
  
  return(re)
} 
 node.eval <- function(x){
  
  trt        <- Treat[,x]
  imb.in.node <-   sum(trt*In.nodes)/2- sum((1-trt)*In.nodes)/2
  imb.ot.node <-   sum(trt*Ot.nodes)/2- sum((1-trt)*Ot.nodes)/2
  imb.in.edge <-   sum(trt* In.ave.edges)/2- sum((1-trt)*In.ave.edges)/2   
  
  T.cluster <- trt*c(1:n.cluster)
  T.cluster <- T.cluster[ which(T.cluster!= 0) ]
  C.cluster <- (1-trt)*c(1:n.cluster)
  C.cluster <- C.cluster[ which(C.cluster!= 0) ]
  T.use.node <-  lapply(1:length(T.cluster), function(x) Inclusion[[T.cluster[x]]][,T.cluster[-x]]  ) 
  C.use.node <-  lapply(1:length(C.cluster), function(x) Inclusion[[C.cluster[x]]][,C.cluster[-x]]  ) 
  T.N.U.node <- lapply(T.use.node, sum) %>% unlist()
  C.N.U.node <- lapply(C.use.node, sum) %>% unlist()
  imb.u.nodes <- sum(T.N.U.node) /2- sum(C.N.U.node)/2  
  
  
  Cal.edges <-  function(X ){
    
    if(length(X)==0 ) {
      re = 0
    }else if(is.null(ncol(X))){
      re =   sum(X)
    }else {
      re = mean(apply(X,1, sum) )
    }
    return(re)
  }
  
  T.u.ave.edge<- sapply(1:length(T.cluster), function(x)  Cal.edges(CL.adj[[T.cluster[x]]][ which( T.use.node[[x]]==1),]))
  C.u.ave.edge<- sapply(1:length(C.cluster), function(x)  Cal.edges(CL.adj[[C.cluster[x]]][ which( C.use.node[[x]]==1),]))
  
  imb.u.ave.edge<-sum(T.u.ave.edge)/2 - sum(C.u.ave.edge)/2
  re<- c(  imb.in.node, imb.ot.node, imb.in.edge,imb.u.nodes ,imb.u.ave.edge)
  return(re)
}

RE.eval<-sapply(1:ncol(Treat),Treat.eval )

colnames(RE.eval) <- c("(1,1,0,0)","(0,0,1,1)", "(1,0,1,0)", "(0,1,0,1)", "(1,0,0,1)","(0,1,1,0)")


Re.mean     <- apply(RE.eval, 1, mean)
Re.var <- apply(RE.eval, 1, function(X)  sum(X^2)/6 - (sum(X)/6 )^2  )
RE.summary<-cbind( bias =  Re.mean -Tau   ,var = Re.var )
rownames(RE.summary) <- rownames(RE.eval)


Imb<-  sapply(1:ncol(Treat), node.eval )

rownames(Imb)<- c("Inner nodes","All outer nodes" , "Ave Inner edges","Used outer nodes","Ave used outer edges")
colnames(Imb)<-  c("(1,1,0,0)","(0,0,1,1)", "(1,0,1,0)", "(0,1,0,1)", "(1,0,0,1)","(0,1,1,0)")

rbind( RE.eval[,c(2,4,6)], Imb[c(1,3,4,5),c(2,4,6)] )%>%xtable(digits = 3)
