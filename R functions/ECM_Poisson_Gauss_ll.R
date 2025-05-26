#Coded by Ricardo Carrizo V
#Version at 25/02/2025

ECM_Poisson_Gauss_ll <- function(  q , 
                           lambda , 
                           P , 
                           eps = 1e-6){
  m <- sapply( q , length)
  n.d <- sum(m) 
  P.vec <- unlist(P$'1') #lambda=1 mean vector
  P.mat <- matrix( 0 , nrow = n.d , ncol = n.d  ) #lambda=1 covariance matrix (initialization)
  for(k1 in 1:n.t){
    P.mat[ sum(m[0:(k1-1)]) + 1:m[k1] , 
           sum(m[0:(k1-1)]) + 1:m[k1] ] <- base::diag(P.vec[sum(m[0:(k1-1)]) + 1:m[k1]], 
                                                      ncol = m[k1], 
                                                      nrow = m[k1])
    if(k1 < n.t){
      for(k2 in (k1+1):n.t){
        P.mat[ sum(m[0:(k1-1)]) + 1:m[k1] , 
               sum(m[0:(k2-1)]) + 1:m[k2] ] <- P$'2'[[k1]][[k2]]
        
      }
    }
  }
  #symetrization
  P.mat[lower.tri(P.mat)] = t(P.mat)[lower.tri(P.mat)]
  
  #Regularization
  diag(P.mat) <- diag(P.mat) + eps
  
  #Gaussian log-likelihood
  ll <- dmnorm( x = (unlist(q) - lambda*P.vec)/sqrt(lambda) , 
                varcov = P.mat  , 
                log = T)
  return(ll)
}