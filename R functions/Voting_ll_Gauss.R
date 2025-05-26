#Coded by Ricardo CarrizoV
library(matrixStats)

Voting_ll_Gauss <- function( P2.1 , Q1 , Q2 ){
  #Number of districts
  n.dist <- nrow(Q1)
  
  #Determinants of the covariance matrix at each district
  det.Covs <- (Q1%*%(P2.1[1 , ]-P2.1[1 , ]^2))*(Q1%*%(P2.1[2 , ]-P2.1[2 , ]^2))- 
    ( Q1%*%( P2.1[1 , ]*P2.1[2 , ]   ) )^2
  
  #Inverse covariance matrices at each district
  inv.Covs <- array( 0 , dim = c(n.dist , 2 , 2)   )
  for(j in 1:n.dist){
    inv.Covs[j , , ] <- rbind( c(sum(Q1[j , ]*P2.1[2, ]*(1-P2.1[2, ])) , sum(Q1[j , ]*P2.1[1, ]*P2.1[2, ])) , 
                               c( 0 , sum(Q1[j , ]*P2.1[1, ]*(1-P2.1[1, ])))  )
    inv.Covs[j , 2 , 1 ] <- inv.Covs [j , 1 , 2 ]
  }
  inv.Covs <- inv.Covs/as.vector(det.Covs)
  
  #Centred values at each district
  Z <- Q2 - Q1%*%t(P2.1)
  
  #Log-likelihoods at each district
  lls <- rep( 0 , n.dist)
  for(j in 1:n.dist){
    lls[j] <- -0.5*sum( (inv.Covs[j , , ]%*%Z[j , ])*Z[j , ] )  
  }
  lls <- lls - 0.5*log(det.Covs)
  ll <- sum(lls) - n.dist*log(2*base::pi)
  print(ll)
  return(ll)
}

Voting_ll_Gauss_Gradient <- function(P2.1 , Q1 , Q2){
  #Number of districts
  n.dist <- nrow(Q1)
  
  #Determinants of the covariance matrix at each district
  det.Covs <- (Q1%*%(P2.1[1 , ]-P2.1[1 , ]^2))*(Q1%*%(P2.1[2 , ]-P2.1[2 , ]^2))- 
    ( Q1%*%( P2.1[1 , ]*P2.1[2 , ]   ) )^2
  
  #Inverse covariance matrices at each district
  inv.Covs <- array( 0 , dim = c(n.dist , 2 , 2)   )
  for(j in 1:n.dist){
    inv.Covs[j , , ] <- rbind( c(sum(Q1[j , ]*P2.1[2, ]*(1-P2.1[2, ])) , sum(Q1[j , ]*P2.1[1, ]*P2.1[2, ])) , 
                               c( 0 , sum(Q1[j , ]*P2.1[1, ]*(1-P2.1[1, ])))  )
    inv.Covs[j , 2 , 1 ] <- inv.Covs [j , 1 , 2 ]
  }
  inv.Covs <- inv.Covs/as.vector(det.Covs)
  
  m <- ncol(P2.1)
  Grads <- array( 0 , dim = c(n.dist , 2 , m  )  )
  for(j in 1:n.dist){
    z <- inv.Covs[j, , ]%*%(Q2[j,] - as.vector(P2.1%*%Q1[j,]))
    for(l in 1:m){
      der.mean <- c( Q1[j,l] , 0 )
      der.Cov <- Q1[j,l]*rbind(  c( 1-2*P2.1[1,l] , -P2.1[2,l] ) , 
                                 c( -P2.1[2,l] , 0 ) )
      Grads[j,1,l] <- 0.5*sum( (der.Cov%*%z)*z ) + sum(der.mean*z) - 0.5*sum(diag(inv.Covs[j, , ]%*%der.Cov))
      
      der.mean <- c( 0 , Q1[j,l] )
      der.Cov <- Q1[j,l]*rbind(  c( 0 , -P2.1[1,l] ) , 
                                 c( -P2.1[1,l] , 1-2*P2.1[2,l] ) )
      Grads[j,2,l] <- 0.5*sum( (der.Cov%*%z)*z ) + sum(der.mean*z) - 0.5*sum(diag(inv.Covs[j, , ]%*%der.Cov))
    }
  }
  grad <- apply( Grads , c(2,3) , sum )
  return(grad)
}

inv.softmax1 <- function(p , lb = NULL ){
  if(!is.null(lb)){
    p <- (p - lb)/(1-lb)
  }
  return(  log(p) - log1p(-sum(p))  )
}


softmax1 <- function( x , lb = NULL){
  s <- exp(x-logSumExp(c(x,0)))
  if(!is.null(lb)){
    return(lb + s*(1-lb) )  
  }
  return(s)         
}


Jac.Trans.Softmax1 <- function( x , 
                                lb = NULL){
  #To transform from 'x' to 'p'
  p.scale <- exp(x-logSumExp(c(x,0)))
  if(!is.null(lb)){
    p <- p.scale*(1-lb) + lb
    Jac <- t(matrix( lb - p , length(x) , length(x)))
  }else{
    p <- p.scale
    Jac <- t(matrix( - p , length(x) , length(x)))
  }
  diag(Jac) <- 1 - p
  Jac <- Jac*p.scale 
  return(Jac)
}

