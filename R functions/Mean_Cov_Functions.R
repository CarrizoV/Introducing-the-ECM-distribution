#Coded by Ricardo CarrizoV
#Version at 03-01-2024

#Mean and covariance functions for Gaussian trajectories


#OU stationary process (exponential...) with independent components

#Constant mean function at z
f.mean.OU.stat <- function( sigma , theta , z  , d){
  f <- function(t){
    return(matrix( rep( z , each = length(t)) , 
                   nrow = length(t) , ncol = d ))
  }
  return(f)
}
#Exponential covariance function
f.cov.OU.stat <- function( sigma , theta , z , d ){
  f <- function(t,s){
    C.each.dim <- (sigma^2/(2*theta))*exp(-abs(outer( theta*t , theta*s , 
                                                      '-' ) ) )
    return(  array( rep( C.each.dim , d) , dim = c(  length(t) , length(s) , d )  )   )
  }
  return(f)
}


