#Coded by Ricardo Carrizo V
#Version at 06-11-2024
library(matrixStats)
library(abind)
library(matrixStats)

Sim_Brown <- function( N , 
                       d = 2 , 
                       t0 = 0 , 
                       par = c(sigma = 1 , adv = rep(0,d) )  , 
                       x0  = NULL , 
                       sim_type = "times" , 
                       sim_sett = NULL,
                       return.t0 = FALSE){
  #
  #Function to simulate a Brownian motion (with independent components, for now).
  #
  # INPUTS:
  #
  # N:        Number of trajectories to simulate
  #
  # d:        Spatial dimension.
  #
  # t0:       Initial time.
  #
  # par:      Parameters of the Brownian motion movement. It must have one of the following formats:
  #               -if it is a double:           par is the standard deviation of the Brownian movement.
  #               -if it is a double vector:    par must be of length 1 + d, with
  #                                                   par[1]:       the standard deviation of the Brownian movement.
  #                                                   par[1+1:d]:   a (constant) advection vector.
  #               -if it is a list:             par must be of length (at least) 2, with
  #                                                   par[[1]]:     a double, the standard deviation of the Brownian movement.
  #                                                   par[[2]]:     a vectorized function of time indicating an average deterministic behaviour. 
  #                                                                 This can be used to model a non-stationary advection or even to impose an initial position (make the function constant, not very optimal though). 
  #           Any other format of "par" will either crash or provide a NULL result.
  #
  #
  # sim_type: String with the type of simulation required:
  #             "times"   -> if the simulation is over a given time grid
  #             "KL"      -> if the simulation is through the Karhunen-Loève expansion over a compact time interval.
  #
  # sim_sett: Vector containing the main settings of the simulation:
  #             if sim_type = "times" -> sim_sett must be the time grid to simulate over.           
  #             if sim_type = "KL"    -> sim_sett must be a 3-length vector, with,
  #                                         sim_sett[1]:  n_KL, the truncation limit of KL approximations
  #                                         sim_sett[2]:  left limit of the compact interval to simulate over
  #                                         sim_sett[3]:  right limit of the compact interval to simulate over
  #
  # return.t0:  Boolean indicating if the initial positions must be also returned. By default FALSE. Not applicable if sim_type = "KL".
  
  # OUTPUT:
  #
  # If sim_type = "times":  array X of dimensions (N , length(sim_sett) , d), with:
  #                             X[j,k,l] = value of the lth coordinate of the position of individual j at time sim_sett[k]        
  #                         if return.t0 = TRUE, a matrix with the initial positions would be added as first temporal values in the array X.
  #
  # If sim_type = "KL":     a list KL of size 2 containing the information of the KL expansion:
  #                             KL[[1]]:   matrix of dimensions (N , n_KL) containing the normalised random coefficients of the KL expansion:
  #                                             KL[[1]][ j , k ] = kth coefficient of the KL expansion of the trajectory of individual j.
  #                             KL[[2]]:   vector of length n_KL containing the standard deviations of the non-normalized random coefficients of the KL expansion.
  #                   
  #                             KL[[3]]:   a (ideally vectorized) function of time, for which at every t, KL[[3]][t] is the n_KL-length vector of values of the orthonormal functions in the KL expansion.
  #
  #(for now, the KL expansion is assumed to be with respect to L2)
  if(is.double(par)){
    sigma <- par[1]
    if(length(par) == 1 + d){
      adv <- par[1 + 1:d]
    }else{
      adv <- NULL
    }
  }
  if(is.list(par)){
    sigma <- par[[1]]
    adv <- par[[2]]
  }
  
  if(sim_type == "times"){
    t <- sim_sett
    n.t <- length(t)
    i <- 1:n.t
    if(is.unsorted(t)){
      i <- order(t)
    }
    diff.t <- diff( c(t0,t[i]) )
    
    X <- array( 0 , dim = c(N , n.t , d) )
    if(is.null(adv)){
      for(l in 1:d){
        X[ , , l] <- t(colCumsums( matrix( rnorm( N*n.t ,
                                                  sd = sigma*sqrt(diff.t) ) , 
                                           n.t , N ) ))
      }
    }else{ 
      #Adding the advection
      if(!is.vector(adv) & !is.function(adv)){
        print("ERROR: the advection must be specified as a vector, a function, or NULL.")
        return(NULL)
      }
      b <- matrix( 0 , n.t ,  d  )
      if(is.vector(adv)){
        b <- outer( diff.t , adv )
      }else{#Here the advection is a function
        if(any( adv(t0) != 0 )){
          print("ERROR: the advection function must be null at t0.")
          return(NULL)
        }
        b <- colDiffs(adv( c(t0 , t[i] ) ))
      }
      for( l in 1:d ){
        X[ , , l] <- t(colCumsums(matrix( rnorm( N*n.t , 
                                                 mean =  b[ , l] , 
                                                 sd = sigma*sqrt(diff.t)) ,
                                          n.t , N)))
      }
    }
    #Adding the initial positions
    if(is.null(x0)){
      if(return.t0){
        return( abind( matrix( 0 , N , d)  , X[ , i , ] , 
                       along = 2   )    )
      }
      return(X[ , i , , drop = F ])
    }
    if(is.vector(x0)){
      if(length(x0) != d){
        print("ERROR: x0 must be a d-length vector.")
        return(NULL)
      }
      for(l in 1:d){
        X[ , , l ] <- X[ , , l] + x0[l]
      }
      if(return.t0){
        return(abind( matrix( rep(x0 , each = N) , N , d      )       , 
                      X[ , i , , drop = F ] , along = 2)   )
      }
      return(X[ , i , , drop = F])
    }
    if(is.matrix(x0)){
      if(any( dim(x0) !=  c(N,d)) ){
        print("ERROR: initial positions matrix must be of dimensions (N,d).")
      }
      for(l in 1:d){
        X[ , , l ] <- X[ , , l] + x0[ ,l]
      }
      if(return.t0){
        return(abind( x0 , 
                      X[ , i , , drop = F] , along = 2)   )
      }
      return(X[ , i ,  , drop = F])
    }
  }
  return(NULL)
}

