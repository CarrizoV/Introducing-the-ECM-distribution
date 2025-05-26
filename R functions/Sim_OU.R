#Coded by Ricardo Carrizo V
#Version at 04-01-2024
library(matrixStats)
library(abind)
library(purrr)


Sim_OU <- function(N , 
                   d = 2 , 
                   t0 = 0,
                   par = c(sigma = 1 , theta = 1 , z = rep(0,d) ) ,
                   x0 = NULL , 
                   sim_type = "times" ,
                   sim_sett = NULL , 
                   method = "recursive" , 
                   return.t0 = FALSE){
  # INPUTS:
  #
  # N:        Number of trajectories to simulate
  # d:        Spatial dimension.
  # t0:       Initial time.
  
  # par:      Parameters of the OU-process movement.
  #           Either a vector or a list.
  #           If a vector:
    #           It must be a double vector of length 2 or 2 + d:
      #               -length(par) == 2:          par[1] is the standard deviation of the stochastic integral part.
      #                                           par[2] is the rate parameter of the movement.
      #               -if length(par) == 2+d:     par[1] and par[2] are as before.
      #                                           par[2+1:d]: coordinates of a long-term mean position. By default the origin.
  #           If a list:
  #             It must be a list of length 2.
  #             The first element must be a vector of length 2:
  #               par[[1]][1] is the standard deviation of the stochastic integral part.
  #               par[[2]][2] is the rate parameter of the movement.
  #             The second element must be a matrix of size (N,d) whose rows are the centers of each individual's OU process.
  #           Any other format of "par" will either crash or provide a NULL result.
  
  # x0:       Initial position(s) of the OU-process movements.
  #           It can be NULL, a vector or length d, or a matrix of dimensions (N,d)
  #               - if x0 is NULL, all individuals start by default at the origin.
  #               - if x0 is a vector, it must be of length d, and it represents the initial position of all the individuals.
  #               - if x0 is a matrix, is must be of dimensions (N,d), each row containing the initial position of the corresponding individual.
  
  # sim_type: String with the type of simulation required:
  #             "times"   -> if the simulation is over a given time grid
  #             "KL"      -> if the simulation is through the Karhunen-Loève expansion over a compact time interval.
  # sim_sett: Vector containing the main settings of the simulation:
  #             if sim_type = "times" -> sim_sett must be the time grid to simulate over.           
  #             if sim_type = "KL"    -> sim_sett must be a 3-length vector, with,
  #                                         sim_sett[1]:  n_KL, the truncation limit of KL approximations
  #                                         sim_sett[2]:  left limit of the compact interval to simulate over
  #                                         sim_sett[3]:  right limit of the compact interval to simulate over
  
  # method: String with the method to generate the simulation when sim_type = "times".
  #           if method = "Doob", the function will perform a Doob transformation.
  #           if method = "recursive", the function will compute the times evaluation through the corresponding recursive formula (Markov property).
  #
  
  # OUTPUT:
  #
  # If sim_type = "times":  the array X of dimensions (N , length(sim_sett) , d), with:
  #                             X[j,k,l] = value of the lth coordinate of the position of individual j at time sim_sett[k]        
  #
  # If sim_type = "KL":     a list KL of size 2 containing the information of the KL expansion:
  #                             KL[[1]]:   matrix of dimensions (N , n_KL) containing the normalised random coefficients of the KL expansion:
  #                                             KL[[1]][ j , k ] = kth coefficient of the KL expansion of the trajectory of individual j.
  #                             KL[[2]]:   vector of length n_KL containing the standard deviations of the non-normalized random coefficients of the KL expansion.
  #                   
  #                             KL[[3]]:   a (ideally vectorized) function of time, for which at every t, KL[[3]][t] is the n_KL-length vector of values of the orthonormal functions in the KL expansion.
  #
  #(for now, the KL expansion is assumed to be with respect to L2)
  
  #Diverse verifications
  if(is.list(par)){
    #List mode for vectorized version
    #Currently: the center z can be vectorized, but not the other parameters
    if(length(par) != 2){
      print("ERROR: if 'par' is a list it must be of length 2.")
      return(NULL)
    }
    if(!is.vector(par[[1]])){
      print("ERROR: if 'par' is a list, 'par[[1]]' must be a vector.")
      return(NULL)
    }
    if(length(par[[1]]) != 2){
      print("ERROR: if 'par' is a list, 'par[[1]]' must be a vector of length 2.")
      return(NULL)
    }
    sigma <- par[[1]][1]
    theta <- par[[1]][2]
    if(!is.matrix(par[[2]])){
      print("ERROR: if 'par' is a list, 'par[[2]]' must be a matrix.")
      return(NULL)
    }
    if(  any( dim(par[[2]]) != c(N,d) )  ){
      print("ERROR: if 'par' is a list, 'par[[2]]' must be a matrix of dimensions N x d.")
      return(NULL)
    }
    z <- par[[2]]
  }else{
    if(is.vector(par)){
      if(length(par) == 2){ 
        #Case with trivial center and initial condition
        sigma <- par[1]
        theta <- par[2]
        z <- rep(0,d)
      }else if(length(par) == 2 + d){
        sigma <- par[1]
        theta <- par[2]
        z <- par[2+1:d]
      }else{
        print("ERROR: 'par' must be a 2-length vector, a 2+d-length vector, or a list.")
        return(NULL)
      }
    }else{
      print("ERROR: 'par' must be a 2-length vector, a 2+d-length vector, or a list.")
      return(NULL)
    }
  } 

  if(is.null(x0)){
    x0 <- rep( 0 , d)
  }
  if(is.vector(x0)){
    if(length(x0) != d){
      print("ERROR: if 'x0' is a vector it must have length d.")
      return(NULL)
    }
  }else if(is.matrix(x0)){
    if(!all(  dim(x0) == c(N,d) )){
      print("ERROR: if 'x0' is a matrix it must have dimensions (N,d).")
      return(NULL)
    }
  }else{
    print("ERROR: 'x0' must be either a d-length vector or a matrix of dimensions (N,d).")
    return(NULL)
  }
  #End of verifications.
  
  if(sim_type == "times"){
    t <- sim_sett
    n.t <- length(t)
    i <- 1:n.t
    if(is.unsorted(t)){
      i <- order(t)
    }
    
    if(method == "Doob"){
      #Doob transformation method.
      scaled.t <-  exp(  (2*theta)*(t[i]-t0) )- 1
      diff.scaled.t <- diff(   c(t0, scaled.t  )   )
      
      
      X <- array( 0 , dim = c(n.t , N , d)  )
      re.scale.exp <- exp( -theta*(t-t0)  ) 
      #Brownian part
      for( l in 1:d){
        X[ , , l] <- colCumsums( matrix( rnorm( N*n.t , 
                                                sd = sqrt(diff.scaled.t) ) , 
                                         n.t , N) )
      }
      #Transformation
      if(is.matrix(x0)){
        for(l in 1:d){
          if(all(z == 0)){
            X[ , , l] <- X[ , , l]*((sigma/sqrt(2*theta))*re.scale.exp ) + outer(re.scale.exp , x0[ ,l]) 
          }else{
            if(is.matrix(z)){
              X[ , , l] <- X[ , , l]*((sigma/sqrt(2*theta))*re.scale.exp ) + outer(1-re.scale.exp , z[ ,l]) + outer(re.scale.exp , x0[ ,l]) 
            }else{
              X[ , , l] <- X[ , , l]*((sigma/sqrt(2*theta))*re.scale.exp ) + z[l]*(1-re.scale.exp) + outer(re.scale.exp , x0[ ,l]) 
            }
          }
        }
        if(return.t0){
          return(abind( x0 , 
                        aperm(X , c(2,1,3))[ , i ,  , drop = F] , along = 2)   )
        }
        return(aperm(X , c(2,1,3))[ , i ,  , drop = F]  ) 
      }
      
      for(l in 1:d){
        if(all(z==0)){
          if(all(x0==0)){
            X[ , , l] <- X[ , , l]*( (sigma/sqrt(2*theta))*re.scale.exp )
          }else{
            X[ , , l] <- X[ , , l]*( (sigma/sqrt(2*theta))*re.scale.exp ) + x0[l]*re.scale.exp
          }
        }else{
          if(is.matrix(z)){
            X[ , , l] <-   X[ , , l]*( (sigma/sqrt(2*theta))*re.scale.exp ) +  x0[l]*re.scale.exp + outer(1-re.scale.exp , z[ ,l])
          }else{
            X[ , , l] <-   X[ , , l]*( (sigma/sqrt(2*theta))*re.scale.exp ) + ((x0[l] - z[l])*re.scale.exp + z[l])
          }
        }
      }
      if(return.t0){
        return(abind( matrix( rep(x0 , each = N) , N , d      )       , 
                      aperm(X , c(2,1,3))[ , i ,  , drop = F] , along = 2)   )
      }
      return(aperm(X , c(2,1,3))[ , i ,  , drop = F]  ) 
    }
    
    if(method == "recursive"){
      diff.t <- diff( c(t0,t[i]))
      exp.diff.t <- exp( -theta*diff.t  )
      
      X <- array( 0 , dim = c(n.t , N , d )  )
      for(l in 1:d){
        #Noisy part
        if(is.matrix(z)){
          W <- matrix( rnorm(  N*n.t , 
                               mean = outer(1-exp.diff.t , z[,l]) , 
                               sd = sigma*sqrt(  ( 1-exp.diff.t^2 )/(2*theta) ) ) , 
                       n.t , N )
        }else{
          W <- matrix( rnorm( N*n.t , 
                              mean = z[l]*(1-exp.diff.t)  ,
                              sd = sigma*sqrt(  ( 1-exp.diff.t^2 )/(2*theta) )  ) , 
                       n.t , N )
        }
        #Recursive accumulation
        j <- 1
        
        if(is.matrix(x0)){
          L2<- accumulate( lapply(seq_len(n.t), function(i){W[i,]}) ,  
                           function(a,b){ x <- exp.diff.t[j]*a + b 
                           j <<- j + 1 
                           return(x)} , 
                           .init = x0[ ,l] , 
                           .simplify = F )
        }else{
          L2<- accumulate( lapply(seq_len(n.t), function(i){W[i,]}) ,  
                           function(a,b){ x <- exp.diff.t[j]*a + b 
                           j <<- j + 1 
                           return(x)} , 
                           .init = rep(x0[l] , N) , 
                           .simplify = F )
        }
        X[ , , l] <- matrix(unlist(L2), byrow=TRUE, nrow=(n.t+1) )[ 2:(n.t+1) , ]
      }
      if(return.t0){
        if(is.matrix(x0)){
          return(abind( x0 , 
                        aperm(X , c(2,1,3))[ , i , , drop = F] , along = 2)   )
        }
        return(abind( matrix( rep(x0 , each = N) , N , d      )       , 
                      aperm(X , c(2,1,3))[ , i , , drop = F ] , along = 2)   )
      }
      return(aperm(X , c(2,1,3))[ , i ,  , drop = F])
      return(X)
    }
  }
  
  print("ERROR: only the 'times' option is currently available for sim_type")
  return(NULL)
}


