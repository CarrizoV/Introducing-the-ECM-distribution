#Coded by Ricardo CarrizoV
#Version at 08-01-2024
library(matrixStats)
library(mc2d)
library(data.table)
#dBivPois function must have been called!

ECM_Poisson_Pairwise_ll <- function( q , 
                                     lambda , 
                                     P , 
                                     method = "AP" , 
                                     weights = NULL){
  #Format:
  #   q <-  List of vectors to evaluate the log-likelihood at.
  #   lambda <-  rate of the Poisson number of individuals.
  #   P <-  The path-probabilities given in a list format with level >= 2.
  
  #   method <- method of computing a pairwise-based composite likelihood.
  #             Currently available:
  #               -"All-pairs-CL":    CL considering all the pairs.

  #   IMPORTANT: EITHER YOU ASK THE USER TO ADD THE LAST COMPONENTS OF EACH MULTINOMIAL VECTOR (N-sum(q[[k]])), OR YOU SPECIFY A BOOLEAN FOR THE USER TO PREFER THE FORMAT AND YOU DO IT INTERNALLY 
  #              THE CURRENT FORMAT ASSUMES THE COUNTS ARE NOT EXHAUSTIVE!
  #Is up to the user to verify that everything is compatible.
  
  n.t <- length(q) #Number of time-steps
  
  m <- rep(0 , n.t) #Number of categories at each time-step
  for(k in 1:n.t){
    m[k] <- length(q[[k]])
  }
  n.d <- sum(m)  #Number of scalar data values
  
  
  if(n.t == 1){#YOU SHOULD DO SOMETHING PARTICULAR IF n.t EQUALS 1!!!!!
    print("ERROR: ECM_Pairwise_ll function not yet coded for only one time-step.")
    return(NULL)
  }
  
  if(method == "AP"){
    ll <- 0
    ll1 <- 0
    #One time pairs
    for(k in 1:n.t){
      if(m[k] > 1){
        ll1<- ll1 + sum( dpois( combn(q[[k]],2) , 
               combn( lambda*P$'1'[[k]],2) , 
               log = T ) )
      }
    }
    
    #print(paste("ll one-time pairs part:" , ll1))
    #ll <- ll + ll1
    ll2 <- 0
    #two-times pairs
    for(k1 in 1:(n.t-1)){
      for(k2 in (k1+1):n.t){ 
        
        q.pairs <- as.matrix(CJ(q[[k1]],q[[k2]], 
                                sorted = F))
        n.pairs <- m[k1]*m[k2]
        
        p.pairs.1 <- rep( 0 , n.pairs ) 
        p.pairs.2 <- rep( 0 , n.pairs )
        p.pairs.12 <- rep( 0 , n.pairs )
        
        for(l1 in 1:m[k1]){
          p.pairs.1[(l1-1)*m[k2]+1:m[k2]] <- P$'1'[[k1]][l1]-P$'2'[[k1]][[k2]][l1, ]
          p.pairs.2[(l1-1)*m[k2]+1:m[k2]] <- P$'1'[[k2]] - P$'2'[[k1]][[k2]][l1, ]
          p.pairs.12[(l1-1)*m[k2]+1:m[k2]] <- P$'2'[[k1]][[k2]][l1,]
        }
        
        ll2 <- ll2 + sum( dBivPois(  q.pairs , 
                                     q2 = NULL , 
                                     lambda = lambda*cbind( p.pairs.1 , p.pairs.2 , p.pairs.12) ,
                                     log = T) )
      }
    }
    #print(paste("ll two-time pairs part:" , ll2))
    ll <- ll1 + ll2
    return(ll)
  }
  return(NULL)
  
  
  
}
#

