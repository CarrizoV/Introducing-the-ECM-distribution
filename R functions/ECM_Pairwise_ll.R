#Coded by Ricardo CarrizoV
#Version at 23/02/2025
library(matrixStats)
library(mc2d)
library(data.table)

ECM_Pairwise_ll <- function( q ,  
                             N , 
                             P ,
                             exhaustive = F){
  #Format:
  #   q <-  List of vectors to evaluate the composite log-likelihood at.
  #   N <-  Number of individuals.
  #   P <-  The path-probabilities given in a list format with level >= 2.
  
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
    print("ERROR: ECM_Pairwise_ll function not yet coded for only one time-step. We suggest to use a multinomial log-likelihood.")
    return(NULL)
  }
  
  ll <- 0
  ll1 <- 0
  #One time pairs
  for(k in 1:n.t){
    if(m[k] > 1){
      n.inner.pairs <- choose(m[k] , 2)
      l.index.mat <- matrix( 0 ,  n.inner.pairs , 2 )
      
      ind <- 1
      for(l1 in 1:(m[k]-1)){
        lim <- ind + m[k]-l1-1
        l.index.mat[ ind:lim , 1] <- rep( l1 , m[k]-l1 )
        l.index.mat[ ind:lim , 2] <- (l1+1):m[k]
        ind <- lim + 1
      }
      x <- matrix( q[[k]][ l.index.mat  ] , 
                   nrow = n.inner.pairs , 
                   ncol = 2)
      x <- cbind( x , N - rowSums(x) )
      probs <- matrix( P$'1'[[k]][ l.index.mat  ] , 
                       nrow = n.inner.pairs , 
                       ncol = 2)
      probs <- cbind( probs , 1 - rowSums(probs) )
      #Apparently a regularization is required
      probs[probs < 1e-128] <- 1e-128
      
      ll1 <- ll1 + sum(dmultinomial( x = x , 
                                     size = N , 
                                     prob = probs , 
                                     log = T))
    }
  }
  #print(paste("ll one-time pairs part:" , ll1))
  ll <- ll + ll1
  #two-times pairs
  ll2 <- 0
  for(k1 in 1:(n.t-1)){
    for(k2 in (k1+1):n.t){ 
      
      q.pairs <- as.matrix(CJ(q[[k1]],q[[k2]], 
                              sorted = F))
      
      p1.p2.pairs <-as.matrix(CJ(P$'1'[[k1]],P$'1'[[k2]], 
                                 sorted = F))
      p12.pairs <- as.vector(t(P$'2'[[k1]][[k2]]))
      p.pairs <- cbind( p12.pairs ,  
                        p1.p2.pairs - p12.pairs , 
                        1 - (rowSums(p1.p2.pairs) - p12.pairs) )
      
      j.max <- rowMins(q.pairs) 
      ind.j.max.0 <- which( j.max == 0  )
      
      if(length(ind.j.max.0)>0){
        #Adicionar los valores de los pares en donde j.max es 0. 
        #En tal caso no es necesario llamar a la función dmultinomial.
        q.pairs.0 <- q.pairs[ind.j.max.0 , ]
        p.pairs.0 <- p.pairs[ind.j.max.0 , ]
        
        #IMPORTANT CURRENT PROBLEM
        #Because of informatic details, sometimes some of the probabilities in p.pairs.0 are informatically 0, therefore allowing the possibility to explode.
        #Easy and fast way out: those are changed to be a small number here.
        p.pairs.0[p.pairs.0 < 1e-128] <- 1e-128
        
        ll2 <- ll2 + sum( lfactorial(N)-
                            lfactorial(N-rowSums(q.pairs.0))-
                            lfactorial(q.pairs.0[,1])-
                            lfactorial(q.pairs.0[,2])+
                            q.pairs.0[,1]*log(p.pairs.0[,2])+ 
                            q.pairs.0[,2]*log(p.pairs.0[,3])+
                            (N - rowSums( q.pairs.0 ))*log( p.pairs.0[,4] ) )
        #Esto puede simplificarse si separamos los casos en que ambos son nulos y distinguir cuál de los dos es nulo.
        #Es probable que el tiempo ganado sea despreciable con respecto a otras partes, así que se deja así por ahora.
      }
      
      if(length(j.max) - length(ind.j.max.0)>0){
        #Now the complicated ones
        if(length(ind.j.max.0) > 0){
          q.pairs <- q.pairs[ -ind.j.max.0,  , drop = F ]
          p.pairs <- p.pairs[ -ind.j.max.0,  , drop = F ]
          j.max <- j.max[-ind.j.max.0] 
        }
        
        n.pairs <- nrow(q.pairs)
        j.min <- pmax( 0 ,  rowSums(q.pairs) - N )
        
        j.lengths <- (j.max-j.min)+1
        n.total <- sum(j.lengths)
    
        j.total <- rep( 0 , n.total )
        for(i in 1:n.pairs){
          j.total[(sum(j.lengths[1:i])-j.lengths[i])+1:j.lengths[i]] <- j.min[i]:j.max[i]
        }
        
        q.total <-cbind(  j.total , 
                          rep(q.pairs[,1] , times=j.lengths)-j.total , 
                          rep(q.pairs[,2] , times=j.lengths)-j.total , 
                          rep(N-rowSums(q.pairs) , times=j.lengths)+j.total )
        p.total <- cbind( rep( p.pairs[,1] , times = j.lengths ) , 
                          rep( p.pairs[,2] , times = j.lengths ) , 
                          rep( p.pairs[,3] , times = j.lengths ) , 
                          rep( p.pairs[,4] , times = j.lengths ))
        #Regularization to avoid all-zero or negative probabilities scenarios (the all-zeros should never happen!).
        p.total[p.total < 1e-128] <- 1e-128
        
        log.mult <- dmultinomial( x = q.total , 
                                  size = N , 
                                  prob = p.total , 
                                  log = T) 
        
        lses <- rep(0 , n.pairs)
        for(i in 1:n.pairs){
          lses[i] <- logSumExp(log.mult[(sum(j.lengths[1:i])-j.lengths[i])+1:j.lengths[i]])
        }
        ll2 <- ll2 + sum(lses)
      }
    }
  }
  #print(paste("ll two-time pairs part:" , ll2))
  ll <- ll + ll2
  return(ll)
}