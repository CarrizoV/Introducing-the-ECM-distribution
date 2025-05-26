#Coded by Ricardo CarrizoV
#Version at 22-01-2024

#Update 22-01-2024: transformed the mixture function and put a simple for loop. Apparently I was using mapply wrongly.
#Update 21-01-2024: added correction for avoiding abs(r) > 1 in P_Snapshot

library(pbivnorm)
library(mnormt)

P_Snapshot <- function( level = 1 , 
                        p , 
                        trajectory = "Brown" , theta.X = c(1) , 
                        t0 = 0 , t  , 
                        d , B.L , dx = NULL , Lengths = NULL , 
                        reg = 0 , 
                        unit.correction = F , reg.unit = 0 , 
                        consistency.correction = T){
  
  n.t <- length(t)
  if(level < 1 || level > n.t){
    print("ERROR: level of computed P-arrays must be at least 1 and maximum length(t).")
    return(NULL)
  }
  if(level > 2){
    print("ERROR: P_Snapshot function currently available for level <= 2.")
  }
  if(trajectory != "Brown" & trajectory != "OU" & trajectory != "Gauss"){
    print("P_Snapshot function is currently available only for 'Brown', 'OU' and 'Gauss' trajectories.")
    return(NULL)
  }
  
  m <- rep( 0 , n.t) #vector containing the number of count regions at each time.
  for(k in 1:n.t){
    m[k] <- nrow(B.L[[k]])
  }
  
  #Fixing the parameters
  if(trajectory == "Brown"){
    if(length(theta.X) < 1){
      print("ERROR: Brownian motion process requires at least the parameter 'sigma'.")
      return(NULL)
    }
    sigma <- theta.X[1]
    if(length(theta.X) < 1 + d ){
      v <- rep( 0 , d)
    }else{
      v <- theta.X[1 + 1:d]
    }
    if(length(theta.X) < 1 + 2*d){
      x0 <- rep(0,d)  
      #print("WARNING: P_Snapshot function with 'Brown' trajectory not yet adapted to not-null initial positions (assumed null in this instance).")
    }else{
      x0 <- theta.X[1 + d + 1:d]
    }
    t.trans <- t - t0 #Translation of t
  }
  if(trajectory == "OU"){
    if(length(theta.X) < 2){
      print("ERROR: OU process requires at least the parameters 'sigma' and 'theta'.")
      return(NULL)
    }
    sigma <- theta.X[1]
    theta <- theta.X[2]
    if(length(theta.X) < 2 + d ){
      z <- rep( 0 , d)
    }else{
      z <- theta.X[2 + 1:d]
    } 
    if(length(theta.X) < 2 + 2*d){
      x0 <- rep(0,d)
    }else{
      x0 <- theta.X[2 + d + 1:d]
    }
    t.trans <- t - t0 #Translation of t
  }
  if(trajectory == "Gauss"){
    Mean <- theta.X$mean(t) #(matrix length(t) x d)
    Cov <- theta.X$cov(t,t) #(array  length(t) x length(t) x d x d  )
    
    #There should be a convention for d = 1, watch out!
    
    # if(theta.X$indep.components){
    #   Mean <- theta.X$mean(t)
    #   Cov <- theta.X$cov(t,t)
    # }
    # else{
    #   print("ERROR: P_Snapshot function with 'Gauss' trajectory currently accepts independent components only.")
    #   return(NULL)
    # }
  }
  
  P <- vector("list" , length = level)
  names(P) <- 1:level
  if(level == 1){
    P$'1' <- vector("list" , length = n.t)
    if(trajectory == "Brown"){
      for(k in 1:n.t){
        pr <- rep( p , m[k] )
        for(l in 1:d){
          pr <- pr*( pnorm( B.L[[k]][,l] + dx , mean = x0[l]+v[l]*t.trans[k] , sd = sigma*sqrt(t.trans[k])  )-
                       pnorm( B.L[[k]][,l] , mean = x0[l]+v[l]*t.trans[k] , sd = sigma*sqrt(t.trans[k]) ) )
        }
        P$'1'[[k]] <- pr
      }
    }
    if(trajectory == "OU"){
      for(k in 1:n.t){
        pr <- rep( p , m[k] )
        for(l in 1:d){
          pr <- pr*(pnorm(B.L[[k]][ , l] + dx , 
                          mean = z[l] + (x0[l]-z[l])*exp(-theta*(t.trans[k]))  , 
                          sd = sqrt((1-exp(-2*theta*(t.trans[k])) )*(sigma^2)/(2*theta) ) )-
                      pnorm(B.L[[k]][ , l] , 
                            mean = z[l] + (x0[l]-z[l])*exp(-theta*(t.trans[k]))  , 
                            sd = sqrt((1-exp(-2*theta*(t.trans[k])) )*(sigma^2)/(2*theta) ) ) )
        }
        P$'1'[[k]] <- pr
      }
    }
    if(trajectory == "Gauss"){
      
      if(theta.X$indep.components || d == 1 ){
        for(k in 1:n.t){
          pr <- rep( p , m[k] )
          #Important: we are assuming a very particular form of delivering the mean and covariance functions.
          #The mean is intuitive.
          #The current code only works if the components of the trajectory are independent.
          #       (if not, then one must use pbivnorm or another one depending on the dimension)
          for(l in 1:d){
            pr <- pr*(pnorm(B.L[[k]][ , l] + dx , 
                            mean = Mean[k,l] , 
                            sd = sqrt(Cov[k,k,l]) )-
                        pnorm(B.L[[k]][ , l] , 
                              mean = Mean[k,l]  , 
                              sd = sqrt(Cov[k,k,l]) ) )
          }
          P$'1'[[k]] <- pr
        }
      }else{
        if(d == 2){
          for(k in 1:n.t){
            pr <- rep( p , m[k] )  
            sd.x <- sqrt(Cov[k,k,1,1])
            sd.y <- sqrt(Cov[k,k,2,2])
            r <- Cov[k,k,1,2]/(sd.x*sd.y)
            mean.x <- Mean[k,1]
            mean.y <- Mean[k,2]
            
            pr <- pr*(pbivnorm( (B.L[[k]][ , 1] + dx - m.x)/sd.x , 
                                (B.L[[k]][ , 2] + dx - m.y)/sd.y ,    
                                rho = r ) - 
                        pbivnorm( (B.L[[k]][ , 1] - m.x)/sd.x , 
                                  (B.L[[k]][ , 2] + dx - m.y)/sd.y ,    
                                  rho = r ) - 
                        pbivnorm( (B.L[[k]][ , 1] + dx - m.x)/sd.x , 
                                  (B.L[[k]][ , 2] - m.y)/sd.y ,    
                                  rho = r ) + 
                        pbivnorm( (B.L[[k]][ , 1] - m.x)/sd.x , 
                                  (B.L[[k]][ , 2] - m.y)/sd.y ,    
                                  rho = r ))
            P$'1'[[k]] <- pr
          }
        }
        if(d > 2){
          for(k in 1:n.t){
            for(l in 1:m[k]){
              P$'1'[[k]][l] <- p*sadmvn(  lower = B.L[[k]][l , ] , 
                                          upper = B.L[[k]][l , ] + dx , 
                                          mean = Mean[k , ] , 
                                          varcov = Cov[k,k, , ] )
            }
          }
        }
      }
    }
    #The zero regularization
    if(reg!=0){
      for(k in 1:n.t){
        ind.reg <- which( P$'1'[[k]] < reg )
        if(any(ind.reg)){
          P$'1'[[k]][ind.reg] <- reg
        }
      }
    }
    #The unitarity correction: everything must sum up to 1 - reg.unit maximum.
    if(unit.correction){
      for(k in 1:n.t){
        s <- sum(P$'1'[[k]])
        if(s > 1 - reg.unit){
          P$'1'[[k]] <- P$'1'[[k]]*( (1-reg.unit)/s )
        }
      }
    }
    return(P)
  }
  if(level == 2){
    P$'1' <- vector("list" , length = n.t)
    P$'2' <- vector("list" , length = n.t-1)
    for(k in 1:(n.t-1)){
      P$'2'[[k]] <- vector("list" , length = n.t)
    }
    
    if(trajectory == "Brown"){
      for(k1 in 1:(n.t-1)){
        for(k2 in (k1+1):n.t){
          #Computing the probabilities within the matrix
          pr.mat <- matrix( p^2 , nrow =  m[k1] , ncol = m[k2])
          r <- sqrt(t.trans[ k1 ]/t.trans[ k2 ])
          if( abs(r) > 1){#For corrections
            r <- sign(r)
          }
          for(l in 1:d){
            low.b.x.nn <-  rep( B.L[[ k1 ]][,l] - x0[l] ,  m[k2]  )
            low.b.y.nn <-  rep( B.L[[ k2 ]][,l] - x0[l] , each = m[k1]  )
            low.b.x <- ( low.b.x.nn - v[l]*t.trans[ k1 ] )/( sigma*sqrt(t.trans[ k1 ]) )
            up.b.x <- low.b.x + dx/( sigma*sqrt(t.trans[ k1 ]) )
            low.b.y <- ( low.b.y.nn - v[l]*t.trans[ k2 ] )/( sigma*sqrt(t.trans[ k2 ]) )
            up.b.y <- low.b.y + dx/( sigma*sqrt(t.trans[ k2 ]) )
            
            pr.mat <- pr.mat*matrix( 
              pbivnorm(  up.b.x , up.b.y , rho =  r   ) - 
                pbivnorm(  low.b.x , up.b.y , rho =  r   ) -
                pbivnorm(  up.b.x , low.b.y , rho =  r   ) +
                pbivnorm(  low.b.x , low.b.y , rho = r   ), 
              m[k1] , 
              m[k2])  
          }
          #regularising if required
          if(reg != 0){
            erratic.too.small <- which(pr.mat < reg)
            erratic.nan <- which(is.nan(pr.mat))
            erratic <- union(erratic.too.small , erratic.nan)
            if( length(erratic) > 0 ){
              pr.mat[erratic] <- reg
            }
          }
          #saving
          P$'2'[[ k1 ]][[ k2 ]] <- pr.mat
        }
      }  
      
      #Computing the one-time vectors
      for(k in 1:n.t){
        pr.vec <- rep( p , nrow(B.L[[k]])  )
        for(l in 1:d){
          pr.vec <- pr.vec*( pnorm( B.L[[k]][,l] + dx , mean = x0[l] + v[l]*t.trans[k] , sd = sigma*sqrt(t.trans[k])  )-
                               pnorm( B.L[[k]][,l] , mean = x0[l] + v[l]*t.trans[k] , sd = sigma*sqrt(t.trans[k]) ) )
        }
        if(consistency.correction){
          #We verify the inequalities these vectors should satisfy with respect to two-times matrices
          for(k2 in 1:n.t){
            if(k2 < k){
              pr.vec <- pmax( pr.vec , colSums( P$'2'[[k2]][[k]] )  )
            }
            if(k2 > k){
              pr.vec <- pmax( pr.vec , rowSums( P$'2'[[k]][[k2]] )  )
            }
          }
        }
        #saving
        P$'1'[[k]] <- pr.vec
      }
    }
    if(trajectory == "OU"){
      #Computing the two-times matrices
      for(k1 in 1:(n.t-1)){
        for(k2 in (k1+1):n.t){
          pr.mat <- matrix( p^2 , nrow =  m[k1] , ncol = m[k2])
          
          r <- sqrt( exp(-2*theta*(t.trans[k2] - t.trans[k1])) - exp(-2*theta*t.trans[k2] ) )
          r <- r/sqrt( 1 - exp(-2*theta*t.trans[k2] ) )
          if(abs(r)>1){
            r <- sign(r)
          }
          for(l in 1:d){
            low.b.x.nn <-  rep( B.L[[ k1 ]][ ,l] ,  m[k2]  )
            low.b.y.nn <-  rep( B.L[[ k2 ]][ ,l] , each = m[k1]  )
            low.b.x <- ( low.b.x.nn - z[l] - (x0[l]-z[l])*exp(-theta*t.trans[k1]) )/( sqrt((1-exp(-2*theta*t.trans[k1]) )*(sigma^2)/(2*theta)) )
            up.b.x <- low.b.x + dx/( sqrt((1-exp(-2*theta*t.trans[k1]) )*(sigma^2)/(2*theta)) )
            low.b.y <- ( low.b.y.nn - z[l] - (x0[l]-z[l])*exp(-theta*t.trans[k2]) )/( sqrt((1-exp(-2*theta*t.trans[k2]) )*(sigma^2)/(2*theta)) )
            up.b.y <- low.b.y + dx/( sqrt((1-exp(-2*theta*t.trans[k2]) )*(sigma^2)/(2*theta)) )
            
            pr.mat <- pr.mat*matrix( 
              pbivnorm(  up.b.x , up.b.y , rho =  r   ) - 
                pbivnorm(  low.b.x , up.b.y , rho =  r   ) -
                pbivnorm(  up.b.x , low.b.y , rho =  r   ) +
                pbivnorm(  low.b.x , low.b.y , rho = r   ), 
              m[k1] , 
              m[k2] )  
          }
          #Regularizing if required
          if(reg != 0){
            erratic.too.small <- which(pr.mat < reg)
            erratic.nan <- which(is.nan(pr.mat))
            erratic <- union(erratic.too.small , erratic.nan)
            if( length(erratic) > 0 ){
              pr.mat[erratic] <- reg
            }
          }
          #saving
          P$'2'[[ k1 ]][[ k2 ]] <- pr.mat
        }
      }
      
      #Computing the one-time vectors
      for(k in 1:n.t){
        pr.vec <- rep( p , m[k] )
        for(l in 1:d){
          pr.vec <- pr.vec*(pnorm(B.L[[k]][ , l] + dx , 
                                  mean = z[l] + (x0[l]-z[l])*exp(-theta*t.trans[k])  , 
                                  sd = sqrt((1-exp(-2*theta*t.trans[k]) )*(sigma^2)/(2*theta) ) )-
                              pnorm(B.L[[k]][ , l] , 
                                    mean = z[l] + (x0[l]-z[l])*exp(-theta*t.trans[k])  , 
                                    sd = sqrt((1-exp(-2*theta*t.trans[k]) )*(sigma^2)/(2*theta) ) ) )
        }
        if(consistency.correction){
          #We verify the inequalities the vector should satisfy with respect to two-times matrices
          for(k2 in 1:n.t){
            if(k2 < k){
              pr.vec <- pmax( pr.vec , colSums( P$'2'[[k2]][[k]] )  )
            }
            if(k2 > k){
              pr.vec <- pmax( pr.vec , rowSums( P$'2'[[k]][[k2]] )  )
            }
          }
        }
        P$'1'[[k]] <- pr.vec
      }
    }
    if(trajectory == "Gauss"){
      if(theta.X$indep.components || d == 1 ){
        #Important: we are assuming a very particular form of delivering the mean and covariance functions.
        #The mean is intuitive.
        #The current code only works if the components of the trajectory are independent.
        #       (if not, then one must use pbivnorm or another one depending on the dimension)
        #Computing the two-times matrices
        for(k1 in 1:(n.t-1)){
          for(k2 in (k1+1):n.t){
            pr.mat <- matrix( p^2 , nrow =  m[k1] , ncol = m[k2])
            
            for(l in 1:d){
              r.cov <- Cov[k1,k2,l]
              if(r.cov == 0){
                r <- 0
                #print(r)
              }else{
                r <- sign(r.cov)*exp( log(abs(r.cov)) - 0.5*log(Cov[k1,k1,l]) -0.5*log(Cov[k2,k2,l]))
              }
              #print(r)
              #r <- Cov[k1,k2,l]/( sqrt(Cov[k1,k1,l]*Cov[k2,k2,l]) )
              # print(r)
              # if(abs(r)>1){
              #   r <- sign(r)
              # }
              
              low.b.x.nn <-  rep( B.L[[ k1 ]][ ,l] ,  m[k2]  )
              low.b.y.nn <-  rep( B.L[[ k2 ]][ ,l] , each = m[k1]  )
              low.b.x <- (low.b.x.nn - Mean[k1,l])/sqrt( Cov[k1,k1,l] )
              up.b.x <- (low.b.x.nn + dx - Mean[k1,l])/sqrt( Cov[k1,k1,l] )
              low.b.y <- (low.b.y.nn - Mean[k2,l])/(sqrt( Cov[k2,k2,l] ))
              up.b.y <- (low.b.y.nn + dx - Mean[k2,l])/(sqrt( Cov[k2,k2,l] ))
              
              pr.mat <- pr.mat*matrix( 
                pbivnorm(  up.b.x , up.b.y , rho =  r   ) - 
                  pbivnorm(  low.b.x , up.b.y , rho =  r   ) -
                  pbivnorm(  up.b.x , low.b.y , rho =  r   ) +
                  pbivnorm(  low.b.x , low.b.y , rho = r   ), 
                m[k1] , 
                m[k2] )  
            }
            #Regularizing if required
            if(reg != 0){
              erratic.too.small <- which(pr.mat < reg)
              erratic.nan <- which(is.nan(pr.mat))
              erratic <- union(erratic.too.small , erratic.nan)
              if( length(erratic) > 0 ){
                pr.mat[erratic] <- reg
              }
            }
            #saving
            P$'2'[[ k1 ]][[ k2 ]] <- pr.mat
          }
        }
        #Computing the one-time vectors
        for(k in 1:n.t){
          pr.vec <- rep( p , m[k]  )
          for(l in 1:d){
            pr.vec <- pr.vec*(pnorm(B.L[[k]][ , l] + dx , 
                                    mean = Mean[k,l] , 
                                    sd = sqrt( Cov[k,k,l] ) )-
                                pnorm(B.L[[k]][ , l] , 
                                      mean = Mean[k,l]  , 
                                      sd = sqrt( Cov[k,k,l]  ) ) )
          }
          if(consistency.correction){
            #We verify the inequalities the vector should satisfy with respect to the two-times matrices
            for(k2 in 1:n.t){
              if(k2 < k){
                pr.vec <- pmax( pr.vec , colSums( P$'2'[[k2]][[k]] )  )
              }
              if(k2 > k){
                pr.vec <- pmax( pr.vec , rowSums( P$'2'[[k]][[k2]] )  )
              }
            }
          }
          P$'1'[[k]] <- pr.vec
        }
      }
      else{
        for(k1 in 1:(n.t-1)){
          for(k2 in (k1+1):n.t){
            pr.mat <- matrix( p^2 , nrow =  m[k1] , ncol = m[k2])
            cov.mat <- rbind( cbind( Cov[k1,k1, , ] , Cov[k1,k2, , ]),
                              cbind( t(Cov[k1,k2, , ]) , Cov[k2,k2, , ]))
            mean.vec <- c( Mean[k1, ] , Mean[k2, ] ) 
            for(l1 in 1:m[k1]){
              for(l2 in 1:m[k2]){
                #Al parecer aquí hay un problema.
                pr.mat[l1,l2] <- pr.mat[l1,l2]*sadmvn( lower = c( B.L[[k1]][l1, ] , B.L[[k2]][l2, ] ) ,
                                         upper = c( B.L[[k1]][l1, ] , B.L[[k2]][l2, ] ) + dx , 
                                         mean  = mean.vec , 
                                         varcov = cov.mat )
              }
            }
            #Regularizing if required
            if(reg != 0){
              erratic.too.small <- which(pr.mat < reg)
              erratic.nan <- which(is.nan(pr.mat))
              erratic <- union(erratic.too.small , erratic.nan)
              if( length(erratic) > 0 ){
                pr.mat[erratic] <- reg
              }
            }
            #saving
            P$'2'[[ k1 ]][[ k2 ]] <- pr.mat
          }
        }
        
        if(d == 2){
          #Computing the one-time vectors
          for(k in 1:n.t){
            pr.vec <- rep( p , m[k] )  
            sd.x <- sqrt(Cov[k,k,1,1])
            sd.y <- sqrt(Cov[k,k,2,2])
            r <- Cov[k,k,1,2]/(sd.x*sd.y)
            if(abs(r)>1){
              r <- sign(r)
            }
            m.x <- Mean[k,1]
            m.y <- Mean[k,2]
            
            pr.vec <- pr.vec*(pbivnorm( (B.L[[k]][ , 1] + dx - m.x)/sd.x , 
                                        (B.L[[k]][ , 2] + dx - m.y)/sd.y ,    
                                        rho = r ) - 
                                pbivnorm( (B.L[[k]][ , 1] - m.x)/sd.x , 
                                          (B.L[[k]][ , 2] + dx - m.y)/sd.y ,    
                                          rho = r ) - 
                                pbivnorm( (B.L[[k]][ , 1] + dx - m.x)/sd.x , 
                                          (B.L[[k]][ , 2] - m.y)/sd.y ,    
                                          rho = r ) + 
                                pbivnorm( (B.L[[k]][ , 1] - m.x)/sd.x , 
                                          (B.L[[k]][ , 2] - m.y)/sd.y ,    
                                          rho = r ))
            if(consistency.correction){
              #We verify the inequalities the vector should satisfy with respect to the two-times matrices
              for(k2 in 1:n.t){
                if(k2 < k){
                  pr.vec <- pmax( pr.vec , colSums( P$'2'[[k2]][[k]] )  )
                }
                if(k2 > k){
                  pr.vec <- pmax( pr.vec , rowSums( P$'2'[[k]][[k2]] )  )
                }
              }
            }
            P$'1'[[k]] <- pr.vec
          }
        }
        if(d > 2){
          #Computing the one-time vectors
          for(k in 1:n.t){
            pr.vec <- rep(0,m[k])
            for(l in 1:m[k]){
              pr.vec[l] <- p*sadmvn(  lower = B.L[[k]][l , ] , 
                                          upper = B.L[[k]][l , ] + dx , 
                                          mean = Mean[k , ] , 
                                          varcov = Cov[k,k, , ] )
            }
            
            if(consistency.correction){
              #We verify the inequalities the vector should satisfy with respect to the two-times matrices
              for(k2 in 1:n.t){
                if(k2 < k){
                  pr.vec <- pmax( pr.vec , colSums( P$'2'[[k2]][[k]] )  )
                }
                if(k2 > k){
                  pr.vec <- pmax( pr.vec , rowSums( P$'2'[[k]][[k2]] )  )
                }
              }
            }
            P$'1'[[k]] <- pr.vec
          }
        }
        
      }
      
      
      
       
    }
    
    if(unit.correction){
      print("WARNING: No unitarity correction is yet implemented for the P_Snapshot function at level = 2.")
    }
    return(P)
  }
  
  print("ERROR: P_Snapshot currently available up to level 2.")
  return(NULL)
}



P_Snapshot_Mixture <- function(alpha = 0.5 , 
                               level = 1 , 
                               p , 
                               traj.1 , theta.X.1 ,
                               traj.2 , theta.X.2 , 
                               t0 = 0 , t  , 
                               d , B.L , dx = NULL , Lengths = NULL , 
                               reg = 0 , 
                               unit.correction = F , reg.unit = 0 , 
                               consistency.correction = T){
  P.X.1 <- P_Snapshot(level = level , 
                      p = p , 
                      trajectory = traj.1 , theta.X = theta.X.1 ,
                      t0 = t0 , t = t  , 
                      d = d , B.L = B.L , dx = dx , Lengths = Lengths , 
                      reg = reg , 
                      unit.correction = unit.correction , reg.unit = reg.unit , 
                      consistency.correction = consistency.correction)
  #return(P.X.1)
  P.X.2 <- P_Snapshot(level = level , 
                      p = p , 
                      trajectory = traj.2 , theta.X = theta.X.2 ,
                      t0 = t0 , t = t  , 
                      d = d , B.L = B.L , dx = dx , Lengths = Lengths , 
                      reg = reg , 
                      unit.correction = unit.correction , reg.unit = reg.unit , 
                      consistency.correction = consistency.correction)
  
  
  P <- vector("list" , length = level)
  names(P) <- 1:level
  
  n.t <- length(t)
  
  #Let's do the for loop in order to have an aidea of what should happen
  if(level >= 1){
    P$'1' <- vector("list" , length = n.t)
    for(k in 1:n.t){
      P$'1'[[k]]  <- alpha*P.X.1$'1'[[k]] + (1-alpha)*P.X.2$'1'[[k]]
    }
  }
  if(level >= 2){
    P$'2' <- vector( "list" , length = n.t)
    for(k1 in 1:(n.t-1)){
      P$'2'[[k1]] <- vector( "list" , length = n.t)
      for(k2 in (k1+1):n.t){
        P$'2'[[k1]][[k2]] <- alpha*P.X.1$'2'[[k1]][[k2]] + (1-alpha)*P.X.2$'2'[[k1]][[k2]]
      }
    }
  }
  
  # if(level >= 1){
  #   P$'1' <- vector("list" , length = n.t)
  #   P$'1' <- mapply( function(x,y){alpha*x + (1-alpha)*y} , 
  #                    P.X.1$'1' , P.X.2$'1' , 
  #                    SIMPLIFY = F)
  # }
  # if(level >= 2){
  #   P$'2' <- vector( "list" , length = n.t)
  #   for(k1 in 1:(n.t-1)){
  #     P$'2'[[k1]] <- vector( "list" , length = n.t)
  #     P$'2'[[k1]][(k1+1):n.t] <- mapply( function(x,y){alpha*x + (1-alpha)*y} , 
  #                                        P.X.1$'2'[[k1]][(k1+1):n.t] , P.X.1$'2'[[k1]][(k1+1):n.t] ,
  #                                        SIMPLIFY = F)
  #   }
  # }
  
  if(level <= 2){
    return(P)
  }
  
  print("ERROR: P_Snapshot currently available up to level 2.")
  return(NULL)
}




