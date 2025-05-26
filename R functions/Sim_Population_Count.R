#Coded by Ricardo CarrizoV
#Version at 18-02-2025
#library(mnormt)

#Update 18-02-2025: Added rm(X), rm(BEL) and gc() to avoid memory issues.
#Update 30-01-2025: Added the option of return the trajectories (in the times setting)

Sim_Population_Count_list <- function(N0 , 
                                      X.traj = list() , x0 = NULL , 
                                      t0 = 0 , t , 
                                      d , B.L , dx = NULL , Lengths = NULL , 
                                      det.error = NULL , 
                                      death = NULL, 
                                      output = "list-of-vectors" , 
                                      return.trajectories = F){
  
  
  n.t <- length(t)
  m <- rep(0 , n.t)
  for(k in 1:n.t){
    m[k] <- nrow( B.L[[k]] )
  }
  #n.d <- sum(m)
  
  #Simulating the trajectories of the elders
  X <- X.traj$sim_func(  N = N0 , 
                         d = d , 
                         t0 = t0 , 
                         par = X.traj$par , 
                         x0 = x0 , 
                         sim_type = X.traj$sim_type , 
                         sim_sett = t  )
  
  #Abundance list
  Q <- vector("list" , length = n.t)
  #Counting the presence
  for(k in 1:n.t){
    #s.locations <- st.locations[ st.locations[,1] == t[k] ,                                2:(1+2*d)]
    BEL <- matrix( TRUE , nrow = m[k] , ncol = dim(X)[1]  )   
    for(l in 1:d){
      BEL <- BEL&outer( B.L[[k]][ ,l] , X[ ,k,l] , '<'  )
      BEL <- BEL&outer( B.L[[k]][ ,l] + dx , X[ ,k,l] , '>=')
    }
    Q[[k]] <- rowSums(BEL)
  }
  rm(BEL)
  if(!is.null(det.error)){
    if(det.error$method == "Snapshot"){
      p <- det.error$par
      for(k in 1:n.t){
        Q[[k]] <- rbinom( m[k] ,  Q[[k]] , p )
      }
      return(Q)
    }
    print("WARNING: Only 'Snapshot' detection error method is currently available.")
  }
  
  if(return.trajectories){
    l <- list(Q = Q , X = X  )
    return(l)
  }
  rm(X)
  gc()
  return(Q)
}

