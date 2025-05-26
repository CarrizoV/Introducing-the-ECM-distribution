# rm(list = ls(all = T))
# gc()
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#####----- Simulation Study for OU (stationary) Process  -----#####

## log-scale for theta and tau, center unknown.
## Toy setting with d = 2
## parallel code
## Unknown N (ECM-Poisson count)

library(parallel)

detectCores()
n.cores <- 95
cl <- makeCluster(n.cores)

print("OU steady-state simulation study. Unknown center. Unknown N. ")

#IMPORTANT: The case MGLE lambda=10^2 should be done separately with other configuration. Read the comments.


#Population size loop
for(theta in c(0.001)){
  for(lambda in 10^c(3:4)){ 
    #Don't do lambda=10^2 with MGLE method unless the correct configuration has been coded.
    #For lambda=10^2 with MCLE, all the objects and functions regarding MGLE must be commented, or, failing that, the MGLE results should be ignored. 

    clusterExport(cl , "lambda")
    print("lambda = ")
    print(lambda)
    clusterExport(cl , "theta")
    print("theta = ")
    print(theta)
    print("IMPORTANT: THE PARAMETRIZATION IS ON (tau,sigma)")
    print("method:")
    print("AP")
    print("MGLE")
    
    clusterEvalQ(cl , {
      library(mnormt)
      library(numDeriv)
      library(gtools)
      
      #Diverse computing functions
      source( file.path( functions.path , "P_Snapshot.R"  )  )
      source( file.path( functions.path , "dBivPois.R" ) )
      source( file.path( functions.path , "ECM_Poisson_Pairwise_ll.R" ) )
      source( file.path( functions.path , "ECM_Poisson_Gauss_ll.R" ) )
      source( file.path( functions.path , "Sim_OU.R" ) )
      source( file.path( functions.path , "Sim_Population_Count.R" ) )
      source( file.path( functions.path , "Mean_Cov_Functions.R") )
      source( "ToySetting2D.R" )
      
      #Setting the parameters
      
      #log population rate
      log.lambda <- log(lambda)
      
      #center
      z <- c(-0.2 , 0.1)
      
      #Individual movement parameters
      tau <- 0.4
      log.tau <- log(tau)
      log.theta <- log(theta)
      
      #sigma parameter
      sigma <- tau*sqrt(2*theta) 
      log.sigma <- log(sigma)
      
      
      
      par <- c(log.tau , log.sigma , z , log.lambda )
      n.par <- length(par)
      
      
      #Number of simulation-estimations per core
      N.sim.per.core <- 11 #Set up to 22 for MGLE lambda = 10^2 case
      
      #optimization domain limits
      lower.par <- c(-8 , -8 , -1 , -1 , log(lambda/10) )  #Use log(lambda/5) for MGLE lambda=10^2 case.
      upper.par <- c(6  , 10  , 1 ,  1 , log(lambda*10) )  #Use log(lambda*6) for MGLE lambda=10^2 case.
      
      #parscale and ndeps
      parscale <- c(1 , 1 , 1 , 1 , 1)
      ndeps <- rep(1e-5 , n.par)
      
    }
    )
    
    seeds <- 13*(  2*(1:n.cores) - 1) - 2*( 0:(n.cores-1) )
    
    start <- Sys.time()
    L <- parLapply(cl , seeds , 
                   function(s){
                     #Number of simulations (per core)
                     N.sim <- N.sim.per.core
                     
                     #List with the simulated values
                     Q.sim <- vector( "list" , length = N.sim)
                     
                     #Each method has its associated objects:
                     # -List of optim outputs
                     # -Matrices with the obtained point estimates
                     # -Vectors with the final obtained (pseudo) log-likelihood value
                     # -Matrices with the gradients at the optimum
                     # -Arrays with the hessians at the optimum
                     
                     #MGLE objects
                     MGLE <- vector( "list" , length = N.sim)
                     MGLE.par <- matrix(0 , N.sim , n.par)
                     MGLE.ll <- rep( 0 , N.sim)
                     MGLE.Grad <- matrix( 0 , N.sim , n.par)
                     MGLE.Hess <- array( 0 , c(N.sim , n.par , n.par))
                     
                     #AP objects
                     AP <- vector( "list" , length = N.sim)
                     AP.par <- matrix(0 , N.sim , n.par)
                     AP.ll <- rep( 0 , N.sim)
                     AP.Grad <- matrix( 0 , N.sim , n.par)
                     AP.Hess <- array( 0 , c(N.sim , n.par , n.par))

                     
                     #Negative log-likelihood functions
                     MGLE.n.ll <- function(par,q){
                       #print(par)
                     P <- P_Snapshot(level = 2 ,
                                     p = 1 ,
                                     trajectory = "Gauss" ,
                                     theta.X = list( mean = f.mean.OU.stat( sigma = exp(par[2]) ,
                                                                            theta = 0.5*exp( 2*(par[2]-par[1]) ),
                                                                            z = par[2 + 1:d],d) ,
                                                     cov = f.cov.OU.stat( sigma = exp(par[2]) ,
                                                                          theta = 0.5*exp( 2*(par[2]-par[1]) ),
                                                                          z = par[2 + 1:d],d) ,
                                                     indep.components = T)  ,
                                     t0 = t0 , t = t ,
                                     d = d , B.L = B.L , dx = dx ,
                                     reg = 1e-32 ,
                                     consistency.correction = T)
                       ll <- ECM_Poisson_Gauss_ll( q = q ,
                                                   lambda = exp(par[5]) ,
                                                   P = P ,
                                                   eps = 1e-6)
                       #print(ll)
                       return(-ll)
                     }
                     AP.n.ll <- function(par,q){
                       #print(par)
                       P <- P_Snapshot(level = 2 ,
                                       p = 1 ,
                                       trajectory = "Gauss" ,
                                       theta.X = list( mean = f.mean.OU.stat( sigma = exp(par[2]) ,
                                                                              theta = 0.5*exp( 2*(par[2]-par[1]) ),
                                                                              z = par[2 + 1:d],d) ,
                                                       cov = f.cov.OU.stat( sigma = exp(par[2]) ,
                                                                            theta = 0.5*exp( 2*(par[2]-par[1]) ),
                                                                            z = par[2 + 1:d],d) ,
                                                       indep.components = T)  ,
                                       t0 = t0 , t = t ,
                                       d = d , B.L = B.L , dx = dx ,
                                       reg = 1e-32 ,
                                       consistency.correction = T)
                       ll <- ECM_Poisson_Pairwise_ll(q = q ,
                                                     lambda = exp(par[5]) ,
                                                     P = P)
                       #print(ll)
                       return(-ll)
                     }

                     #The simulation loop
                     set.seed(s)
                     for(j in 1:N.sim){
                       print(paste("---- Simulation" , j , "----"))
                       ##Simulation
                       N <- rpois( 1 , lambda  )

                       #initial positions.
                       x0 <- t(matrix( rnorm( d*N , mean = z , sd = tau  ) , 
                                       nrow = d , 
                                       ncol = N ))
                       
                       #watch out with the case N = 0...
                       Q.sim[[j]] <- Sim_Population_Count_list(N0 = N , 
                                                               X.traj = list( sim_func = Sim_OU , 
                                                                              par = c(sigma , theta , z ) , 
                                                                              sim_type = "times") , 
                                                               x0 = x0 , 
                                                               t0 = t0 , t = t , 
                                                               d = d , B.L = B.L , dx = dx , 
                                                               det.error = NULL , 
                                                               return.trajectories = F)
                       
                       gc()
                       
                       #Proposed initial values for the optimization
                       lambda.0 <- lambda/2
                       log.lambda.0 <- log(lambda.0)
                       
                       z.0 <- c(0,0)
                       
                       tau.0 <- tau/2 
                       log.tau.0 <- log(tau.0)
                       
                       #Try.out with different thetas as initial value
                       thetas.0 <- 0.5*c( theta/10 , theta , 10*theta , 100*theta )
                       
                       #We retain the best obtained from such initial values
                       ll.best.mgle <- -Inf
                       ll.best.ap <- -Inf

                       
                       for(k in 1:length(thetas.0)){
                         
                         theta.0 <- thetas.0[k]
                         log.theta.0 <- log(theta.0)
                         sigma.0 <- tau.0*sqrt(2*theta.0)
                         log.sigma.0 <- log(sigma.0)
                         
                         par.0 <- c(log.tau.0 , log.sigma.0 ,  z.0 , log.lambda.0)
                         
                         #MGLE
                         #start <- Sys.time()
                         mgle <- optim( par.0 ,
                                        fn = function(x){MGLE.n.ll(x,Q.sim[[j]])},
                                        method = "L-BFGS-B" ,
                                        lower = lower.par ,
                                        upper = upper.par ,
                                        control = list(parscale = parscale ,
                                                       ndeps = ndeps) ,
                                        hessian = F)

                         if(ll.best.mgle < -mgle$value){
                           MGLE[[j]] <- mgle
                           ll.best.mgle <- -mgle$value
                         }
                         #print("MGLE finished. Required time:")
                         #print(Sys.time()-start)
                         
                         #All Pairs composite likelihood
                         #start <- Sys.time()
                         ap <- optim( par.0 ,
                                      function(x){AP.n.ll(x,Q.sim[[j]])} ,
                                      method = "L-BFGS-B" ,
                                      lower = lower.par ,
                                      upper = upper.par ,
                                      control = list(parscale = parscale ,
                                                     ndeps = ndeps) ,
                                      hessian = F)
                         if( ll.best.ap < -ap$value ){
                           AP[[j]] <- ap
                           ll.best.ap <- -ap$value
                         }
                       }
                       
                       # #print("AP finished. Required time:")
                       # #print(Sys.time()-start)
                       # 
                       
                       AP.par[j , ] <- AP[[j]]$par
                       AP.ll[j] <- -AP[[j]]$value
                       AP.Grad[j , ] <- grad(function(x){AP.n.ll(x,Q.sim[[j]])} ,
                                             AP[[j]]$par  )
                       AP.Hess[j , ,  ] <- hessian(function(x){AP.n.ll(x,Q.sim[[j]])} ,
                                                   AP[[j]]$par  )
                       
                       MGLE.par[j , ] <- MGLE[[j]]$par
                       MGLE.ll[j] <- -MGLE[[j]]$value
                       MGLE.Grad[j , ] <- grad(func = function(x){MGLE.n.ll(x,Q.sim[[j]])} ,
                                               MGLE[[j]]$par  )
                       MGLE.Hess[j , ,  ] <- hessian(func = function(x){MGLE.n.ll(x,Q.sim[[j]])} ,
                                                     MGLE[[j]]$par)
                     }
                     
                     
                     l <- list(Q.sim = Q.sim, 
                               MGLE = MGLE , MGLE.par = MGLE.par , MGLE.ll = MGLE.ll , MGLE.Grad = MGLE.Grad , MGLE.Hess = MGLE.Hess ,
                               AP = AP , AP.par = AP.par , AP.ll = AP.ll , AP.Grad = AP.Grad , AP.Hess = AP.Hess ,
                               lambda = lambda , par = par , N.sim = N.sim )
                     gc()
                     return(l)
                   } )
    #stopCluster(cl)
    gc()
    print(Sys.time()-start)
    
    #Retrieving some key values
    lambda <- L[[1]]$lambda
    N.sim.per.core <- L[[1]]$N.sim
    par <- L[[1]]$par
    n.par <- length(par)
    
    tau <- exp(par[1])
    sigma <- exp(par[2])
    theta <- 0.5*exp( 2*(par[2]-par[1]) )
    z <- par[3:4]

    
    #Merging the results into large lists/arrays.
    Q.sim <- vector( "list" , length = N.sim.per.core*n.cores)

    
    #MGLE objects
    MGLE <- vector( "list" , length = N.sim.per.core*n.cores  )
    MGLE.par <- matrix( 0 , N.sim.per.core*n.cores , n.par )
    MGLE.ll <- rep( 0 , N.sim.per.core*n.cores)
    MGLE.Grad <- matrix( 0 , N.sim.per.core*n.cores , n.par)
    MGLE.Hess <- array( 0 , c(N.sim.per.core*n.cores , n.par , n.par) )
    
    #AP objects
    AP <- vector( "list" , length = N.sim.per.core*n.cores )
    AP.par <- matrix( 0 , N.sim.per.core*n.cores , n.par )
    AP.ll <- rep( 0 , N.sim.per.core*n.cores )
    AP.Grad <- matrix( 0 , N.sim.per.core*n.cores , n.par)
    AP.Hess <- array( 0 , c(N.sim.per.core*n.cores , n.par , n.par) )

    
    #Merging
    for(i in 1:n.cores){
      for(j in 1:N.sim.per.core){
        Q.sim[[ (i-1)*N.sim.per.core + j  ]] <- L[[i]]$Q.sim[[j]]
        MGLE[[ (i-1)*N.sim.per.core + j  ]] <- L[[i]]$MGLE[[j]]
        AP[[ (i-1)*N.sim.per.core + j  ]] <- L[[i]]$AP[[j]]
      }
      
      MGLE.par[(i-1)*N.sim.per.core + 1:N.sim.per.core ,  ] <- L[[i]]$MGLE.par
      MGLE.ll[(i-1)*N.sim.per.core + 1:N.sim.per.core ] <- L[[i]]$MGLE.ll
      MGLE.Grad[(i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$MGLE.Grad
      MGLE.Hess[(i-1)*N.sim.per.core + 1:N.sim.per.core , , ] <- L[[i]]$MGLE.Hess
      
      AP.par[(i-1)*N.sim.per.core + 1:N.sim.per.core ,  ] <- L[[i]]$AP.par
      AP.ll[(i-1)*N.sim.per.core + 1:N.sim.per.core ] <- L[[i]]$AP.ll
      AP.Grad[(i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$AP.Grad
      AP.Hess[(i-1)*N.sim.per.core + 1:N.sim.per.core , , ] <- L[[i]]$AP.Hess
    }
    
    #Simulations
    saveRDS( Q.sim , file = paste("output/Q.sim", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                  ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "")      )
    
    #MGLE results
    saveRDS( MGLE , file = paste("output/MGLE", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                 ", lambda=" , lambda, ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( MGLE.par , file = paste("output/MGLE.par", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                     ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( MGLE.ll , file = paste("output/MGLE.ll", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                    ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( MGLE.Grad , file = paste("output/MGLE.Grad", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                      ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( MGLE.Hess , file = paste("output/MGLE.Hess", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                      ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )

    #AP results
    saveRDS( AP , file = paste("output/AP", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                               ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( AP.par , file = paste("output/AP.par", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                   ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( AP.ll , file = paste("output/AP.ll", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                  ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( AP.Grad , file = paste("output/AP.Grad", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                    ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    saveRDS( AP.Hess , file = paste("output/AP.Hess", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                    ", lambda=" , lambda , ", N.sim=" , N.sim.per.core*n.cores, ".rds",  sep = "") )
    }
}
