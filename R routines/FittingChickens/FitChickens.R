# rm(list = ls(all = T) )
# gc()
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


Q <- readRDS(file.path(getwd() , "Q") )
B.L <- readRDS(file.path(getwd() , "B.L"))
dx <- readRDS(file.path(getwd() , "dx"))
t <- readRDS(file.path(getwd() ,  "t"))

m <- sapply(B.L , nrow)
n.d <- sum(m)
n.t <- length(t)
d <- ncol(B.L[[1]])

#Population size
N <- 93

# library(mnormt)
library(numDeriv)
library(gtools)

functions.path <- file.path( dirname(dirname(getwd()))  , "R functions" )

#Diverse computing functions
source( file.path( functions.path , "P_Snapshot.R"  )  )
source( file.path( functions.path , "ECM_Pairwise_ll.R" ) )
source( file.path( functions.path , "Sim_OU.R" ) )
source( file.path( functions.path , "Sim_Brown.R" ) )
source( file.path( functions.path , "Sim_Population_Count.R" ) )
source( file.path( functions.path , "Mean_Cov_Functions.R") )



####---- Fitting mixture Brownian motion + OU process (non-stationary) ----

print("Fitting (convenient) mixture Brownian motion + OU process (non-stationary). Unknown release time.")

#Here we must also estimate the initial time.

#Parametrization for the OU process in the form (tau,sigma).

#AP negative composite log-likelihood function
AP.Mix.Brown.OU.n.ll <- function(par,q){
  #print(par)
  P.Brown <- P_Snapshot(level = 2 ,
                        p = 1 ,
                        trajectory = "Brown" ,
                        theta.X = exp(par[2])  ,
                        t0 = -exp(par[3]) , t = t ,
                        d = d , B.L = B.L , dx = dx ,
                        reg = 1e-32 ,
                        consistency.correction = T)
  P.OU <- P_Snapshot(level = 2 ,
                     p = 1 ,
                     trajectory = "OU" ,
                     theta.X = c(  sigma = exp(par[2]) ,
                                   theta = 0.5*exp( 2*(par[2]-par[1])) ) ,
                     t0 = -exp(par[3]) , t = t ,
                     d = d , B.L = B.L , dx = dx ,
                     reg = 1e-32 ,
                     consistency.correction = T)
  #Mixture
  level <- 2
  P <- vector("list" , length = level)
  names(P) <- 1:level
  
  n.t <- length(t)
  P$'1' <- vector("list" , length = n.t)
  for(k in 1:n.t){
    P$'1'[[k]]  <- inv.logit(par[4])*P.Brown$'1'[[k]] + (1-inv.logit(par[4]))*P.OU$'1'[[k]]
  }
  P$'2' <- vector( "list" , length = n.t)
  for(k1 in 1:(n.t-1)){
    P$'2'[[k1]] <- vector( "list" , length = n.t)
    for(k2 in (k1+1):n.t){
      P$'2'[[k1]][[k2]] <- inv.logit(par[4])*P.Brown$'2'[[k1]][[k2]] + (1-inv.logit(par[4]))*P.OU$'2'[[k1]][[k2]]
    }
  }
  
  ll <- ECM_Pairwise_ll(q = q ,
                        N = N ,
                        P = P)
  #print(ll)
  return(-ll)
}


#Initial values of the parameters
sigma.0 <- 0.1  #Brownian sd and infinitesimal quadratic variation of OU (at steady state)
t0.0 <- -1 #Initial common release time
alpha.0 <- 0.5  #Proportion of Brownian individuals


log.sigma.0 <- log(sigma.0)
log.minus.t0.0 <- log(-t0.0)
logit.alpha.0 <- logit(alpha.0)


n.par <- 4
#optimization domain limits
lower.par <- c(-8 , -8 , -8 , -7 )
upper.par <- c(6  , 10 , 10 , 7 )

#parscale and ndeps
parscale <- c( 1 , 1 , 1 , 1)
#parscale <- c(1,1)
ndeps <- rep(1e-5 , n.par)

#The optimization is sensitive to the initial value.
#Apparently the parameters which intervenes the most in that aspect is tau.0 (or theta.0)
#We are going to try out with 5 values of theta.0 in different scales, and we will retain the best of those.
thetas.0 <- 10^((-3):3)
best.L <- -Inf
AP.mix.Brown.OU <- NULL
est.par <- NULL
best.theta.0 <- NULL
for(i in 1:length(thetas.0)){
  theta.0 <- thetas.0[i]
  tau.0 <- sigma.0/sqrt(2*theta.0) #Radius of the OU process
  
  log.theta.0 <- log(theta.0)
  log.tau.0 <- log(tau.0)
  
  par.0 <- c( log.tau.0 , log.sigma.0 , log.minus.t0.0 , logit.alpha.0)  
  
  #Optimization
  start <- Sys.time()
  try.optim <- optim(  par.0 ,
                       fn = function(par){AP.Mix.Brown.OU.n.ll(par,Q)} ,
                       # gr = function(par){grad( function(x){AP.Mix.Brown.OU.n.ll(x,Q)} ,
                       #                          x = par)  } ,
                       method = "L-BFGS-B" ,
                       lower = lower.par ,
                       upper = upper.par ,
                       control = list( parscale = parscale ,
                                       ndeps = ndeps ,
                                       maxit = 1e3 ,
                                       factr = 1e5) ,
                       hessian = F )
  print(Sys.time()-start)
  
  if( best.L < -try.optim$value ){
    best.L <- -try.optim$value
    AP.mix.Brown.OU <- try.optim
    est.par <- try.optim$par
    best.theta.0 <- thetas.0[i]
  }
}

###

#est.par

#-AP.Mix.Brown.OU.n.ll( est.par , Q )


G <- grad( function(x){AP.Mix.Brown.OU.n.ll( x , Q )} ,
           x = est.par)
H <- hessian( function(x){AP.Mix.Brown.OU.n.ll( x , Q )} ,
              x = est.par)

print("-----RESULT------")

print("Point estimate:")
est.par
print("Composite log-likelihood:")
-AP.Mix.Brown.OU.n.ll( est.par , Q )
print("Gradient at the estimate:")
G
print("Eigenvalues of the Hessian at the estimate:")
eigen(H)$values
print("---------")
print("Point estimates in original scale:")
tau.est <- exp(est.par[1])
sigma.est <- exp(est.par[2])
theta.est <- 0.5*(sigma.est/tau.est)^2
t0.est <- -exp(est.par[3])
alpha.est <- inv.logit( est.par[4] )
print("Estimate of tau:")
tau.est
print("Estimate of sigma:")
sigma.est
print("Estimate of theta (for information):")
theta.est
print("Estimate of t0:")
t0.est
print("Estimate of alpha:")
alpha.est
print("-----")


#Saving results
saveRDS(AP.mix.Brown.OU , "output/AP.Mix.Brown.OU.rds")


print("Doing parametric bootstrap.")

##Parametric Bootstrap for uncertainty quantification
library(parallel)

n.cores <- 95
cl <- makeCluster(n.cores)

clusterExport(cl , "Q")
clusterExport(cl , "B.L")
clusterExport(cl , "dx")
clusterExport(cl , "t")
clusterExport(cl , "m")
clusterExport(cl , "n.d")
clusterExport(cl , "n.t")
clusterExport(cl , "d")
clusterExport(cl , "N")
clusterExport(cl , "AP.mix.Brown.OU")
clusterExport(cl , "n.par")
clusterExport(cl , "functions.path")


clusterEvalQ(cl , {
  library(numDeriv)
  library(gtools)
  
  #Diverse computing functions
  source( file.path( functions.path , "P_Snapshot.R"  )  )
  source( file.path( functions.path , "ECM_Pairwise_ll.R" ) )
  source( file.path( functions.path , "Sim_OU.R" ) )
  source( file.path( functions.path , "Sim_Brown.R" ) )
  source( file.path( functions.path , "Sim_Population_Count.R" ) )
  source( file.path( functions.path , "Mean_Cov_Functions.R") )
  
  AP.Mix.Brown.OU.n.ll <- function(par,q){
    #print(par)
    P.Brown <- P_Snapshot(level = 2 ,
                          p = 1 ,
                          trajectory = "Brown" ,
                          theta.X = exp(par[2])  ,
                          t0 = -exp(par[3]) , t = t ,
                          d = d , B.L = B.L , dx = dx ,
                          reg = 1e-32 ,
                          consistency.correction = T)
    P.OU <- P_Snapshot(level = 2 ,
                       p = 1 ,
                       trajectory = "OU" ,
                       theta.X = c(  sigma = exp(par[2]) ,
                                     theta = 0.5*exp( 2*(par[2]-par[1])) ) ,
                       t0 = -exp(par[3]) , t = t ,
                       d = d , B.L = B.L , dx = dx ,
                       reg = 1e-32 ,
                       consistency.correction = T)
    level <- 2
    P <- vector("list" , length = level)
    names(P) <- 1:level
    
    n.t <- length(t)
    P$'1' <- vector("list" , length = n.t)
    for(k in 1:n.t){
      P$'1'[[k]]  <- inv.logit(par[4])*P.Brown$'1'[[k]] + (1-inv.logit(par[4]))*P.OU$'1'[[k]]
    }
    P$'2' <- vector( "list" , length = n.t)
    for(k1 in 1:(n.t-1)){
      P$'2'[[k1]] <- vector( "list" , length = n.t)
      for(k2 in (k1+1):n.t){
        P$'2'[[k1]][[k2]] <- inv.logit(par[4])*P.Brown$'2'[[k1]][[k2]] + (1-inv.logit(par[4]))*P.OU$'2'[[k1]][[k2]]
      }
    }
    
    ll <- ECM_Pairwise_ll(q = q ,
                          N = N ,
                          P = P)
    #print(ll)
    return(-ll)
  }
  
  
  #Setting the parameters to the estimates
  tau <- exp(AP.mix.Brown.OU$par[1])
  sigma <- exp(AP.mix.Brown.OU$par[2])
  theta <- 0.5*(sigma/tau)^2
  t0 <- -exp(AP.mix.Brown.OU$par[3])
  alpha <- inv.logit( AP.mix.Brown.OU$par[4] )
  
  
  #Initial values of the parameters
  sigma.0 <- 0.1  #Brownian sd and infinitesimal quadratic variation of OU (at steady state)
  t0.0 <- -1 #Initial common release time
  alpha.0 <- 0.5  #Proportion of Brownian individuals
  
  
  log.sigma.0 <- log(sigma.0)
  log.minus.t0.0 <- log(-t0.0)
  logit.alpha.0 <- logit(alpha.0)
  
  
  n.par <- 4
  #optimization domain limits
  lower.par <- c(-8 , -8 , -8 , -7 )
  upper.par <- c(6  , 10 , 10 , 7 )
  
  #parscale and ndeps
  parscale <- c( 1 , 1 , 1 , 1)
  #parscale <- c(1,1)
  ndeps <- rep(1e-5 , n.par)
  
  
  N.boot.per.core <- 11
}
)

seeds <- 13*(  2*(1:n.cores) - 1) - 2*( 0:(n.cores-1) )


start <- Sys.time()
L <- parLapply(cl , seeds ,
               function(s){
                 #Number of simulations (per core)
                 N.boot <- N.boot.per.core
                 
                 #Returning objects
                 est.par.boot <- matrix( 0 ,  N.boot , n.par)
                 ll.boot <- rep(  0 , N.boot )
                 grad.boot <- matrix(  0 , N.boot , n.par  )
                 hess.boot <- array( 0 , dim = c(N.boot , n.par , n.par)  )
                 
                 
                 #The simulation-fitting loop
                 set.seed(s)
                 for(i in 1:N.boot){
                   #Simulation with the estimated parameters
                   N.Brown <- rbinom( 1 , N , alpha  )
                   N.OU <- N - N.Brown
                   
                   #count of the Brownian
                   q.Brown <- Sim_Population_Count_list(N0 = N.Brown ,
                                                        X.traj = list( sim_func = Sim_Brown ,
                                                                       par = c(sigma) ,
                                                                       sim_type = "times") ,
                                                        x0 = rep(0,d) ,
                                                        t0 = t0 , t = t ,
                                                        d = d , B.L = B.L , dx = dx ,
                                                        det.error = NULL ,
                                                        return.trajectories = F)
                   gc()
                   #count of the OU
                   q.OU <- Sim_Population_Count_list(N0 = N.OU ,
                                                     X.traj = list( sim_func = Sim_OU ,
                                                                    par = c(sigma , theta ) ,
                                                                    sim_type = "times") ,
                                                     x0 = rep(0,d) ,
                                                     t0 = t0 , t = t ,
                                                     d = d , B.L = B.L , dx = dx ,
                                                     det.error = NULL ,
                                                     return.trajectories = F)
                   gc()
                   #Addition of the counts
                   q <- mapply(  '+' , q.Brown , q.OU  )
                   gc()
                   
                   
                   
                   #Estimation
                   #With many par.0 try-outs
                   thetas.0 <- 10^( (-3):3 )
                   best.L <- -Inf
                   for(k in 1:length(thetas.0) ){
                     theta.0 <- thetas.0[k]
                     
                     tau.0 <- sigma.0/sqrt(2*theta.0) #Radius of the OU process
                     
                     log.theta.0 <- log(theta.0)
                     log.tau.0 <- log(tau.0)
                     
                     par.0 <- c( log.tau.0 , log.sigma.0 , log.minus.t0.0 , logit.alpha.0)  
                     
                     #Optimization
                     start <- Sys.time()
                     try.ap.boot.mix <- optim(  par.0 ,
                                                fn = function(par){AP.Mix.Brown.OU.n.ll(par,q)} ,
                                                method = "L-BFGS-B" ,
                                                lower = lower.par ,
                                                upper = upper.par ,
                                                control = list( parscale = parscale ,
                                                                ndeps = ndeps ,
                                                                maxit = 1e3 ,
                                                                factr = 1e5) ,
                                                hessian = F )
                     print(Sys.time()-start)
                     
                     if( best.L < -try.ap.boot.mix$value ){
                       best.L <- -try.ap.boot.mix$value
                       ap.boot.mix <- try.ap.boot.mix
                       #est.par <- try.ap.boot.mix$par
                     }
                   }
                   
                   G <- grad( function(x){AP.Mix.Brown.OU.n.ll( x , q )} ,
                              x = ap.boot.mix$par)
                   H <- hessian( function(x){AP.Mix.Brown.OU.n.ll( x , q )} ,
                                 x = ap.boot.mix$par) 
                   
                   est.par.boot[i , ] <- ap.boot.mix$par
                   ll.boot[i] <- -ap.boot.mix$value
                   grad.boot[i , ] <- G
                   hess.boot[i , , ] <- H
                   gc()
                 }
                 
                 l <- list( #q.2.boot = q.2.boot ,
                   est.par.boot = est.par.boot , ll.boot = ll.boot , grad.boot = grad.boot , hess.boot = hess.boot ,
                   N.boot=N.boot)
                 gc()
                 return(l)
               } )
print(Sys.time()-start)
stopCluster(cl)
gc()



print("Saving results of parametric bootstrap")

saveRDS(L, "L.bootstrap.AP.Mix.Brown.OU.rds")

N.boot.per.core <- L[[1]]$N.boot

#Merging the results into large lists/arrays.
Est.par.boot <- matrix( 0 , N.boot.per.core*n.cores , n.par )
Boot.ll <- rep( 0 , N.boot.per.core*n.cores)
Grad.boot <- matrix( 0 , N.boot.per.core*n.cores , n.par  )
Hess.boot <- array( 0 , dim = c(N.boot.per.core*n.cores , n.par , n.par ) )

#Merging
for(i in 1:n.cores){
  Est.par.boot[(i-1)*N.boot.per.core + 1:N.boot.per.core, ] <- L[[i]]$est.par.boot
  Boot.ll[(i-1)*N.boot.per.core + 1:N.boot.per.core] <- L[[i]]$ll.boot[i]
  Grad.boot[(i-1)*N.boot.per.core + 1:N.boot.per.core, ] <- L[[i]]$grad.boot
  Hess.boot[(i-1)*N.boot.per.core + 1:N.boot.per.core, , ] <- L[[i]]$hess.boot
}

saveRDS(Est.par.boot , "Est.par.boot AP.Mix.Brown.OU.rds")
saveRDS(Boot.ll , "Boot.ll AP.Mix.Brown.OU.rds")
saveRDS(Grad.boot , "Grad.boot AP.Mix.Brown.OU.rds")
saveRDS(Hess.boot , "Hess.boot AP.Mix.Brown.OU.rds")