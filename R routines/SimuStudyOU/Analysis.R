####---- Script to analyse the results of simulation studies ----
rm(list = ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(gtools)
library(data.table)
library(matrixStats)

#Setting the trajectory model scenario to analyse
traj.model <- "OU stat (abundance, unknown center)"

#Setting theoretical parameters (except N and lambda)

#Proposed theoretical parameters
tau <- 0.4
theta <- 0.001 #(for information)

sigma <- tau*sqrt(2*theta)

z <- c(-0.2,0.1)

par <- c( log(tau) , log(sigma) , z   )
n.par <- length(par)


#Proposed (gross) number of simulations
N.sim <- 1045


#Proposed Ns and lambdas
N.s <- 10^(2:4)
lambda.s <- 10^(2:4)


N.sim.wanted <- 1000


EST.PAR.ECM.AP <- array( 0 , dim = c( N.sim.wanted , n.par , length(N.s) ) )   
EST.PAR.ECM.MGLE <- array( 0 , dim = c( N.sim.wanted , n.par , length(N.s) ) )  

#Loading results for ECM
for(method in c("AP","MGLE")){
  for(i.N in 1:length(N.s)){
    print("--------Loading results for ECM estimation-------")
    print(paste("--------N: " , N.s[i.N] ) )
    print(paste("--------method: " , method ) )
    
    #Proposed parameters and size population
    N <- N.s[i.N]
    
    #Change file directory to the one where the results are
    cur.dir <- getwd()
    setwd( file.path(  getwd() , "output" )  )
    
    
    #Loading the simulations
    Q.sim <- readRDS(paste( "Q.sim", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                            ", N=" , N ,   ", N.sim=" , N.sim , ".rds",  sep = ""))
    
    #Loading the estimation key results
    optim.result <- readRDS(paste( method ,   ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                   ", N=" , N , ", N.sim=" , N.sim, ".rds",  sep = ""))
    est.par <- readRDS(paste( method , ".par" ,    ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                              ", N=" , N , ", N.sim=" , N.sim, ".rds",  sep = ""))
    ll <- readRDS(paste( method , ".ll" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                         ", N=" , N , ", N.sim=" , N.sim, ".rds",  sep = ""))
    grads <- readRDS(paste( method , ".grad" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                            ", N=" , N , ", N.sim=" , N.sim, ".rds",  sep = ""))
    hessians <- readRDS(paste( method , ".Hess" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                               ", N=" , N , ", N.sim=" , N.sim, ".rds",  sep = ""))
    ##Re-changing the working directory
    setwd(cur.dir)
    
    
    #We filter the "erratic" estimations, that is, the cases where the Hessian is not positive definite.
    n.erratic <- 0
    ind.erratic <- numeric(0)
    for(i in 1:N.sim){
      #Important: as positive-definiteness criterion, we require all the eigenvalues to be bigger than 1e-12. 
      if( !all( eigen(hessians[i, , ])$values > 1e-12 ) ){
        print(paste("Warning: not positive variances at simulation " , i))
        print("Eigenvalues of the Hessian:")
        print(eigen(hessians[i, , ])$values)
        n.erratic <- n.erratic + 1
        ind.erratic <- c(ind.erratic , i)
      }
    }
    print(paste("Number of erratic estimations: " , n.erratic))
    
    #Only non-erratic estimations are considered.
    if( n.erratic > 0 ){
      Q.sim <- Q.sim[-ind.erratic]
      est.par <- est.par[ -ind.erratic , ]
      ll <- ll[ -ind.erratic ]
      grads <- grads[-ind.erratic , ]
      hessians <- hessians[ -ind.erratic , , ]
    }
    
    
    N.sim.considered <- N.sim.wanted #Number of simulations we want to consider
    if(N.sim - n.erratic < N.sim.considered){
      print("WARNING: You are trying to analize more estimations than available and valid.")
      print(paste("Number of desired estimations: " , N.sim.considered))
      print(paste("Number of available estimations: " , N.sim))
      print(paste("Number of erratic estimations: " , n.erratic))
    }
    #Real final number of analysed simulations
    N.sim.considered <- min(N.sim.considered , N.sim - n.erratic) 
    print(paste("Final number of effectively analysed simulations: " , N.sim.considered))
    
    #Retaining only the considered ones.
    Q.sim <- Q.sim[1:N.sim.considered]
    est.par <- est.par[ 1:N.sim.considered , ]
    ll <- ll[ 1:N.sim.considered ]
    grads <- grads[ 1:N.sim.considered , ]
    hessians <- hessians[ 1:N.sim.considered , , ]
    
    
    if( method == "AP" ){
      EST.PAR.ECM.AP[  , , i.N ] <- est.par
    }
    if( method == "MGLE" ){
      EST.PAR.ECM.MGLE[  , , i.N ] <- est.par
    }
  }
}


EST.PAR.ECM.POISSON.AP <- array( 0 , dim = c( N.sim.wanted , n.par+1 , length(lambda.s) ) )  
EST.PAR.ECM.POISSON.MGLE <- array( 0 , dim = c( N.sim.wanted , n.par+1 , length(lambda.s) ) )  

#Loading results for ECM-Poisson
for(method in c( "AP" , "MGLE")){
  for(i.lambda in 1:length(lambda.s)){ #1:length(lambda.s)
    print("--------Loading results for ECM-Poisson estimation-------")
    print(paste("--------lambda: " , lambda.s[i.lambda] ) )
    print(paste("--------method: " , method ) )
    
    #Proposed size rate population
    lambda <- lambda.s[i.lambda]
    
    
    if(method == "MGLE" & i.lambda==1){
      N.sim <- 2090
      N.sim.wanted <- 800
    }else{
      N.sim <- 1045
      N.sim.wanted <- 1000
    }
    
    
    #Change file directory to the one where the results are
    cur.dir <- getwd()
    setwd( file.path(  getwd() , "output" )  )
    
    #Loading the simulations
    Q.sim <- readRDS(paste( "Q.sim", ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                            ", lambda=" , lambda ,   ", N.sim=" , N.sim , ".rds",  sep = ""))
    
    #Loading the estimation key results
    optim.result <- readRDS(paste( method ,   ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                                   ", lambda=" , lambda , ", N.sim=" , N.sim, ".rds",  sep = ""))
    est.par <- readRDS(paste( method , ".par" ,    ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                              ", lambda=" , lambda , ", N.sim=" , N.sim, ".rds",  sep = ""))
    ll <- readRDS(paste( method , ".ll" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                         ", lambda=" , lambda , ", N.sim=" , N.sim, ".rds",  sep = ""))
    grads <- readRDS(paste( method , ".grad" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                            ", lambda=" , lambda , ", N.sim=" , N.sim, ".rds",  sep = ""))
    hessians <- readRDS(paste( method , ".Hess" , ", par=(log(",tau , "), log(" , sigma , "), (" , z[1] , " , " , z[2] , "))" ,
                               ", lambda=" , lambda , ", N.sim=" , N.sim, ".rds",  sep = ""))
    ##Re-changing the working directory
    setwd(cur.dir)
    
    
    #We filter the "erratic" estimations, that is, the cases where the Hessian is not positive definite.
    n.erratic <- 0
    ind.erratic <- numeric(0)
    for(i in 1:N.sim){
      #Important: as positive-definiteness criterion, we require all the eigenvalues to be bigger than 1e-12. 
      if( !all( eigen(hessians[i, , ])$values > 1e-12 ) ){
        print(paste("Warning: not positive variances at simulation " , i))
        print("Eigenvalues of the Hessian:")
        print(eigen(hessians[i, , ])$values)
        n.erratic <- n.erratic + 1
        ind.erratic <- c(ind.erratic , i)
      }
    }
    print(paste("Number of erratic estimations: " , n.erratic))
    
    
    N.sim.considered <- N.sim.wanted #Number of simulations we want to consider
    
    if( n.erratic > 0 ){
      Q.sim <- Q.sim[-ind.erratic]
      est.par <- est.par[ -ind.erratic , ]
      ll <- ll[ -ind.erratic ]
      grads <- grads[-ind.erratic , ]
      hessians <- hessians[ -ind.erratic , , ]
    }
    
    if(N.sim - n.erratic < N.sim.considered){
      print("WARNING: You are trying to analize more estimations than available and valid.")
      print(paste("Number of desired estimations: " , N.sim.considered))
      print(paste("Number of available estimations: " , N.sim))
      print(paste("Number of erratic estimations: " , n.erratic))
    }
    #Real final number of analysed simulations
    N.sim.considered <- min(N.sim.considered , N.sim - n.erratic) 
    print(paste("Final number of effectively analysed simulations: " , N.sim.considered))
    
    
    #Retaining only the considered ones.
    Q.sim <- Q.sim[1:N.sim.considered]
    est.par <- est.par[ 1:N.sim.considered , ]
    ll <- ll[ 1:N.sim.considered ]
    grads <- grads[ 1:N.sim.considered , ]
    hessians <- hessians[ 1:N.sim.considered , , ]
    
    
    if( method == "AP" ){
      EST.PAR.ECM.POISSON.AP[  , , i.lambda ] <- est.par
    }
    if( method == "MGLE" ){
      EST.PAR.ECM.POISSON.MGLE[ 1:N.sim.considered ,  , i.lambda ] <- est.par
    }
  }
}
#The optimization is really poor for the case ECM-Poisson MGLE, lambda = 10^2.


####---- Boxplots generator ----  


#Positions in the boxplots
positions.boxplots <- c(1:4 , 6:9 , 11:14)

#Tau parameter
tau
log(tau)
par.index <- 1

png("BoxplotTau.png", width = 3500, height = 1500, res = 300)

#Parameters for boxplots
par(mar = c(4.5, 2.5, 2, 1) + 0.1)


values.for.box <- list(
  EST.PAR.ECM.MGLE[ ,par.index,1] , EST.PAR.ECM.AP[ ,par.index,1] , EST.PAR.ECM.POISSON.MGLE[1:800 ,par.index,1] , EST.PAR.ECM.POISSON.AP[ ,par.index,1],
  EST.PAR.ECM.MGLE[ ,par.index,2] , EST.PAR.ECM.AP[ ,par.index,2] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,2] , EST.PAR.ECM.POISSON.AP[ ,par.index,2],
  EST.PAR.ECM.MGLE[ ,par.index,3] , EST.PAR.ECM.AP[ ,par.index,3] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,3] , EST.PAR.ECM.POISSON.AP[ ,par.index,3]
)

#Boxplot without names
boxplot(values.for.box, at = positions.boxplots,
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 3),
        las = 1,
        #main = expression("Estimation of " * log(tau) == log(0.4) ~ "\u2248" ~ -0.916) ,
        main = expression(log(tau) )
)
abline( h = log(tau) , lty = 6 , col = "gold" , 
        lwd = 1.8 )

labels.box <- expression(
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^2)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^2))), atop("MGLE*" ,  displaystyle(atop("ECM-P" , lambda == 10^2)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^2))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^3)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^3))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^3)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^3))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^4)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^4))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^4)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^4)))
)

axis(side = 1, at = positions.boxplots, 
     labels = labels.box, tick = FALSE, line = 2)
dev.off()




#Sigma parameter
sigma
log(sigma)
par.index <- 2

png("BoxplotSigma.png", width = 3500, height = 1500, res = 300)

#Parameters for boxplots
par(mar = c(4.5, 2.5, 2, 1) + 0.1)


values.for.box <- list(
  EST.PAR.ECM.MGLE[ ,par.index,1] , EST.PAR.ECM.AP[ ,par.index,1] , EST.PAR.ECM.POISSON.MGLE[ 1:800 ,par.index,1] , EST.PAR.ECM.POISSON.AP[ ,par.index,1],
  EST.PAR.ECM.MGLE[ ,par.index,2] , EST.PAR.ECM.AP[ ,par.index,2] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,2] , EST.PAR.ECM.POISSON.AP[ ,par.index,2],
  EST.PAR.ECM.MGLE[ ,par.index,3] , EST.PAR.ECM.AP[ ,par.index,3] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,3] , EST.PAR.ECM.POISSON.AP[ ,par.index,3]
)

#Boxplot without names
boxplot(values.for.box, at = positions.boxplots,
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 3),
        las = 1,
        #main = expression("Estimation of " * log(sigma) ~ "\u2248" ~ log(0.0179) ~ "\u2248" ~ -4.024),
        main = expression(log(sigma))
)
abline( h = log(sigma) , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^2)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^2))), atop("MGLE*" ,  displaystyle(atop("ECM-P" , lambda == 10^2)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^2))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^3)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^3))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^3)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^3))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^4)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^4))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^4)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^4)))
)


axis(side = 1, at = positions.boxplots, 
     labels = labels.box, tick = FALSE, line = 2)
dev.off()



#Z1 parameter
z[1]
par.index <- 3

png("BoxplotZ1.png", width = 3500, height = 1500, res = 300)

#Parameters for boxplots
par(mar = c(4.5, 2.5, 2, 1) + 0.1)


values.for.box <- list(
  EST.PAR.ECM.MGLE[ ,par.index,1] , EST.PAR.ECM.AP[ ,par.index,1] , EST.PAR.ECM.POISSON.MGLE[1:800 ,par.index,1] , EST.PAR.ECM.POISSON.AP[ ,par.index,1],
  EST.PAR.ECM.MGLE[ ,par.index,2] , EST.PAR.ECM.AP[ ,par.index,2] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,2] , EST.PAR.ECM.POISSON.AP[ ,par.index,2],
  EST.PAR.ECM.MGLE[ ,par.index,3] , EST.PAR.ECM.AP[ ,par.index,3] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,3] , EST.PAR.ECM.POISSON.AP[ ,par.index,3]
)

#Boxplot without names
boxplot(values.for.box, at = positions.boxplots,
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 3),
        las = 1,
        #main = expression("Estimation of " * z[1] == -0.2  ),
        main = expression(z[1])
)
abline( h = z[1] , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^2)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^2))), atop("MGLE*" ,  displaystyle(atop("ECM-P" , lambda == 10^2)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^2))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^3)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^3))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^3)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^3))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^4)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^4))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^4)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^4)))
)


axis(side = 1, at = positions.boxplots, 
     labels = labels.box, tick = FALSE, line = 2)
dev.off()


#Z2 parameter
z[2]
par.index <- 4

png("BoxplotZ2.png", width = 3500, height = 1500, res = 300)

#Parameters for boxplots
par(mar = c(4.5, 2.5, 2, 1) + 0.1)


values.for.box <- list(
  EST.PAR.ECM.MGLE[ ,par.index,1] , EST.PAR.ECM.AP[ ,par.index,1] , EST.PAR.ECM.POISSON.MGLE[1:800 ,par.index,1] , EST.PAR.ECM.POISSON.AP[ ,par.index,1],
  EST.PAR.ECM.MGLE[ ,par.index,2] , EST.PAR.ECM.AP[ ,par.index,2] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,2] , EST.PAR.ECM.POISSON.AP[ ,par.index,2],
  EST.PAR.ECM.MGLE[ ,par.index,3] , EST.PAR.ECM.AP[ ,par.index,3] , EST.PAR.ECM.POISSON.MGLE[ ,par.index,3] , EST.PAR.ECM.POISSON.AP[ ,par.index,3]
)

#Boxplot without names
boxplot(values.for.box, at = positions.boxplots,
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 3),
        las = 1,
        #main = expression("Estimation of " * z[2] == 0.1  ),
        main = expression(z[2])
)
abline( h = z[2] , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^2)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^2))), atop("MGLE*" ,  displaystyle(atop("ECM-P" , lambda == 10^2)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^2))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^3)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^3))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^3)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^3))),
  atop("MGLE" , displaystyle(atop("ECM" , N == 10^4)) ), atop("MCLE", displaystyle(atop("ECM" , N == 10^4))), atop("MGLE" ,  displaystyle(atop("ECM-P" , lambda == 10^4)) ), atop("MCLE",  displaystyle(atop("ECM-P" , lambda == 10^4)))
)


axis(side = 1, at = positions.boxplots, 
     labels = labels.box, tick = FALSE, line = 2)
dev.off()



#lambda parameter
par.index <- 5
log(lambda.s)


png("BoxplotLambda.png" , height = 1200 , width = 2000 , res = 300)

par(mfrow = c(1,3))
#Parameters for boxplots
par(mar = c(2, 2, 2, 1) + 0.1)

#lambda = 10^2
values.for.box <- list(
  EST.PAR.ECM.POISSON.MGLE[ 1:800 , par.index,1] , EST.PAR.ECM.POISSON.AP[ , par.index,1]
)

#Boxplot without names
boxplot(values.for.box, 
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 1),
        las = 1,
        ylim = c(4,9.5),
        #main = expression("Estimation of " * log(lambda) ~ "\u2248" ~ log(10^2) ~ "\u2248" ~ 4.605)
        main = expression(log(lambda) == log(10^2)),
        cex.main = 1.5
)
abline( h = log(lambda.s[1]) , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  "MGLE*", "MCLE"
)
axis(side = 1, at = c(1,2) , cex.axis = 1.5,
     labels = labels.box, tick = FALSE)
#dev.off()

#lambda = 10^3
values.for.box <- list(
  EST.PAR.ECM.POISSON.MGLE[ , par.index,2] , EST.PAR.ECM.POISSON.AP[ , par.index,2]
)


boxplot(values.for.box, 
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 1),
        las = 1,
        ylim = c(4,9.5),
        #main = expression("Estimation of " * log(lambda) ~ "\u2248" ~ log(10^3) ~ "\u2248" ~ 6.908),
        main = expression(log(lambda) == log(10^3)),
        cex.main = 1.5
)
abline( h = log(lambda.s[2]) , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  "MGLE", "MCLE"
)
axis(side = 1, at = c(1,2) , cex.axis = 1.5,
     labels = labels.box, tick = FALSE)


#lambda = 10^3
values.for.box <- list(
  EST.PAR.ECM.POISSON.MGLE[ , par.index,3] , EST.PAR.ECM.POISSON.AP[ , par.index,3]
)

#Boxplot without names
boxplot(values.for.box, 
        xaxt = "n", 
        lwd = 0.5 , 
        col = rep(c(2, 4), 1),
        las = 1,
        ylim = c(4,9.5),
        cex.main = 1.5,
        #main = expression("Estimation of " * log(lambda) ~ "\u2248" ~ log(10^4) ~ "\u2248" ~ 9.210) , 
        main = expression(log(lambda) == log(10^4) )
)
abline( h = log(lambda.s[3]) , lty = 6 , col = "gold" , 
        lwd = 1.8 )


labels.box <- expression(
  "MGLE", "MCLE"
)
axis(side = 1, at = c(1,2) , cex.axis = 1.5,
     labels = labels.box, tick = FALSE)


dev.off()



#####---- Correlograms generator ----


library(GGally)
library(ggplot2)
library(rlang)


###For ECM:
for(method in c("AP","MGLE")){
  for(i.N in 1:length(N.s)){
    
    if( method == "AP" ){
      est.par <- EST.PAR.ECM.AP[  , , i.N ]
    }
    if( method == "MGLE" ){
      est.par <- EST.PAR.ECM.MGLE[  , , i.N ]
    }
    
    DF.est.par <- as.data.frame(est.par)
    names(DF.est.par) <- c( "log_tau" , "log_sigma" , "z1" , "z2" )
    
    #par
    # Medias teóricas con los nuevos nombres
    medias_teoricas <- c(log_tau = par[1], 
                         log_sigma = par[2], 
                         z1 = par[3], 
                         z2 = par[4])
    
    # Etiquetas como expresiones matemáticas
    labels <- c(
      log_tau = "log(tau)",
      log_sigma = "log(sigma)",
      z1 = "z[1]",
      z2 = "z[2]"
    )
    
    
    # Función para histogramas personalizados
    custom_hist <- function(data, mapping, ...) {
      variable <- as_label(mapping$x)
      if (!variable %in% names(medias_teoricas)) return(ggplot())
      
      mu_teorico <- medias_teoricas[[variable]]
      mu_muestral <- mean(data[[variable]], na.rm = TRUE)
      
      tmp_plot <- ggplot(data, mapping) +
        geom_histogram(bins = 25, fill = "skyblue", color = "white")
      y_max <- max(ggplot_build(tmp_plot)$data[[1]]$count)
      
      step <- y_max / 20
      segments_rojos <- seq(0, y_max, by = 2 * step)
      segments_azules <- segments_rojos + step
      
      p <- ggplot(data, mapping) +
        geom_histogram(bins = 25, fill = "skyblue", color = "white") +
        lapply(segments_rojos, function(y0) {
          annotate("segment",
                   x = mu_teorico, xend = mu_teorico,
                   y = y0, yend = y0 + step,
                   color = "red", linewidth = 0.8)
        }) +
        lapply(segments_azules, function(y0) {
          annotate("segment",
                   x = mu_muestral, xend = mu_muestral,
                   y = y0, yend = y0 + step,
                   color = "blue", linewidth = 0.8)
        })
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      
      return(p)
    }
    # Función para puntos personalizados
    custom_points <- function(data, mapping, ...) {
      x_var <- as_label(mapping$x)
      y_var <- as_label(mapping$y)
      
      mu_x_teorico <- medias_teoricas[[x_var]]
      mu_y_teorico <- medias_teoricas[[y_var]]
      
      mu_x_muestral <- mean(data[[x_var]])
      mu_y_muestral <- mean(data[[y_var]])
      
      p <- ggplot(data, mapping) +
        geom_point(alpha = 0.5, ...) +
        annotate("point", x = mu_x_teorico, y = mu_y_teorico,
                 color = "red", size = 2) +
        annotate("point", x = mu_x_muestral, y = mu_y_muestral,
                 color = "blue", fill = NA, shape = 21, size = 2, stroke = 1)
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      
      return(p)
    }
    
    theme_custom <- theme(strip.text = element_text(size = 11)) 
    
    #For png generating
    png.name <- paste("Correlogram" ,"OU" , method , "N=" , N.s[i.N] , ".png" , sep="" )
    #png( png.name, width = 2500, height = 2400, res = 360)
    
    # Correlograma final
    correlo <- ggpairs(DF.est.par,
                       diag = list(continuous = custom_hist),
                       lower = list(continuous = custom_points),
                       upper = list(continuous = "cor"),
                       labeller = label_parsed,
                       columnLabels = labels)+ theme_custom
    ggsave(png.name , 
           plot = correlo ,
           width = 2500, height = 2400 , units = "px")
    
    #dev.off()
  }
}


##For ECM-Poisson
for(method in c( "AP" , "MGLE")){
  for(i.lambda in 1:length(lambda.s)){
    
    lambda <- lambda.s[i.lambda]
    
    if(method == "MGLE" & i.lambda==1){
      est.par <- EST.PAR.ECM.POISSON.MGLE[ 1:800 ,  , i.lambda ]
    }else{
      if( method == "AP" ){
        est.par <- EST.PAR.ECM.POISSON.AP[  , , i.lambda ] 
      }
      if( method == "MGLE" ){
        est.par <- EST.PAR.ECM.POISSON.MGLE[ ,  , i.lambda ]
      }
    }
    
    
    
    DF.est.par <- as.data.frame(est.par)
    names(DF.est.par) <- c( "log_tau" , "log_sigma" , "z1" , "z2" , "log_lambda")
    
    #par
    # Medias teóricas con los nuevos nombres
    medias_teoricas <- c(log_tau = par[1], 
                         log_sigma = par[2], 
                         z1 = par[3], 
                         z2 = par[4] , 
                         log_lambda = log(lambda))
    
    # Etiquetas como expresiones matemáticas
    labels <- c(
      log_tau = "log(tau)",
      log_sigma = "log(sigma)",
      z1 = "z[1]",
      z2 = "z[2]",
      log_lambda = "log(lambda)"
    )
    
    # Función para histogramas personalizados
    custom_hist <- function(data, mapping, ...) {
      variable <- as_label(mapping$x)
      if (!variable %in% names(medias_teoricas)) return(ggplot())
      
      mu_teorico <- medias_teoricas[[variable]]
      mu_muestral <- mean(data[[variable]], na.rm = TRUE)
      
      tmp_plot <- ggplot(data, mapping) +
        geom_histogram(bins = 25, fill = "skyblue", color = "white")
      y_max <- max(ggplot_build(tmp_plot)$data[[1]]$count)
      
      step <- y_max / 20
      segments_rojos <- seq(0, y_max, by = 2 * step)
      segments_azules <- segments_rojos + step
      
      p <- ggplot(data, mapping) +
        geom_histogram(bins = 25, fill = "skyblue", color = "white") +
        lapply(segments_rojos, function(y0) {
          annotate("segment",
                   x = mu_teorico, xend = mu_teorico,
                   y = y0, yend = y0 + step,
                   color = "red", linewidth = 0.8)
        }) +
        lapply(segments_azules, function(y0) {
          annotate("segment",
                   x = mu_muestral, xend = mu_muestral,
                   y = y0, yend = y0 + step,
                   color = "blue", linewidth = 0.8)
        })
      
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      
      return(p)
    }
    # Función para puntos personalizados
    custom_points <- function(data, mapping, ...) {
      x_var <- as_label(mapping$x)
      y_var <- as_label(mapping$y)
      
      mu_x_teorico <- medias_teoricas[[x_var]]
      mu_y_teorico <- medias_teoricas[[y_var]]
      
      mu_x_muestral <- mean(data[[x_var]])
      mu_y_muestral <- mean(data[[y_var]])
      
      p <- ggplot(data, mapping) +
        geom_point(alpha = 0.5, ...) +
        annotate("point", x = mu_x_teorico, y = mu_y_teorico,
                 color = "red", size = 2) +
        annotate("point", x = mu_x_muestral, y = mu_y_muestral,
                 color = "blue", fill = NA, shape = 21, size = 2, stroke = 1)
      
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      return(p)
    }
    
    theme_custom <- theme(strip.text = element_text(size = 11)) 
    
    
    #For png generating
    png.name <- paste("Correlogram" ,"OU" , method , "lambda=" , lambda.s[i.lambda] , ".png" , sep="" )
    #png( png.name, width = 2500, height = 2400, res = 360)
    
    # Correlograma final
    correlo <- ggpairs(DF.est.par,
                       diag = list(continuous = custom_hist),
                       lower = list(continuous = custom_points),
                       upper = list(continuous = "cor"),
                       labeller = label_parsed,
                       columnLabels = labels)+ theme_custom
    ggsave(png.name , 
           plot = correlo ,
           width = 2500, height = 2400 , units = "px")
    
    #dev.off()
  }
}




#####---- Bias analysis (for Table in the paper)  ----


#Averages errors
#ECM, MGLE
t(colMeans(EST.PAR.ECM.MGLE) - par)
#ECM, MCLE
t(colMeans(EST.PAR.ECM.AP) - par)

#ECM-Poisson, MGLE
t(colMeans(EST.PAR.ECM.POISSON.MGLE) - 
    c( par , log(10^2) , par , log(10^3) , par , log(10^4) )
)
#ECM-Poisson, MCLE
t(colMeans(EST.PAR.ECM.POISSON.AP) - 
    c( par , log(10^2) , par , log(10^3) , par , log(10^4) )
)



#For the special case ECM-Poisson, MGLE, lambda = 1e2
colMeans(EST.PAR.ECM.POISSON.MGLE[1:800 , , 1])

t(colMeans(EST.PAR.ECM.POISSON.MGLE[1:800, , ]) - 
    c( par , log(10^2) , par , log(10^3) , par , log(10^4) )
)[1,]


#square-root of  mean quadratic distance distance error
#ECM, MCLE
DIFF <- apply( EST.PAR.ECM.AP , 
               c(1,3) , 
               function(x){ (x - par)^2 } )
t( sqrt( apply( DIFF , c(1,3) , mean  ) ) )


#ECM-Poisson, MCLE
DIFF <- EST.PAR.ECM.POISSON.AP
for(i in 1:1000){
  DIFF[i , , ] <- DIFF[i, , ] - rbind( matrix( par , 4 , 3) , 
                                       log(10^(2:4)) )
}
t( sqrt(colMeans(DIFF^2)) )


#ECM-Poisson, MGLE  
DIFF <- EST.PAR.ECM.POISSON.MGLE
for(i in 1:1000){
  DIFF[i , , ] <- DIFF[i, , ] - rbind( matrix( par , 4 , 3) , 
                                       log(10^(2:4)) )
}
round( t( sqrt(colMeans(DIFF^2)) ) , 3 )

  

#For the special case:
#ECM-Poisson, MGLE
DIFF <- EST.PAR.ECM.POISSON.MGLE[1:1000, , ]
for(i in 1:1000){
  DIFF[i , , ] <- DIFF[i, , ] - rbind( matrix( par , 4 , 3) , 
                                       log(10^(2:4)) )
}
round( t( sqrt(colMeans(DIFF^2)) ) , 3 )


#Use the value here for ECM-Poisson, MGLE lambda = 1e2 (special case)
DIFF <- EST.PAR.ECM.POISSON.MGLE[1:800, , ]
for(i in 1:800){
  DIFF[i , , ] <- DIFF[i, , ] - rbind( matrix( par , 4 , 3) , 
                                       log(10^(2:4)) )
}
round( t( sqrt(colMeans(DIFF^2)) ) , 3 )

