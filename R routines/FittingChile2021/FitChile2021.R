rm(list = ls(all = T) )
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readxl)
library(mc2d)

functions.path <- file.path(  dirname( dirname(getwd()) ) , "R functions" )
source( file.path(  functions.path ,  "Voting_ll_Gauss.R") )

#Number of candidates
m <- 7
#Importing first round results
Chile2021.1 <- read_xlsx(path = "Resultados-Ofciales-Eleccion-Presidencial-21-de-noviembre-de-2021.xlsx" , 
                       sheet = 19)
#Take out last line which contains unnecesary information
Chile2021.1 <- Chile2021.1[-nrow(Chile2021.1) ,  ]



#Set to "EXTRANJERO" all the districts abroad
Chile2021.1$Comuna[which(is.na(Chile2021.1$Comuna))] <- "EXTRANJERO"

#Number of districts (comunas)
n.dist <- length(unique(Chile2021.1$Comuna))
names.dist <- unique(Chile2021.1$Comuna) 

#Number of voters in each district (initialization)
Nd <- rep( 0 , n.dist)

#Matrix of results in first round (initialization)
Q.1 <- matrix( 0, nrow = n.dist , ncol = m+1  )
colnames(Q.1) <- c("KAST" , 
                   "BORIC" , 
                   "PARISI" , 
                   "SICHEL" , 
                   "PROVOSTE" , 
                   "ME-O" , 
                   "ARTES" , 
                   "ABSTENCION")
row.names(Q.1) <- names.dist

order <- c(17,16,22,19,18,21,20) 

#Puting the total counts per district into the matrix.
#If there is any NA, it is transformed into 0
for(j in 1:n.dist){
  Data.dist <- as.matrix(Chile2021.1[ Chile2021.1$Comuna == names.dist[j] , c(order , 26)  ])
  Data.dist[ which(is.na(Data.dist)) ] <- 0
  Q.1[j , 1:m] <- colSums(Data.dist[ , 1:m , drop = F]) 
  Nd[j] <- sum(Data.dist[ , m+1])
  Q.1[j , m+1] <- Nd[j] - sum(Data.dist[ , 1:m])
}

#Total number of voters (abstention included)
N <- sum(Nd)


#Importing second round results
Chile2021.2 <- read_xlsx(path = "Resultados-definitivos-eleccion-Presidencial-2021.xlsx" , 
                         sheet = 19)
#Take out last line which contains unnecesary information
Chile2021.2 <- Chile2021.2[-nrow(Chile2021.2) ,  ]


#Set to "EXTRANJERO" all the districts abroad
Chile2021.2$Comuna[which(is.na(Chile2021.2$Comuna))] <- "EXTRANJERO"



#Matrix of results in first round (initialization)
Q.2 <- matrix( 0, nrow = n.dist , ncol = 3  )
colnames(Q.2) <- c("BORIC" , 
                   "KAST" , 
                   "ABSTENCION")
row.names(Q.2) <- names.dist

#Puting the total counts per district into the matrix.
#If there is any NA, it is transformed into 0
for(j in 1:n.dist){
  Data.dist <- as.matrix(Chile2021.2[ Chile2021.2$Comuna == names.dist[j] , c(16,17,21)  ])
  Data.dist[ which(is.na(Data.dist)) ] <- 0
  Q.2[j , 1:2] <- colSums(Data.dist[ , 1:2 , drop = F]) 
  Q.2[j , 3] <- Nd[j] - sum(Data.dist[ , 1:2])
}

#The data is ready.



#Results (quantities):
colSums(Q.1)
colSums(Q.2)


#Results (percentages):
#1st round
round( colSums(Q.1)/N , 4 )
round( colSums(Q.1[,1:m])/sum(Q.1[,1:m]) , 4 ) #(valid votes)

#2nd round
round( colSums(Q.2)/N , 4 )
round( colSums(Q.2[,1:2])/sum(Q.2[,1:2]) , 4 ) #(valid votes)

#Number of districts
length(Nd)




##Fitting the model

#First round estimates (just for information)
P.1.est <- Q.1/Nd

#For the second round, we have a Poisson-multinomial with Gaussian approximation.

#Number of parameters to estimate
n.par <- 2*(m+1)

# #Let us use a uniform scenario as an initial vector
P2.1.initial <- rbind( c(1/3 ,  1/3 , 1/3 , 1/3 , 1/3 , 1/3 , 1/3 , 1/3  ) ,
                       c(1/3 ,  1/3 , 1/3 , 1/3 , 1/3 , 1/3 , 1/3 , 1/3  ))


P2.1.initial<- rbind( P2.1.initial , 
                          1 - colSums(P2.1.initial))
colnames(P2.1.initial) <- colnames(Q.1)
rownames(P2.1.initial) <- colnames(Q.2)


#Reference likelihood and gradient at the initial value
P2.1.initial
Voting_ll_Gauss( P2.1.initial[1:2, ] , Q.1 , Q.2[ , 1:2]   )
Voting_ll_Gauss_Gradient( P2.1.initial[1:2, ] , Q.1 , Q.2[ , 1:2]   )

#par.0 <- as.vector( P2.1.initial[1:2 , ] )

par.0.p.scale <- as.vector(P2.1.initial[1:2, ])


#Go to the softmax scale

#Initial values
par.0 <- rep( 0 , 2*(m+1))
#lb <- c(0 , 0.95 , 0.95 , rep(0 , 2*(m+1) - 3 )  )
for(l in 1:(m+1)){
  par.0[ (2*l-1):(2*l) ] <- inv.softmax1( par.0.p.scale[(2*l-1):(2*l)] 
                                          #, lb = lb[(2*l-1):(2*l) ] 
                                          )
}


par.0


start <- Sys.time()
mgle <- optim(  par.0 , 
                fn = function(par){
                  soft.par <- rep( 0 , 2*(m+1))
                  for(l in 1:(m+1)){
                    soft.par[ (2*l-1):(2*l) ] <- softmax1( par[(2*l-1):(2*l)]
                                                           #,lb = lb[(2*l-1):(2*l) ] 
                    )
                  }
                  P2.1 <- matrix( soft.par , 2 , m+1 )
                  return( -Voting_ll_Gauss( P2.1 , Q.1 , Q.2[ , 1:2] ) )} , 
                gr = function(par){
                  soft.par <- rep( 0 , 2*(m+1))
                  for(l in 1:(m+1)){
                    soft.par[ (2*l-1):(2*l) ] <- softmax1( par[(2*l-1):(2*l)]
                                                           #,lb = lb[(2*l-1):(2*l) ] 
                    )
                  }
                  P2.1 <- matrix( soft.par , 2 , m+1 )
                  G <- -Voting_ll_Gauss_Gradient(P2.1 , Q.1 , Q.2[ , 1:2] )
                  for(l in 1:(m+1)){
                    G[ , l] <- Jac.Trans.Softmax1( par[(2*l-1):(2*l)]
                                                   #,lb = lb[(2*l-1):(2*l) ] 
                    )%*%G[,l]
                  }
                  return(as.vector(G))
                } , 
                method = "L-BFGS-B" , 
                lower = rep(  -15 , 2*(m+1) ) , 
                upper = rep( 15 , 2*(m+1)) , 
                control = list( maxit = 5e3 ,
                                factr = 1e4) , 
                hessian = T )
print(Sys.time()-start)



mgle
mgle$value
mgle$par

eigen(mgle$hessian)$values
#optim provided hessian is positive definite: good point.


mgle.p.scale <- mgle$par
for(l in 1:(m+1)){
  mgle.p.scale[ (2*l-1):(2*l) ] <- softmax1( mgle$par[(2*l-1):(2*l)]
                                             #,lb = lb[(2*l-1):(2*l) ] 
                                             )
}


P2.1.est <- matrix( mgle.p.scale , 2 , m+1)

P2.1.est <- rbind( P2.1.est , 
                   1 -  colSums(P2.1.est)  )

colnames(P2.1.est) <- colnames(Q.1)
rownames(P2.1.est) <- colnames(Q.2)

P2.1.initial
P2.1.est
round(P2.1.est, 5)


#### ----Bootstrapping

library(parallel)

n.cores <- 3
cl <- makeCluster(n.cores)

clusterExport(cl , "Nd")
clusterExport(cl , "Q.1")
clusterExport(cl , "Q.2")
clusterExport(cl , "P2.1.est")
clusterExport(cl , "par.0")
clusterExport(cl , "m")
clusterExport(cl , "n.dist")
clusterExport(cl , "n.par")


clusterEvalQ(cl , {
  library(readxl)
  library(mc2d)
  
  source( "Voting_ll_Gauss.R")
  
  N.boot.per.core <- 335
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
                 ll.boot <- matrix(  0 , N.boot )
                 q.2.boot <- array( 0 , dim = c(N.boot , n.dist , 3) )
                 
                 #The simulation loop
                 set.seed(s)
                 for(i in 1:N.boot){

                   #Simulation of second round results
                   for(j in 1:n.dist){
                     q.2.boot[i , j , ] <- colSums(rmultinomial( n = m+1 , 
                                                                 size = Q.1[j , ] , 
                                                                 prob = t(P2.1.est) ))
                   }
                   
                   #Optimization in the softmax scale 
                   mgle.boot <- optim(  par.0 , 
                                        fn = function(par){
                                          soft.par <- rep( 0 , 2*(m+1))
                                          for(l in 1:(m+1)){
                                            soft.par[ (2*l-1):(2*l) ] <- softmax1( par[(2*l-1):(2*l)]
                                                                                   #,lb = lb[(2*l-1):(2*l) ] 
                                            )
                                          }
                                          P2.1 <- matrix( soft.par , 2 , m+1 )
                                          return( -Voting_ll_Gauss( P2.1 , Q.1 , q.2.boot[i ,  ,1:2] ) )} , 
                                        gr = function(par){
                                          soft.par <- rep( 0 , 2*(m+1))
                                          for(l in 1:(m+1)){
                                            soft.par[ (2*l-1):(2*l) ] <- softmax1( par[(2*l-1):(2*l)]
                                                                                   #,lb = lb[(2*l-1):(2*l) ] 
                                            )
                                          }
                                          P2.1 <- matrix( soft.par , 2 , m+1 )
                                          G <- -Voting_ll_Gauss_Gradient(P2.1 , Q.1 , q.2.boot[i, ,1:2] )
                                          for(l in 1:(m+1)){
                                            G[ , l] <- Jac.Trans.Softmax1( par[(2*l-1):(2*l)]
                                                                           #,lb = lb[(2*l-1):(2*l) ] 
                                            )%*%G[,l]
                                          }
                                          return(as.vector(G))
                                        } , 
                                        method = "L-BFGS-B" , 
                                        lower = rep(-15 , 2*(m+1)) , 
                                        upper = rep( 15 , 2*(m+1)) , 
                                        control = list( maxit = 5e3 ,
                                                        factr = 1e4) , 
                                        hessian = F )
                   
                   
                   est.par.boot[i , ] <- mgle.boot$par
                   ll.boot[i , ] <- -mgle.boot$value
                   gc()
                   
                 }
                 
                 l <- list( q.2.boot = q.2.boot , 
                            est.par.boot = est.par.boot , ll.boot = ll.boot , 
                            N.boot=N.boot)
                 gc()
                 return(l)
               } )
#stopCluster(cl)
gc()
print(Sys.time()-start)


saveRDS(L,"L.bootstrap")
L <- readRDS("L.bootstrap")

N.boot.per.core <- L[[1]]$N.boot

#Merging the results into large lists/arrays.
Q.2.boot <- vector( "list" , length = N.boot.per.core*n.cores)
Est.par.boot <- matrix( 0 , N.boot.per.core*n.cores , n.par )
Boot.ll <- rep( 0 , N.boot.per.core*n.cores)


#Merging
for(i in 1:n.cores){
  for(j in 1:N.boot.per.core){
    Q.2.boot[[ (i-1)*N.boot.per.core + j  ]] <- L[[i]]$q.2.boot
  }
  Est.par.boot[(i-1)*N.boot.per.core + 1:N.boot.per.core, ] <- L[[i]]$est.par.boot
  Boot.ll[(i-1)*N.boot.per.core + 1:N.boot.per.core] <- L[[i]]$ll.boot[i]
}



# hist( Est.par.boot[,1] , breaks = 20)
# hist( Est.par.boot[,2] , breaks = 20)
# hist( Est.par.boot[,3] , breaks = 20)
# hist( Est.par.boot[,4] , breaks = 20)
# hist( Est.par.boot[,5] , breaks = 20)
# hist( Est.par.boot[,6] , breaks = 20)
# hist( Est.par.boot[,7] , breaks = 20)
# hist( Est.par.boot[,8] , breaks = 20)
# hist( Est.par.boot[,9] , breaks = 20)
# hist( Est.par.boot[,10] , breaks = 20)
# hist( Est.par.boot[,11] , breaks = 20)
# hist( Est.par.boot[,12] , breaks = 20)
# hist( Est.par.boot[,13] , breaks = 20)
# hist( Est.par.boot[,14] , breaks = 20)
# hist( Est.par.boot[,15] , breaks = 20)
# hist( Est.par.boot[,16] , breaks = 20)

#

#Transform to p-scale
est.boot.p.scale <- matrix( 0 , 
                            nrow = nrow(Est.par.boot) , 
                            ncol = ncol(Est.par.boot) )
for(i in 1:nrow(Est.par.boot)){                            
  for(l in 1:(m+1)){
    est.boot.p.scale[i , (2*l-1):(2*l) ] <-  softmax1( Est.par.boot[i , (2*l-1):(2*l)]
                                                       #,lb = lb[(2*l-1):(2*l) ] 
    )
  }
}

#Obtaining 95% confidence intervals

#Centred at the means.

means <- colMeans(est.boot.p.scale)

CI <- matrix( 0 , n.par , 2)
confi <- 0.95

n.sample <- ceiling(  nrow(est.boot.p.scale)*confi  )
for(i in 1:n.par){
  dists <- abs(est.boot.p.scale[,i] - means[i])
  ord <- order(dists)
  sub.sample <- est.boot.p.scale[ ord[1:n.sample] ,i]
  CI[i  ,  ] <- c(min(sub.sample) , max(sub.sample))
}

CI
round(CI , 5)

# hist(est.boot.p.scale[,1] , breaks = 50)
# hist(est.boot.p.scale[,2] , breaks = 50)
# hist(est.boot.p.scale[,3] , breaks = 50)
# hist(est.boot.p.scale[,4] , breaks = 50)
# hist(est.boot.p.scale[,5] , breaks = 50)
# hist(est.boot.p.scale[,6] , breaks = 50)
# hist(est.boot.p.scale[,7] , breaks = 50)
# hist(est.boot.p.scale[,8] , breaks = 50)
# hist(est.boot.p.scale[,9] , breaks = 50)
# hist(est.boot.p.scale[,10] , breaks = 50)
# hist(est.boot.p.scale[,11] , breaks = 50)
# hist(est.boot.p.scale[,12] , breaks = 50)
# hist(est.boot.p.scale[,13] , breaks = 50)
# hist(est.boot.p.scale[,14] , breaks = 50)
# hist(est.boot.p.scale[,15] , breaks = 50)
# hist(est.boot.p.scale[,16] , breaks = 50)

#Obtaining the estimates of abstention
est.abstention.boot <- t(apply(  est.boot.p.scale   , 1 , function(x){  1-x[ 2*(1:(m+1)) ]-x[ 2*(1:(m+1))-1]}))


means.abs <- colMeans( est.abstention.boot )

#means.abs

CI.abs <- matrix( 0 , n.par/2 , 2)
confi <- 0.95

n.sample.abs <- ceiling(  nrow(est.abstention.boot)*confi  )
for(i in 1:(n.par/2)){
  dists <- abs(est.abstention.boot[,i] - means.abs[i])
  ord <- order(dists)
  sub.sample <- est.abstention.boot[ ord[1:n.sample] ,i]
  CI.abs[i  ,  ] <- c(min(sub.sample) , max(sub.sample))
}


round( CI.abs , 5 )


# hist(est.abstention.boot[ , 1] , breaks = 20)
# hist(est.abstention.boot[ , 2] , breaks = 20)
# hist(est.abstention.boot[ , 3] , breaks = 20) 
# hist(est.abstention.boot[ , 4] , breaks = 20) 
# hist(est.abstention.boot[ , 5] , breaks = 20) 
# hist(est.abstention.boot[ , 6] , breaks = 20) 
# hist(est.abstention.boot[ , 7] , breaks = 20) 
# hist(est.abstention.boot[ , 8] , breaks = 20) 