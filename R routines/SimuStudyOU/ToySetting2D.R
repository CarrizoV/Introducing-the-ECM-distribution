##Toy setting 2D
library(data.table)

set.seed(131313)
t0 <- 0
n.t <- 10
t <- sort(runif(n.t ,  0 , 10 ))

d <- 2
B.L <- vector("list", length = n.t)

m <- ceiling(runif(n.t , 10 , 50)) #Between 10 and 50 windows per time
dx <- 0.1  #Unit length

#Domain: [-1 , 1]x[-1 , 1] +- dx/2
b.l.grid.coords <- as.matrix(CJ( (-10:10)*dx - dx/2 , (-10:10)*dx) - dx/2 )

#the associated regular grid has 21*21 = 441 centers (or square regions).
for(k in 1:n.t){
  B.L[[k]] <- b.l.grid.coords[ sample(1:441,m[k]) ,  ] 
}

n.d <- sum(m)
rm(b.l.grid.coords)
gc()


