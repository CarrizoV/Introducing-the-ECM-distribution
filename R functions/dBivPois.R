#Coded by Ricardo Carrizo V
#Version at the end of January

dBivPois <- function( q1 , q2 =NULL, lambda , log = T ){
  #q1 assumed matrix (n x 2)
  if(is.matrix(lambda)){
    return(apply(  cbind( q1 , lambda  ) , 1 , 
                 function(x){
                   j <- 0:min(x[1],x[2])
                   return(logSumExp( 
                     colSums( dpois(  x = rbind( x[1]-j , x[2] - j , j) , 
                                      lambda = c( x[3]  , x[4] , x[5] ) , 
                                      log = T  ) )))
                 } ) )
  }
}





