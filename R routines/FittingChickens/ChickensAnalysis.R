#rm(list = ls(all = T) )
#gc()
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Q <- readRDS("Q" )
B.L <- readRDS("B.L")
dx <- readRDS("dx")
t <- readRDS( "t")

m <- sapply(B.L , nrow)
n.d <- sum(m)
n.t <- length(t)
d <- ncol(B.L[[1]])


# 
# #Code for ploting
# 
# #png("CountChickens.png" , height = 800 , width = 1500 )
# 
# # Crear una grilla de layout: 2 filas x 5 columnas
# layout(matrix(1:10, nrow = 2, byrow = TRUE))
# 
# # Márgenes exteriores comunes
# par(oma = c(2, 2, 2, 1))  # bottom, left, top, right
# 
# # Iterar sobre los tiempos
# for(k in 1:n.t){
#   show_y_axis <- (k %% 5 == 1)
#   
#   # Márgenes internos por gráfico
#   par(mar = c(2.5, ifelse(show_y_axis, 0.1, 0.1), 2.2, 0.1))  # bottom ajustado
#   
#   plot(NULL, xlim = c(-0.7, 0.7), ylim = c(-0.7, 0.7),
#        xlab = "",
#        ylab = ifelse(show_y_axis, "", ""), 
#        xaxt = "s", yaxt = ifelse(show_y_axis, "s", "n"),
#        main = bquote(t[.(k)] == .(round(t[[k]], 2)) ~ "days")
#   )
#   
#   rect(B.L[[k]][,1], B.L[[k]][,2], 
#        B.L[[k]][,1] + dx, B.L[[k]][,2] + dx)
#   
#   text(B.L[[k]][,1] + dx/2, B.L[[k]][,2] + dx/2, 
#        labels = Q[[k]], cex = 0.9)
# }
# 
# #dev.off()



#Population size
N <- 93

library(numDeriv)
library(gtools)
library(GGally)

#Analysis of the fitting
L <- readRDS("L.bootstrap.AP.Mix.Brown.OU.rds")
AP.mix.Brown.OU <- readRDS("AP.Mix.Brown.OU.rds")
Est.par.boot <- readRDS("Est.par.boot AP.Mix.Brown.OU.rds")
Boot.ll <- readRDS("Boot.ll AP.Mix.Brown.OU.rds")
Grad.boot <- readRDS("Grad.boot AP.Mix.Brown.OU.rds")
Hess.boot <- readRDS("Hess.boot AP.Mix.Brown.OU.rds")


#Taking out the erratic ones
ind.erratic <- NULL
for(k in 1:dim(Hess.boot)[1]){
  if( any(eigen(Hess.boot[k , , ])$values < 1e-12) ){
    ind.erratic <- c(ind.erratic , k)
  }
}
ind.erratic
n.erratic <- length(ind.erratic)
n.erratic

#We take out the erratic bootstrap samples 
Est.par.boot <- Est.par.boot[-ind.erratic , ]


#We fix a desired quantity of bootstrap samples
N.desired <- 1000
if(nrow(Est.par.boot) > N.desired){
  Est.par.boot <- Est.par.boot[ 1:N.desired , ]
}



AP.mix.Brown.OU
AP.mix.Brown.OU$par
colMeans(Est.par.boot)
cov(Est.par.boot)
cor(Est.par.boot)




##Correlogram in log and logit scales
library(GGally)
library(ggplot2)
library(rlang)

DF.est.par <- as.data.frame(Est.par.boot)
names(DF.est.par) <- c( "log_tau" , "log_sigma" , "log_negt0" , "logit_alpha" )

#par
# Medias teóricas con los nuevos nombres
medias_teoricas <- c(log_tau = AP.mix.Brown.OU$par[1], 
                     log_sigma = AP.mix.Brown.OU$par[2], 
                     log_negt0 = AP.mix.Brown.OU$par[3], 
                     logit_alpha = AP.mix.Brown.OU$par[4])

# Etiquetas como expresiones matemáticas
labels <- c(
  log_tau = "log(tau)",
  log_sigma = "log(sigma)",
  log_negt0 = "log(-t[0])",
  logit_alpha = "logit(alpha)"
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
png.name <- paste("CorrelogramBootstrapChickens",".png" , sep="" )
png( png.name, width = 2500, height = 2400, res = 360)

# Correlograma final
ggpairs(DF.est.par,
        diag = list(continuous = custom_hist),
        lower = list(continuous = custom_points),
        upper = list(continuous = "cor"),
        labeller = label_parsed,
        columnLabels = labels) + theme_custom


dev.off()




#Transforming to the original scales
Est.par.boot.original <- cbind( exp(Est.par.boot[,1:2]), 
                                -exp(Est.par.boot[,3]) ,
                                inv.logit(Est.par.boot[,4]))

est.par.original <- c( exp(AP.mix.Brown.OU$par[1:2]) , 
                       -exp(AP.mix.Brown.OU$par[3])  ,
                       inv.logit(AP.mix.Brown.OU$par[4]) )

est.par.original


#Obtention of 95% bootstrap confidence intervals
N.boot <- nrow(Est.par.boot.original)
confi <- 0.95

n.sample <- ceiling(confi*N.boot)

n.par <- ncol(Est.par.boot.original)


#Distances to the estimate
CI <- matrix( 0 , 4 , 2)
for(i in 1:n.par){
  ord <- order( abs(Est.par.boot.original[,i]-est.par.original[i]) )
  CI[i,1] <- min( Est.par.boot.original[ord[1:n.sample],i] )
  CI[i,2] <- max( Est.par.boot.original[ord[1:n.sample],i] )
}
est.par.original
CI


png("HistogramBootstrapChickens.png" , height = 500 , width=2000 , res = 300)

#Histograms in original scales
par(mfrow=c(1,4) , 
    mar = c(2,2,2,2))

hist(Est.par.boot.original[,1] , breaks = 20 , 
     xlim = c(0,0.1) , main = expression(tau) , xlab = "" , ylab="",
     cex.main = 1.5)
abline( v = est.par.original[1] , col = "red" , lwd = 1.5 )
abline( v = colMeans(Est.par.boot.original)[1] , col = "blue" , lwd = 1.5  )
abline( v = CI[1,1] , lty = 2 , lwd = 2 )
abline( v = CI[1,2] , lty = 2 , lwd = 2 )



hist(Est.par.boot.original[,2] , breaks = 20 , 
     main = expression(sigma) , xlab = "" , ylab = "" , 
     cex.main = 1.5)
abline( v = est.par.original[2] , col = "red" , lwd = 1.5  )
abline( v = colMeans(Est.par.boot.original)[2] , col = "blue" , lwd = 1.5  )
abline( v = CI[2,1] , lty = 2 , lwd = 2 )
abline( v = CI[2,2] , lty = 2 , lwd = 2 )



hist(Est.par.boot.original[,3] , breaks = 20 , 
     main = expression(t[0]) , xlab = "" , ylab = "",
     cex.main = 1.5)
abline( v = est.par.original[3] , col = "red" , lwd = 1.5   )
abline( v = colMeans(Est.par.boot.original)[3] , col = "blue" , lwd = 1.5   )
abline( v = CI[3,1] , lty = 2 , lwd = 2 )
abline( v = CI[3,2] , lty = 2 , lwd = 2 )



hist(Est.par.boot.original[,4] , breaks = 20, xlim = c(0,1) ,
     main = expression(alpha) , xlab = "" , ylab ="",
     cex.main = 1.5)
abline( v = est.par.original[4] , col = "red" , lwd = 1.5  )
abline( v = colMeans(Est.par.boot.original)[4] , col = "blue" , lwd = 1.5   )
abline( v = CI[4,1] , lty = 2 , lwd = 2 )
abline( v = CI[4,2] , lty = 2 , lwd = 2 )


dev.off()
