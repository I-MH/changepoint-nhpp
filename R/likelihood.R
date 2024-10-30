
library(fda)

loglambda <- function(s,w,basis){
  # log intensity function in terms of a b-spline basis function
  # args:
  #   s: continuous parameter 
  #   w: vector of weight with length nbasis
  #   basis: basis object of the type bspline
  #
  # values: numeric 
  #   the value of the logintensity function at s

  w.fd <- fd(w,basis)
  eval <- as.numeric( eval.fd(s,w.fd) )
  return(eval)
}

lambda <- function(s,w,basis){
  # Intensity function
  # args:
  #   s: continuous parameter 
  #   w: vector of weigts with length nbasis
  #   basis: basis object of the type bspline
  #
  # values: numeric 
  #   the value of the intensity function at s
  exp( loglambda(s,w=w, basis=basis) )
}


neg.log.likelihood <- function(w, data, basis, L=0, U=24) {
  # negative log-likelihood with replicates
  # args:
  #   w: coeff vector to be estimated
  #   data: data frame with variable 
  #         time: vector at where data is observed for each period
  #         n: a vector indicator for periods.  
  #         this could be a list (not implemented yet)
  #   basis: basis object of the type bspline
  #   (L,U): time interval where data is observed
  #
  # values: numeric 
  #   evaluation of the negative log-likelihood
  
  n.id <- unique(data$n) # period's index
  n <- length( n.id ) # total number of periods
  value1 <- n*integrate(lambda, L, U, w, basis ,subdivisions=100)$value # first term
  value2 <- 0
  for (j in n.id) {
    value2 <- value2 + sum( loglambda(s=data$time[data$n==j], w=w, basis = basis) )
  }
value1 -value2
}

#"Nelder-Mead"
#'BFGS'  (mejor)
#"SANN" 
# 'L-BFGS-B'

mle <- function(basis, data, L=0, U=24, iterations=1000) {
  # it optimizes the neg.log.likelihood 
  # it contains the cost function
  w0 <- rep.int(0,basis$nbasis)
  res <- optim(par=w0, fn=neg.log.likelihood, NULL, method = 'Nelder-Mead',
               data, basis, L, U,
               control=list(maxit=iterations))
  
  # If method = 'BFGS' this returns Error in many cases Error in integrate(lambda, L, U, w, basis, subdivisions = 100) : 
  # non-finite function value. So better use the default 'Nelder-Mead' method.
  return(res)
}



# Example----------------------------------
#
if(FALSE){
  source('R/sim.fd.ihpp.R')
  source('R/sim.ihpp.R')
  
  set.seed(3128)
  
  lamb2 <- function(s){
    
    z <- 20*dnorm(s, mean = 5, sd = 2) +1
  }
  curve(lamb2,  from = 0, to =24, ylim=c(0,6))
  toy.data <- sim.fd.ihpp(N=20, xlim = c(0,24), lambda = lamb2 )
  
  toy.data1 <- toy.data[toy.data$n==1,] # to test with only one replication
  toy.data2 <- toy.data                 # to test with replications
  plot(toy.data1$time, toy.data1$count, type = 's')
  
  # Estimate the basis coeff
  # define knots 
  knots <- c(0,3:8, 16, 24)
  norder <- 4
  basis <- create.bspline.basis(rangeval = c(-.5,24.5), 
                                breaks=knots, 
                                norder = norder)# spline basis
  result1 <- mle(basis = basis, 
                 data=toy.data1, 
                 L=0, U=24, 
                 iterations=1000) # one replicate 
  result2 <- mle(basis = basis, 
                 data=toy.data2, 
                 L=0, U=24, 
                 iterations=1000) # several replicates
  
  w.hat1 <- round(result1$par, 3)
  w.hat2 <- round(result2$par, 3)
  s <- seq(0,24, length= 100)
  plot(s,lambda(s=s,w=w.hat1, basis = basis), type = 'l' , 
       lwd=2, ylab='',lty=1, ylim = c(0,6))
  curve(lamb2,  from = 0, to =24, add = TRUE, col=4, lwd=2, lty=2)
  lines(s,lambda(s=s,w=w.hat2, basis = basis), col=2, lwd=2, lty=1 )
  legend("topright", legend = c('1 rep','n rep', 'true'), col = c(1,2,4),
         ncol = 3, cex = 1, lwd = 2, bty='n')
  
  # try PELT
  source('PELT.tpp.R')
  
  PELT.tpp(toy.data2,basis=basis) # 23:17-23:21, 23:25
  
}
