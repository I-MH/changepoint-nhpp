## it simulates one trajectory of an inhomogeneous Poisson processes 

sim.ihpp <- function(Tmax=24, lambda)
{
  # Simulates an inhomogeneous Poisson process with intensity lambda
  #
  # args :
  #   Tmax: maximum observation time 
  #   lambda: a function, f(t), corresponding to the intensity function.
  #
  # returns: a data frame with 
  #   time = a vector at where events occur 
  #   count = a vector representing the cumulative events
  
  X <- 0 
  t.aux <- 0
  lamb.max <- max( lambda(seq(0,Tmax, length.out = 256)) ) # max value (approx) of lambda 
  while(X[length(X)]<Tmax) {
    e <- -log( runif(1) )/lamb.max
    t.aux <- t.aux + e
    if(t.aux<Tmax){
      if( runif(1)<= lambda(t.aux)/lamb.max ) X=c(X,t.aux)  
    }else{
      break
    }
  }
  count = 1:(length(X)-1) 
  return( data.frame(time=X[-1], count=count) )
}


# Example -------------------------------------------
if(FALSE){
  #log intensity 1
  lamb1 <- function(s){
      z <- 20*dnorm(s, mean = 15, sd = 2) +1
  }
  #log intensity 2
  lamb2 <- function(s){
    z <- 20*dnorm(s, mean = 5, sd = 2) +1
  }
    
  curve(lamb1,  from = 0, to =24)
  curve(lamb2,  from = 0, to =24)
  
  ex.sim1 <- sim.ihpp(Tmax = 24, lambda = lamb1 ) 
  ex.sim2 <- sim.ihpp(Tmax = 24, lambda = lamb2 ) 
  plot(ex.sim1$time, ex.sim1$count, type='s')
  plot(ex.sim2$time, ex.sim2$count, type='s')
}
# End example -------------------------------------------

