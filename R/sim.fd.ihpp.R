# it simulates a sequence of inhomogeneous Poisson processes

# sim.ihpp.R contains the function sim.ihpp() 
# which simulates one trayectory of an inhomogeneous PP
source("R/sim.ihpp.R")

sim.fd.ihpp <- function(N=50, 
                        xlim=c(0,24), 
                        lambda){
  # it simulates a sequence of ihpp using function sim.ihpp
  # args:
  #   N: sample size of the functional data
  #   xlim: time interval where ihpp is simulated
  #   lambda: intensity function
  #
  # values: a data frame with
  #   time: the time at where events are observed
  #   n: indicator of periods (index of functional data)
  
  # this conditional allow us to simulate Cox processes.
  all.lambda <- lambda
  if(!is.list(lambda)){
    all.lambda <- vector('list', N) # each 
    for (j in 1:N) {
      all.lambda[[j]] <- lambda
    }
  }else{
    if( is.list(lambda) & length(lambda)!=N) 
      stop("length of lambda must be the same as 'N' ")
  }
  
  time <- NULL # time at where events are observed
  n.indx <- NULL # index of functional data
  count <- NULL # cumulative events
  
  for (j in 1:N) {
    data <- sim.ihpp(Tmax = xlim[2], 
                     lambda = all.lambda[[j]]) 
    time <- c(time,data$time)
    n.indx <- c(n.indx, rep(j,length(data$time)) )
    count <- c(count, data$count )
  }
  
  return( data.frame(time=time, count=count, n=n.indx) )
}


# Example ---------------------------------
if(FALSE){
  lamb2 <- function(s){
    z <- 20*dnorm(s, mean = 5, sd = 2) +1
  }
  
  ex.simfd <- sim.fd.ihpp(N=20, xlim = c(0,24), lambda = lamb2 ) 
  library(ggplot2)
  pl=ggplot(data = ex.simfd)
  pl+ geom_path( aes(time, count, group=n, color = as.factor(n) ) )
}

#end example -----------------------------


