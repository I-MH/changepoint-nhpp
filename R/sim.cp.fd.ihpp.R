
# It simulates a sequence of ihpp with multiple changepoints across periods.

# sim.fd.ihpp.R contains the function sim.fd.ihpp() 
# which simulates a sequence of inhomogeneous PP (with no changepoints)
source("R/sim.fd.ihpp.R")

sim.cp.fd.ihpp <- function(N=200, tau = 0.5, changepoints=NULL, xlim=c(0,24), lambda){
  # it simulates a sequence of ihpp with change points at tau,
  # args:
  #   N: sample size of the functional ihpp
  #   tau: location of the change
  #   xlim:time interval where CP is simulated
  #   lambda: a list with intensity functions 
  #
  # values: a data frame with
  #   time: the time at where events are observed
  #   n: indicator of periods
  #   group: the group to which the ihpp belong. 
  
  ntau <- length(tau) #number of changepoints
  if(is.null(changepoints)) changepoints=floor(N*tau) #changepoints
  seg.at <- c(0, changepoints, N)
  sample.seg <- diff(seg.at)
  
  if( length(lambda) != length(sample.seg) ){
    stop(paste0('lambda must be of length=',length(sample.seg),'.') ) } 
  
  data <- NULL
  group <- NULL
  for (j in 1:length(sample.seg) ) {
    sim.data <- sim.fd.ihpp(N=sample.seg[j] , 
                            xlim=xlim, 
                            lambda=lambda[[j]])
    sim.data$n <- sim.data$n + seg.at[j] # correction of index n
    sim.data$group <- rep(j,length(sim.data$time)) # add ind of segmentation
    data=rbind(data, sim.data)
  }
  
  return(data)
}


# example ---------------------------------------------
if(FALSE){
  #log intensity 1
  lamb1 <- function(s){
    z <- 20*dnorm(s, mean = 15, sd = 2) +.5
  }
  #log intensity 2
  lamb2 <- function(s){
    z <- 20*dnorm(s, mean = 5, sd = 2) +1
  }
  
  N=20
  tau=0.5
  all.lambda <- c(lamb1, lamb2)
  data <- sim.cp.fd.ihpp(N=N, tau = tau, xlim=c(0,24), lambda=all.lambda)
  
  library(ggplot2)
  pl=ggplot(data = data)
  pl+ geom_step( aes(time, count, group=n, color = as.factor(group) ) )
}
#end example ---------------------------------------------
