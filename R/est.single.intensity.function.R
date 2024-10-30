est.single.intensity.function <- function(data, basis, L=0, U=24, iterations=1000){
  # it estimates intensity functions for a single period (or sequence)
  # - does not assume replicates -
  # args:
  #   data: data frame with variable 
  #         time: vector at where data is observed for each period
  #         n: a vector indicator for periods.  
  #         this could be a list (not implemented yet)
  #   basis: basis object of the type bspline
  #   (L,U): time interval where data is observed
  #
  # values: a list 
  #   n: the index of the periods
  #   w.hat: a matrix containing the coeff estimated
  #   
  
  t0 <- Sys.time()
  
  n.id <- unique(data$n) # period's index
  n <- length( n.id ) # total number of periods
  
  w.hat <- matrix(NA, basis$nbasis, n) # coeff estimated
  cost.value <- rep(NA,n) 
  for (j in 1:n ) {
    # cat('estimatign intensity for period ', j, '\n',sep="")
    data_n <- data[which(data$n == n.id[j]), ]
    
    # use mle instead find.optimal.parameters
    result <- mle(basis = basis, 
                  data=data_n, 
                  L=L, U=U, 
                  iterations=iterations)
    
    w.hat[,j] <- result$par 
    cost.value[j] <- result$value 
  }
  
  sgrid <- seq(L,U, length=64)
  intensity.all <- exp( eval.fd(sgrid, fd( w.hat,basis)) )
  
  t1 <- Sys.time()
  #cat('============================================== \n',sep="")
  #cat('cpu: ', t1-t0, ' \n', sep="")
  
  return(list(n=n.id,
              w.hat = w.hat,
              intensity.mat =intensity.all,
              sgrid=sgrid,
              log.intensity.fda = fd( w.hat,basis),
              cost.value=cost.value,
              cpu=t1-t0))
  
}

