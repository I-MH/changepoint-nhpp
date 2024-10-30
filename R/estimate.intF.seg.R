
#source('~/Dropbox/personal/LancasterSRA/Project1/RCode/R/likelihood.R')

estimate.intF.seg <- function(data, basis, cpts, L, U, iterations=100){
  # it computes intensity functions given a segmentation
  # args:
  #   data: data frame with variable 
  #         time: vector at where data is observed for each period
  #         n: a vector indicator for periods.  
  #   basis: basis object of the type bspline
  #   cpts: a vector with changepoints
  #   (L,U): time interval where data is observed
  #   iterations=1000
  # values: a list  
  # w.hat: matrix containing the coeff estimated for the intensity function
  #   intensity.mat: matrix with descrite values of each intensity function 
  #   sgrid: grid use to evaluate intensity.mat
  #   log.intensity.fda: fda object 
  #   cost.value
  #   data.seg: list with segmented data as elements 

  n <-  length(unique(data$n))
  w.hat <- matrix(NA, basis$nbasis, length(cpts) +1 ) # coeff estimated
  cost.value <- rep(NA,length(cpts) +1) 
  data.segm <- vector('list', length(cpts) +1) # segmented data
  
  segm <- c(0, cpts, n)
  for ( k in 1:(length(segm)-1) ){
    datak <- data[(data$n>segm[k])&(data$n<=segm[k+1] ),]
    estimation <- mle(basis = basis, 
                                      data=datak, 
                                      L=L, U=U, 
                                      iterations=iterations)
    w.hat[,k] <- estimation$par 
    cost.value[k] <- estimation$value 
    data.segm[[k]] <- datak
  }
  names(data.segm) <- paste0('seg',1:(length(segm)-1) )
  sgrid <- seq(L,U, length=96) # to evaluate the cont version of the int functions
  intensity.all <- exp( eval.fd(sgrid, fd( w.hat,basis)) )
  
  return(list(w.hat = w.hat,
              intensity.mat =intensity.all,
              cpts=cpts,
              sgrid=sgrid,
              log.intensity.fda = fd( w.hat,basis),
              cost.value=cost.value,
              data.seg= data.segm))
}



