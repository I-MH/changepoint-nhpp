
# it uses fda package to create bspline basis functions.


library(fda) 

spline.basis.quantile <- function(data ,n.knots=5, 
                                  rangeval=NULL, 
                                  norder=4){
  # n.knots: total number of interior knots
  all.times <- sort(unique(data$time))
  if(is.null(rangeval)){
    a = min(all.times)        	# lower endpoint
    b = max(all.times)  
  }else{
    a = rangeval[1]
    b = rangeval[2] 
  }
  int.knots = quantile(all.times,
                       seq(0,1,length= (n.knots+2))[-c(1,(n.knots+2))]) 	# interior knots
  basis = create.bspline.basis(c(a,b),
                               norder = norder,
                               breaks=c(a,int.knots,b))	
  return(list(basis=basis, L=a, U=b, intKnots=int.knots))
}

spline.basis.eq <- function(data, n.knots=5, 
                            rangeval=NULL,
                            norder=4){
  all.times <- sort(unique(data$time))
  if(is.null(rangeval)){
    a = min(all.times)        	# lower endpoint
    b = max(all.times)  
  }else{
    a = rangeval[1]
    b = rangeval[2] 
  }
  int.knots = seq(a,b,length=n.knots+2)[-c(1,n.knots+2)]	# interior knots
  basis = create.bspline.basis(c(a,b),
                               norder = norder,
                               breaks=c(a,int.knots,b))	
  return(list(basis=basis, L=a, U=b, intKnots=int.knots))
}

select.optimal.basis.quan <- function(data, n.int.knots=5:40,
                                      L=0, U=24,
                                      iterations=100){
  # returns the best basis function (quantile method) based on AIC
  # args:
  #   data: data frame containing data, this can be a segment of data
  #   n.int.knots: a vector containing the total interior knots to be tested
  #   (L,U): domain of intensity function
  # val: a list
  #   basis: fda object - optimal basis function
  #   optimal.int.knots: a number representing the optimal total 
  #                      number of interior knots.
  #   AIC: a vector of AIC values for different values of n.int.knots
  #   AICn: a vector of AIC (finite sample version) values 
  #         for different values of n.int.knots
  
  nk <- length(n.int.knots)
  npoints <- dim(data)[1] # sample size
  AIC <- AICn <- rep(0, nk)
  
  for (j in 1:nk){
    bsf <- spline.basis.quantile(data = data, 
                                 rangeval= c(-1,25),
                                 n.knots=n.int.knots[j]) 
    npar <- bsf$basis$nbasis
    lh <-  mle(data=data, 
               basis=bsf$basis, 
               L=L, U=U,
               iterations=iterations)
    
    AIC[j] <- 2*npar - 2*lh$value
    AICn[j] <- AIC[j] + ( 2*(npar^2) + 2*npar )/(npoints-npar-1)
  }
  
  best.int.knots <- n.int.knots[which(AIC==min(AIC))] 
  best.basis <- spline.basis.quantile(data = data, 
                                      rangeval= c(-1,25),
                                      n.knots=best.int.knots) 
  
  return(list(basis=best.basis, 
              opt.nknots= best.int.knots,
              AIC= AIC,
              AICn= AICn,
              nknots.tested= n.int.knots))
}

select.optimal.basis.eq <- function(data, n.int.knots=5:40,
                                      L=0, U=24,
                                      iterations=100){
  # returns the best basis function (eq spaced method) based on AIC
  # args:
  #   data: data frame containing data, this can be a segment of data
  #   n.int.knots: a vector containing the total interior knots to be tested
  #   (L,U): domain of intensity function
  # val: a list
  #   basis: fda object - optimal basis function
  #   optimal.int.knots: a number representing the optimal total 
  #                      number of interior knots.
  #   AIC: a vector of AIC values for different values of n.int.knots
  #   AICn: a vector of AIC (finite sample version) values 
  #         for different values of n.int.knots
  
  nk <- length(n.int.knots)
  npoints <- dim(data)[1] # sample size
  AIC <- AICn <- rep(0, nk)
  
  for (j in 1:nk){
    bsf <- spline.basis.eq(data = data, 
                           rangeval= c(-1,25),
                           n.knots=n.int.knots[j]) 
    npar <- bsf$basis$nbasis
    lh <-  mle(data=data, 
               basis=bsf$basis, 
               L=L, U=U,
               iterations=iterations)
    
    AIC[j] <- 2*npar - 2*lh$value
    AICn[j] <- AIC[j] + ( 2*(npar^2) + 2*npar )/(npoints-npar-1)
  }
  
  best.int.knots <- n.int.knots[which(AIC==min(AIC))] 
  best.basis <- spline.basis.eq(data = data, 
                                rangeval= c(-1,25),
                                n.knots=best.int.knots) 
  
  return(list(basis=best.basis, 
              opt.nknots= best.int.knots,
              AIC= AIC,
              AICn= AICn,
              nknots.tested= n.int.knots))
}

if(FALSE){
  load("../Israel-Jess/HowzData_10sites2.RData")
  Data <- HowzData_10sites2
  ndays <- length( unique(Data$n) ) # total of days (periods)
  sites <- unique(Data$hid) # sites
  ns=6
  data.one.site <- Data[Data$hid==sites[ns],]
  data <- data.one.site[which(data.one.site$n %in% 1:3), ]
  best.basis <- select.optimal.basis.quan(data=data, 
                                          n.int.knots=5:40,
                                          L=0, U=24,
                                          iterations=100)
  plot(best.basis$basis$basis)
  plot(best.basis$AIC)
  abline(v=best.basis$opt.nknots, col=2)
}
