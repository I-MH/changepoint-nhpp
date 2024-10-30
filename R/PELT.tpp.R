#### Function to Run PELT to identify changes in #####
#### the intensity function of a inhomogeneous   #####
#### Poisson process                             #####

# likelihood.R contains the function find.optimal.parameters()
# which provides the optimal parameters and likelihood for a
# segment

#source('R/likelihood.R')

# here data is assumed to be a three column data.frame with columns: time, count
# and n.  time is the within-day indexing, count is the incidence count for the
# day and n is index of the day.  We want to work on the changes across days so
# we see the column n as containing the index for segmentation.

PELT.tpp=function(data,basis,
                  pen=(basis$nbasis+1)*log(length(unique(data$n))),
                  minseglen=1,L=0,U=24,iterations=100){
  # PELT for temporal point processes
  cost=mle
  
  lastchangecpts=lastchangelike=NULL
  lastchangelike[1]=-pen # no cpt scenario
  lastchangecpts[1]=0

  checklist=0

  for(i in 1:(2*minseglen-1)){
    lastchangelike[i+1]=cost(basis = basis, 
                             data=data[data$n<=i,], 
                             L=L, U=U, 
                             iterations=iterations)$value
    lastchangecpts[i+1]=0
  }

  for(tstar in (2*minseglen):max(data$n)){
    tmplike=NULL
    tmpt=c(checklist, tstar-minseglen)
    for(t in 1:length(tmpt)){
      tmplike[t]=lastchangelike[tmpt[t]+1]+cost(basis = basis, 
                                          data=data[(data$n>tmpt[t])&(data$n<=tstar),], 
                                          L=L, U=U, 
                                          iterations=iterations)$value+pen
    }
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=tmpt[tmplike==lastchangelike[tstar+1]][1]
    checklist=tmpt[tmplike<=lastchangelike[tstar+1]+pen] # step 4 in PELT
  }
  fcpt=NULL
  last=max(data$n)
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(cpt=sort(fcpt))
}

# PELT with best basis function


PELT.tpp.estbasis=function(data,
                           basis.minknots=6, 
                           basis.maxknots=48, 
                           basis.type=c("quantile","classic"),
                           pen=basis.maxknots*log(length(unique(data$n))),
                           minseglen=1,L=0,U=24,iterations=1000){
  # PELT for temporal point processes
  if(basis.type!="classic" & basis.type!="quantile"){stop("Unknown basis.type, must be 'quantile' or 'classic'")}
  if(basis.minknots>basis.maxknots){stop("basis.minknots must be less than or equal to basis.maxknots")}
  
  mle.tpp=function(data,
                   basis.minknots=basis.minknots,
                   basis.maxknots=basis.maxknots,
                   basis.type=basis.type,
                   L=L,U=U, iterations=100){
    if(basis.type=="classic"){
      basis=select.optimal.basis.eq(data=data, 
                                    n.int.knots=basis.minknots:basis.maxknots,
                                    L=L, U=U,iterations=100)
    }
    else{
      basis=select.optimal.basis.quan(data=data, 
                                      n.int.knots=basis.minknots:basis.maxknots,
                                      L=L, U=U,iterations=100)
    }
    return(mle(basis$basis$basis, data=data, 
               L=L, U=U, iterations=iterations)$value)
  }
  
  cost=mle.tpp
  
  lastchangecpts=lastchangelike=NULL
  lastchangelike[1]=-pen # no cpt scenario
  lastchangecpts[1]=0
  
  checklist=0
  
  for(i in 1:(2*minseglen-1)){
    lastchangelike[i+1]=cost(data=data[data$n<=i,] )
    lastchangecpts[i+1]=0
  }
  
  cost(basis = basis, 
       data=data[data$n<=i,], 
       L=L, U=U, 
       iterations=iterations)$value
  
  for(tstar in (2*minseglen):max(data$n)){
    tmplike=NULL
    tmpt=c(checklist, tstar-minseglen)
    for(t in 1:length(tmpt)){
      tmplike[t]=lastchangelike[tmpt[t]+1]+cost(basis = basis, 
                                                data=data[(data$n>tmpt[t])&(data$n<=tstar),], 
                                                L=L, U=U, 
                                                iterations=iterations)$value+pen
    }
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=tmpt[tmplike==lastchangelike[tstar+1]][1]
    checklist=tmpt[tmplike<=lastchangelike[tstar+1]+pen]
  }
  fcpt=NULL
  last=max(data$n)
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(cpt=sort(fcpt))
}



