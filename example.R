
# an example to detect a multiple changepoint ------------------------

# required packages:
library(fda); library(ggplot2)
# sim.cp.fd.ihpp.R contains function sim.cp.fd.ihpp() which is used 
# to simulate seq of ipp with changepoints
source("R/sim.cp.fd.ihpp.R")
source("R/likelihood.R")
#source("R/estimate.single.change.R")
source('R/spline.basis.R')
#source('R/choose_best_basis.R')
source('R/PELT.tpp.R')

############################################################################
# Example 1
############################################################################

#simulate data with 3 diff intensity functions
lamb1 <- function(s){
  z <- 20*dnorm(s, mean = 6, sd = 3) +1
}
lamb2 <- function(s){
  z <- 20*dnorm(s, mean = 12, sd = 3) +1
}
lamb3 <- function(s){
  z <- 20*dnorm(s, mean = 18, sd = 3) +1
}
curve(lamb1, 0,24, col=1, ylim=c(0,5), lwd=2, xlab='s', ylab = '',
      main='Intensity functions used to simulate data')
curve(lamb2, 0,24, col=2, lwd=2, add = TRUE)
curve(lamb3, 0,24, col=3, lwd=2, add = TRUE)

# set parameters---------------------------------------------------------------
N=20 # sample size (total number of periods) 
tau=c(0.35,0.65) # changepoints
changepoints <- floor(N*tau) # location of the changepoints
sample.seg <- diff(c(0, changepoints, N)) # sample size for each segment
s.range <- c(0,24) # domain of the process
all.lambda <- c(lamb1, lamb2, lamb3 )

# simulate data----------------------------------------------------------------
set.seed(2949)
toy.data <- sim.cp.fd.ihpp(N=N, tau = tau, xlim=s.range, lambda=all.lambda)

# plot data simulated --------------------------------------------------------
pl=ggplot(data = toy.data)
pl + geom_step( aes(time, count, group=n, color = as.factor(group) ) ) +
  scale_colour_discrete("segment")+
  theme(plot.title = element_text(lineheight =.8,face = "bold",size = 16,hjust = 0.5 ),
        legend.position="top",
        axis.text=element_text(size=9),
        legend.text = element_text(size = 12),
        axis.title.y=element_blank() )

# estimating changepoints -----------------------------------------------------

#define a basis functions

basis0 <- spline.basis.eq(data = toy.data, 
                          n.knots=2, 
                          rangeval=c(-1,25),
                          norder = 4)
t0 <- Sys.time()
result <-  PELT.tpp(data=toy.data[,c(1,2,3)],
                    basis=basis0$basis, 
                    L=0,U=24,
                    iterations=100) 
Sys.time()-t0

result
changepoints 

cat('sample zise for each segment:', sample.seg)

# visulaize the intensity function estimated for each segment

source('R/estimate.intF.seg.R')
hat.intF <- estimate.intF.seg(data= toy.data, 
                              basis=basis0$basis, 
                              cpts=result[-1], 
                              L=0, U=24, 
                              iterations=200)
matplot(hat.intF$sgrid, hat.intF$intensity.mat, type = 'l', lwd=2, lty=2,
        ylim = c(0,4))
curve(lamb1, 0,24, col=1, ylim=c(0,5), lwd=2, xlab='s', ylab = '',
      main='Intensity functions used to simulate data', add = TRUE)
curve(lamb2, 0,24, col=2, lwd=2, add = TRUE)
curve(lamb3, 0,24, col=3, lwd=2, add = TRUE)

