##################################
## HMMS FOR POPULATION DYNAMICS ##
##        ROBERT PATON          ##
##################################

rm(list=ls(all=T))

# The following packages are required to run this code.

require(MASS);require(ggplot2);require(gridExtra);require(grid);require(boot);require(lattice);

# NB: this code requires several function included in the zip file. To run, the working directory must be 
# changed to your own, and include all the stated files. See the files for details in functionality. 

setwd("/Users/Rob/Google Drive/Masters/Hidden Markov Models/")
source("estimation_vole_yearly.R") 
#source("simulate_tpm_vole_yearly.R")
source("viterbi_vole_yearly.R")
source("extrapolation_vole_yearly.R")


## DATA INPUT ##

dat<-read.csv("yearly_vole_data.csv", header=T)
obs<-dat$Count

## MODEL FITTING ##

# Fit 'repeats' number of models with stochastic starting values, then select the model with the smallest
# log likelihood

repeats<-100
models<-list()
the.liks<-c()

for(each.model in 1:repeats)
{
  # The number of states
  N <- 3
  
  # Random mu's
  mu0<-c(runif(1, 0.3, 0.5), 1, runif(1, 1.5, 2.5))
  
  # Random params for the left hand column
  beta10<-matrix(rnorm(9,c(-1,-1,-1,rep(0, 6)),
                       c(0.5,0.5,0.5,rep(0.1, 6))),byrow=F, nrow=3)
  
  # Random params for the right hand column
  beta30<-matrix(rnorm(9,c(1,1,1,rep(0, 6)),
                       c(0.5,0.5,0.5,rep(0.1, 6))),byrow=F, nrow=3)
  
  # Initial distribution
  delta<-sample(c(1,0,0), prob=c(rep(1/3,3)), size=3)
  
  # Estimator
  mod<- mle(obs[-1], 
            mu0=mu0, 
            beta10=beta10,
            beta30=beta30,
            delta=delta,
            N=N, hessian=T, print.level=0, t=t, two.lag=obs[-length(obs)])
  
  the.liks<-c(the.liks, mod$mllk)
  models[[each.model]]<-mod
  print(each.model)
  
}

mod<-models[[which.min(the.liks)]]
plot(the.liks)

source("viterbi_vole_yearly.R")

states<-viterbi.algorithm(obs=obs[-1], N=N, beta1=mod$beta1, beta3=mod$beta3, 
                          mu=mod$mu, delta=mod$delta, t=t, two.lag=obs[-length(obs)])[[1]]

source("extrapolation_vole_yearly.R")

n<-25
repeats<-10000
store<-matrix(NA, ncol=n, nrow=repeats)

for(i in 1:repeats){
  store[i,]<-extrapolate(n=n, N=N, 
                  x.init=tail(obs, 1), s.init=tail(states, 1),two.year=tail(obs,3),
                  random_params=c(mod$mu,mod$beta1, mod$beta3))[[1]]
}

means<-apply(store, 2, mean)
lower<-apply(store, 2, sort)[round(repeats*0.025),]
upper<-apply(store, 2, sort)[round(repeats*0.975),]

plot(x=1:(length(obs)+n), c(obs, rep(NA, n)), type='l', xlab='Years', ylab='N', ylim=c(0,max(upper)))
lines(x=(length(obs)):(length(obs)-1+n), means, lty=2)
lines(x=(length(obs)):(length(obs)-1+n), upper, col="grey", lty=2)
lines(x=(length(obs)):(length(obs)-1+n), lower, col="grey", lty=2)

