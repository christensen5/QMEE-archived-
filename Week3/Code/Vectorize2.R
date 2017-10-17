# Runs the stochastic (with gaussian fluctuations) Ricker Eqn .
rm(list=ls())

stochrick<-function(p0=runif(1000,.5,1.5),r=1.2,K=1,sigma=0.2,numyears=100)
{
  #initialize
  N<-matrix(NA,numyears,length(p0))
  N[1,]<-p0
  
  for (pop in 1:length(p0)) #loop through the populations
  {
    for (yr in 2:numyears) #for each pop, loop through the years
    {
      N[yr,pop]<-N[yr-1,pop]*exp(r*(1-N[yr-1,pop]/K) + rnorm(1,0,sigma) )
    }
  }
  return(N)
}

print("Stochastic Ricker takes:")
print(system.time(res1<-stochrick()))

# Now write another code called stochrickvect that vectorizes the above 
# to the extent possible, with improved performance: 
  
stochrickvect <- function(p0 = runif(1000,.5,1.5), r = 1.2, K = 1, sigma = 0.2, numyears = 100) {
  rands = replicate(length(p0),rnorm(numyears,1,sigma)) #pre-generate random fluctuations for each generation of each population.
  N <- matrix(NA, numyears, length(p0)) #pre-allocate output matrix N for efficiency.
  N[1,] <- p0 #store initial populations.
  # Loop through the years, vectorise-computing the new value for each population.
  for (yr in 2:numyears) {
    N[yr,] <- N[yr-1,]*exp(r*(1-(N[yr-1,]/K)) + rands[yr-1,]) #essentially apply the Ricker eqn to the last computed row of N and then store the result as the next row.
  }
  return(N)
}
print("Vectorized Stochastic Ricker takes:")
print(system.time(res2<-stochrickvect()))