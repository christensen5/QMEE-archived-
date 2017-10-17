M <- matrix (runif(1000000),1000,1000)

SumAllElements <- function(M){
  Dimensions <- dim(M)
  Tot <- 0
  for (i in 1:Dimensions[1]){
    for (j in 1:Dimensions[2]){
      Tot <- Tot + M[i,j]
    }
  }
  return(Tot)
}

## This on Samraat's computer takes about 1 seconds.
print(system.time(SumAllElements(M)))

## While this takes about 0.01 seconds.
print(system.time(sum(M)))