## run a simulation that involves sampling from a population

x <- rnorm(50)
doit <- function(x){
  x <- sample(x,replace=TRUE) 
  if(length(unique(x)) > 30) { #only take mean if sample sufficient
    print(paste("Mean fo this sample was:", as.character(mean(x))))
  } else {
    stop("Couldn't calculate mean: too few unique points!")
  }
}

#run 100 iterations using vectorization
result <- lapply(1:100, function(i) try(doit(x),FALSE))