MyFunction <- function(Arg1, Arg2){
  
  # Statements involving Arg1, Arg2
  print(paste("Argument", as.character(Arg1), "is a", class(Arg1))) # print Arg1's type
  print(paste("Argument", as.character(Arg2), "is a", class(Arg2))) # likewise Arg2
  
  return (c(Arg1, Arg2)) #this is optional but very useful
}

MyFunction(1,2) #test the function
MyFunction("Riki","Tiki") #a different test