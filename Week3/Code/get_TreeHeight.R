# This function calculates heights of trees in a user-supplied csv file from the angle of
# elevation and the distance from the base using the trigonometric
# formula height = distance * tan(radians)
#
# INPUTS:
# csv file containing angle of elevation and distance from base for each desired tree
#
# ARGUMENTS:
# degrees
# distance
# The angle of elevation
# The distance from base
# 
# OUTPUT:
# csv file containing the input data as well as the height of each tree.

# Define TreeHeight to compute height of a tree given angle of elevation and base distance
TreeHeight <- function(degrees, distance){
  radians <- degrees * pi / 180
  height <- distance * tan(radians)
  print(paste("Tree height is:", height))
  
  return(height)
}



# Define the input data file.
args = commandArgs(trailingOnly=TRUE)
Input.Dir = paste("../Data/",args[1], sep="")
Input.File = read.csv(Input.Dir, header = TRUE)

# Run TreeHeight on each entry in the input file, saving the outputs as a column vector.
Tree.Height.m = c();
for (i in 1:dim(Input.File)[1]) {
  Tree.Height.m = rbind(Tree.Height.m, TreeHeight(Input.File[i,3],Input.File[i,2]))
}
Output.Frame = cbind(Input.File, Tree.Height.m) #Append the vector of tree heights to the input data.
Output.Dir = paste("../Results/",substr(args[1],1,nchar(args[1])-4),"_treeheights.csv", sep="") #Generate the output filename
write.csv(Output.Frame, Output.Dir) #Save the data frame to disk.
