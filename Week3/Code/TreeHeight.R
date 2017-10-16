# This function calculates heights of trees in a given csv file from the angle of
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
# The height of the tree, same units as "distance"

# Define TreeHeight to compute height of a tree given angle of elevation and base distance
TreeHeight <- function(degrees, distance){
  radians <- degrees * pi / 180
  height <- distance * tan(radians)
  print(paste("Tree height is:", height))
  
  return(height)
}

# TreeHeight(37, 40)

# Define the input data file.
Input.File = read.csv("../Data/trees.csv", header = TRUE)

# Run TreeHeight on each entry in the input file, saving the outputs as a column vector.
Tree.Height.m = c();
for (i in 1:dim(Input.File)[1]) {
  Tree.Height.m = rbind(Tree.Height.m, TreeHeight(Input.File[i,3],Input.File[i,2]))
}
Output.Frame = cbind(Input.File, Tree.Height.m) #Appent the vector of tree heights to the input data.
write.csv(Output.Frame, "../Results/TreeHts.csv") #Save the data frame to disk.
