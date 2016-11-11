# This script applies TreeHeight Function to file from command line.
# This function calculates heights of trees
# from the angle of elevation and the distance
# from the base using the trigonometric formula
# height = distance * tan(radians)
# Arguments:
#   degrees: The angle of elevation
#   distance: The distance from base
# Output: The height of the tree, same units as "distance"
TreeHeight <- function(degrees, distance)
  {radians <- degrees * pi / 180
  height <- distance * tan(radians)
  return (height)
}

TreeHeight(37, 40) # tests functions out on two numers 

#takes file from command line
arg1 <- commandArgs(TRUE)
Trees <- read.csv(paste("../Data/", arg1[1],".csv", sep=""))

Tree.Height.m <- TreeHeight(distance = Trees[,2], degrees = Trees[,3]) # applies 
# function to column data and saves it as vector
Trees[,"Tree.Height.m"] <- Tree.Height.m # add vector to data frame

#Output new data frame as a csv file.
write.csv(Trees, paste("../Results/",arg1[1],"_treeheigts.csv", sep = ""),row.names=FALSE) # append
