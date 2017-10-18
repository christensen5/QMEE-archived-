## This function will accept processed data on global species distribution and plot it on a global map using the map() R package.

## IMPORTS
load("../Data/GPDDFiltered.RData") # Load default dataset
library(maps) # load the maps package

## DEFINE WORKHORSE FUNCTION
GPDD_Plotter <- function(gpdd_data=gpdd) {
  graphics.off() #clear previous plots
  
  
  png('../Results/GPDD_plot.png',width = 1440, height = 1440, units = "px") #Open plot device
  
  fig = map('world', resolution=0) #plot base world map
  points(gpdd_data[,3], gpdd_data[,2], pch = 20, col='red') #For each species in dataset, plot its location over the base map, as a solid red circle.
  
  dev.off() #close plot device
}

GPDD_Plotter() #for testing purposes, call the function with default argument when the file is sourced.
  
## Biases would stem principally from there being no data for the vast majority of the planet, besides the UK, Europe and NA. UK is particularly overrepresented.