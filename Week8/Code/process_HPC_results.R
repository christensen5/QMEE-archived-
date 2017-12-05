#====================================================================================================================================
#====================================================================================================================================
### This script will process and plot the output of the HPC simulation in model_functions_HPC.R
### Author: Alexander Kier Christensen
### QMEE CDT, High-Performance Computing week.
#====================================================================================================================================
#====================================================================================================================================

# Initial stuff
rm(list=ls())
graphics.off()
library(ggplot2)
library(grid)
#library(gridExtra)

#====================================================================================================================================

# Load functions

sum_vect <- function(x,y) {
  # This function takes two vectors, extends the shorter one to the same length as the longer by appending zeros, and then sums them.
  x_len = length(x)
  y_len = length(y)
  diff = x_len - y_len
  
  # Deal with empty vectors first
  # if (x_len == 0) {
  #   return(y)
  # }
  # if (y_len == 0) {
  #   return(x)
  # }
  
  if (diff > 0) {
    y = c(y, numeric(diff))
  } else if (diff < 0) {
    x = c(x, numeric(-diff))
  }
  
  return(x+y)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout)+1, ncol(layout), heights = unit(c(1,4,4), "null")) ))
    grid.text("Mean Species Abundances", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row+1,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#====================================================================================================================================

## Perform the analysis & processing

setwd("/home/alexander/Documents/QMEE/Week8/Results/hpcrun")
setwd("/home/alexander/Documents/QMEE/Week8/Results/hpcrun")
filenames = list.files(getwd(), pattern = "AKC_NTSresult_iter[:digit:]*")
numResults = length(filenames)

octave_mean_500 = numeric()
octave_mean_1000 = numeric()
octave_mean_2500 = numeric()
octave_mean_5000 = numeric()

# Community size 500
for (i in seq(1,numResults,4)) {
  loadstr = paste0(getwd() , "/" , filenames[i])
  load(loadstr)
  tmp_mean = numeric()
  burned_in_index = burn_in_generations / interval_oct  # index of abundances containing last burn-in generation
  for (j in seq(burned_in_index + 1, length(abundances))) {
    tmp_mean = sum_vect(tmp_mean, abundances[[j]])
  }
  tmp_mean = tmp_mean / (length(abundances) - burned_in_index)
  octave_mean_500 = sum_vect(octave_mean_500 , tmp_mean)
  numSims = floor(i/4) + 1  # save #sims so that at the end we can use it to divide sum to get mean
}
octave_mean_500 = octave_mean_500 / numSims
burn_in_500 = burn_in_generations

# Community size 1000
for (i in seq(2,numResults,4)) {
  loadstr = paste0(getwd() , "/" , filenames[i])
  load(loadstr)
  tmp_mean = numeric()
  burned_in_index = burn_in_generations / interval_oct  # index of abundances containing last burn-in generation
  for (j in seq(burned_in_index + 1, length(abundances))) {
    tmp_mean = sum_vect(tmp_mean, abundances[[j]])
  }
  tmp_mean = tmp_mean / (length(abundances) - burned_in_index)
  octave_mean_1000 = sum_vect(octave_mean_1000, tmp_mean)
  numSims = floor(i/4) + 1  # save #sims so that at the end we can use it to divide sum to get mean
}
octave_mean_1000 = octave_mean_1000 / numSims
burn_in_1000 = burn_in_generations

# Community size 2500
for (i in seq(3,numResults,4)) {
  loadstr = paste0(getwd() , "/" , filenames[i])
  load(loadstr)
  tmp_mean = numeric()
  burned_in_index = burn_in_generations / interval_oct  # index of abundances containing last burn-in generation
  for (j in seq(burned_in_index + 1, length(abundances))) {
    tmp_mean = sum_vect(tmp_mean, abundances[[j]])
  }
  tmp_mean = tmp_mean / (length(abundances) - burned_in_index)
  octave_mean_2500 = sum_vect(octave_mean_2500, tmp_mean)
  numSims = floor(i/4) + 1  # save #sims so that at the end we can use it to divide sum to get mean
}
octave_mean_2500 = octave_mean_2500 / numSims
burn_in_2500 = burn_in_generations

# Community size 5000
for (i in seq(4,numResults,4)) {
  loadstr = paste0(getwd(), "/" ,filenames[i])
  load(loadstr)
  tmp_mean = numeric()
  burned_in_index = burn_in_generations / interval_oct  # index of abundances containing last burn-in generation
  for (j in seq(burned_in_index + 1, length(abundances))) {
    tmp_mean = sum_vect(tmp_mean, abundances[[j]])
  }
  tmp_mean = tmp_mean / (length(abundances) - burned_in_index)
  octave_mean_5000 = sum_vect(octave_mean_5000, tmp_mean)
  numSims = floor(i/4) + 1  # save #sims so that at the end we can use it to divide sum to get mean
}
octave_mean_5000 = octave_mean_5000 / numSims
burn_in_5000 = burn_in_generations

# Clean up workspace
rm("abundances", "burn_in_generations", "community", "i", "interval_oct", "interval_rich", "j", "numGenerations", "output_file_name", "richness", "runtime_mins", "size", "speciation_rate", "wall_time")

# Plot
names = list(list(), list(), list(), list())
plots = list()
octave_means = list(octave_mean_500, octave_mean_1000, octave_mean_2500, octave_mean_5000)
titles = c("Initial population = 500", "Initial population = 1000", "Initial population = 2500", "Initial population = 5000")

for (k in c(1,2,3,4)) { 
  names[[k]] = c(names[[k]], parse(text="2^0 - 2^1"))
  for (l in seq(2,length(octave_means[[k]]))) {
    names[[k]] = c(names[[k]], parse(text=paste0("2^",toString(l-1)," - 2^",toString(l))))
  }
  plot_df = data.frame(bins = seq(1,length(octave_means[[k]])), counts = octave_means[[k]])
  plots[[k]] = ggplot(data = plot_df, aes(x = bins, y = counts)) +
    geom_bar(stat="identity", fill="navy") +
    labs(title = titles[k], x = "Abundance Octaves", y = "Counts") +
    scale_x_continuous(breaks = 1:length(octave_means[[k]]), labels = names[[k]]) +
    theme(plot.title = element_text(size=10, hjust=0.5))
}
multiplot(plotlist = plots, cols = 2)
