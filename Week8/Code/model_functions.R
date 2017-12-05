# This script defines all the functions required for our population model.

# Initial stuff
rm(list=ls())
graphics.off()
#====================================================================================================================================
# 1
species_richness <- function(community) {
  
  # This function accepts a vector wherein each entry is the identity of the species in the relevant position in the environment.
  return(length(unique(community)))
}

#====================================================================================================================================
# 2
initialise_max <- function(size) {
  
  # This function creates an initial state for the simulation community with the maximal possible number of species for the community size.
  return(seq(size))
}

#====================================================================================================================================
# 3
initialise_min <- function(size){
  
  # This function creates an initial state for the simulation community with the minimum possible number of species for the community size.
  return(rep(1,size))
}

#====================================================================================================================================
# 4
choose_two <- function(x) {
  
  # Generate a random pair of DIFFERENT integers between 1 and x (inclusive).
  return(sample(seq(x),2))
}

#====================================================================================================================================
# 5
neutral_step <- function(community) {
  
  # This function performs a single step of a simple neutral model simulation, without speciation, on the community vector.
  fate = choose_two(length(community))  # Choose the indexes of the individual to die and to reproduce.
  community[fate[1]] <- community[fate[2]] # Kill the individual in entry fate[1] and replace with another individual from species in fate[2]
  return(community)
}  

#====================================================================================================================================
# 6
neutral_generation <- function(community) {
  
  # This function performs neutral_step() on a community for a single generation's worth of timesteps (= communitysize/2)
  generationTime = length(community)/2
  i = 1
  while (i < generationTime) {
    community = neutral_step(community)  # Perform the simulation step.
    i = i+1
  }
  return(community)
}

#====================================================================================================================================
# 7
neutral_time_series <- function(initial, duration) {
  
  # This function produces a neutral theory simulation and returns a time series of species richness in the system.
  time_series = numeric(duration+1)  # Initialise time series
  time_series[1] = species_richness(initial)
  
  Community_old = initial
  i=1
  while (i < duration+1) {
    Community_new = neutral_generation(Community_old)
    Community_old = Community_new
    time_series[i+1] = species_richness(Community_new)
    i = i+1
  }
  return(time_series)
  
}

#====================================================================================================================================
# 8 
question_8 <- function() {
  
  # This function produces the output for question 8, which requires us to plot a time series graph of neutral_time_series
  # with a maximally diverses initial community of 100 individuals over 200 generations.
  library(ggplot2)
  time_series = neutral_time_series(initialise_max(100),200)
  plotframe = data.frame(generations = seq(length(time_series)), richness = time_series)
  
  # Plot
  ggplot(plotframe, aes(generations)) +
    geom_line(aes(y=richness), colour="Navy") +
    labs(title = "Species Diversity in a Neutral Model Simulation (without speciation)", x = "Generation", y = "Species Richness", colour = "")
}


#====================================================================================================================================
# 9
neutral_step_speciation <- function(community, v) {
  
  # This performs a step of a neutral model with speciation (with speciation probability v).
  roll = runif(1)
  if (roll > v) {
    return(neutral_step(community))  # reduce to non-speciation case if speciation doesnt occur this step.
  } else {
    fate = sample(length(community),1)  # Choose the index of the individual to speciate.
    community[fate] = max(community) + 1  # Give the new species a new label.
    # It doesn't matter if this number has been used as a label before and died out, as long as its not currently in use (since we only care about species richness, not species history).
    return(community)
  }
  
}


#====================================================================================================================================
# 10
neutral_generation_speciation <- function(community,v) {
  
  # This function performs neutral_step_speciation() on a community for a single generation's worth of timesteps (= communitysize/2)
  generationTime = length(community)/2
  i = 1
  while (i < generationTime) {
    community = neutral_step_speciation(community,v)  # Perform the simulation step.
    i = i+1
  }
  return(community)
}


#====================================================================================================================================
# 11
neutral_time_series_speciation <- function(initial, duration, v) {
  
  # This function produces a neutral theory simulation and returns a time series of species richness in the system.
  time_series = numeric(duration+1)  # Initialise time series
  time_series[1] = species_richness(initial)
  
  Community_old = initial
  i=1
  while (i < duration+1) {
    Community_new = neutral_generation_speciation(Community_old,v)
    Community_old = Community_new
    time_series[i+1] = species_richness(Community_new)
    i = i+1
  }
  return(time_series)
  
}


#====================================================================================================================================
# 12 
question_12 <- function() {
  
  # This function produces the output for question 12, which requires us to plot a time series graph of neutral_time_series_speciation
  # with both a maximally diverse and a minimally diverse initial community of 100 individuals over 200 generations, with speciation rate 0.1
  library(ggplot2)
  v = 0.1
  initialSize = 100
  numGen = 200
  time_series_max = neutral_time_series_speciation(initialise_max(initialSize), numGen, v)
  time_series_min = neutral_time_series_speciation(initialise_min(initialSize), numGen, v)
  
  # Plot
  plotframe <- data.frame(generations = seq(length(time_series_max)), richness = cbind(as.matrix(time_series_max), as.matrix(time_series_min)))
  ggplot(plotframe, aes(generations)) +
    geom_line(aes(y=richness.1, colour = "Maximally diverse start")) + 
    geom_line(aes(y=richness.2, colour = "Minimally diverse start")) +
    labs(title = "Species Diversity in a Neutral Model Simulation (with speciation)", x = "Generation", y = "Species Richness", colour = "") +
    theme(plot.title = element_text(hjust=0.5))
}


#====================================================================================================================================
# 13
species_abundance <- function(community) {
   
  # This function computes the species abundances in the community.
  return(as.numeric(sort(table(community), decreasing=TRUE)))
  
}


#====================================================================================================================================
# 14
octaves <- function(abundances) {
  
  # This function bins the abundances of species into octave classes (based on powers of 2).
  log_abundances = log2(abundances)
  binCeiling = floor(log2(max(abundances)))+1  # Find the power of two relating to the final bin.
  binning = numeric(binCeiling)  # Initialise vector to hold result.
  i=1
  while (i<=binCeiling) {
    binning[i] = sum(as.numeric(i-1 == floor(log_abundances)))
    i = i+1
  }
  
  return(binning)
  
  
}


#====================================================================================================================================
# 15
sum_vect <- function(x,y) {
  
  # This function takes two vectors, extends the shorter one to the same length as the longer by appending zeros, and then sums them.
  x_len = length(x)
  y_len = length(y)
  diff = x_len - y_len
  
  if (diff > 0) {
    y = c(y, numeric(diff))
  } else if (diff < 0) {
    x = c(x, numeric(-diff))
  }
  
  return(x+y)

    
}


#====================================================================================================================================
# 16
question_16 <- function() {
  
  # This function produces the output for question 16, which requires us to run a full Neutral simulation with speciation, including
  # burn in period, periodically recording species abundances, and plotting the results.
  library(ggplot2)
  v = 0.1
  initialSize = 100
 
  burn_in = 200
  numGen = 2000
  recPeriod = 20
  numRecs = 1 + (numGen/20)  # 1 initial record + records every recPeriod generations for all of numGen
  
  # Burn-in run & then record abundances.
  burn_in_sim = neutral_generation_speciation(initialise_max(initialSize), v)
  for (i in 2:burn_in) {
    burn_in_sim = neutral_generation_speciation(burn_in_sim, v)
  }
  octave_mat = as.matrix(octaves(species_abundance(burn_in_sim)))  # typecast to matrix to avoid issues with dim() later.
  
  # Run the post-burn-in simulations.
  simulations = neutral_generation_speciation(burn_in_sim, v)
  for (i in 2:numGen) {
    simulations = neutral_generation_speciation(simulations, v)
    
    # Compute and record abundances every recPeriod generations
    if (i%%recPeriod == 0) {
      octave_vect = octaves(species_abundance(simulations))
      
      # Edit dimensions of octave_mat or octave_vect so we can cbind them.
      diff = length(octave_vect) - dim(octave_mat)[1]
      if (diff == 0) {
        octave_mat = cbind(octave_mat, octave_vect)
      } else if (diff > 0) {
        octave_mat = rbind(octave_mat, matrix(0,diff,dim(octave_mat)[2]))
        octave_mat = cbind(octave_mat, octave_vect)
      } else {
        octave_vect = c(octave_vect, numeric(-diff))
        octave_mat = cbind(octave_mat, octave_vect)
      }
    }
  }
  
  octave_means = rowMeans(octave_mat)
  
  # Plot
  title = paste("Mean Species Abundances after", numGen, "Generations with", burn_in, "Gen burn-in")
  barplot(octave_means,
          ylim = c(0,2*ceiling(max(octave_means)/2)),
          main = title, 
          col = "navy",
          xlab = "Abundance Octaves",
          ylab = "Species Counts",
          names.arg=c(parse(text=("2^0 - 2^1")), parse(text=("2^1 - 2^2")), parse(text=("2^2 - 2^3")),
                      parse(text=("2^3 - 2^4")), parse(text=("2^4 - 2^5")), parse(text=("2^5 - 2^6")))
  )
  #print(ceiling(max(octave_means)))
  
  
}
