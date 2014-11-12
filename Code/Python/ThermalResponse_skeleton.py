#!/usr/bin/env python

# This code performs non-linear least squares fitting of different
# unimodal functions to experimental thermal response curves.
# Scroll down to the section called "MAIN CODE" first.

from math import log, exp, pi
from lmfit import minimize, Parameters
import csv
import numpy
import re
import sys

#############################
# F  U  N  C  T  I  O  N  S #
#############################

def GauGo_eq(temps, E, E_D, T_pk, theta):
	"""Gaussian-Gompertz model, used for calculating trait values at 
	a given temperature"""
	
	global max_trait

	function = max_trait * exp((-E * (temps - T_pk) * (temps \
	 - T_pk) - exp((E_D * (temps - T_pk) - theta))))
	# Note the log, and float64 (for numerical precision):                     
	return numpy.array(map(log, function), dtype=numpy.float64) 

def schoolf_eq(temps, B0, E, E_D, T_pk):
	"""Schoolfield model, used for calculating trait values at 
	a given temperature"""
	global K
	
	function = B0 * exp(-E * ((1/(K*temps)) - (1/(K*283.15)))) \
	/ (1 + (E/(E_D - E)) * exp(E_D / K * (1 / T_pk - 1 / temps)))
	return numpy.array(map(log, function), dtype=numpy.float64)

def GauGo(params, temp, data):
	"""Gaussian - Gompertz model to be called by the gaugo_model() 
	function"""

	# return ??

def schoolf(params, temp, data):
	"""Schoolfield model, to be called by schoolfield_model()"""

	# return ??

def gaugo_model(dataset, E_start, T_pk_start, E_D_start):
	"""NLLS fitting to the Gaussian-Gompertz model; this 
	function will contain the lmfit.minimize calls to the GauGo() 
	function. It will return the results of the fitting"""

	# return ??

def schoolfield_model(dataset, B0_start, E_start, T_pk_start, E_D_start):
	"""NLLS fitting to the Schoolfield model; this function will 
	contain the lmfit.minimize calls to the schoolf() function. This is 
	where you can constrain the parameters."""

# Prepare the parameters and their bounds:
	# params = Parameters()
	# params.add('theta', value=7)
	# params.add('E', value=E_start, min=??, max=??)
	# params.add('E_D', value=E_D_start, min=??, max=??)
	# params.add('T_pk', value=T_pk_start - 273.15, min=??, max=??)
	
	# return ??

def AICrss(n, k, rss):
	"""Calculate the Akaike Information Criterion value, using:

	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""

	return n * log((2 * pi) / n) + n + 2 + n * log(rss) + 2 * k

def BICrss(n, k, rss):
	"""Calculate the Bayesian Information Criterion value, using:
	
	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""

	return n + n * log(2 * pi) + n * log(rss / n) + (log(n)) * (k + 1)
	
############################
# M  A  I  N    C  O  D  E #
############################

def main(argv):
	 
	# Define the Boltzmann constant (units of eV * K^-1).
	global K
	K = 8.617 * 10 ** (-5)
	
	# Raise an error if an input dataset wasn't provided.
	if len(sys.argv) != 2:
		sys.exit("USAGE: " + sys.argv[0] + " input_dataset")
	
	# Read the dataset file into a csv object.
	with open(sys.argv[1]) as csvfile:
		csv_dataset = csv.reader(csvfile)

		# Store the data in a list.
		original_dataset = [row for row in csv_dataset]
		
	counter = 0 #initialize counter

	# Create an output CSV file with parameters estimates for all the 
	# models (you will write to this csv, one thermal curve at a time):
	results = open("../Results/results.csv", 'w')
	results_csv = csv.writer(results, delimiter="\t")
	results_csv.writerow(
				['Species_orig', 'Species_stand', 'Reference', 'Trait', 
				'Latitude', 'Longitude', 'Temp_Vals', 'Trait_Vals', 'B0', 
				'B0_stderr', 'E', 'E_stderr', 'T_pk', 'T_pk_stderr', 'E_D', 
				'E_D_stderr', 'Theta', 'Theta_stderr', 'R_Squared', 'Formula', 
				'Model_name', 'DataPoints_rise', 'DataPoints_fall', 
				'AIC_GauGo', 'AIC_Schoolf', 'BIC_GauGo', 'BIC_Schoolf', 
				'SelectedModel'])

	import ipdb; ipdb.set_trace()	# Move this down as code develops 

	# Extract unique IDs for thermal responses in the dataset:  
			 
	# Now fit models to each thermal response by iterating over unique IDs : 
	for entry in ids:

		# Get the subset of the data for this ID, in appropriate format:
		# tmp_data = ??
		
		# Extract the starting values for the NLLS fitting:
		# (max_trait, T_pk_start, B0_start, E_start, E_D_start) = ??
	
		# Fit the models:
		GauGo_fit = gaugo_model(tmp_data, E_start, T_pk_start, E_D_start)
	   
		Schoolf_fit = schoolfield_model(tmp_data, B0_start, E_start,
											T_pk_start, E_D_start)
   
				# Calculate AIC and BIC using the NLLS fitting results in 
				# GauGo_fit and Schoolf_fit:
					
		# Decide on the best model according to AIC and BIC and add to 
		# the results csv:

		# Write the results for this thermal response to the CSV:

		# results_csv.writerow(??) 
						
		counter += 1 #update counter
		print("Finished thermal response number: " + str(counter))

	return 0

if __name__ == "__main__":
	main(sys.argv)
