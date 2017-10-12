#!/usr/bin/python

"""This script will exemplify the use of list comprehensions. Three
    tasks will be carried out, each using both a list comprehension and
    a loop. The tasks consist of extracting information from the 
    provided dataset 'birds'. """
__author__ = 'Alexander Kier Christensen'
__version__ = '0.0.1'

# Average UK Rainfall (mm) for 1910 by month
# http://www.metoffice.gov.uk/climate/uk/datasets
rainfall = (('JAN', 111.4),
            ('FEB', 126.1),
            ('MAR', 49.9),
            ('APR', 95.3),
            ('MAY', 71.8),
            ('JUN', 70.2),
            ('JUL', 97.1),
            ('AUG', 140.2),
            ('SEP', 27.0),
            ('OCT', 89.4),
            ('NOV', 128.4),
            ('DEC', 142.2),
            )

# (1) Use a list comprehension to create a list of month,rainfall tuples where
# the amount of rain was greater than 100 mm.

# (2) Use a list comprehension to create a list of just month names where the
# amount of rain was less than 50 mm. 

# (3) Now do (1) and (2) using conventional loops (you can choose to do 
# this before 1 and 2 !). 

# ANNOTATE WHAT EVERY BLOCK OR IF NECESSARY, LINE IS DOING! 

# ALSO, PLEASE INCLUDE A DOCSTRING AT THE BEGINNING OF THIS FILE THAT 
# SAYS WHAT THE SCRIPT DOES AND WHO THE AUTHOR IS

# ======================================================================================================================
# LIST COMPREHENSIONS
# 1
# Populate a list directly with tuples from rainfall whose second element exceeds 100.
lc1 = list(entry for entry in rainfall if entry[1] > 100)

# 2
# Populate a list directly with months from rainfall whose respective rainfall is less than 50.
lc2 = list(entry[0] for entry in rainfall if entry[1] < 50)

# ======================================================================================================================
# CONVENTIONAL LOOPS
# 1
# Instantiate an empty list, then run a for loop over each line in rainfall to find months with >100mm rain and save
# each such line in the new list.
loops1 = list()
for entry in rainfall:
    if entry[1] > 100:
        loops1.append(entry)

# 2
# Instantiate an empty list, then run a for loop over each line in rainfall to find months with <50mm rain, and save
# only the month names in the new list.
loops2 = list()
for entry in rainfall:
    if entry[1] < 50:
        loops2.append(entry[0])
