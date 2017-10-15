#!/usr/bin/python

"""This script will exemplify the use of list comprehensions. Three
    tasks will be carried out, each using both a list comprehension and
    a loop. The tasks consist of extracting information from the 
    provided dataset 'birds'. """
__author__ = 'Alexander Kier Christensen'
__version__ = '0.0.1'

## Imports
import sys

birds = ( ('Passerculus sandwichensis','Savannah sparrow',18.7),
          ('Delichon urbica','House martin',19),
          ('Junco phaeonotus','Yellow-eyed junco',19.5),
          ('Junco hyemalis','Dark-eyed junco',19.6),
          ('Tachycineata bicolor','Tree swallow',20.2),
         )

#(1) Write three separate list comprehensions that create three different
# lists containing the latin names, common names and mean body masses for
# each species in birds, respectively. 

# (2) Now do the same using conventional loops (you can shoose to do this 
# before 1 !). 

# ANNOTATE WHAT EVERY BLOCK OR, IF NECESSARY, LINE IS DOING! 

# ALSO, PLEASE INCLUDE A DOCSTRING AT THE BEGINNING OF THIS FILE THAT 
# SAYS WHAT THE SCRIPT DOES AND WHO THE AUTHOR IS


## LIST COMPREHENSIONS
#1
# Create and populate a list directly with latin names.
latin_names_lc = list(tuple[0] for tuple in birds)

#2
# Create and populate a list directly with common names.
common_names_lc = list(tuple[1] for tuple in birds)

#3
# Create and populate a list directly with mean masses.
mean_mass_lc = list(tuple[2] for tuple in birds)

## CONVENTIONAL LOOPS
#1
#Initialise empty list, then populate with latin names with a loop.
latin_names_loop = list()
for tuple in birds:
    latin_names_loop.append(tuple[0])
    
#2
#Initialise empty list, then populate with common names with a loop.
common_names_loop = list()
for tuple in birds:
    common_names_loop.append(tuple[1])
    
#3
#Initialise empty list, then populate with mean masses with a loop.
mean_mass_loop = list()
for tuple in birds:
    mean_mass_loop.append(tuple[2])    










