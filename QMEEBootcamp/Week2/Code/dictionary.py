#!/usr/bin/python

"""This script will take taxa and populate a dictionary with it,
    mapping order names to sets of taxa."""
__author__ = 'Alexander Kier Christensen'
__version__ = '0.0.1'

## Imports
import sys

taxa = [ ('Myotis lucifugus','Chiroptera'),
         ('Gerbillus henleyi','Rodentia',),
         ('Peromyscus crinitus', 'Rodentia'),
         ('Mus domesticus', 'Rodentia'),
         ('Cleithrionomys rutilus', 'Rodentia'),
         ('Microgale dobsoni', 'Afrosoricida'),
         ('Microgale talazaci', 'Afrosoricida'),
         ('Lyacon pictus', 'Carnivora'),
         ('Arctocephalus gazella', 'Carnivora'),
         ('Canis lupus', 'Carnivora'),
        ]

# Write a short python script to populate a dictionary called taxa_dic 
# derived from taxa so that it maps order names to sets of taxa. 
# E.g. 'Chiroptera' : set(['Myotis lucifugus']) etc. 

# ANNOTATE WHAT EVERY BLOCK OR IF NECESSARY, LINE IS DOING! 

# ALSO, PLEASE INCLUDE A DOCSTRING AT THE BEGINNING OF THIS FILE THAT 
# SAYS WHAT THE SCRIPT DOES AND WHO THE AUTHOR IS

# Write your script here:


def ext_order_set(any_taxa=taxa): #default to taxa if no input given
    """Extract the order names and create a set from them."""
    order_set = set()
    for entry in any_taxa:
        order_set.add(entry[1]) 
    print "The orders in the dictionary are " + str(list(order_set))
    
    return order_set
    
    
def main(argv):
    #Initialise taxa_dic with keys defined by ext_order_set and a 
    #DIFFERENT empty set as the value assigned to each key.
    taxa_dic = {k: set() for k in ext_order_set(taxa)}
    
    #Iterate over taxa and populate the taxa_dic keys with the relevant
    #species.
    for entry in taxa:
        taxa_dic[entry[1]].add(entry[0])
    
    print "\nThe dictionary is " + str(taxa_dic)
    return taxa_dic
     
        
if (__name__ == "__main__"):
    main(sys.argv)







