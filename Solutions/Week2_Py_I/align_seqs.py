"""
	This script reads two sequences from a csv file, aligns them
	and returns the best alignment along with alignment score
"""

import csv
import pickle
# These are the two sequences to match
#~seq2 = "ATCGCCGGATTACGGG"
#~seq1 = "CAATTCGGAT"

def align_seqs(InputFile):

    # Open csv file for reading.
    f = open(InputFile,'rb')
    seqs = csv.reader(f)
    
    # Open output file
    g = open('../Results/AlignedSeqs.p','wb')
    
    # Read and store the sequences
    for sequence in seqs:
		seq1 = sequence[0]
		seq2 = sequence[1]
    
    f.close()
    
    # assign the longest sequence 
    # to s1, and the shortest to s2
    # l1 is the length of the longest,
    # l2 that of the shortest
    l1 = len(seq1)
    l2 = len(seq2)
    if l1 >= l2:
        s1 = seq1
        s2 = seq2
    else:
        s1 = seq2
        s2 = seq1
	l1, l2 = l2, l1 # swap the two lengths
	
    # now try to find the best match (highest score)
    my_best_align = None
    my_best_score = -1

    for i in range(l1):
        z = calculate_score(s1, s2, l1, l2, i)
        if z > my_best_score:
            my_best_align = "." * i + s2
            my_best_score = z

    print my_best_align
    print s1
    print "Best score:", my_best_score

    # Create a dictionary to hold the results
    Alignment = {}

    # Now populate it with results
    Alignment["Sequence"] = s1
    Alignment["Best alignment"] = my_best_align
    Alignment["score"] = my_best_score

    # Dump the results
    pickle.dump(Alignment, g)
    g.close()

    return 0

# function that computes a score
# by returning the number of matches 
# starting from arbitrary startpoint
def calculate_score(s1, s2, l1, l2, startpoint):
    """
    This function computes the score of the alignment.
    """
	
    # startpoint is the point at which we want to start
    matched = "" # contains string for alignement
    score = 0
    for i in range(l2):
        if (i + startpoint) < l1:
            # if it's matching the character
            if s1[i + startpoint] == s2[i]:
                matched = matched + "*"
                score = score + 1
            else:
                matched = matched + "-"
    return score
