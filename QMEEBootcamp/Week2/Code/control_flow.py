#!/user/bin/python

"""Some functions exemplifying the use of control statements"""
#docstrings are considered part of the running code (normal comments are
#stripped). Hence you can access your docstrings at run time.
__author__ = 'Samraat Pawar (email here)'
__version__ = '0.0.1'

# imports
import sys

# constants

# functions
def even_or_odd(x=0): #if not specified, x takes value 0
    """Find whether a number is even or odd."""
    if x % 2 == 0:
        return "%d is even!" %x
    return "%d is odd!" %x
	
def largest_divisor_five(x=120):
    """ Find which is the largest divisor of x among 2,3,4,5."""
    largest = 0
    if x % 5 == 0:
        largest = 5
    elif x % 4 == 0:
        largest = 4
    elif x % 3 == 0:
        largest = 3
    elif x % 2 == 0:
        largest = 2
    else: #When all other conditions are not met...
        return "No divisor found for %d!" % x
        
    return "The largest dicisor of %d is %d" % (x, largest)

def is_prime(x=70):
    """Find whether an integer is prime."""
    for i in range(2,x):
        if x % i == 0:
            print "%d is not a prime: %d is a divisor" % (x, i)
            
            return False
    print "%d is a prime!" % x
    return True
    
def find_all_primes(x=22):
    """Find all primes up to x"""
    allprimes = []
    for i in range (2, x + 1):
        if is_prime(i):
            allprimes.append(i)
        print "There are %d primes between 2 and %d" % (len(allprimes), x)
        
def main(argv):
    # sys.exit("don't want to do this right now!")
    print even_or_odd(22)
    print even_or_odd(33)
    print largest_divisor_five(120)
    print largest_divisor_five(121)
    print is_prime(60)
    print is_prime(59)
    print find_all_primes(100)
    return 0
    
if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)



def main(argv):
    print 'This is a boilerplate'
    return 0
	
if (__name__ == "__main__"): #makes sure the "main" function is called from commandline
    status = main(sys.argv)
    sys.exit(status)
