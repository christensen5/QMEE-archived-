    #!/user/bin/python

"""Some functions exemplifying the use of control statements"""
#docstrings are considered part of the running code (normal comments are
#stripped). Hence you can access your docstrings at run time.
__author__ = 'Samraat Pawar'
__version__ = '0.0.1'

# imports
import sys
import doctest

# constants

# functions
def even_or_odd(x=0): #if not specified, x takes value 0
    """Find whether a number is even or odd.
    
    >>> even_or_odd(10)
    '10 is Even!'
    
    >>> even_or_odd(5)
    '5 is Odd!'
    
    whenever a float is provided, then the closest integer is used:
    >>> even_or_odd(3.2)
    '3 is Odd!'
    
    in case of negative numbers, the positive is taken:
    >>> even_or_odd(-2)
    '-2 is Even!'
    
    """
    if x % 2 == 0:
        return "%d is Even!" % x
    return "%d is Odd!" % x
	
#def largest_divisor_five(x=120):
    #""" Find which is the largest divisor of x among 2,3,4,5."""
    #largest = 0
    #if x % 5 == 0:
        #largest = 5
    #elif x % 4 == 0:
        #largest = 4
    #elif x % 3 == 0:
        #largest = 3
    #elif x % 2 == 0:
        #largest = 2
    #else: #When all other conditions are not met...
        #return "No divisor found for %d!" % x
        
    #return "The largest divisor of %d is %d" % (x, largest)

#def is_prime(x=70):
    #"""Find whether an integer is prime."""
    #for i in range(2,x):
        #if x % i == 0:
            #print "%d is not a prime: %d is a divisor" % (x, i)
            
            #return False
    #print "%d is a prime!" % x
    #return True
    
#def find_all_primes(x=22):
    #"""Find all primes up to x"""
    #allprimes = []
    #for i in range (2, x + 1):
        #if is_prime(i):
            #allprimes.append(i)
        #print "There are %d primes between 2 and %d" % (len(allprimes), x)
        
#def main(argv):
    ## sys.exit("don't want to do this right now!")
    #print even_or_odd(22)
    #print even_or_odd(33)
    #print largest_divisor_five(120)
    #print largest_divisor_five(121)
    #print is_prime(60)
    #print is_prime(59)
    #print find_all_primes(100)
    #return 0
    
#if (__name__ == "__main__"):
    #main(sys.argv)


doctest.testmod() # to run with embedded tests!
