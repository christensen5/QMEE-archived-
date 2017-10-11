#!/user/bin/python

"""Some statements printing 'hello', and some functions demonstrating
    loops and recursion"""
__author__ = 'Samraat Pawar, edited by Alexander Christensen'
__version__ = '0.0.2'

## Imports
import sys

## How many times will 'hello' be printed?
## 1)
for i in range(3, 17):
    print 'hello'
    
# 2)
for j in range(12):
    if j % 3 == 0:
        print 'hello'
    
# 3)
for j in range(15):
    if j % 5 == 3:
        print 'hello'
    elif j % 4 == 3:
        print 'hello'
        
# 4)
z = 0
while z != 15:
    print 'hello'
    z = z + 3
    
# 5)
z = 12
while z < 100:
    if z == 31:
        for k in range(7):
            print 'hello'
    elif z == 18:
        print 'hello'
    z = z + 1
    
## What does fooXX do?
def foo1(x=4): 
    """Find the square root of x"""
    print "The square root of x is %d" % (x ** 0.5)
    return x ** 0.5
    
def foo2(x=1, y=2):
    """Find the greater of x and y. If equal return y"""
    if x > y:
        print "%d is greater" % x
        return x
    print "%d is greater" % y
    return y
    
def foo3(x=1, y=2, z=3):
    """Swap x and y if x>y, THEN swap y and z if y>z"""
    if x > y:
        print "Swapping %d and %d" % (x, y)
        tmp = y
        y = x
        x = tmp
    if y > z:
        print "Swapping %d and %d" % (y, z)
        tmp = z
        z = y
        y = tmp
    print "The result is (%d, %d, %d)" % (x, y, z)
    return [x, y, z]
    
def foo4(x=10):
    """Find x! (using a for-loop)"""
    result = 1
    for i in range (1, x + 1):
        result = result * i
    print "%d! = %d" % (x, result)
    return result
    
# This is a recursive function, meaning that the function calls itself.
def foo5(x=10):
    """Find x! (using recursion)"""
    if x == 1:
        return 1
    result = x * foo5(x - 1)
    # Can't print result this time - would print at every recursive call
    return result
    
    
def main(argv):
    # sys.exit("don't want to do this right now!")
    foo1()
    foo2()
    foo3()
    foo4()
    print "10! = %d" % foo5()
    return 0
    
if (__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)
