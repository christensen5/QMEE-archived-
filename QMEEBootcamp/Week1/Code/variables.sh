#!/bin/bash
# shows the use of variables

MyVar='some string'

echo 'the current value of the string is' $MyVar
echo 'Please enter a new string'
read MyVar
echo 'the current value of the string is' $MyVar

##Reading multiple values.
echo 'enter two numbers separated by space(s)'
read a b
echo 'you entered' $a 'and' $b '. Their sum is:'
mysum=`expr $a + $b`
echo $mysum
