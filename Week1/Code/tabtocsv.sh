#!/bin/bash
# Author: Alexander Kier Christensen a.christensen17@imperial.ac.uk
# Script: tabtocsv.sh
# Desc: simple shell script to replace tabs with commas and save output to csv file.
# Arguments: 1, a tab delimited file
# Date: Oct 2017

echo "Creating a comma delimited version of $1..."

FileName=$(echo $1 | head -c-5) #get the filename without the '.txt' on the end. 

cat $1 | tr -s "\t" "," >> $FileName.csv

echo "Done!"

exit
