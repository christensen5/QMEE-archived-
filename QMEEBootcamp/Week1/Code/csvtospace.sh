#!/bin/bash
# Author: Alexander Kier Christensen a.christensen17@imperial.ac.uk
# Script: csvtospace.sh
# Desc: simple shell script to replace commas with spaces and save output to a different file.
# Arguments: 1, a comma delimited file
# Date: Oct 2017

echo "Creating a comma delimited version of $1..."

FileName=$(echo $1 | head -c-5) #get the filename without the '.csv' on the end. 

cat $1 | tr -s "," " " > $FileName.txt

echo -e "Done!"
