#!/bin/bash
cat $1 > $3
cat $2 >> $3
echo "Merged file is"
cat $3
