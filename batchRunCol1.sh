#!/bin/bash
args=("$@")
total=10
echo "#======================queen5_5" 
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/queen5_5.col 4 
done

