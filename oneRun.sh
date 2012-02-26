#!/bin/bash
cp solver test
args=("$@")
total=1
echo "#======================scen" $1
for (( c=1; c<=$total; c++ ))
do
   ./test celar/CELAR/scen$1/DOMt.TXT celar/CELAR/scen$1/VARt.TXT celar/CELAR/scen$1/CTRt.TXT s
done
