#!/bin/bash
args=("$@")
total=5
echo "#======================graph" $1
for (( c=1; c<=$total; c++ ))
do
   ./solver celar/FullRLFAP/GRAPH/graph$1/dom.txt celar/FullRLFAP/GRAPH/graph$1/varm.txt celar/FullRLFAP/GRAPH/graph$1/ctrm.txt
done
