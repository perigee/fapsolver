#!/bin/bash
cp solver test
args=("$@")
total=1
echo "#======================scen" $1
for (( c=1; c<=$total; c++ ))
do
   ./test celar_bk/FullRLFAP/GRAPH/graph$1/dom.txt celar_bk/FullRLFAP/GRAPH/graph$1/varm.txt celar_bk/FullRLFAP/GRAPH/graph$1/ctrm.txt
done
