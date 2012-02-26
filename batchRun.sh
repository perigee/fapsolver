#!/bin/bash
args=("$@")
total=50
echo "#======================scen 01" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen01/DOMt.TXT celar/CELAR/scen01/VARt.TXT celar/CELAR/scen01/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 02" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen02/DOMt.TXT celar/CELAR/scen02/VARt.TXT celar/CELAR/scen02/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 03" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen03/DOMt.TXT celar/CELAR/scen03/VARt.TXT celar/CELAR/scen03/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 04" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen04/DOMt.TXT celar/CELAR/scen04/VARt.TXT celar/CELAR/scen04/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 05" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen05/DOMt.TXT celar/CELAR/scen05/VARt.TXT celar/CELAR/scen05/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 06" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen06/DOMt.TXT celar/CELAR/scen06/VARt.TXT celar/CELAR/scen06/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 07" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen07/DOMt.TXT celar/CELAR/scen07/VARt.TXT celar/CELAR/scen07/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 08" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen08/DOMt.TXT celar/CELAR/scen08/VARt.TXT celar/CELAR/scen08/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 09" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen09/DOMt.TXT celar/CELAR/scen09/VARt.TXT celar/CELAR/scen09/CTRt.TXT s >> resultCELAR_marked.dat
done
echo "#======================scen 11" >> resultCELAR_marked.dat
for (( c=1; c<=$total; c++ ))
do
   ./testCELAR celar/CELAR/scen11/DOMt.TXT celar/CELAR/scen11/VARt.TXT celar/CELAR/scen11/CTRt.TXT s >> resultCELAR_marked.dat
done
