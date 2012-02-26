#!/bin/bash
args=("$@")
total=10
echo "#======================myciel3" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/myciel3.col 3 >> resultCOL.dat
done
echo "#======================le450_5a" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/le450_5a.col 4 >> resultCOL.dat
done
echo "#======================le450_5b" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/le450_5b.col 4 >> resultCOL.dat
done
echo "#======================le450_5c" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/le450_5c.col 4 >> resultCOL.dat
done
echo "#======================le450_5d" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/le450_5d.col 4 >> resultCOL.dat
done
echo "#======================anna" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/anna.col 10 >> resultCOL.dat
done
echo "#======================david" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/david.col 10 >> resultCOL.dat
done
echo "#======================huck" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/huck.col 10 >> resultCOL.dat
done
echo "#======================jean" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/jean.col 9 >> resultCOL.dat
done
echo "#======================DSJC125.1" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/DSJC125.1.col 4 >> resultCOL.dat
done
echo "#======================DSJR500.1" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/DSJR500.1.col 11 >> resultCOL.dat
done
echo "#======================queen5_5" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/queen5_5.col 4 >> resultCOL.dat
done
echo "#======================queen7_7" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/queen7_7.col 6 >> resultCOL.dat
done
echo "#======================queen10_10" >> resultCOL.dat
for (( c=1; c<=$total; c++ ))
do
   ./dimacs benchmark/dimacsGCP/easy/queen10_10.col 9 >> resultCOL.dat
done
