#!/bin/bash
cp $1 chosen.solution.txt
plane=`echo $1 | cut -c 17- | cut -c1-3`
suffix=`tail -n 1 chosen.solution.txt | awk -v p="$plane" '{print "_p"p"_b"$1$2$3}'`
gfortran -o gccm_generate GCCM_GenerateCrystal1.f90
./gccm_generate
sed -n 1,14p data_cell.ref > data.0
sed -n 7,12p newcell.lmp > data.1
sed -n 21,59p data_cell.ref > data.2
sed -n 20,77p newcell.lmp > data.3
sed -n 118,473p data_cell.ref > data.4
cat data.* > data_cell.HMX$suffix
rm data.*
