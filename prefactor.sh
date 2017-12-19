#!/bin/bash
cp $1 chosen.solution.txt
cat chosen.solution.txt | tail -n 1 > chosen.mat
./gccm_generate | grep -A3 Net | tail -n 3 > rotation.mat
matlab -nosplash -nodisplay -r "munday_hmx1;quit" 
