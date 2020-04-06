#!/bin/bash

for ((i=7;i<=100;i++)); do
    #root -l -b -q makeMatrix_isPVGood.C\(0\,0\,0\,$i\)
    root -l -b -q 'makeMatrix_isPVGood.C(0,0,0,'$i')'
    root -l -b -q 'makeMatrix_isPVGood.C(0,0,1,'$i')'
    root -l -b -q 'makeMatrix_isPVGood.C(0,0,2,'$i')'
    root -l -b -q 'makeMatrix_isPVGood.C(0,0,3,'$i')'

done
