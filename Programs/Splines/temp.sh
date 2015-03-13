#!/bin/bash


for penalty in  100 200 300 400 500 700 1000 1400; do
#for penalty in  1500; do
	./BSpline1D TEST/SmoothPhase_pos.SAC.txt $penalty
	fout=test_spline_${penalty}.txt
	mv TEST/SmoothPhase_pos.SAC.txt_spline $fout
	echo $fout
done
