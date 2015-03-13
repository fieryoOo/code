#!/bin/bash

SAC_dump COR_J47A_J37A.SAC.pos.ph COR_J47A_J37A.SAC.pos.ph.txt
SAC_dump COR_J47A_J37A.SAC.neg.ph COR_J47A_J37A.SAC.neg.ph.txt
psout=COR_J47A_J37A.SAC.compare_phase.txt.ps
/usr/bin/gnuplot <<- END
	set term postscript enhanced color
	set out '$psout'
	set multiplot layout 2,1
	set xrange[0.:0.2]
	set origin -0.05,0.5
	set size 1.05,0.4
	plot 'Phase_pos.SAC_rms' w l, 'Phase_neg.SAC_rms' lt 1 lc 3 w l
	set origin 0.,0.1
	set size 1.,0.4
	plot 'COR_J47A_J37A.SAC.pos.ph.txt' w l, 'COR_J47A_J37A.SAC.neg.ph.txt' lt 1 lc 3 w l
	unset multiplot
	set term x11
	set out
	q
END

echo $psout
