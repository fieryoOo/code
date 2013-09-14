#!/bin/csh
python << END
from math import *

h=40000.
dst1=2700.
dst2=3300.
B1=3500.
B2=4500.
u1=dst1*B1**2
u2=dst2*B2**2
for i in range(3):
	ff = open(str(i),'w')
	for c in range(3510,4500,10):
		p=1/float(c)
		w=(i*pi+atan((u2*sqrt(p**2-1/B2**2))/(u1*sqrt(1/B1**2-p**2))))/(h*sqrt(1/B1**2-p**2))
		ff.write("%.0f   %.3f\n" % (c,2*pi/w))
	ff.close()

END

gnuplot << END
set xlabel "Period (sec)"
set ylabel "Phase Velocity (m/sec)"

plot '0' using 2:1 with linespoints title 'fundamental mode'
set term postscript enhanced color
set output 'c_of_T_0.ps'
replot

plot '1' using 2:1 with linespoints title '1st mode'
set term postscript enhanced color
set output 'c_of_T_1.ps'
replot

plot '2' using 2:1 with linespoints title '2nd mode'
set term postscript enhanced color
set output 'c_of_T_2.ps'
replot
set output
set term x11
END
