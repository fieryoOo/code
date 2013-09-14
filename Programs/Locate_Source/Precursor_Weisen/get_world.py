import math;
import string;
import sys;

if len(sys.argv) !=3:
	print "usage: input regional file, average v";
	sys.exit();

for j in range(-90,91):
	for i in range(0,360):
		v = float(sys.argv[2]);
		if j<28 or j >55 or i<100 or i>152:
			print i,j,v;
			continue;
		for line in open(sys.argv[1]) :
			line1 = line.rstrip();
			line2 = line.split();
			lon = float(line2[0]);
			lat = float(line2[1]);
			vel = float(line2[2]);
			if math.fabs(lon-i) < 0.01 and math.fabs(lat-j) < 0.01:
				v = vel;
#				print v;
		print i,j,v;
