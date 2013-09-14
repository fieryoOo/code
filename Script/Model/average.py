import sys;
import math;
import string;

for line in open(sys.argv[1],"r"):
	line1 = line.rstrip();
	line2 = line1.split();
	lon = float(line2[0]);
	lat = float(line2[1]);
	vel1 = float(line2[2]);
	for l2 in open(sys.argv[2],"r"):
		line1 = l2.rstrip();
		line2 = line1.split();
		if float(line2[0]) == lon and float(line2[1]) == lat:
			print lon, lat, (vel1 + float(line2[2]))/2
