import sys;
import math;
import string;

am = [];
for line in open(sys.argv[1],"r"):
	line1 = line.rstrip();
	line2 = line1.split();
	am.append(float(line2[2]));
mam = max(am);
for line in open(sys.argv[1],"r"):
	line1 = line.rstrip();
        line2 = line1.split();
	print line2[0],line2[1],float(line2[2])/mam;
