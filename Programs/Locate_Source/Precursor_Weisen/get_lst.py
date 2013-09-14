import sys;
import math;
import string;

for line in open(sys.argv[1]):
	line1 = line.rstrip();
	line2 = line1.split("_");
	per = int(line2[0]);
	sta1 = line2[2];
	tline = line2[3];
	line3 = tline.split(".");
	sta2 = line3[0];
	print line1,per,sta1,sta2;
