import sys;
import math;
import string;

lon1 = [];
lat1 = [];
value1 = [];
lon2 = [];
lat2 = [];
value2 = [];

for l1 in open(sys.argv[1]):
	line = l1.rstrip();
	line1 = line.split();
	lon1.append(float(line1[0]));
	lat1.append(float(line1[1]));
	value1.append(float(line1[2]));	
i = 0;
for l2 in open(sys.argv[2]):
	line = l2.rstrip();
	line1 = line.split();
	lon2.append(float(line1[0]));
        lat2.append(float(line1[1]));
        value2.append(float(line1[2]));

#		lon2 = line1[0];
#		lat2 = line1[1];
#		value1 = line1[2];
#		if lon1 == float(line1[0]) and lat1 == float(line1[1]):
#			value1= value1*float(line1[2]);
	print lon2[i], lat2[i], value2[i]*value1[i];
	i = i+1;
