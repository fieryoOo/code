import math;
import string;
import sys;

for l1 in open(sys.argv[1],"r"):
	line1 = l1.rstrip();
	ll1 = line1.split();
	for l2 in open(sys.argv[2],"r"):
		line1 = l2.rstrip();
		ll2 = line1.split();
		str1 = "%5d%5d %-6s %-6s %10.5f %10.5f %10.5f %10.5f" % (int(ll1[0]), int(ll2[0]), ll1[1], ll2[1], float(ll1[3]), float(ll1[2]), float(ll2[3]), float(ll2[2]));
		print str1;

# 2i5,1x,a16,1x,a6,4(f10.3)
