#! /usr/bin/env python

import re
import math


f=open("events.raw.txt",'r')
rawe = f.read()
rawe = rawe.strip()
rawe = rawe.replace("ALERT,","ALERT")
edlist = rawe.split('\n')

for i in range(len(edlist)):
	edlist[i] = re.split('\n| \| | \||\| |\||,',edlist[i])
	edlist[i][0:2] = edlist[i][1].split(' ')
	edlist[i] = edlist[i][0].replace('/','-')
edlist.sort()

dlist = edlist[:]
rcount = 0
for i in range(len(dlist)-1):
	if dlist[i] == dlist[i+1]:
		dlist[i] = 'remove'
		rcount += 1
edlist = dlist[:]
if rcount != 0:
	for i in range(rcount):
		edlist.remove('remove')

edstr = ""
for i in range (len(edlist)):
	edstr = edstr + str(edlist[i])+"\n"

outfile = open("daylist.txt",'w')
outfile.write(edstr)
