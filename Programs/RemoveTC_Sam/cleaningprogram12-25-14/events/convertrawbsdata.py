#! /usr/bin/env python

import re
import math

f=open("rawbsdata.txt",'r')
rdata = f.read()
rdata=rdata.strip()
rdlist=rdata.split('\n')

for i in range(len(rdlist)):
	rdlist[i] = rdlist[i].split(" ")
	rdlist[i] = [rdlist[i][0],rdlist[i][15]]

removelist = []
for i in range(len(rdlist)):
	if rdlist[i][1] == 'BXH':
		removelist.append([rdlist[i][0],'BDH'])
	if rdlist[i][1] == 'BX1':
		removelist.append([rdlist[i][0],'BH1'])
	if rdlist[i][1] == 'BX2':
		removelist.append([rdlist[i][0],'BH2'])
	if rdlist[i][1] == 'BXZ':
		removelist.append([rdlist[i][0],'BHZ'])

for i in range(len(removelist)):
	rdlist.remove(removelist[i])
	
rdstr = ""
for i in range (len(rdlist)):
	rdstr = rdstr + str(rdlist[i][0]) + "/" + str(rdlist[i][1]) + "\n"

outfile = open("bslist.txt",'w')
outfile.write(rdstr)
