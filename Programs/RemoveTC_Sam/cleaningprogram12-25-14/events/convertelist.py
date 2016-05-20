#! /usr/bin/env python

import re
import math
import os
import fnmatch

def log(x):
	return math.log(float(x), 10.0)

def ln(x):
	return math.log(float(x))

def readfile(filepath):
	infile = open(filepath, 'r')
	return infile.read()

def writefile(filepath, outstr):
	infile = open(filepath, 'w')
	infile.write(outstr)

def findin(regexstr, strlist):
	regex = re.compile(regexstr)
	return [m.group(0) for lv in strlist for m in [regex.search(lv)] if m][0]
	
def cases(regexstr, strlist):
	regex = re.compile(regexstr)
	return [m.group(0) for lv in strlist for m in [regex.search(lv)] if m]

def latloncoordinatestolineardistance(lat1, lon1, lat2, lon2):
	return 6367.5 * math.acos(math.sin(float(lat1) * math.pi / 180) * math.sin(float(lat2) * math.pi / 180) + math.cos(float(lat1) * math.pi / 180) * math.cos(float(lat2) * math.pi / 180) * math.cos(math.fabs(scaleangle(float(lon2) * math.pi / 180 - float(lon1) * math.pi / 180))))

def scaleangle(rang):
	ang = rang % (2.0 * math.pi)
	if ang > math.pi:
		ang -= 2.0* math.pi
	if ang < -math.pi:
		ang += 2.0 * math.pi
	return ang

edlist = ((readfile('elist.txt').strip()).replace('ALERT,','ALERT')).split('\n')

for i in range(len(edlist)):
	edlist[i]=re.split('\n| \| | \||\| |\||,',edlist[i])
	if(re.match('usc.*', edlist[i][8])):
		edlist[i].remove(edlist[i][8])
	edlist[i][0:2]=edlist[i][1].split(' ')
	dli=edlist[i][0:4]
	dli.append(latloncoordinatestolineardistance(dli[2], dli[3], 47.0, -126.0))
	dli.extend(edlist[i][8:10])
	dli.append(edlist[i][11])
	nstr=dli[1].split(':')
	dli.insert(2,3600*float(nstr[0])+60*float(nstr[1])+float(nstr[2]) + dli[4] / 7.0 - 3.2**float(dli[6]))
	dli.insert(3,3600*float(nstr[0])+60*float(nstr[1])+float(nstr[2]) + dli[5] / 4.0 + 500 +4.6**float(dli[7]))
	print dli
	edlist[i]=dli
	edlist[i][0]=edlist[i][0].replace('/','-')

blist=range(len(edlist))
for i in range(len(edlist)):
	blist[i]=edlist[i][2:4]
blist.sort()

dblist = blist[:]
for i in range(len(blist)-1):
	for j in range(i+1, len(blist)):
		if blist[j][1] < blist[i][1]:
			dblist[j] = 'remove'
blist = dblist[:]
if blist.count('remove') != 0:
	for i in range(blist.count('remove')):
		blist.remove('remove')		

dblist = blist[:]
for i in range(len(blist)-1):
	for j in range(i+1, len(blist)):
		if blist[j][0] < blist[i][1]:
			dblist[i] = 'remove'
			dblist[j][0] = blist[i][0]
blist = dblist[:]
if blist.count('remove') != 0:
	for i in range(blist.count('remove')):
		blist.remove('remove')

print blist

glist = range(len(blist)+1)
glist[0] = [400,blist[0][0]]
glist[-1] = [blist[-1][1],86000]

for i in range(1,len(blist)):
	glist[i] = [blist[i-1][1],blist[i][0]]
	
wlist = []
for i in range(len(glist)):
	width = glist[i][1]-glist[i][0]
	if width >= 2000.0:
		wnum = math.floor(width/2000.0)
		inc = (width-wnum*2000.0)/2.0
		for j in range(int(wnum)):
			wlist.append(int(glist[i][0]+inc+float(j)*2000.0))

writefile('estring.txt', str(len(wlist)) + ' 2000\n' + '\n'.join(map(str, wlist)) + '\n')
