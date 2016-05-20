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

def yeardashmonthdashday2yeardotdaynumber(nstr):
	year = int(nstr[0:4])
	daynum = int(nstr[-2:])
	month = nstr[5:7]
	if year % 400 == 0:
		feblength = 28
	else:
		if year % 4 == 0:
			feblength = 29
		else:
			feblength = 28
	if month == "01":
		day = daynum + 0
	else:
		if month == "02":
			day = daynum + 31	
		else:
			if month == "03":
				day = daynum + (31 + feblength)
			else:
				if month == "04":
					day = daynum + (31 + feblength + 31)
				else:
					if month == "05":
						day = daynum + (31 + feblength + 31 + 30)
					else:
						if month == "06":
							day = daynum + (31 + feblength + 31 + 30 + 31)
						else:
							if month == "07":
								day = daynum + (31 + feblength + 31 + 30 + 31 + 30)
							else:
								if month == "08":
									day = daynum + (31 + feblength + 31 + 30 + 31 + 30 + 31)
								else:
									if month == "09":
										day = daynum + (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31)
									else:
										if month == "10":
											day = daynum + (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30)
										else:
											if month == "11":
												day = daynum + (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31)
											else:
												if month == "12":
													day = daynum + (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30)
	return str(year) + '.' + (3 - len(str(day))) * '0' + str(day)

edlist = ((readfile('events.raw.txt').strip()).replace('ALERT,','ALERT')).split('\n')

for i in range(len(edlist)):
	edlist[i]=re.split('\n| \| | \||\| |\||,',edlist[i])
	edlist[i][0:2]=edlist[i][1].split(' ')
	dli=edlist[i][0:5]
	dli.append(edlist[i][9])
	dli.append(edlist[i][11])
	nstr=dli[1].split(':')
	dli.insert(2, str(3600*float(nstr[0])+60*float(nstr[1])+float(nstr[2])))
	edlist[i]=dli
	edlisti0 = edlist[i][0]
	edlist[i][0] = edlisti0[-5:] + '/' + edlisti0[:4]
	edlist[i] = dli
	dli.insert(7, yeardashmonthdashday2yeardotdaynumber(edlisti0.replace('/','-')))
	edlist[i]=dli
	edlist[i] = ",".join(edlist[i])

writefile('events.txt', '\n'.join(edlist))
