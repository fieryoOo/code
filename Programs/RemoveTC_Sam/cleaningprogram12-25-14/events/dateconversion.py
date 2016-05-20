#! /usr/bin/env python

#This script converts between yyyy-mm-dd and yyyy.ddd date formats
#It reads in a list of yyyy-mm-dd days from dateconvertinput.txt
#It writes a list of yyyy.ddd days to dateconvertoutput.txt
#It can handle days in yyyy/mm/dd format, too

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

#Searches a list of strings strlist and returns the string that matches the regular expression string regexstr.	
def findin(regexstr, strlist):
	regex = re.compile(regexstr)
	return [m.group(0) for lv in strlist for m in [regex.search(lv)] if m][0]

#Converts a yyyy.ddd day to a yyyy-mm-dd day
#This function is not actually used in this script, but it is included so that the script can be altered to convert yyyy.ddd days to yyyy-mm-dd days instead
def yeardotdaynumber2yeardashmonthdashday(nstr):
	year = int(nstr[0:4])
	daynum = int(nstr[5:8])
	if year % 400 == 0:
		feblength = 28
	else:
		if year % 4 == 0:
			feblength = 29
		else:
			feblength = 28
	if daynum <= 31:
		month = "01"
		day = daynum - 0
	else:
		if daynum <= 31 + feblength:
			month = "02"
			day = daynum - 31
		else:
			if daynum <= 31 + feblength + 31:
				month = "03"
				day = daynum - (31 + feblength)
			else:
				if daynum <= 31 + feblength + 31 + 30:
					month = "04"
					day = daynum - (31 + feblength + 31)
				else:
					if daynum <= 31 + feblength + 31 + 30 + 31:
						month = "05"
						day = daynum - (31 + feblength + 31 + 30)
					else:
						if daynum <= 31 + feblength + 31 + 30 + 31 + 30:
							month = "06"
							day = daynum - (31 + feblength + 31 + 30 + 31)
						else:
							if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31:
								month = "07"
								day = daynum - (31 + feblength + 31 + 30 + 31 + 30)
							else:
								if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31 + 31:
									month = "08"
									day = daynum - (31 + feblength + 31 + 30 + 31 + 30 + 31)
								else:
									if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30:
										month = "09"
										day = daynum - (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31)
									else:
										if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31:
											month = "10"
											day = daynum - (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30)
										else:
											if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30:
												month = "11"
												day = daynum - (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31)
											else:
												if daynum <= 31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31:
													month = "12"
													day = daynum - (31 + feblength + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30)
	tempdaystr = str(day)
	if len(tempdaystr) == 2:
		daystr = tempdaystr
	else:
		daystr = "0" + tempdaystr
	rstring = str(year) + "-" + month + "-" + daystr
	return rstring

#Converts a yyyy-mm-dd day to a yyyy.ddd day
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

#Prints out the list of days in yyyy-mm-dd format, reading them in from dateconvertinput.txt
print re.split('\n|\r', ((readfile('dateconvertinput.txt').replace('/', '-')).strip('\n')).strip('\r'))

#Reads in the list of yyyy-mm-dd (or yyyy/mm/dd) days from dateconvertinput.txt
#Writes the list of converted yyyy.ddd days to dateconvertoutput.txt
writefile('dateconvertoutput.txt', '\n'.join([yeardashmonthdashday2yeardotdaynumber(date) for date in re.split('\n|\r', ((readfile('dateconvertinput.txt').replace('/', '-')).strip('\n')).strip('\r'))]))

