#! /usr/bin/env python

#This script automatically files the data in the specified directory fromdir into subdirectories in the main directory for each station.
#It needs to be modified to be able to read the station name from whatever file naming system you're using.

import re
import math
import os
import fnmatch

fromdir = 'events/y2data'

dirfiles = os.listdir(fromdir)
sacfiles = fnmatch.filter(dirfiles, "20??.???.*.SAC")

for i in range(len(sacfiles)):
	sfi = sacfiles[i]
	enddropped = sfi.partition('B.')[0] + 'B'
	station = enddropped.rpartition('.')[2]
	if (os.listdir('.')).count(station) == 0:
		os.mkdir(station)
	print 'mv' + ' ' + fromdir + '/' + sfi + ' ' + station
	os.system('mv' + ' ' + fromdir + '/' + sfi + ' ' + station)
