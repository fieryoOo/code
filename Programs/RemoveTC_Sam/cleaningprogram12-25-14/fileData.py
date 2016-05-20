#! /usr/bin/env python

#This script automatically files the data in the main directory into subdirectories for each station.
#It needs to be modified to be able to read the station name from whatever file naming system you're using.

import re
import math
import os
import fnmatch

dirfiles = os.listdir('.')
sacfiles = fnmatch.filter(dirfiles, "20??.???.*.SAC")

for i in range(len(sacfiles)):
	sfi = sacfiles[i]
	enddropped = sfi.partition('..')[0]
	station = enddropped.rpartition('.')[2]
	if (os.listdir('.')).count(station) == 0:
		os.mkdir(station)
	os.system('mv' + ' ' + sfi + ' ' + station)
