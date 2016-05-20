#! /usr/bin/env python

import time
print(time.time())
import re
import math
import os
import os.path
import fnmatch
#import numpy as np
#import matplotlib
#print(time.time())
#import matplotlib.pyplot as plt
#print(time.time())

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

getvar = lambda searchList, ind: [searchList[j] for j in ind]

def getvalues(list, ilist):
    rlist = [list[ilist[0]]]
    if len(ilist) > 1:
        for i in range(1,len(ilist)):
            rlist.append(list[ilist[i]])
    return rlist
	
def getcolumns(table, clist):
	rtable = list(range(len(table)))
	for i in range(len(table)):
		rtable[i] = getvar(table[i], clist)
	return rtable
	
def getcolumn(table, cnum):
	rtable = list(range(len(table)))
	for i in range(len(table)):
		rtable[i] = table[i][cnum]
	return rtable
	
def floatlist(list):
	return [float(li) for li in list]
	
def floattable(table):
	return [[float(li) for li in lni] for lni in table]

def settoifnegativeorzero(num, setval):
	if num <= 0:
		rval = setval
	else:
		rval = num
	return rval

#Reads in the station depths from a file and formats them as a nested list.
depths = readfile('depths.csv')
depths = re.split('\n|\r', depths)
for i in range(len(depths)):
	depths[i] = depths[i].split(',')

#Returns the depth for a given station.
def depth(station):
	depths0 = list(range(len(depths)))
	for i in range(len(depths)):
		depths0[i] = depths[i][0]
	ind = depths0.index(station)
	return math.fabs(float(depths[ind][1]))
	
#Imports and implements a function giving the compliance edge as a function of station depth.
em = readfile('em.csv').split('\n')
def cedge(d):
	return 10.0 ** (float(em[0]) + float(em[1]) * log(d))

#Loads the raw list of input files
rawd = (readfile('inputfiles.txt').strip()).split('\n')

#Searches the list of input files to get the pressure and vertical displacement filenames.
fz = findin('.*Z\.dec\.SAC', rawd)
fh = findin('.*H\.dec3\.SAC', rawd)
f1 = findin('.*1\.dec\.SAC', rawd)
f2 = findin('.*2\.dec\.SAC', rawd)

#Calculates the output filenames for the noise-removed data.
fdep = fz[:-7] + 'dep.SAC'
fnot = fz[:-7] + 'not.SAC'
fnop = fz[:-7] + 'nop.SAC'

#Find the day number and station
station = fz.partition('/')[0]
dayn = (fz.partition('/')[2])[:8]
print(station)
print(dayn)

#Create subdirectories
os.system('mkdir out/' + station)
os.system('mkdir plots/' + station)
os.system('mkdir plots/' + station + '/' + dayn)

#Initialize the filenames for the transfer functions.
ptfdecfn = 'out/' + station + '/' + station + '.'  + dayn + '.dec.pt'
ptfdepfn = 'out/' + station + '/' + station + '.'  + dayn + '.dep.pt'
ptfnotfn = 'out/' + station + '/' + station + '.'  + dayn + '.not.pt'
ptfnopfn = 'out/' + station + '/' + station + '.'  + dayn + '.nop.pt'
ttfdecfn = 'out/' + station + '/' + station + '.'  + dayn + '.dec.tt'
ttfdepfn = 'out/' + station + '/' + station + '.'  + dayn + '.dep.tt'
ttfnotfn = 'out/' + station + '/' + station + '.'  + dayn + '.not.tt'
ttfnopfn = 'out/' + station + '/' + station + '.'  + dayn + '.nop.tt'

#Reads in the string of windows from a file, based off the day number.  If there is a hand window file, it uses that.  Otherwise, it uses the default file.
if os.path.isfile('events/handwindows/' + dayn + '/estring.' + dayn + '.' + station + '.txt'):
	windowfile = 'events/handwindows/' + dayn + '/estring.' + dayn + '.' + station + '.txt'
	print('\n************************************************************\n   ALERT: Hand windows have been chosen and will be used!   \n************************************************************\n')
else:
	windowfile = 'events/estring.' + dayn + '.txt'
print(windowfile)
windowstr = readfile(windowfile)
print(windowstr)

#Set the number of bands for the fit
nbands = 3

#Writes the file used to give the windows for the pressure removal SAC macro.
writefile('p2z.w.inp', '1\n.005 .06\n1\n.005 ' + str(cedge(depth(station))*0.85) + '\n' + windowstr + '\n')

#Writes the pressure removal SAC macro.
writefile('macros/correctPressure.m', '# remove pressure\n\nbinoperr npts warning delta warning\n\nsc echo ' + fz + '\n\nsetbb dpgfile ' + fh + '\nsetbb notfile ' + fz + '\nsetbb outptrans ../p2zout/transferp2zw.' + dayn + '.' + station + '\nsetbb outpcoeff ../ptcp2zs/' + station + '.' + dayn + '.ptc\nsetbb outrecord ../' + fdep + '\n\n\ncp ../%dpgfile tempd\ncp ../%notfile tempztiltc\n\nsc ./../findP2ZtransferfnsAuto < ../p2z.w.inp\ncp tempP2Ztransfer %outptrans\ncp tempP2Zcoeff %outpcoeff\ncut off\nr tempd\nrmean\ntaper\n# tapering because very long-period tidal components\n# cause discontinuities at ends of record\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoiseAuto\nreadsp amph tempf\nifft\nrmean\nw tempPred\n# enter vertical filename\nr tempztiltc\nrmean\ntaper\nsubf tempPred\nw %outrecord\n\nquit\n\n')

#Executes the pressure removal SAC macro.
os.system('cd macros\necho m correctPressure.m | sac\ncd ..')
print('pressure removed')
#print(time.time())

#Exports the horizontal to vertical transfer function input file.
if cedge(depth(station))*0.85 < 0.055:
	ttffile = '1\n' + ttfdepfn + '\n1 1\n' + str(nbands) + '\n.005 ' + str((cedge(depth(station))*0.85-0.005)/2.0 + 0.005) + '\n' + str((cedge(depth(station))*0.85-0.005)/2.0 + 0.0050001) + ' ' + str(cedge(depth(station))*0.85) + '\n' + str(cedge(depth(station))*0.85 + 0.0000001) + ' .06\n' + f1 + '\n' + f2 + '\n' + fdep + '\n' + windowstr + '\n'
else:
	ttffile = '1\n' + ttfdepfn + '\n1 1\n' + str(nbands) + '\n.005 ' + str((0.06-0.005)/3.0 + 0.005) + '\n' + str((0.06-0.005)/3.0 + 0.0050001) + ' ' + str((0.06-0.005)*2.0/3.0 + 0.005) + '\n' + str((0.06-0.005)*2.0/3.0 + 0.0050001) + ' ' + str(0.06) + '\n' + f1 + '\n' + f2 + '\n' + fdep + '\n' + windowstr + '\n'
writefile('ttf.inp', ttffile)

print(ttffile)

#Calculates the horizontal to vertical transfer function file.
os.system('./findtransferfns < ttf.inp')
print('tilt transfer functions made')
#print(time.time())

#Imports and extracts the data from the tilt transfer function file.
rimp = readfile(ttfdepfn).split('\n')
ptcinfodep = rimp[-1*(1+4*nbands):-1]
for i in range(1,len(rimp)):
	rimp[i] = re.split('  | ',((rimp[i]).strip()))
directiondep = ((rimp[0]).strip()).split('  ')[0]
coherencedep = ((rimp[0]).strip()).split('  ')[1]
tddep = rimp[2:-1*(1+4*nbands)]

print(directiondep)
print(coherencedep)
print(ptcinfodep)

#Writes the polynomial transfer function model coefficients file (.ptc).
ptclinelist = list(range(nbands * 4))
for i in range(nbands):
	ptclinelist.remove(4 * i)
writefile('ptcs/' + station + '.' + dayn + '.dep.ptc', station + '.' + dayn + '\n' + str(nbands) + '\n' + ('\n'.join(getvalues(ptcinfodep, ptclinelist))).strip() + '\n')

#Writes the tilt removal SAC macro.
writefile('macros/correcttiltBYtrans.m', 'cut off\n\nbinoperr npts warning delta warning\n\n\nsc cd ../' + station + '\n\nsc cp ../ptcs/' + station + '.' + dayn + '.dep.ptc ../Platetransfercoeff\n\nsetbb azim ' + directiondep + '\nsetbb bh1fn ../' + f1 + '\nsetbb bh2fn ../' + f2 + '\nsetbb bhzfn ../' + fz + '\nsetbb outfn ../' + fnot + '\n\nevaluate to conv2rad 3.14159 / 180.\nevaluate to a1 %azim * %conv2rad\nevaluate to cos1 cos %a1\nevaluate to sin1 sin %a1\n\n# rotate horizontals to tilt direction \nr %bh1fn\nrmean\nmul %cos1\nw temp1c\nr %bh2fn\nrmean\nmul %sin1\naddf temp1c\nw temph\n\n#predict vertical noise from horizontal record and transferfunction\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoise\nreadsp amph tempf\nifft\nrmean\nw tempPred\n\n# subtract from vertical\nr %bhzfn\nsubf tempPred\nw %outfn\n\nquit\n\n')

#Executes the tilt removal SAC macro.
os.system('cd macros\necho m correcttiltBYtrans.m | sac\ncd ..')
print('tilt removed')
#print(time.time())

#Writes the file used to give the windows for the pressure removal SAC macro.
writefile('p2z.w.inp', '1\n.005 .06\n1\n.005 ' + str(cedge(depth(station))*0.85) + '\n' + windowstr + '\n')

#Writes the pressure removal SAC macro.
writefile('macros/correctPressure.m', '# remove pressure\n\nbinoperr npts warning delta warning\n\nsc echo ' + fnot + '\n\nsetbb dpgfile ' + fh + '\nsetbb notfile ' + fnot + '\nsetbb outptrans ../p2zout/transferp2zw.' + dayn + '.' + station + '\nsetbb outpcoeff ../ptcp2zs/' + station + '.' + dayn + '.ptc\nsetbb outrecord ../' + fnop + '\n\n\ncp ../%dpgfile tempd\ncp ../%notfile tempztiltc\n\nsc ./../findP2ZtransferfnsAuto < ../p2z.w.inp\ncp tempP2Ztransfer %outptrans\ncp tempP2Zcoeff %outpcoeff\ncut off\nr tempd\nrmean\ntaper\n# tapering because very long-period tidal components\n# cause discontinuities at ends of record\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoiseAuto\nreadsp amph tempf\nifft\nrmean\nw tempPred\n# enter vertical filename\nr tempztiltc\nrmean\ntaper\nsubf tempPred\nw %outrecord\n\nquit\n\n')

#Executes the pressure removal SAC macro.
os.system('cd macros\necho m correctPressure.m | sac\ncd ..')
print('pressure removed')
print(time.time())


