#! /usr/bin/env python

import time
#print(time.time())
import re
import math
import os
import fnmatch
import numpy as np
import matplotlib
#print(time.time())
import matplotlib.pyplot as plt
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

#Initialize the functions that plot the pressure to vertical transfer functions:

#Plots the admittance
def aplots(data, extension):
	freqs = floatlist(getcolumn(data, 0))
	logamps = [log(float(li)) for li in getcolumn(data, 1)]
	ampanderrs = floattable(getcolumns(data, [1, 2]))
	logamperr = [[-1 * (log(settoifnegativeorzero(li[0] - li[1], 0.00000000000001)) - log(li[0]))  for li in ampanderrs], [log(li[0] + li[1]) - log(li[0])  for li in ampanderrs]]
	rfig = plt.figure()
	plt.errorbar(freqs, logamps, yerr = logamperr, fmt = 'o')
	plt.axis([0.00, 0.405, min(logamps), max(logamps)])
	plt.title(dayn + '  ' + station, fontsize = 28)
	plt.ylabel('log(Admittance)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.p.aplot.' + extension + '.png')
	plt.axis([0.004, 0.12, min(logamps), max(logamps)])
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.p.aplotZoom.' + extension + '.png')
	plt.close(rfig)

#Plots the coherence
def cplots(data, extension):
	freqs = floatlist(getcolumn(data, 0))
	cohers = floatlist(getcolumn(data, 5))
	rfig = plt.figure()
	plt.plot(freqs, cohers, 'bo')
	plt.axis([0.00, 0.405, 0.0, 1.0])
	plt.title(dayn + '  ' + station, fontsize = 28)
	plt.ylabel('Coherence', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.p.cplot.' + extension + '.pdf')
	plt.axis([0.004, 0.12, 0.0, 1.0])
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.p.cplotZoom.' + extension + '.pdf')
	plt.close(rfig)

#Plots the phase	
def pplot(data, extension):
	freqs = floatlist(getcolumn(data, 0))
	phases = floatlist(getcolumn(data, 3))
	phaseerrelement = floatlist(getcolumn(data, 4))
	rfig = plt.figure()
	plt.errorbar(freqs, phases, yerr = [phaseerrelement, phaseerrelement], fmt = 'o')
	plt.axis([0.00, 0.405, -0.5, 0.5])
	plt.title(dayn + '  ' + station, fontsize = 28)
	plt.ylabel('Phase (cycles)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.p.pplot.' + extension + '.png')
	plt.close(rfig)

#Initialize the functions that plot the horizontal to vertical transfer functions:

#Plots the admittance
def tiltaplots(data, extension, direction, coherence):
	freqs = floatlist(getcolumn(data, 0))
	logamps = [log(float(li)) for li in getcolumn(data, 1)]
	ampanderrs = floattable(getcolumns(data, [1, 2]))
	logamperr = [[-1 * (log(settoifnegativeorzero(li[0] - li[1], 0.00000000000001)) - log(li[0]))  for li in ampanderrs], [log(li[0] + li[1]) - log(li[0])  for li in ampanderrs]]
	rfig = plt.figure()
	plt.errorbar(freqs, logamps, yerr = logamperr, fmt = 'o')
	plt.axis([0.00, 0.405, min(logamps), max(logamps)])
	plt.title(dayn + '  ' + station + '   ' + direction + '   ' + coherence, fontsize = 28)
	plt.ylabel('log(Admittance)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.aplot.' + extension + '.png')
	plt.axis([0.004, 0.12, min(logamps), max(logamps)])
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.aplotZoom.' + extension + '.png')
	plt.close(rfig)


#Plots the coherence
def tiltcplots(data, extension, direction, coherence):
	freqs = floatlist(getcolumn(data, 0))
	cohers = floatlist(getcolumn(data, 5))
	rfig = plt.figure()
	plt.plot(freqs, cohers, 'bo')
	plt.axis([0.00, 0.405, 0.0, 1.0])
	plt.title(dayn + '  ' + station + '   ' + direction + '   ' + coherence, fontsize = 28)
	plt.ylabel('Coherence', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.cplot.' + extension + '.pdf')
	plt.axis([0.004, 0.12, 0.0, 1.0])
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.cplotZoom.' + extension + '.pdf')
	plt.close(rfig)

#Plots the phase	
def tiltpplot(data, extension, direction, coherence):
	freqs = floatlist(getcolumn(data, 0))
	phases = floatlist(getcolumn(data, 3))
	phaseerrelement = floatlist(getcolumn(data, 4))
	rfig = plt.figure()
	plt.errorbar(freqs, phases, yerr = [phaseerrelement, phaseerrelement], fmt = 'o')
	plt.axis([0.00, 0.405, -0.5, 0.5])
	plt.title(dayn + '  ' + station + '   ' + direction + '   ' + coherence, fontsize = 28)
	plt.ylabel('Phase (cycles)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.pplot.' + extension + '.png')
	plt.close(rfig)

#Plots the fitted tilt
def tiltaplotfit(data, extension, direction, coherence, actable, btable):
	freqs = floatlist(getcolumn(data, 0))
	logamps = [log(float(li)) for li in getcolumn(data, 1)]
	ampanderrs = floattable(getcolumns(data, [1, 2]))
	logamperr = [[-1 * (log(settoifnegativeorzero(li[0] - li[1], 0.00000000000001)) - log(li[0]))  for li in ampanderrs], [log(li[0] + li[1]) - log(li[0])  for li in ampanderrs]]
	rfig = plt.figure()
	plt.errorbar(freqs, logamps, yerr = logamperr, fmt = 'o')
	x = np.arange(0.005, btable[-1][1], 0.0001)
	def segfunc(i):
		return lambda x: np.log10(actable[i][0] + actable[i][1] * x + actable[i][2] * x**2)
	y = np.piecewise(x, [(x >= btable[i][0]) & (x <= btable[i][1]) for i in range(len(btable))], [segfunc(j) for j in list(range(len(actable)))])
	plt.plot(x, y)
	plt.axis([0.004, btable[-1][1], min(logamps), max(logamps)])
	plt.title(dayn + '  ' + station + '   ' + direction + '   ' + coherence, fontsize = 28)
	plt.ylabel('log(Admittance)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.aplotFit.' + extension + '.png')
	plt.close(rfig)

#Plots the fitted phase
def tiltpplotfit(data, extension, direction, coherence, pctable, btable):
	freqs = floatlist(getcolumn(data, 0))
	phases = floatlist(getcolumn(data, 3))
	phaseerrelement = floatlist(getcolumn(data, 4))
	rfig = plt.figure()
	plt.errorbar(freqs, phases, yerr = [phaseerrelement, phaseerrelement], fmt = 'o')
	x = np.arange(0.005, btable[-1][1], 0.0001)
	def segfunc(i):
		return lambda x: pctable[i][0] + pctable[i][1] * x + pctable[i][2] * x**2
	y = np.piecewise(x, [(x >= btable[i][0]) & (x <= btable[i][1]) for i in range(len(btable))], [segfunc(j) for j in list(range(len(pctable)))])
	plt.plot(x, y)
	plt.axis([0.004, btable[-1][1], -0.5, 0.5])
	plt.title(dayn + '  ' + station + '   ' + direction + '   ' + coherence, fontsize = 28)
	plt.ylabel('Phase (cycles)', fontsize = 24)
	plt.xlabel('Frequency (Hz)', fontsize = 24)
	rfig.savefig('plots/' + station + '/' + dayn + '/' + station + '.' + dayn + '.t.pplotFit.' + extension + '.png')
	plt.close(rfig)


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
nbands = 1

#Exports the horizontal to vertical transfer function input file.
ttffile = '1\n' + ttfdecfn + '\n1 1\n' + str(nbands) + '\n.005 .06\n' + f1 + '\n' + f2 + '\n' + fz + '\n' + windowstr + '\n'
writefile('ttf.inp', ttffile)

print(ttffile)

#Calculates the horizontal to vertical transfer function file.
os.system('./findtransferfns < ttf.inp')
print('tilt transfer functions made')
#print(time.time())

#Imports and extracts the data from the tilt transfer function file.
rimp = readfile(ttfdecfn).split('\n')
ptcinfodec = rimp[-1*(1+4*nbands):-1]
for i in range(1,len(rimp)):
	rimp[i] = re.split('  | ',((rimp[i]).strip()))
directiondec = ((rimp[0]).strip()).split('  ')[0]
coherencedec = ((rimp[0]).strip()).split('  ')[1]
tddec = rimp[2:-1*(1+4*nbands)]

actable = getvalues(ptcinfodec, [2 + 4 * i for i in range(nbands)])
for i in range(len(actable)):
	actable[i] = floatlist(re.split('  | ',((actable[i]).strip())))
pctable = getvalues(ptcinfodec, [3 + 4 * i for i in range(nbands)])
for i in range(len(pctable)):
	pctable[i] = floatlist(re.split('  | ',((pctable[i]).strip())))
btable = getvalues(ptcinfodec, [1 + 4 * i for i in range(nbands)])
for i in range(len(btable)):
	btable[i] = floatlist(re.split('  | ',((btable[i]).strip())))

print(actable)
print(pctable)
print(btable)

#Plots the tilt transfer functions.
tiltaplots(tddec, 'dec', directiondec, coherencedec)
tiltcplots(tddec, 'dec', directiondec, coherencedec)
tiltpplot(tddec, 'dec', directiondec, coherencedec)
tiltaplotfit(tddec, 'dec', directiondec, coherencedec, actable, btable)
tiltpplotfit(tddec, 'dec', directiondec, coherencedec, pctable, btable)
print('tilt plots made')
#print(time.time())

#Exports the pressure to vertical transfer function input file.
writefile('ptf.inp', ptfdecfn + '\n1 1\n' + fh + '\n' + fz + '\n' + windowstr + '\n')
#print(time.time())

#Calculates the pressure to vertical transfer function file.
os.system('./transferp2z < ptf.inp')
print('pressure transfer functions made')
#print(time.time())

#Imports and extracts the data from the pressure transfer function file.
pddec = (readfile(ptfdecfn).strip()).split('\n')
for i in range(0,len(pddec)):
	pddec[i] = re.split('  | ',((pddec[i]).strip()))

#Plots the pressure transfer functions.
aplots(pddec, 'dec')
cplots(pddec, 'dec')
pplot(pddec, 'dec')
print('pressure plots made')
#print(time.time())

#Writes the file used to give the windows for the pressure removal SAC macro.
writefile('p2z.w.inp', '1\n.005 .06\n1\n.005 ' + str(cedge(depth(station))*0.85) + '\n' + windowstr + '\n')

#Writes the pressure removal SAC macro.
writefile('macros/correctPressure.m', '# remove pressure\n\nbinoperr npts warning delta warning\n\nsc echo ' + fz + '\n\nsetbb dpgfile ' + fh + '\nsetbb notfile ' + fz + '\nsetbb outptrans ../p2zout/transferp2zw.' + dayn + '.' + station + '\nsetbb outpcoeff ../ptcp2zs/' + station + '.' + dayn + '.ptc\nsetbb outrecord ../' + fdep + '\n\n\ncp ../%dpgfile tempd\ncp ../%notfile tempztiltc\n\nsc ./../findP2ZtransferfnsAuto < ../p2z.w.inp\ncp tempP2Ztransfer %outptrans\ncp tempP2Zcoeff %outpcoeff\ncut off\nr tempd\nrmean\ntaper\n# tapering because very long-period tidal components\n# cause discontinuities at ends of record\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoiseAuto\nreadsp amph tempf\nifft\nrmean\nw tempPred\n# enter vertical filename\nr tempztiltc\nrmean\ntaper\nsubf tempPred\nw %outrecord\n\nquit\n\n')

#Executes the pressure removal SAC macro.
os.system('cd macros\necho m correctPressure.m | sac\ncd ..')
print('pressure removed')
#print(time.time())

#Exports the pressure to vertical transfer function input file.
writefile('ptf.inp', ptfdepfn + '\n1 1\n' + fh + '\n' + fdep + '\n' + windowstr + '\n')
#print(time.time())

#Calculates the pressure to vertical transfer function file.
os.system('./transferp2z < ptf.inp')
print('pressure transfer functions made')
#print(time.time())

#Imports and extracts the data from the pressure transfer function file.
pddep = (readfile(ptfdepfn).strip()).split('\n')
for i in range(0,len(pddep)):
	pddep[i] = re.split('  | ',((pddep[i]).strip()))

#Plots the pressure transfer functions.
aplots(pddep, 'dep')
cplots(pddep, 'dep')
pplot(pddep, 'dep')
print('pressure plots made')
#print(time.time())

#Exports the horizontal to vertical transfer function input file.
ttffile = '1\n' + ttfdepfn + '\n1 1\n' + str(nbands) + '\n.005 .06\n' + f1 + '\n' + f2 + '\n' + fz + '\n' + windowstr + '\n'
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

#Plots the tilt transfer functions.
tiltaplots(tddep, 'dep', directiondep, coherencedep)
tiltcplots(tddep, 'dep', directiondep, coherencedep)
tiltpplot(tddep, 'dep', directiondep, coherencedep)
print('tilt plots made')
#print(time.time())

#Writes the polynomial transfer function model coefficients file (.ptc).
ptclinelist = list(range(nbands * 4))
for i in range(nbands):
	ptclinelist.remove(4 * i)
writefile('ptcs/' + station + '.' + dayn + '.dec.ptc', station + '.' + dayn + '\n' + str(nbands) + '\n' + ('\n'.join(getvalues(ptcinfodec, ptclinelist))).strip() + '\n')

#Writes the tilt removal SAC macro.
writefile('macros/correcttiltBYtrans.m', 'cut off\n\nbinoperr npts warning delta warning\n\n\nsc cd ../' + station + '\n\nsc cp ../ptcs/' + station + '.' + dayn + '.dec.ptc ../Platetransfercoeff\n\nsetbb azim ' + directiondep + '\nsetbb bh1fn ../' + f1 + '\nsetbb bh2fn ../' + f2 + '\nsetbb bhzfn ../' + fz + '\nsetbb outfn ../' + fnot + '\n\nevaluate to conv2rad 3.14159 / 180.\nevaluate to a1 %azim * %conv2rad\nevaluate to cos1 cos %a1\nevaluate to sin1 sin %a1\n\n# rotate horizontals to tilt direction \nr %bh1fn\nrmean\nmul %cos1\nw temp1c\nr %bh2fn\nrmean\nmul %sin1\naddf temp1c\nw temph\n\n#predict vertical noise from horizontal record and transferfunction\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoise\nreadsp amph tempf\nifft\nrmean\nw tempPred\n\n# subtract from vertical\nr %bhzfn\nsubf tempPred\nw %outfn\n\nquit\n\n')

#Executes the tilt removal SAC macro.
os.system('cd macros\necho m correcttiltBYtrans.m | sac\ncd ..')
print('tilt removed')
#print(time.time())

#Writes the pressure to vertical transfer function input file.
writefile('ptf.inp', ptfnotfn + '\n1 1\n' + fh + '\n' + fnot + '\n' + windowstr + '\n')

#Calculates the pressure to vertical transfer function file.
os.system('./transferp2z < ptf.inp')
print('pressure transfer functions made')
#print(time.time())

#Imports and extracts the data from the pressure transfer function file.
pdnot = (readfile(ptfnotfn).strip()).split('\n')
for i in range(0,len(pdnot)):
	pdnot[i] = re.split('  | ',((pdnot[i]).strip()))

#Plots the pressure transfer functions.
aplots(pdnot, 'not')
cplots(pdnot, 'not')
pplot(pdnot, 'not')
print('pressure plots made')
#print(time.time())

#Writes the horizontal to vertical transfer function input file.
writefile('ttf.inp', '1\n' + ttfnotfn + '\n1 1\n1\n.005 .06\n' + f1 + '\n' + f2 + '\n' + fnot + '\n' + windowstr + '\n')

#Calculates the horizontal to vertical transfer function file.
os.system('./findtransferfns < ttf.inp')
print('tilt transfer functions made')
#print(time.time())

#Imports and extracts the data from the tilt transfer function file.
rimp = readfile(ttfnotfn).split('\n')
ptcinfonot = rimp[-5:-1]
for i in range(1,len(rimp)):
	rimp[i] = re.split('  | ',((rimp[i]).strip()))
directionnot = ((rimp[0]).strip()).split('  ')[0]
coherencenot = ((rimp[0]).strip()).split('  ')[1]
tdnot = rimp[2:-5]

#Plots the tilt transfer functions.
tiltaplots(tdnot, 'not', directionnot, coherencenot)
tiltcplots(tdnot, 'not', directionnot, coherencenot)
tiltpplot(tdnot, 'not', directionnot, coherencenot)
print('tilt plots made')
#print(time.time())

#Writes the file used to give the windows for the pressure removal SAC macro.
writefile('p2z.w.inp', '1\n.005 .06\n1\n.005 ' + str(cedge(depth(station))*0.85) + '\n' + windowstr + '\n')

#Writes the pressure removal SAC macro.
writefile('macros/correctPressure.m', '# remove pressure\n\nbinoperr npts warning delta warning\n\nsc echo ' + fnot + '\n\nsetbb dpgfile ' + fh + '\nsetbb notfile ' + fnot + '\nsetbb outptrans ../p2zout/transferp2zw.' + dayn + '.' + station + '\nsetbb outpcoeff ../ptcp2zs/' + station + '.' + dayn + '.ptc\nsetbb outrecord ../' + fnop + '\n\n\ncp ../%dpgfile tempd\ncp ../%notfile tempztiltc\n\nsc ./../findP2ZtransferfnsAuto < ../p2z.w.inp\ncp tempP2Ztransfer %outptrans\ncp tempP2Zcoeff %outpcoeff\ncut off\nr tempd\nrmean\ntaper\n# tapering because very long-period tidal components\n# cause discontinuities at ends of record\nfft amph\nwritesp amph temp\nsc ./../remveWaterNoiseAuto\nreadsp amph tempf\nifft\nrmean\nw tempPred\n# enter vertical filename\nr tempztiltc\nrmean\ntaper\nsubf tempPred\nw %outrecord\n\nquit\n\n')

#Executes the pressure removal SAC macro.
os.system('cd macros\necho m correctPressure.m | sac\ncd ..')
print('pressure removed')
#print(time.time())

#Writes the pressure to vertical transfer function input file.
writefile('ptf.inp', ptfnopfn + '\n1 1\n' + fh + '\n' + fnop + '\n' + windowstr + '\n')

#Calculates the pressure to vertical transfer function file.
os.system('./transferp2z < ptf.inp')
print('pressure transfer functions made')
#print(time.time())

#Imports and extracts the data from the pressure transfer function file.
pdnop = (readfile(ptfnopfn).strip()).split('\n')
for i in range(0,len(pdnop)):
	pdnop[i] = re.split('  | ',((pdnop[i]).strip()))

#Plots the pressure transfer functions.
aplots(pdnop, 'nop')
cplots(pdnop, 'nop')
pplot(pdnop, 'nop')
print('pressure plots made')
#print(time.time())

#Writes the horizontal to vertical transfer function input file.
writefile('ttf.inp', '1\n' + ttfnopfn + '\n1 1\n1\n.005 .06\n' + f1 + '\n' + f2 + '\n' + fnop + '\n' + windowstr + '\n')

#Calculates the horizontal to vertical transfer function file.
os.system('./findtransferfns < ttf.inp')
print('tilt transfer functions made')
#print(time.time())

#Imports and extracts the data from the tilt transfer function file.
rimp = readfile(ttfnopfn).split('\n')
ptcinfonop = rimp[-5:-1]
for i in range(1,len(rimp)):
	rimp[i] = re.split('  | ',((rimp[i]).strip()))
directionnop = ((rimp[0]).strip()).split('  ')[0]
coherencenop = ((rimp[0]).strip()).split('  ')[1]
tdnop = rimp[2:-5]

#Plots the tilt transfer functions.
tiltaplots(tdnop, 'nop', directionnop, coherencenop)
tiltcplots(tdnop, 'nop', directionnop, coherencenop)
tiltpplot(tdnop, 'nop', directionnop, coherencenop)
print('tilt plots made')
#print(time.time())

