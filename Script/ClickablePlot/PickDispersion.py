#!/curc/tools/x86_64/rh6/software/python/2.7.10/gcc/5.1.0/bin/python

def onclick(event):
	s = '%f %f' % (event.xdata, event.ydata)
	print s
	fout.write(s+'\n')


### main ###
# ask for file name
import sys
if len(sys.argv) != 6:
	print "Usage:",sys.argv[0],"[inname1] [col1] [inname2] [col2] [outname]"
	sys.exit()
inname1=sys.argv[1]
col1=int(sys.argv[2])
inname2=sys.argv[3]
col2=int(sys.argv[4])
outname=sys.argv[5]

# load x and y columns from input file
x1 = []; y1 = []
with open(inname1) as fin:
	for line in fin:
		x1.append( line.split()[2] )
		y1.append( line.split()[col1-1] )

x2 = []; y2 = []
with open(inname2) as fin:
	for line in fin:
		x2.append( line.split()[2] )
		y2.append( line.split()[col2-1] )

# plot input disp
import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x1, y1, 'bo', x2, y2, 'ro', ms=1)
ax.plot(x1, y1, 'bo', ms=0.1, mec='b')
ax.plot(x2, y2, 'ro', ms=0.1, mec='r')
#ax.scatter(x1, y1, marker='o', c='b')
ax.axis([0,45,0.,5.])

# open outfile
fout = open(outname, 'w')
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

fout.close()
