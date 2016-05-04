#!/curc/tools/x86_64/rh6/software/python/2.7.10/gcc/5.1.0/bin/python

def onclick(event):
	s = '%f %f' % (event.xdata, event.ydata)
	print s
	fout.write(s+'\n')


### main ###
# ask for file name
import sys
if len(sys.argv) != 4:
	print "Usage:",sys.argv[0],"[inname1] [inname2] [outname]"
	sys.exit()
inname1=sys.argv[1]
inname2=sys.argv[2]
outname=sys.argv[3]

# load x and y columns from input file
x1 = []; y1 = []
with open(inname1) as fin:
	for line in fin:
		x1.append( line.split()[2] )
		y1.append( line.split()[3] )

x2 = []; y2 = []
with open(inname2) as fin:
	for line in fin:
		x2.append( line.split()[2] )
		y2.append( line.split()[3] )

# plot input disp
import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x1, y1, 'bo', x2, y2, 'ro')
ax.axis([0,50,0.,5.])

# open outfile
fout = open(outname, 'w')
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

fout.close()
