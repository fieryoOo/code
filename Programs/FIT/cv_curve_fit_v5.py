#/usr/bin/env pytho:$n

#python version of cv_curve_fit.
#using linear inversion algorithm.

import math;
import string;
import sys;
import numpy as np;
import scipy as sp;

def cv_B_spline(nBs, degBs, zmin_Bs, zmax_Bs, disfacBs, npts, Bs_basis, depth1):
# nBs = 5
# normally disfacBs = 2
# degBs = 4 # cubic B-splie
# return B-spline basis functions in Bs_basis
	m = nBs-1+degBs;
	t = [];
	for i in range (m+1):
		t.append(0);
	for i in range (degBs):
		temp = zmin_Bs + i*(zmax_Bs-zmin_Bs)/1000.;
		t[i] = temp;
	for i in range (degBs,m+1-degBs):
		n_temp = m+1-degBs-degBs+1;
		if (disfacBs !=1):
			temp = (zmax_Bs-zmin_Bs)*(disfacBs-1)/(math.pow(disfacBs,n_temp)-1);
		else:
			temp = (i-degBs+1)*(zmax_Bs-zmin_Bs)/n_temp;
#		print temp,i,math.pow(disfacBs,(i-degBs));
		t[i] = temp*math.pow(disfacBs,(i-degBs)) + zmin_Bs;
	for i in range (m+1-degBs,m+1):
		t[i] = zmax_Bs-(zmax_Bs-zmin_Bs)/1000.*(m-i);

	print t;

	step = (zmax_Bs-zmin_Bs)/(npts-1);
	obasis = [];
	nbasis = [];
	for i in range (m):
		temp = [];
		obasis.append(temp);
		nbasis.append(temp);
	for i in range (m):
		for j in range (npts):
			nbasis[i].append(0);		
			obasis[i].append(0);
	depth = [];


	for i in range (npts):
		depth.append(i);
	for i in range (npts):
		depth[i] = i*step + zmin_Bs;
		depth1.append(i*step+zmin_Bs);
	for i in range (m):
		for j in range (npts):
			if (depth[j] >=t[i] and depth[j]<t[i+1]):
				obasis[i][j] = 1;
			else:
				obasis[i][j] = 0;
	for pp in range (1,degBs):
		for i in range (m-pp):
			for j in range (npts):
				temp = (depth[j]-t[i])/(t[i+pp]-t[i])*obasis[i][j]+(t[i+pp+1]-depth[j])/(t[i+pp+1]-t[i+1])*obasis[i+1][j];
				nbasis[i][j] = temp;
		for i in range (m-pp):
			for j in range (npts):
				obasis[i][j] = nbasis[i][j];
	nbasis[0][0] = 1;
	nbasis[nBs-1][npts-1] = 1;


	for i in range (m-pp):
		for j in range (npts):
			temp = nbasis[i][j];
			Bs_basis.append(temp-0.);


#if __name__ == '__main__':
#if (len(sys.argv) != 5):
#	print "input [top] [botom] [nBs] [velocity_model]";
#	sys.exit();

def invert_Bs ( top, botom, nBB, indep, invalue, disfacBs):
	zmin_Bs = top;
	zmax_Bs = botom;
	nBs = int(nBB);
	degBs = 4; # cubic B-splines;
	npts = 50;
#	disfacBs = 2.;
	Bs_basis = [];
	depth = [];
#	print "disfacBs, ", disfacBs;
	cv_B_spline(nBs, degBs, zmin_Bs, zmax_Bs, disfacBs, npts, Bs_basis, depth);

	ff = open("Bs.dat","w");
	for i in range (npts):
		outBs = [];
		ff.write("%g " % depth[i])
		for j in range (nBs):
			tn = Bs_basis[j*npts+i];
			outBs.append(tn);
			ff.write("%g " % tn);
		ff.write("\n");
	ff.close();

	indepth = [];
	invel = [];
	flag = 0;
#	for line in open(sys.argv[4]):
#		line.rstrip();
#		line1 = line.split();
	for i in range (len(indep)):
		if indep[i] >= zmin_Bs and indep[i] <= zmax_Bs:
			if flag == 0 and indep[i] >zmin_Bs:
				indepth.append(zmin_Bs);
				invel.append(float(invalue[i]) - 0.);
				flag = 1;
			indepth.append(indep[i] - 0.);
			invel.append(invalue[i] - 0.);
			ttvel = invalue[i] - 0.;
			if (flag == 0):
				flag = 1;
	indepth.append(zmax_Bs);
	invel.append(ttvel);
	ndata = len(invel);
	iid  = [];
	G = [];
	for i in range(ndata):
#	print indepth[i],invel[i];
		tempG = [];
		tflag = 0;
		for j in range (npts-1):
			if (depth[j] <= indepth[i] and depth[j+1] >= indepth[i]):	
				iid.append(j);
				tflag = 1;
				break;
		if (tflag == 0):
			iid.append(npts-1);
		for j in range (nBs):
#			print iid;
#			print i,j,tn,j*npts + iid[i],iid[i],"!!";
			tn = Bs_basis[j*npts + iid[i]];
#			print i,j,tn,j*npts + iid[i],iid[i];
			tempG.append(tn);
#		print indepth[i],invel[i],iid[i],depth[iid[i]];
		G.append(tempG);
	
	pBs = [];

#	pBs = np.dot(np.linalg.pinv( np.dot(np.transpose(G),G),0.00001), np.dot(np.transpose(G),invel));
#	print pBs;

	pBs = np.linalg.lstsq(G,invel);
#        print pBs;

#	pBs = numpy.linalg.pinv()i
        ff = open("in_data.dat","w");
        for i in range (len(indepth)):
                ff.write("%g %g\n" % (indepth[i],invel[i]));
        ff.close();
	
	ff = open("out_Bs.dat","w");
	for i in range (len(pBs[0])):
		ff.write("%g\n" % pBs[0][i]);
	ff.close();
	
	ff = open("out_fit.dat","w");
	for i in range(npts):
		tn = 0;
		for j in range (nBs):
			tn = tn + pBs[0][j]*Bs_basis[j*npts + i];
		ff.write(" %g %g \n" % (depth[i],tn));
	ff.close();
	return pBs[0];
