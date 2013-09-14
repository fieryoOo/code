import math;
import string;
import sys;
import os;
import numpy as np;

def invert_1 ( inbaz, indat, inun):
        m = len(inbaz); # data space;
        n = 5; # model space;
        G = [];
        A0 = 0;
        A1 = 0;
        A2 = 0;
        A3 = 0;
        fi1 = 0;
        fi2 = 0;
        fi3 = 0;
#       print m,len(inun);
        U = np.zeros((m,m));
        for i in range (m):
                U[i][i] = 1./inun[i];
        for i in range (m):
                tG = [];
                tbaz = inbaz[i]/180*(math.pi);
                tG.append(1);
                tG.append(math.sin(tbaz));
                tG.append(math.cos(tbaz));
                tG.append(math.sin(tbaz*2));
                tG.append(math.cos(tbaz*2));
                tG.append(math.sin(tbaz*4));
                tG.append(math.cos(tbaz*4));
                G.append(tG);
        G1 = np.array(G);
        G1 = np.dot(U,G1);
#       print G1;
        d = np.array(indat);
        d = d.T;
        d = np.dot(U,d);
#       print d;
        model = np.linalg.lstsq(G1,d)[0];
        resid = np.linalg.lstsq(G1,d)[1];
#       print result,resid;
#       sys.exit();
        A0 = model[0];
        A1 = math.sqrt(model[1]**2 + model[2]**2);
        fi1 = math.atan2(model[2],model[1]);
        A2 = math.sqrt(model[3]**2 + model[4]**2); 
        fi2 = math.atan2(model[4],model[3]);
        A3 = math.sqrt(model[5]**2 + model[6]**2);
        fi3 = math.atan2(model[6],model[5])

	return(A0,A1,A2,A3,fi1,fi2,fi3);

inb = [];
inv = [];
inu = []; 
outn = sys.argv[1] + ".fit"
outf = open(outn,"w");
for l1 in open(sys.argv[1]):
	l2 = l1.rstrip().split();
	inb.append(float(l2[0]));
	inv.append(float(l2[1]));
	inu.append(float(l2[2]));

(A0,A1,A2,A3,fi1,fi2,fi3) = invert_1 ( inb,inv,inu);
fi1 = fi1*180./math.pi;
fi2 = fi2*180./math.pi;
fi3 = fi3*180./math.pi;
print A0,A1,A2,A3,fi1,fi2,fi3;
for i in range (361):
#	v0 = A0 + A1*math.sin(i + fi1) + A2*math.sin(2*i + fi2);
	v1 = A0;
	v2 = A1*math.sin(i*math.pi/180. + fi1*math.pi/180.);
	v3 = A2*math.sin(2*i*math.pi/180. + fi2*math.pi/180.);
        v4 = A3*math.sin(4*i*math.pi/180. + fi3*math.pi/180.);
	v0 =v1 + v2 + v3 + v4;
	outf.write("%g %g %g %g %g %g\n" % (i,v0,v1,v2,v3,v4));
outf.close();
	




