#!/usr/bin/python
#this is used to generate the mhr input file, pathfile
import string
import sys

#fin = "/media/JAPAN/jiayi/WC_cv/staCEA_IRIS_90to110.lst"
fin = "./station.lst"
fine = "./event.loc"
stanm=[];stla=[];stlo=[]
for line in open(fin):
	l=line.rstrip().split()
	stanm.append(l[0])
	stlo.append(float(l[1]))
	stla.append(float(l[2]))
N=len(stanm)
evnm=[];evla=[];evlo=[]
for line in open(fine):
        l=line.rstrip().split()
        evnm.append(l[0])
        evlo.append(float(l[1]))
        evla.append(float(l[2]))

Ne=len(evnm)
fout=open("pathfile","w")
for i in range(Ne):
	ev=evnm[i]
	for j in range(N):
		sta=stanm[j]
		fout.write("%5d%5d%16s%8s%10.4f%10.4f%10.4f%10.4f\n"%(i+1,j+1,ev,sta,evla[i],evlo[i],stla[j],stlo[j]))	

fout.close()
