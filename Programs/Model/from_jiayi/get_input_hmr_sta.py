#!/usr/bin/python
#this is used to generate the mhr input file, pathfile
import string
import sys

#fin = "/home/jiayi/work/SC/Info/station.lst"
#fin = "/media/CHINA/jiayi/WesternChina/sta.lst"
#fin = "/media/CHINA/jiayi/WesternChina/staSmallReg.lst"
fine = ""
fin = "/media/JAPAN/jiayi/WC_cv/staCEA_IRIS_90to110.lst"
stanm=[];stla=[];stlo=[]
for line in open(fin):
	l=line.rstrip().split()
	stanm.append(l[0])
	stlo.append(float(l[1]))
	stla.append(float(l[2]))

N=len(stanm)
fout=open("pathfile","w")
for i in range(N):
	sta1=stanm[i]
	for j in range(i+1,N):
		sta2=stanm[j]
		fout.write("%5d%5d %10s %10s%10.3f%10.3f%10.3f%10.3f\n"%(i+1,j+1,sta1,sta2,stla[i],stlo[i],stla[j],stlo[j]))	

fout.close()
