#here is the script I used to ask for Receiver Function data. You only need the sentence 'taup ...' which is highlighted in Bold
#only write events whose distance satisfeies some criteria
#write out two file: one is the email for requesting data; Another one is used to change the SAC header, by default, the data is named after its beginning time.
#STATIONS HAVE DIFFERENT NUMBER OF EVENTS, EVENTS WITHIN THE RUNING PERIOD OF EACH STATION
import time
import sys
import string
import datetime
import os
sys.path.append('/home/jiayi/progs/jy/python')
import get_dist
from copy import deepcopy
if(len(sys.argv)!=4):
    print "Usage python xx.py 1]event_iris.lst 2]station.lst 3]title"
    sys.exit()
#====parameters========
Ncri=10
tinter=datetime.timedelta(seconds=400)
mag_criSm = 5.4 #  mag_criSm<=mag<=mag_criBg
mag_criBg = 5.7
deg_min=30.
deg_max=120.
Toedge = 20.
c_stla=33.
c_stlo=88.
tob=60 #time before t0
toe=90 #time after t0
cmp="BH?"
samplesac="/home/jiayi/Tool/sac1.sac"
fevt=sys.argv[1]
fsta=sys.argv[2]
title=sys.argv[3]
stlo=[];stla=[];ntnm=[];stnm=[];st_b=[];st_e=[];
Nst=0;
#======station_info.lst, in iris format: NETWORK STATION STARTTIME ENDTIME  LAT LON SITE =====================
t_early=datetime.datetime(2012,12,12)
t_late=datetime.datetime(1000,12,12)
for line in open(fsta): 
    l=line.rstrip().split()
    ntnm.append(l[0])
    stnm.append(l[1])
    stla.append(float(l[6]))
    stlo.append(float(l[7]))
    date1=l[2].split("-")
     date2=l[4].split("-")
    time1=l[3].split(":")
    time2=l[5].split(":")
    tb=datetime.datetime(int(date1[0]),int(date1[1]),int(date1[2]),int(time1[0]),int(time1[1]),int(time1[2]))
    te=datetime.datetime(int(date2[0]),int(date2[1]),int(date2[2]),int(time2[0]),int(time2[1]),int(time2[2]))
    st_b.append(tb)
    st_e.append(te)
    if tb<t_early and tb > datetime.datetime(1990,1,1):
        t_early=tb;
    if te> t_late:
        t_late=te;
    Nst=Nst+1
print t_early, "~", t_late
#=====read event file, iris format============ DEFAULT: The event time is decending!!!!!
def not_number(s):
   try:
      float(s)
      return False
   except ValueError:
      return True
evla=[];evlo=[];evdep=[];evmag=[];evtime=[];
Nev=0;iindex=0;
indexnm=[];indexloc=[]; #used to search the location of a given year.minth
for line in open(fevt):
    l1=line.rstrip().split( )
    if(len(l1)!=10): #problem in length of line
#    print line;
    continue
    if not_number(l1[5][:-1]): #lack of depth information
#    print line,l1[5]
    continue;
    tla=float(l1[3][:-1]);tlo=float(l1[4][:-1]);
    tdist=get_dist.get_dist(tlo,tla,c_stlo,c_stla)
    tmag = float(l1[9])
    if tmag < mag_criSm or tmag > mag_criBg or tdist < 111*(deg_min-Toedge) or tdist > 111*(deg_max+Toedge): #throw some evt according to its dist to station centra
    continue;
    ya=l1[1][0:4]
    mon=l1[1][5:7]
    da=l1[1][8:10]
    hr=l1[2][0:2]
    min=l1[2][3:5]
    sec=l1[2][6:8]
    evla.append(tla);evlo.append(tlo);
    evdep.append(float(l1[5][:-1]))
    evmag.append(tmag)
    evtime.append(datetime.datetime(int(ya),int(mon),int(da),int(hr),int(min),int(sec)))
    #===== keep the index======
    yamo=int(str(ya)+str(mon))
    if yamo not in indexnm:
      indexnm.append(yamo)
    indexloc.append(Nev)
    iindex=iindex+1
    Nev=Nev+1;
print "Total events:%d stations:%d\n"%(Nev,Nst)
#print indexnm,"\n",indexloc
#print len(indexnm)
#====================================================
def find_indexb(timeb,tindexnm,tindexloc,nindex):
    if timeb<=tindexnm[nindex-1]:
        loc_b=tindexloc[nindex-1]
        return loc_b,nindex-1;
    for i in range(nindex):
        if tindexnm[i]<timeb:
        loc_b=tindexloc[i]-1;
        return loc_b,i
    return -1,-1;
def find_indexe(timee,tindexnm,tindexloc,nindex):   
    if timee>=tindexnm[0]:
        loc_e=tindexloc[0]
        return loc_e,0
    for i in range(nindex-1,-1,-1): #nindex-1~0
        if tindexnm[i]>timee:
        loc_e=tindexloc[i+1];
        return loc_e,i
    return -1,-1;

os.system("mkdir sacInfo\n")
for ist in range(Nst):
    styamob=int("%04d%02d"%(int(str(st_b[ist].year)),int(str(st_b[ist].month)))) #here the b and e are named in favor of time, not order in the sequence.
    styamoe=int("%04d%02d"%(int(str(st_e[ist].year)),int(str(st_e[ist].month))))
    Ne,ne=find_indexe(styamoe,indexnm,indexloc,iindex)
    Nb,nb=find_indexb(styamob,indexnm,indexloc,iindex)
    if Nb<0 or Ne <0:
        print "####problem! Index can't be found for station %s### %d %d"%(stnm[ist],styamob,styamoe)
        print Nb,Ne
        sys.exit()
#    print st_b[ist],"vs",st_e[ist]
#    print styamob,indexnm[nb],indexloc[nb],Nb,evtime[Nb],evtime[Nb+1],evtime[Nb-1],"\t",styamoe,indexnm[ne],indexloc[ne],Ne,evtime[Ne],evtime[Ne+1],evtime[Ne-1]
#    print evtime[Ne+4]
    print "wrkong on station %s %d, %d~%d,%d~%d,   %s~%s "%(stnm[ist],ist,styamob,styamoe,Ne,Nb,evtime[Nb],evtime[Ne])
    if (Nb-Ne+1<Ncri):
        print "###lack of evts, skip!!"
        continue;
   
    #=======use Taup to get t0(first arrival)==
    os.system("rm -f temp_taup.txt\n")
    cout=0;
    for iev in range(Ne,Nb+1):
        if evtime[iev]<st_b[ist] or evtime[iev]+tinter>st_e[ist]: #400s--travel time from evt to sta (30 deg)
            continue;
        dist=get_dist.get_dist(stlo[ist],stla[ist],evlo[iev],evla[iev])
        if dist < deg_min*111 or dist > deg_max*111:
            continue;
        string = "taup_time -mod prem -h %f -ph P,Pdiff -km %f | awk '{if(NR==6){printf\"%d %d\";print $0}}' >> temp_taup.txt"%(evdep[iev],dist,iev,ist);
        os.system(string);
        cout=cout+1;
    if cout<10:
        print "####lack of events after time/dist criteria. skip!!!"
        continue;
    #==============write mail and ch_header macro file
    mailnm="sacInfo/email_to_iris_%s_%s.txt"%(ntnm[ist],stnm[ist])
    fmail=open(mailnm,"w")
#    fhead=open("sacInfo/SAC_ch_%s_%s.txt"%(ntnm[ist],stnm[ist]),"w")   
    fhead=open("sacInfo/SAC_ch_%s_%s.txt"%(ntnm[ist],stnm[ist]),"a") # remember to remove old files first!!   
    fmail.write(".NAME Jiayi Xie\n.INST University of Colorado\n.MAIL Boulder,CO\n.EMAIL jiayi.seis@gmail.com\n.PHONE\n.FAX\n.MEDIA: Electronic(FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL %s_%s_%s\n.QUALITY B\n.END\n"%(title,ntnm[ist],stnm[ist]))   
    for line in open("temp_taup.txt"):
        l=line.rstrip().split()
        evid=int(l[0])
        stid=int(l[1])
        p=float(l[6]) # ray parameter??
        tstnm=stnm[stid]
        tntnm=ntnm[stid]
        t_o=evtime[evid]
        dt=datetime.timedelta(seconds=int(float(l[5]))) #ingore the time smaller than 1 sec
        t0=t_o+dt #the first p arrival time(GMT time)
        dt1=datetime.timedelta(seconds=tob)
        dt2=datetime.timedelta(seconds=toe)
        b=t0-dt1;
        e=t0+dt2;
        fmail.write("%-5s %5s %4d %02d %02d %02d %02d %02d.00 %4d %02d %02d %02d %02d %02d.00 1 %s\n"%(tstnm,tntnm,b.year,b.month,b.day,b.hour,b.minute,b.second,e.year,e.month,e.day,e.hour,e.minute,e.second,cmp))
        fhead.write("r %s\nr %4d.%03d.%02d.%02d.%02d*\n ch O GMT %4d %03d %02d %02d %02d 000 evlo %10.4f evla %10.4f evdp %f mag %f t0 %d user0 %f\n wh over\n"%(samplesac,b.year,b.timetuple().tm_yday,b.hour,b.minute,b.second,t_o.year,t_o.timetuple().tm_yday,t_o.hour,t_o.minute,t_o.second,evlo[evid],evla[evid],evdep[evid],evmag[evid],tob,p))
    fmail.close()
    fhead.close()
    print time.ctime()
    #=====send out the email ==========
    os.system("mail BREQ_FAST@iris.washington.edu < %s\n"%(mailnm))
    #==================
