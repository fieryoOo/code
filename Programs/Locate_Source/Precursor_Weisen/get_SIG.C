#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "/home/weisen/progs/NOISE_CODA/HEADFILE/mysac.h"

using namespace std;

#define N 200000


SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)

{
 FILE *fsac;
/*..........................................................................*/
        if((fsac = fopen(fname, "rb")) == NULL) {
          printf("could not open sac file %s to write\n",fname);
          return NULL;
        }

        if ( !fsac )
        {
          /*fprintf(stderr,"file %s not find\n", fname);*/
         return NULL;
        }

//        if ( !SHD ) SHD = &SAC_HEADER;

         fread(SHD,sizeof(SAC_HD),1,fsac);

         if ( SHD->npts > nmax )
         {
          fprintf(stderr,
           "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);

          SHD->npts = nmax;
         }

         fread(sig,sizeof(float),(int)(SHD->npts),fsac);

        fclose (fsac);

   /*-------------  calcule de t0  ----------------*/
   {
        int eh, em ,i;
        float fes;
        char koo[9];

        for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
        koo[8] = 0;

        SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
         SHD->nzsec + SHD->nzmsec*.001;

        sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

        SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

        return SHD;
}

        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        if((fsac = fopen(fname, "wb"))==NULL) {
           printf("could not open sac file to write\n");
           exit(1);
        }

//        if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];

   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

         fwrite(SHD,sizeof(SAC_HD),1,fsac);

         fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


        fclose (fsac);
}

////////////////////////////////////////////////////////////////////

//DEFINE N 20000 

double get_dist2(double lat1, double lon1, double lat2, double lon2, double * afa1, double * afa2)
{
double pi;
pi = 4.0*atan(1.0);
double cva = 6378.137;
double cvb = 6356.7523142;
double f = 1/298.257223563;
double tempafa1, tempafa2;
double L = 0.00;
L = lon2-lon1;
if (L > 180.000)  L =360.000000 - L;
if (L < -180.000) L =  360.000 - abs(L);
L = fabs(L);
double U1 = 0;
U1 = atan((1-f)*tan(lat1/180*pi));
double U2 = 0;
U2 = atan((1-f)*tan(lat2/180*pi));
double cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
L = L*pi/180;
double numda = L;
numda1 = numda;
do {
  numda = numda1;
  cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))+ (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))*(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) );
  cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda);
  cv = atan2(cv1,cv2);
  cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv);
  cv4 = 1 - cv3*cv3;
  cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4;
  cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
  numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
 } while (fabs(numda - numda1) > 0.0000000001);
double mius, cvA, cvB, deltacv,s;
mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb);
cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));

cvB = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)));
deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));


s = cvb * cvA *(cv - deltacv);

tempafa1 = atan2(cos(U2)*sin(numda1), cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(numda));
tempafa2 = atan2(cos(U1)*sin(numda1), cos(U1)*sin(U2)*cos(numda) - sin(U1)*cos(U2));

if (lon1-lon2> 0 && lon1 -lon2<180) { tempafa1 = 2*pi - tempafa1; tempafa2 = pi-tempafa2; }
else { tempafa2 = pi+tempafa2; }

* afa1 = tempafa1/pi*180;
* afa2 = tempafa2/pi*180;
return s;
}


int main ( int argn , char *arg[] ) {

 int i,j,k;
 SAC_HD shd;
 float sig[N],sig1[N],lat0,lon0,lat1,lon1,lat2,lon2,lat[2000],lon[2000];
 float ddist,dtime1,dtime2,mag1,mag2,snr;
 int flag1,flag2,nst,t1,t2,per;
 double dist,azi,bazi;
 FILE *ff,*ff1;
 char sta1[6],sta2[6],tempst[6],filename[300],filename1[300],stations[2000][6];

 flag1 = 0;
 flag2 = 0;


 if (argn != 3) {
   cout<<"usage: envelope.lst station.lst"<<endl;
   return 0;
   }
 
 if ((ff = fopen(arg[2],"r"))==NULL) {
   fprintf(stderr,"cannot open file %s\n",arg[2]);
   return 0;
   }
 for (i=0;;i++) {
   if (fscanf(ff,"%s %g %g",stations[i],&(lon[i]), &(lat[i]))==EOF) break;
   }
 fclose(ff);
 nst = i;


 if ( (ff1 = fopen(arg[1],"r"))==NULL) {
   fprintf(stderr,"cannot open file %s\n",arg[1]);
   return 0;
   }
ff = fopen("cv_temp.result","w");

for (k=0;;k++) {
 if (fscanf(ff1, "%s %s %s", filename, sta1, sta2)==EOF) break;
 sprintf(filename1,"cut_%s\0",filename);
 if (read_sac(filename,sig,&shd,N)==0) {
   fprintf(stderr,"cannot read sac file %s\n",arg[1]);
   return 0;
   }
 flag1 = 0;
 flag2 = 0;
 for (i=0;i<nst;i++) {
   if (strcmp(sta1,stations[i])==0) {
	lat1 = lat[i];
	lon1 = lon[i];
	flag1 = 1;
	}
   if (strcmp(sta2,stations[i])==0) {
        lat2 = lat[i];
        lon2 = lon[i];
        flag2 = 1;
        }
   if (flag1 == 1 && flag2 == 1 ) break;
   }
 if (flag1 ==0 || flag2 ==0 ) {
   fprintf(stderr,"no such stations!!!! %s %s\n",sta1,sta2);
   break;
   }
 lat0 = 33;
 lon0 = 131.5;
 ddist = get_dist2( (double)lat0, (double)lon0, (double)lat1, (double)lon1, &azi, &bazi ) - get_dist2( (double)lat0, (double)lon0, (double)lat2, (double)lon2, &azi, &bazi );
 dtime1 = -ddist/2.2;
 dtime2 = -ddist/3.8;
 if ((int)dtime1>(int)dtime2)  {
    t1 = (int)dtime2;
    t2 = (int)dtime1;
    }
 else {
    t1 = (int)dtime1;
    t2 = (int)dtime2;
    }
 mag1 = 0;
 if (   (int)dtime1 < 5000 && (int)dtime1 > -5000    ) {
//	mag = sig[(int)dtime+(shd.npts-1)/2];
	for (i=(int)t1+(shd.npts-1)/2;i<=(int)t2+(shd.npts-1)/2;i++)
//		if (mag1 < sig[i] ) mag1 = sig[i];           // mag1 is the largest number;
		sig1[i-((int)t1+(shd.npts-1)/2)] = sig[i];
	}
// shd.npts = i-((int)t1+(shd.npts-1)/2);
/* else 
	mag1 = 0;

 mag2 = 0;
 for (i=0;i<1000;i++) {
	mag2 = mag2 + sig[i];
	}
 for (i=9000;i<10000;i++) {
	mag2 = mag2 + sig[i];
	}
 mag2 = mag2/2000;
 if (mag2 == 0) mag2 = 0.01; 
 snr = mag1/mag2;
 fprintf(ff,"%s %d %s %s %g %g %g\n",filename, per, sta1, sta2, snr,mag1,mag2);
*/
 shd.b = (int)t1;
 shd.e = (int)t2;
// shd.npts = i-((int)t1+(shd.npts-1)/2);
 shd.npts = (int)t2- (int)t1;
 write_sac(filename1,sig1,&shd);
}

fclose(ff);
return 0;
}
