#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <sys/time.h>
#include "/home/tianye/code/Programs/head/64_mysac.h"

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
if (L < -180.000) L =  360.000 - fabs(L);
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

 int i,j,k,ii;
 SAC_HD shd;
 float sig[N],lat,lon,lat1,lon1,lat2,lon2,templat[3000],templon[3000];
 float ddist,statime,dtime,noise_p,noise_n,vel0,vel1,vel2,mag,tempper,tempvel,veles[10000][2000];
 int per,ista1,ista2,nsta;
 double dist,azi,bazi;
 FILE *ff,*ff1,*ffout,*fvel;
 char st1[6],st2[6],tempst[6],filename[30], stations[3000][10],outname[100],event[30],veldir[300],velfile[100];
 char tmp[300];
 time_t t_start, t_end;
 double t_diff;
 timespec t_s;

 if (argn != 5) {
   cout<<"usage: envelope.lst station.lst event.dat veldir"<<endl;
   return 0;
   }
strcpy (veldir, arg[4]);
fprintf(stderr,"%s\n",veldir);
 if ((ff = fopen(arg[2],"r"))==NULL) {
   fprintf(stderr,"cannot open file %s\n",arg[2]);
   return 0;
   }
 for (i=0;;i++) {
	if (fscanf(ff,"%s %g %g",stations[i],&(templon[i]), &(templat[i]))==EOF) break;
	}
 nsta = i;
 fclose(ff);
 if ((ff=fopen(arg[3],"r"))==NULL) {
	fprintf(stderr,"cannot open file %s\n",arg[3]);
	return 0;
	}
 if ((ff1=fopen(arg[1],"r"))==NULL) {
   fprintf(stderr,"cannot open file %s\n",arg[1]);
   return 0;
   }
 fscanf(ff1,"%s %d %s %s",filename,&per,st1,st2);
 rewind(ff1);
//fprintf(stderr,"     !!! %s\n",arg[1]);

 for (i=0;;i++) {
    if (fscanf(ff,"%s %g %g",event,&lon,&lat)==EOF) break;
    printf("Reading velocities for the %dth event...\n",i);
    for (j=0;j<nsta;j++){
        sprintf(velfile,"%s/%s/%s_%s.dat\0",veldir,event,event,stations[j]);
//printf("j: %d  file: %s\n",j,velfile);
        if (! (fvel = fopen(velfile,"r"))) {
                fprintf(stderr,"no file: %s\n",velfile);
                veles[i][j] = 2.6;
                return 0;//continue;
                }
        for (k=0;;k++) {
                if (fscanf(fvel,"%g %g",&tempper,&tempvel)==EOF) {
                        fprintf(stderr,"no per: %d %s\n",per,velfile);
                        veles[i][j] = 2.6;
                        break;
                        }
                if (tempper == float(per)) {
                        veles[i][j] = tempvel;
                        break;
                        }
                }
           fclose(fvel);
          }
     }
//rewind(ff);

for (k=0;;k++) {
 if (fscanf(ff1,"%s %d %s %s",filename,&per,st1,st2)==EOF) break;

// cout<<st1<<" "<<st2<<endl;
// cout<<filename;
 for (i=0;i<nsta;i++)
   if (strcmp(st1,stations[i])==0) {
	lat1 = templat[i];
	lon1 = templon[i];
	break;
 	}
 ista1 = i;
 if (ista1 == nsta) {
   fprintf(stderr,"Can not find station %s in station.lst!\n",st1);
   continue;
   }
 
 for (i=0;i<nsta;i++)
   if (strcmp(st2,stations[i])==0) {
        lat2 = templat[i];
        lon2 = templon[i];
	break;
        }
 ista2 = i;
 if (ista2 == nsta) {
   fprintf(stderr,"Can not find station %s in station.lst!\n",st2);
   continue;
   } 

 if (read_sac(filename,sig,&shd,N)==0) {
   fprintf(stderr,"cannot read sac file %s\n",filename);
   return 0;
   }
 sprintf(outname,"%d_%s_%s.result",per,st1,st2);
 ffout = fopen(outname,"w");

 sprintf(velfile,"%s/%s/%s_%s.dat\0",veldir,st1,st1,st2);
 if ((fvel = fopen(velfile,"r"))==NULL) {
         fprintf(stderr,"no file: %s\n",velfile);
         vel0 = 2.6;
         }
 else {
//    fgets(tmp,300,fvel);
    for (j=0;;j++) {
         if (fscanf(fvel,"%g %g",&tempper,&tempvel)==EOF) {
                 fprintf(stderr,"no per: %d %s\n",per,velfile);
                 vel0 = 2.6;
                 break;
                 }
         if (tempper == float(per)) {
                 vel0 = tempvel;
                 break;
                 }
         }
     fclose(fvel);
    }
/* time(&t_start);
 t_s.tv_sec = 0;
 t_s.tv_nsec = 0;
 clock_settime(CLOCK_PROCESS_CPUTIME_ID, &t_s);
// clock_settime(CLOCK_REALTIME, &t_s); */
// ff=fopen(arg[3],"r");
 rewind(ff);
 for (i=0;;i++) {
	if (fscanf(ff,"%s %g %g",event,&lon,&lat)==EOF) break;

//printf("AAA  i: %d\n",i);
	vel1 = veles[i][ista1];
	vel2 = veles[i][ista2];
/*        time (&t_end);
	t_diff = difftime(t_end, t_start);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t_s);
	printf("Real Time taken is: %.2f sec\n",t_diff);
        printf("CPU Time taken is: %d.%d sec\n", t_s.tv_sec, t_s.tv_nsec);
//printf("BBB  i: %d\n",i);      */
	statime = get_dist2( (double)lat1, (double)lon1, (double)lat2, (double)lon2, &azi, &bazi )/vel0;
	dtime = get_dist2( (double)lat, (double)lon, (double)lat2, (double)lon2, &azi, &bazi )/vel2 - get_dist2( (double)lat, (double)lon, (double)lat1, (double)lon1, &azi, &bazi )/vel1;
	noise_p=0; noise_n=0;
	for(ii=2000;ii<2500;ii++) { noise_p += sig[ii+(shd.npts-1)/2]; noise_n += sig[-ii+(shd.npts-1)/2]; }
        noise_p=noise_p/500; noise_n=noise_n/500;
	if (   (int)dtime < -3000 || (int)dtime > 3000 )
		mag = (noise_p+noise_n)/2.;
	else if ( ((int)dtime < statime*1.1+80 && (int)dtime > statime*0.9-80) || ((int)dtime > -statime*1.1-80 && (int)dtime < -statime*0.9+80) )
		mag = (noise_p+noise_n)/2.;
	else if ( (int)dtime > 0 ) {
		mag = sig[(int)dtime+(shd.npts-1)/2];
                if ( mag < 5*noise_p ) mag = noise_p;
		}
	else {
		mag = sig[(int)dtime+(shd.npts-1)/2];
                if ( mag < 5*noise_n ) mag = noise_n;
		}
	fprintf(ffout,"%g %g %g\n",lon, lat, mag);
	}
// fclose(ff);
 fprintf(stderr,"%dth pair is ok!!\n",k);
 fclose(ffout);
 } //k
fclose(ff);
fclose(ff1);
return 0;
}
