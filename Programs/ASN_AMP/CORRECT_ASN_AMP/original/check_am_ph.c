//this program check all the *am and *ph files in the directory $month/$month/5to150/$day
//it will creat a directory rec/ and put the result in it
//the file it creats is SAC file, which shows the data in 1 year.
//If there is am file and ph file in $d day, the $d point in SAC file would be 1,otherwise 0
//it calls to enter the year you're trying to process


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "/home/weisen/progs/NOISE_CODA/HEADFILE/sac_db.h"
#include "/home/weisen/progs/NOISE_CODA/HEADFILE/mysac.h"

//-----------------------------------------------

        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
 * ----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        //if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (int)ITIME;
        SHD->leven = (int)TRUE;

        SHD->lovrok = (int)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];

   for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

   fwrite(SHD,sizeof(SAC_HD),1,fsac);

   fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


   fclose (fsac);
}



/*------------------------------------------------------------------------*/
 int jday ( int y, int m, int d )
/*------------------------------------------------------------------------*/
{
  int jd = 0;
  int i;
 
  for ( i = 1; i < m; i++ )
    {
      if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
      else if (i==2)
        {
          if ( y == 4*(y/4) ) jd += 29;
          else jd += 28;
        }
      else jd += 30;
    }
 
  return jd + d;
}


//------------------------------------------------
int main(int na, char *arg[]) {

//int y;
if (na != 4)
{fprintf(stderr,"check_am_ph [year] [in_path] [output_path]\n");
exit(1);
}

int y;
char out_path[100],in_path[100];
sscanf(arg[1],"%d",&y);
sscanf(arg[2],"%s",in_path);
sscanf(arg[3],"%s",out_path);

printf("%d  %s\n",y,in_path);

FILE *ff,*TEMP;
char buff[300],staname[10],month[4],location[300],filename[300],buff1[300];
int n,mo,day,jd,nday,j;
float lon,lan;
float sig_am[366],sig_ph[366],sig[368];
SAC_HD shd;
//system("mkdir rec");
for (n=0;n<366;n++)
{
sig_am[n]=0;
sig_ph[n]=0;
sig[n]=0;
}
sig[366]=0;
sig[367]=0;
char sta_lst[100];
sprintf(sta_lst,"%s/station.lst",in_path);
ff=fopen(sta_lst,"r");

for (n=0;;n++) // read the first colume and the first station information
{    
  if (!fgets(buff,40,ff))
  break;
  
  for (j=0;j<366;j++)
  {
  sig_am[j]=0;
  sig_ph[j]=0;
  sig[j]=0;
  }
  sig[366]=0;
  sig[367]=0;
  
//fprintf(stderr,"station.number:%d  ",n);

  sscanf(&(buff[0]),"%s",staname);
  fprintf(stderr,"station.name:%s\n",staname);
  sscanf(&(buff[7]),"%g%g",&lon,&lan);
  fprintf(stderr,"lon:%g,lan:%g\n",lon,lan);

    for (mo=0;mo<12;mo++)
     {
     if (mo==0) {strcpy(month,"JAN"); nday=31; }
     if (mo==1) {
           strcpy(month,"FEB");
          //  fprintf(stderr,"%-4d,%d",y,4*(y/4)); 
           if (4*(y/4)!=y)
                nday=28;
           else nday=29; }
     if (mo==2) {strcpy(month,"MAR"); nday=31; }
     if (mo==3) {strcpy(month,"APR"); nday=30; }
     if (mo==4) {strcpy(month,"MAY"); nday=31; }
     if (mo==5) {strcpy(month,"JUN"); nday=30; }
     if (mo==6) {strcpy(month,"JUL"); nday=31; }
     if (mo==7) {strcpy(month,"AUG"); nday=31; }
     if (mo==8) {strcpy(month,"SEP"); nday=30; }
     if (mo==9) {strcpy(month,"OCT"); nday=31; }
     if (mo==10) {strcpy(month,"NOV"); nday=30; }
     if (mo==11) {strcpy(month,"DEC"); nday=31; }
     
     //fprintf(stderr,"month:%s,nday:%d",month,nday);
     
     for (day=0;day<nday;day++)
        {
        int dd=day;
        //fprintf(stderr,"dd:%d\n",dd);
        jd=jday(y,mo+1,dd+1);
//        sprintf(location,"%3s/5to150/%4d_%d_%d_0_0_0/ft_%s.LHZ.SAC.am",month,y,mo+1,dd+1,staname);
        sprintf(location,"%s/%4d_%d_%d_0_0_0/ft_%s.LHZ.SAC.am",in_path,y,mo+1,dd+1,staname);
        //fprintf(stderr,"location: %s\n",location);
        //scanf("%*1d");
        if ((TEMP = fopen( location,"rb"))!=NULL) {sig_am[jd-1]=1; fclose(TEMP);}
        //for (j=0;j<30;j++) location[j]='\0'; 
//        sprintf(location,"%3s/5to150/%4d_%d_%d_0_0_0/ft_%s.LHZ.SAC.ph",month,y,mo+1,dd+1,staname);
        sprintf(location,"%s/%4d_%d_%d_0_0_0/ft_%s.LHZ.SAC.ph",in_path,y,mo+1,dd+1,staname);
        //fprintf(stderr,"location: %s\n",location);
        if ((TEMP = fopen( location,"rb"))!=NULL) {sig_ph[jd-1]=1; fclose(TEMP);}
        else sig_ph[jd-1]=0;      
        //fprintf(stderr,"day:%d,am:%g,ph:%g\n", jd,sig_am[jd-1], sig_ph[jd-1]);
        //strcpy(location,"\0\0");
        //if (sig_am[jd-1]) scanf("%*1d");        
        sig[jd]=sig_am[jd-1]*sig_ph[jd-1];
     /*   if (!strcmp(staname,"SC14"))
             if (jd==19){
                fprintf(stderr,"%s\n",staname);
                scanf("%*1d");
                fprintf(stderr,"location: %s\n, %g,%g\n",location,sig_am[jd-1],sig_ph[jd-1]);
                scanf("%*1d");
                
                }*/
        } //day
     }//month 
  shd.npts=368;
  shd.delta=1;
  shd.stla=lan;
  shd.stlo=lon;
  shd.nzyear=y;
  shd.nzjday=0;
  shd.nzhour=0;
  //shd.nzmonth=1;
  strcpy(shd.kstnm,staname);
  sprintf(filename,"%s/%s_rec.SAC",out_path,staname);
  write_sac(filename,&(sig[0]),&shd);
//  system("mv *SAC rec/");
}  //for

fclose(ff);

}  //main



