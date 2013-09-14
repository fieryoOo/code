#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int calc_avg_daynum(char *fname,int dayrange)
{
FILE *f1,*f2;
char out_name[300];
int N=3000,NN=800;
int i,ir,Nr,ip[NN],Np,maxdaynum,daynum[N];
float day_num;
double noisesqr[N],noiserms[N],sigamp[N],normdamp[N];
maxdaynum=0;
if((f2 = fopen(fname, "rb")) == NULL) {
    printf("could not open data file to read: %s \n", fname);
    fclose(f2);
    return 0;
  }
for(i=0;i<N;i++)
  {
   if(fscanf(f2,"%lf %d %lf %lf %lf",&noisesqr[i],&daynum[i],&noiserms[i],&sigamp[i],&normdamp[i])==EOF)
   break;
   if(daynum[i]>maxdaynum)
   maxdaynum=daynum[i];
  }
fclose(f2);
Np=i;
Nr=(int)(maxdaynum-1)/dayrange+1;
double noise_sqr[NN],noise_rms[NN],sig_amp[NN],normd_amp[NN];
for(i=0;i<Nr;i++)
  {
   noise_sqr[i]=0;
   noise_rms[i]=0;
   sig_amp[i]=0;
   normd_amp[i]=0;
   ip[i]=0;
  }
for(i=0;i<Np;i++)
  {
   ir=(int)(daynum[i]-1)/dayrange;
   noise_sqr[ir]+=noisesqr[i];
   noise_rms[ir]+=noiserms[i];
   sig_amp[ir]+=sigamp[i];
   normd_amp[ir]+=normdamp[i];
   ip[ir]+=1;
  }
sprintf(out_name,"avgd_%s",fname);
f1=fopen(out_name,"w");
for(i=0;i<Nr;i++)
  {
   noise_sqr[i]=noise_sqr[i]/ip[i];
   noise_rms[i]=noise_rms[i]/ip[i];
   sig_amp[i]=sig_amp[i]/ip[i];
   normd_amp[i]=normd_amp[i]/ip[i];
   day_num=(i+0.5)*dayrange;
   fprintf(f1,"%lf %f %lf %lf %lf\n",noise_sqr[i],day_num,noise_rms[i],sig_amp[i],normd_amp[i]);
  }
fclose(f1);
return 1;
}

int main (int argc, char *argv[])
{
FILE *ff;
char filename[300];
int daynum_range;
if(argc != 3)
  {
   printf("Usage: calc_avg_daynum_range.C [avg_daynum_range] [input_file_lst]\n");
   exit(-1);
  }
if((ff = fopen(argv[2], "r"))==NULL) {
   printf("Cannot open input file list %s.\n",argv[2]);
   exit(1);
  }
daynum_range=(int)atof(argv[1]);
for(;;) {
   if(fscanf(ff,"%s",&(filename[0]))==EOF)
   break;
   calc_avg_daynum(filename,daynum_range);
  }
fclose(ff);
return 1;
}
