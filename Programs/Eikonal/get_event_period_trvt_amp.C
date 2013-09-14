#define MAIN
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/mysac64.h"
using namespace std;
/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
  /*----------------------------------------------------------------------------
    ----------------------------------------------------------------------------*/
{
  FILE *fsac;
  /*..........................................................................*/
  if((fsac=fopen(fname, "rb")) == NULL) return NULL;
  
  if ( !SHD ) SHD = &SAC_HEADER;
  
  fread(SHD,sizeof(SAC_HD),1,fsac);
  
  if ( SHD->npts > nmax ) {
    fprintf(stderr,
	    "ATTENTION !!! in the file %s npts exceeds limit  %d",fname,nmax);
    SHD->npts = nmax;
  }
  
  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  
  fclose (fsac);
  
  /*-------------  calculate from t0  ----------------*/
  {
    int eh, em ,i;
    float fes;
    char koo[9];
    
    for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
    koo[8] = '\0';
    SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
      SHD->nzsec + SHD->nzmsec*.001;
    
    sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);
    
    SHD->o  -= (eh*3600. + em*60. + fes);
    /*-------------------------------------------*/}
  
  return SHD;
}

SAC_HD shd_cor;

void get_phvel(char filename[300],double per,double *phvel_out,double *amp_temp2)
{
  *phvel_out=0;
  FILE *f1;
  double cp,ap_1,ap_2,phvel_1,phvel_2,snr,wd,gv_1,gv_2,tri;
  double amp1,amp2;
  int k,i;
  if((f1=fopen(filename,"r"))==NULL)
    {
      return;
    }
  ap_1=0;
  for(;;)
    {
      if(fscanf(f1,"%d %lf %lf %lf %lf %lf %lf %lf",&k,&cp,&ap_2,&gv_2,&phvel_2,&amp2,&tri,&wd)==EOF)
	break;
      if(ap_2>per)
	{
	  if(ap_1==0) 
	    {
	      fclose(f1);
	      return;
	    }
	  *phvel_out=(phvel_2-phvel_1)/(ap_2-ap_1)*(per-ap_1)+phvel_1;
	  *amp_temp2=(amp2-amp1)/(ap_2-ap_1)*(per-ap_1)+amp1;
	  //*phvel_out=(gv_2-gv_1)/(ap_2-ap_1)*(per-ap_1)+gv_1;
	  break;
	}
      ap_1=ap_2;
      phvel_1=phvel_2;
      amp1=amp2;
      //gv_1=gv_2;
      i++;
    }
  fclose(f1);
}



int main (int na,char *arg[])
{
  
  FILE *ff,*f2,*f3;
  char sta1[10],sta2[10],channel[10];
  char filename[300],filename_snr[300],path_month[300],month_filename_snr[300],disp_filename[300];
  double dis,per,cri_snr;
  double snr_per,snr,snr_temp,per_temp;
  double amp, amp_temp,amp_temp2;
  double phvel_out,phvel,time;
  double snr_temp2,snr_trash;
  int i;
  
  double dist_cri;
  if(na!=6)
    {
      cout<<"usage:get_event_period_trvt_amp period event_station_dist.lst cri_snr dist_cri(km) path"<<endl;
      return 1;
    }
  per=atoi(arg[1]);
  cri_snr=atof(arg[3]);
  dist_cri=atof(arg[4]);
  if((ff = fopen(arg[2], "r"))==NULL) {
    printf("cannot open %s\n",arg[2]);
    exit(1);
  }
  cout<<"test1"<<endl;
  long event;
  char station[10];
  char buff1[300];
  double lon,lat;
  //  sprintf(buff1,"data_%ss_snr_%s_dist_%s.txt",arg[1],arg[4],arg[5]);
  //f3=fopen(buff1,"w");
  sprintf(buff1,"if [ ! -d %ssec_%ssnr_%sdis ]; then mkdir %ssec_%ssnr_%sdis; fi",arg[1],arg[3],arg[4],arg[1],arg[3],arg[4]);
  system(buff1);
  long event_old=0;
  for(;;)
    {    
      if(fscanf(ff,"%ld %s %lf %lf %lf",&event,&station[0],&dis,&lon,&lat)==EOF)
	{
	  if(event_old!=0)
	    fclose(f3);
	  break;
	}
      cout<<event<<" "<<station<<endl;
      if(dis<dist_cri)
	continue;
      
      if(event!=event_old)
	{
	  if(event_old!=0)
	    fclose(f3);

	  sprintf(buff1,"%ssec_%ssnr_%sdis/%ld.ph.txt",arg[1],arg[3],arg[4],event);
	  cout<<buff1<<endl;
	  f3=fopen(buff1,"w");
	  event_old=event;
	}

      phvel=0;
      sprintf(filename,"%s%ld/%ld.%s.LHZ.sac",arg[5],event,event,station);
      strcpy(filename_snr,filename);
      strcat(filename_snr,"_snr_rms.txt");
      strcpy(disp_filename,filename);
      strcat(disp_filename,"_1_DISP.0");
      
      if((f2=fopen(filename_snr,"r"))==NULL)
	{
	  cout<<" snr_err "<<filename_snr<<endl;
	  continue;
	}
      snr_per=0;
      snr=0;
      snr_temp=0;
      for(;;)
	{
	  if(fscanf(f2,"%lf %lf %lf",&snr_per,&snr,&amp)==EOF)
	    {
	      cout<<" snr_err"<<endl;
	      break;		    
	    }
	  if(snr_per>per)
	    {
	      snr_temp2=(snr-snr_temp)/(snr_per-per_temp)*(per-per_temp)+snr_temp;
	      if(snr_temp2>cri_snr)
		{
		  
		  get_phvel(disp_filename,per,&phvel_out,&amp_temp2);
		  if(phvel_out==0)
		    cout<<" phvel_err";
		  else
		    {
		      //		      amp_temp2=(amp-amp_temp)/(snr_per-per_temp)*(per-per_temp)+amp_temp;
		      cout<<" "<<phvel_out;
		      phvel=phvel_out;
		      time=dis/phvel_out;
		    }
		}
	      else
		cout<<" snr_small";
	      break;
	    }
	  snr_temp=snr;
	  per_temp=snr_per;
	  //amp_temp=amp;
	}
      fclose(f2);
      cout<<endl;
      if(phvel!=0)
	fprintf(f3,"%lf %lf %lf %lf %lf\n",lon,lat,time,phvel,amp_temp2);
    }
  

}
