#define MAIN
#include "/home/tianye/code/Programs/head/koftan.h"
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64.h"
#include <iostream>

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

int main(int na, char *arg[])
{ 
  if(na!=3)
    {
      cout<<"usage:get_dist event_lst station.lst"<<endl;
      return 0;
    }
  FILE *ff,*f1,*f2;
  char file_name[400];
  float sig[100000];
  
  if((f1 = fopen(arg[1], "r"))==NULL) {
    printf("cannot open %s\n",arg[1]);
    exit(1);
  }
  if((f2 = fopen(arg[2], "r"))==NULL) {
    printf("cannot open %s\n",arg[2]);
    exit(1);
  } 
  
  ff=fopen("sta_dist.lst","w");
  int N=3000;
  char event[N][20];
  char name[N][6];
  double lat[N],lon[N];
  int i,j,k,ii,jj,iev,ist;
  for(i=0;i<N;i++)
    {
      if(fscanf(f1,"%s",event[i])==EOF) break;
      cout<<event[i]<<endl;
    }
  fclose(f1);
  cout<<"number of events read "<<i<<endl;
  iev=i;
  for(i=0;i<N;i++)
    {
      if(fscanf(f2,"%s %lf %lf",name[i],&lon[i],&lat[i])==EOF) break;
      cout<<name[i]<<endl;
    }
  fclose(f2);
  ist=i;
  cout<<"number of stations read "<<i<<endl;
  for(j=0;j<iev;j++)
    {
      for(k=0;k<ist;k++)
        {
	  sprintf(file_name,"/home/tianye/data_Eikonal/SAC_XR/%s/%s.%s.LHZ.sac",event[j],event[j],name[k]);
	  
	  if(!read_sac (file_name, sig, &shd_cor, 100000 ))
            {
              cout<<file_name<<" not found!!"<<endl;
              continue;
            }
	  //	  cout<<name[j]<<" "<<name[k]<<" "<<shd_cor.dist<<endl;
	  fprintf(ff,"%s %s %lf %lf %lf\n",event[j],name[k],shd_cor.dist,shd_cor.stlo+360,shd_cor.stla);
	}
    }
  fclose(ff);
  return 0;
}
