#define MAIN

#include <stdio.h>
#include "../inc/mysac64.h"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Function prototypes */

void dcommon_(int *len, float *amp,float *phase);
void dmultifft_(int *len,float *amp,float *phase, int *lag,float *seis_out, int *ns);

//void read_sac(char *name,char *stnam,float *stlat,float *stlon,int *n,float *sei);
void swapn(unsigned char *b, int N, int n);
//void write_cor(char *nam_f1,char *nam_corr,int *lag,float *corr,
//     char *stnam,float *stlat,float *stlon);


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




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
        fsac = fopen(fname, "wb");

        if ( !SHD ) SHD = &SAC_HEADER;


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




/*c/////////////////////////////////////////////////////////////////////////*/
 float sig_old[100000];
 SAC_HD shdamp1, shdph1, shdamp2, shdph2, shd_cor;

/*--------------------------------------------------------------------------*/
int do_cor( int lag, int iii, int kkk, int ns1, int ns2 )
/*--------------------------------------------------------------------------*/
{
  FILE *fst,*fev;
  if((fst=fopen("stalst","r"))==NULL)
    return 0;
  if((fev=fopen("evtlst","r"))==NULL)
    return 0;
  char sta_name[2000][10];
  //  double sta_lat[2000],sta_lon[2000];
  char ev_name[200][20];
  int ist,iev,ine, jsta1, jsta2, k,count;
  ist=0;
  for(;;)
    {
      if(fscanf(fst,"%s",sta_name[ist])==EOF)
	break;
      ist++;
    }
  fclose(fst);
  iev=0;
  for(;;)
    {
      if(fscanf(fev,"%s",ev_name[iev])==EOF)
	break;
      iev++;
    }
  fclose(fev);
  if(iev>200)
    {
      fprintf(stderr,"too many events %d\n",iev);
    }
  if(ns1==0 && ns2==0)
    {
      ns1=1;
      ns2=ist;
    }
  if(ns1<1)
    {
      ns1=1;
    }
  if(ns2>ist)
    {
      ns2=ist;
    }
  ns1=ns1-1;
  ns2=ns2-1;

  int len,ns,i; 
  // int iii,jjj,kkk;
  char name1_N[200],name1_E[200];
  char name2_N[200],name2_E[200]; 
  char filename[200];
  char amp_sac[200], phase_sac[200];
  float amp1[200][100000], phase1[200][100000];
  int yesno[200],len1[200];
  float amp[100000], phase[100000];
  float cor[100000];
  float seis_out[100000];
  double lat1,lat2,lon1,lon2,cmpaz;
  FILE *ff;
  char cp_am_ph[200];
  char buff[300];
  sprintf(buff,"if [ ! -d COR ]; then mkdir COR; fi");
  system(buff);

  //fprintf(stderr,"%d %d %s %s\n",ns1+1,ns2+1,sta_name[ns1],sta_name[ns2]);
  //for( jsta1 = 0; jsta1 < ist; jsta1++ ) 
  for( jsta1 = ns1; jsta1 <= ns2; jsta1++ ) 
    { 
      sprintf(buff,"if [ ! -d COR/%s ]; then mkdir COR/%s; fi",sta_name[jsta1],sta_name[jsta1]);
      system(buff);
    }
  //for( jsta1 = 0; jsta1 < ist; jsta1++ ) 
  for( jsta1 = ns1; jsta1 <= ns2; jsta1++ ) 
    { 
      for( ine = 0; ine < iev; ine++ ) 
	{
	  if(iii==0)
	    {
	      sprintf( amp_sac,   "%s/wt%s.LHZ.SAC.am", ev_name[ine],sta_name[jsta1] );
	      sprintf( phase_sac, "%s/wt%s.LHZ.SAC.ph", ev_name[ine],sta_name[jsta1] ); 
	    }
	  else if(iii==1)
	    {
	      sprintf( amp_sac,   "%s/wt%s.LHN.SAC.am", ev_name[ine],sta_name[jsta1] );
	      sprintf( phase_sac, "%s/wt%s.LHN.SAC.ph", ev_name[ine],sta_name[jsta1] );
	    }
	  else
	    {
	      sprintf( amp_sac,   "%s/wt%s.LHE.SAC.am", ev_name[ine],sta_name[jsta1] );
	      sprintf( phase_sac, "%s/wt%s.LHE.SAC.ph", ev_name[ine],sta_name[jsta1] );
	    }
	  yesno[ine]=0;
	  if( read_sac(amp_sac, amp1[ine], &shdamp1, 100000 )==NULL )
	    {
	      yesno[ine]=1;
	      continue;
	    }
	  if( read_sac(phase_sac, phase1[ine], &shdph1, 100000)== NULL )
	    {
	      //		  fprintf( stderr,"file %s did not found\n", phase_sac );
	      yesno[ine]=1;
	      continue;
	    }
	  len1[ine] = shdamp1.npts;
          printf("%s %s %d\n",amp_sac,phase_sac,iev);
	}
      for( jsta2 = (jsta1+0); jsta2 < ist; jsta2++ )
	{	
	  count=0;
	  if(iii==0)
	    {
	      sprintf(filename, "COR/%s/COR_%s_%s.SAC_Z",
		      sta_name[jsta1], sta_name[jsta1], sta_name[jsta2]);
	    }
	  else if(iii==1)
	    {
	      sprintf(filename, "COR/%s/COR_%s_%s.SAC_N",
		      sta_name[jsta1], sta_name[jsta1], sta_name[jsta2]);
	    }
	  else if(iii==2)
	    {
	      sprintf(filename, "COR/%s/COR_%s_%s.SAC_E",
		      sta_name[jsta1], sta_name[jsta1], sta_name[jsta2]);
	    }
	  if(kkk==0)
	    strcat(filename,"Z");
	  else if(kkk==1)
	    strcat(filename,"N");
	  else
	    strcat(filename,"E");
	  if(access(filename, F_OK) == 0)
	    {
	      fprintf(stderr,"%s exist!! skip\n",filename);
	      continue;
	    }
	  fprintf(stderr,"%s %d %s %d\n",sta_name[jsta1],jsta1,sta_name[jsta2],jsta2);
	  for( ine = 0; ine < iev; ine++ ) 
	    {
	      if(yesno[ine]!=0)
		continue;

	      // read amp and phase files and read into common memory
	      //if ( read_sac(amp_sac, amp, &shdamp1, 1000000 )==NULL )
	      //{
		  //		  fprintf( stderr,"file %s did not found\n", amp_sac );
	      //  continue;
		  //return 0;
	      //}
	      //if ( read_sac(phase_sac, phase, &shdph1, 1000000)== NULL )
	      //{
		  //		  fprintf( stderr,"file %s did not found\n", phase_sac );
	      //  continue;
		  //return 0;
	      //}
	      //	      len = len1[ine];
	      
	      //dcommon_( &len, amp1[ine], phase1[ine] ); // reads amp and phase files into common memory
	      
	      // compute correlation
	      
	      if(kkk==0)
		{
		  sprintf( amp_sac,   "%s/wt%s.LHZ.SAC.am", ev_name[ine],sta_name[jsta2] );
		  sprintf( phase_sac, "%s/wt%s.LHZ.SAC.ph", ev_name[ine],sta_name[jsta2] );
		}
	      else if(kkk==1)
		{
		  sprintf( amp_sac,   "%s/wt%s.LHN.SAC.am", ev_name[ine],sta_name[jsta2] );
		  sprintf( phase_sac, "%s/wt%s.LHN.SAC.ph", ev_name[ine],sta_name[jsta2] );
		}
	      else
		{
		  sprintf( amp_sac,   "%s/wt%s.LHE.SAC.am", ev_name[ine],sta_name[jsta2] );
		  sprintf( phase_sac, "%s/wt%s.LHE.SAC.ph", ev_name[ine],sta_name[jsta2] );
		}
              printf("%s %s %d\n",amp_sac,phase_sac,iev);
	      
	      
	      // get array of floats for amp and phase of first signal
	      
	      if ( read_sac(amp_sac, amp, &shdamp2, 100000) ==NULL ) 
		{
		  //		  fprintf(stderr,"file %s did not found\n", amp_sac );
		  continue;
		  //return 0;
		}
	      
	      if ( read_sac(phase_sac, phase, &shdph2, 100000)==NULL ) 
		{
		  //		  fprintf(stderr,"file %s did not found\n", phase_sac );
		  continue;
		  //		  return 0;
		}
	      len = len1[ine];
	      
	      dcommon_( &len, amp1[ine], phase1[ine] ); // reads amp and phase files into common memory
	      len = shdamp2.npts;
		      
	      //if(!check_info(sdb_N, sdb_E, ine, jsta1, jsta2,iii,kkk )) 
	      //		{
	      //  fprintf(stderr,"files incompatible\n");
	      //  return 0;
	      //}
	      
	      //else
	      //{
			  
	      dmultifft_(&len, amp, phase, &lag, seis_out,&ns);
	      cor[lag] = seis_out[0];
	      for( i = 1; i< (lag+1); i++)
		{ 
		  cor[lag-i] =  seis_out[i];
		  cor[lag+i] =  seis_out[ns-i];
		}
	      count++;
	      if(count!=1)
		{
		  for(k = 0; k < (2*lag+1); k++) sig_old[k] += cor[k];
		}
	      else
		{
		  //		  shdamp1.delta = 1;
		  lat1 = shdamp1.stla; 
		  lon1 = shdamp1.stlo;
		  lat2 = shdamp2.stla;
		  lon2 = shdamp2.stlo;
		  cmpaz= shdamp1.cmpaz;
		  for(k = 0; k < (2*lag+1); k++) sig_old[k] = cor[k];
		}
	      
	    }
	  if(count>0)
	    {
	      shdamp2.delta = 1;
	      shdamp2.user8=  count;
	      shdamp2.user9=  cmpaz;
	      //	      shdamp1.delta = sdb->rec[ine][jsta1].dt;
	      shdamp2.evla =  lat1;
	      shdamp2.evlo =  lon1;
	      shdamp2.stla =  lat2;
	      shdamp2.stlo =  lon2;
	      shdamp2.npts =  2*lag+1;
	      shdamp2.b    = -(lag)*shdamp2.delta;
	      for(k = 0; k < (2*lag+1); k++) sig_old[k] = sig_old[k]/count;
	      write_sac (filename, sig_old, &shdamp2);
	      //	      write_sac (filename, sig, &shd_cor );
		  
	    }
	      
	}  //loop over jsta2
    }  //loop over jsta1

  return 0;
}


//SAC_DB3 sdb;
//SAC_DB sdb_N,sdb_E;

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
int main (int na, char *arg[])
/*--------------------------------------------------------------------------*/
{
  FILE *ff;
  int iii,kkk;
  int ns1 = 0, ns2 = 0,lag;
  char str[600], filename[200], ch1[5],ch2[5];

  if ( na != 6 )
    {
      fprintf(stderr,"usage: justCOR lag iii(Z=0;N=1;E=2) kkk(Z=0;N=1;E=2) ns1 ns2\n");
      exit(1);
    }

  sscanf(arg[1],"%d", &lag );
  iii=atoi(arg[2]);
  kkk=atoi(arg[3]);
  ns1=atoi(arg[4]);
  ns2=atoi(arg[5]);
  //  if((iii!=1 && iii!=2) || (kkk!=1 && kkk!=2) )
  // {
  //  fprintf(stderr,"iii or kkk wrong!!\n");
  //  return;
  //}
  if((ff = fopen("stalst","r"))==NULL)
    {
      fprintf(stderr,"stalst not exist!!\n");
      return;
    }
  
  fclose(ff);

  //for(iii=0;iii<=0;iii++)
  //for(kkk=0;kkk<=0;kkk++)
  do_cor(lag,iii,kkk,ns1,ns2);  


  fprintf(stderr, "finished correlations\n");
 
  return 0;
}
