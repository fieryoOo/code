// this programs is to extract seismc data from seed files.
// Each seed file constains seimic data for one event
// recorded at a number of stations. 
#define MAIN

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <cstdlib>
#include "/home/tianye/code/Programs/head/mysac64.h"
#include "/home/tianye/code/Programs/head/sac_db64_1ev.h"
using namespace std;

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        int isign(double f)
/*--------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
/*..........................................................................*/
        if (f < 0.)     return -1;
        else            return 1;
}
/*c/////////////////////////////////////////////////////////////////////////*/
/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
        int nint(double f)
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
{
        int i;
        double df;
/*..........................................................................*/
                           i=(int)f;
        df=f-(double)i;
        if (fabs(df) > .5) i=i+isign(df);

        return i;
}

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
	if((fsac = fopen(fname, "rb"))==NULL) {
           cout<<"Can't open file "<<fname<<endl;
           return NULL;
        }
	if ( !SHD ) SHD = &SAC_HEADER;

	 fread(SHD,sizeof(SAC_HD),1,fsac);
	 if ( SHD->npts > nmax )
	 {
	   /*fprintf(stderr,
	     "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);*/

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
	koo[8] = NULL;

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
 
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

	 fwrite(SHD,sizeof(SAC_HD),1,fsac);

	 fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


	fclose (fsac);
}


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	void sac_db_write_to_asc ( SAC_DB *sdb, char *fname )
/*--------------------------------------------------------------------------*/
{
  int ie, is;
  FILE *fi, *ff;
  static SAC_HD shd;

  ff = fopen(fname,"w");

  for ( ie = 0; ie < sdb->nev; ie++ ) for ( is = 0; is < sdb->nst; is++ )
    {
      fprintf(ff,"%s  %s  ", sdb->ev[ie].name, sdb->st[is].name );

      if ( sdb->rec[ie][is].n <= 0 ) fprintf(ff,"NO DATA\n");

      else 
	{
	  fi = fopen(sdb->rec[ie][is].fname,"rb");
	  fread(&shd, sizeof(SAC_HD), 1, fi );
	  fclose(fi);

	  fprintf(ff,"%s  t0: %d/%d:%d:%d:%g  %g s of record\n", sdb->rec[ie][is].fname, 
		   shd.nzyear, shd.nzjday, shd.nzhour, shd.nzmin, 
		   (shd.nzsec + 0.001*shd.nzmsec), shd.delta*shd.npts );
	}
    } 

  fclose(ff);
}


/*////////////////////////////////////////////////////////////////////////*/
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


/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 double abs_time ( int yy, int jday, int hh, int mm, int ss, int ms )
/*--------------------------------------------------------------------------
     computes time in s relative to 1900
--------------------------------------------------------------------------*/
{
  int nyday = 0, i;

  for ( i = 1901; i < yy; i++ )
    {
      if ( 4*(i/4) == i ) nyday += 366;
      else nyday += 365;
    }

  return 24.*3600.*(nyday+jday) + 3600.*hh + 60.*mm + ss + 0.001*ms;
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
float av_sig (float *sig, int i, int N, int nwin )
/*------------------------------------------------------------------------*/
{
  int n1, n2, j, nav = 0;
  float av = 0.;

  if ( nwin > N ) nwin = N;

  n1 = i - nwin/2;

  if ( n1 < 0 ) n1 = 0;

  n2 = n1 + nwin - 1;

  if ( n2 > N-1 ) n2 = N-1;

  n1 = n2 - nwin + 1;

  for ( j = n1; j <= n2; j++ ) if ( sig[j] < 1.e29 )
    {
      av += sig[j];
      nav++;
    }

  if ( nav < 1 ) av = 1.e30;

  else av = av/(float)nav;

  return av;
}

/*/////////////////////////////////////////////////////////////////////////*/


char str[300];

/*/////////////////////////////////////////////////////////////////////////*/
char fname[2000][100];
double t1[2000], t2[2000], T1, T2;
SAC_HD sd[2000], s0;
int nf;

#define NPTSMAX 300000
float sig0[NPTSMAX], sig1[NPTSMAX];

/*-------------------------------------------------------------------------*/
 int merge_sac(char *sta, char *chan, double *t0, float *dt, int *nrec, EVENT *ev)
/*-------------------------------------------------------------------------*/
{
  FILE *fi;
  int i, n, j, N, nfirst, Nholes;



  T1 = 1.e25;
  T2 = -100.;

  //  sprintf(str,"ls *%s*%s*SAC > list_sac\0", sta, chan);
  //system(str);

  fi = fopen("list_sac1","r");

  if ( !fi )
    {
      fprintf(stderr,"no list_sac1\n");
      //sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      //system(str);
      //      system("/bin/rm list_sac");
      fclose(fi);
      return 0;
    }

  if ( fscanf(fi,"%s", &(fname[0]) ) == EOF )
    {
      fprintf(stderr,"void list_sac1\n" );

      //      sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      //system(str);
      //      system("/bin/rm list_sac");

      fclose(fi);
      return 0;
    }

  fclose(fi);

  fprintf(stderr,"inside1  merge_sac\n" );

  fi = fopen("list_sac1","r");

  for ( i = 0; ; )
    {
      if ( fscanf(fi,"%s", &(fname[i]) ) == EOF ) break;

      if ( !read_sac ( fname[i], sig1, &(sd[i]), NPTSMAX) )
	{
	  fprintf(stderr,"file %s not found\n", fname[i] );
	  continue;
	}

      t1[i] = abs_time (sd[i].nzyear, sd[i].nzjday, sd[i].nzhour, sd[i].nzmin, sd[i].nzsec, sd[i].nzmsec );
      t2[i] = t1[i] + sd[i].npts*sd[i].delta;

      if ( t1[i] < T1 )
	{
	  T1 = t1[i];
	  nfirst = i;
	}

      if ( t2[i] > T2 ) T2 = t2[i];

      i++;
    }

   fclose(fi);

  fprintf(stderr,"inside2  merge_sac\n" );
  memcpy(&s0, &(sd[nfirst]), sizeof(SAC_HD) );

  
  if( T1 > ev->t0 - 900 ) T1 = ev->t0 - 900;
  T2 = ev->t0 + 8100; 

  
  fprintf(stderr,"T1 %f T2  %f\n",T1,T2 );


  N = nint((T2-T1)/s0.delta);
  if ( N > NPTSMAX ) N = NPTSMAX;

  s0.npts = N;

  //  *t0 = T1;

  // *dt = s0.delta;

  //   *nrec = s0.npts;

  for ( j = 0; j < N; j++ ) sig0[j] = 1.e30;

  fi = fopen("list_sac1","r");

  for ( i = 0; ; )
    {
      int nb;
      double ti;

      if ( fscanf(fi,"%s", &(fname[i]) ) == EOF ) break;

      if ( !read_sac ( fname[i], sig1, &(sd[i]), NPTSMAX) )
	{
	  fprintf(stderr,"file %s not found\n", fname[i] );
	  continue;
	}

      if ( fabs(sd[i].delta-s0.delta) > .0001 )
	{
	  fprintf(stderr,"incompatible dt in file file %s\n", fname[i] );
	  continue;
	}

      ti = abs_time (sd[i].nzyear, sd[i].nzjday, sd[i].nzhour, sd[i].nzmin, sd[i].nzsec, sd[i].nzmsec );

      nb = nint((ti-T1)/s0.delta);

      for ( j = 0; j < sd[i].npts; j++ )
	{
	  int jj = nb+j;

	  if ( sig0[jj] > 1.e29 ) sig0[jj] = sig1[j];
	}

      i++;
    }

  fclose(fi);


  Nholes = 0;

  for ( j = 0; j < N; j++ ) if ( sig0[j] > 1.e29 ) Nholes++;

  if ( (float)Nholes/(float)N > 0.7 )
    {
      fprintf(stderr,"too many holes\n");

      //      sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      //system(str);
      //      system("/bin/rm list_sac");

      return 0;
    }

        fprintf(stderr,"holes numer %f\n",(float)Nholes/(float)N );

  for ( j = 0; j < N; j++ ) if ( sig0[j] > 1.e29 )
    {
      float av;
      int npart = 16;

      for ( ;;) 
	{
	  av = av_sig (sig0, j, N, N/npart );

	  /*fprintf(stderr,"av %g  npart %d", av, npart );
	    scanf("%*1d");*/

	  if ( av < 1.e29 ) break;

	  if ( npart = 1 )
	    {
	      av = 0.;
	      break;
	    }

	  npart = npart/2;
	}

      sig0[j] = av;
    }


  printf(" T1 %f  T2 %f ev->t0 %f\n", T1,T2, ev->t0);

  if ( T1 > (ev->t0) || T2 < (ev->t0 + 8000) ) 
    { printf("%s_%s_sac is not enough int\n",sta,chan);
      //sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      //system(str);
      //      system("/bin/rm list_sac");      
      return 0;
    }
  
  int diff = nint((ev->t0 - 900 - T1)/s0.delta);
 
  s0.nzyear = ev->yy;
  s0.nzjday = ev->jday;
  s0.nzhour = ev->h;
  s0.nzmin  = ev->m;
  s0.nzsec  = ev->s;
/*
  if(s0.nzmin>=15) s0.nzmin-=15;
  else if(s0.nzhour>0) {s0.nzmin+=45; s0.nzhour-=1;}
  else if(s0.nzjday>1) {s0.nzmin+=45; s0.nzhour+=23; s0.nzjday-=1;}
  else {s0.nzmin+=45; s0.nzhour+=23; s0.nzyear-=1;
        if(((s0.nzyear-1)%4==0 && (s0.nzyear-1)%100!=0) ||(s0.nzyear-1)%400==0) s0.nzjday=366;
        else s0.nzjday=365;}
*/
  s0.evlo   = ev->lon;
  s0.evla   = ev->lat;
  s0.npts  = int(8000/s0.delta)+1;
  s0.b = -900;
  s0.e = 7100;
 *nrec = s0.npts;

  write_sac ("merged.sac", &sig0[diff], &s0);


  //sprintf(str,"/bin/rm  *%s*%s*SAC\0", sta, chan);
  //system(str);
  //  system("/bin/rm list_sac");
  printf("it is good here\n");
  return 1;
}



/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void mk_one_rec (SAC_DB *sdb, int ne, int ns, char *nseed, char *ch)
/*------------------------------------------------------------------------*/
{
  int change_sac_file (char *buff1);
  FILE *ff;
  static char resp_name[150];

  float fl1, fl2, fl3, fl4;

  fl1= 0.001;		
  fl2= 0.0015;
  fl3= 10;
  fl4= 15;

  sdb->rec[ne][ns].n = 0;

  FILE *fi;
/*---------------change every sac file -----------------------------*/
  sprintf(str,"ls *.%s.*.%s.*.SAC>list_sac1\0",sdb->st[ns].name,ch);
  system(str);
  int i=0;
  fi=fopen("list_sac1","r");
  
  if(fi)
     for (i=0;;)
      {
      if(fscanf(fi,"%s",fname[i])!=EOF)
        {
        puts(fname[i]);
        change_sac_file (fname[i]);
    
        i++;
        
        }
       else break;
       }
  fclose(fi);
  if(i==0)
    return;
/*---------------------------over----------------------------------*/


//  fprintf(stderr,"before merge_sac\n" );
  int iflag =  merge_sac(sdb->st[ns].name, ch, &(sdb->rec[ne][ns].t0), &(sdb->rec[ne][ns].dt), &(sdb->rec[ne][ns].n),&(sdb->ev[ne]));

  if ( iflag == 0 )
    { printf( " not good here");
      sdb->rec[ne][ns].n = 0;
      return;
    }

 /*---------- response file -----------*/
  //  sprintf(str,"ls /data/minos-2/linf/Earthquake/all_seed_files/all_BHZ_RESP/%d/%s/RESP*%s*%s* > list_resp\0", sdb->ev[ne].yy, sdb->ev[ne].name, sdb->st[ns].name,  ch);
  sprintf(str,"ls RESP*.%s.*.%s*  > list_resp\0",sdb->st[ns].name,  ch);
  system(str);


  ff = fopen("list_resp","r");
  if ( fscanf(ff,"%s", resp_name ) == EOF )
    {
      sdb->rec[ne][ns].n = 0;
      return;
    }
  fclose(ff);
  
  //  fprintf(stderr,"before evalresp\n" );
/* remove instrument response  GENERATE TEMPORARY FILES */
//  sprintf(str,"cp /data/minos-2/linf/Earthquake/all_seed_files/all_BHZ_RESP/%d/%s/RESP*%s*%s*  resp1\0",sdb->ev[ne].yy, sdb->ev[ne].name,  sdb->st[ns].name,  ch);
  sprintf(str," mv %s  resp1\0",resp_name);
  system(str);
    ff = fopen("sac_bp_respcor","w");
  fprintf(ff,"/home/tianye/Software/sac/bin/sac << END\n");
  fprintf(ff,"r merged.sac\n");
  fprintf(ff,"rmean\n");
  fprintf(ff,"rtrend\n");
  fprintf(ff,"taper\n");
  fprintf(ff,"transfer from evalresp fname resp1 to vel freqlimits %f %f %f %f\n", fl1,fl2,fl3,fl4 );
  fprintf(ff,"w merged_v1.sac\n");
  fprintf(ff,"cut 0 5000\n");
  fprintf(ff,"r merged_v1.sac\n");
  fprintf(ff,"w over\n");
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");
  fclose(ff);
  system("sh sac_bp_respcor");

  // fprintf(stderr,"before evalresp 0.5\n" );
/*
  sprintf(str,"/home/tianye/Software/evalresp/evalresp-3.3.3/evalresp %s %s %d %d 0.001 10 1000 -f %s\0",sdb->st[ns].name, ch, sdb->ev[ne].yy, sdb->ev[ne].jday,resp_name);
  system(str);

  // fprintf(stderr,"done evalresp\n" );
  
  sprintf(str,"/bin/mv AMP*%s*%s amp1.txt\0",sdb->st[ns].name, ch);
  system(str);
  sprintf(str,"/bin/mv PHASE*%s*%s phase1.txt\0",sdb->st[ns].name, ch);
  system(str);
  
  // fprintf(stderr,"before rm_instrument_response\n" );
  sprintf(str,"/home/tianye/code/Programs/RM_RESPONSE/rm_instrument_response merged.sac phase1.txt amp1.txt merged_v1.sac\0");
  system(str);
*/
  // fprintf(stderr,"done rm_instrument_response\n" );
  sprintf(sdb->rec[ne][ns].resp_fname,"%s\0",  resp_name);
  //sprintf(str,"/bin/mv %s %s\0", resp_name, sdb->rec[ne][ns].resp_fname);
  //system(str);
   //  system("/bin/rm list_resp");
  // fprintf(stderr,"done rm_instrument_response 0.1\n" );
  sprintf(str,"/bin/rm list_resp list_sac1 resp1 merged.sac sac_bp_respcor\0");
  system(str);
  //  system("/bin/rm list_resp phase1.txt amp1.txt merged.sac");
  //system("/bin/rm AMP*");
  //system("/bin/rm PHASE*");
  // fprintf(stderr,"done rm_instrument_response 0.2\n" );
  printf("/bin/mv merged_v1.sac %s/%s.%s.%s.sac\0", sdb->ev[ne].name, sdb->ev[ne].name,sdb->st[ns].name, ch);
  sprintf(str,"/bin/mv merged_v1.sac %s/%s.%s.%s.sac\0", sdb->ev[ne].name, sdb->ev[ne].name,sdb->st[ns].name, ch);
  system(str);
return;
  //fprintf(stderr,"done rm_instrument_response 0.3\n" );
  sprintf(sdb->rec[ne][ns].fname,   "%s/%s.%s.%s.sac\0", sdb->ev[ne].name, sdb->ev[ne].name,sdb->st[ns].name, ch);
  sprintf(sdb->rec[ne][ns].ft_fname,"%s/%s.%s.%s.norm.sac\0", sdb->ev[ne].name, sdb->ev[ne].name,sdb->st[ns].name, ch);
  sprintf(sdb->rec[ne][ns].chan,"%s\0", ch );
  //fprintf(stderr,"done rm_instrument_response 0.5\n" );
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void  fill_one_sta (STATION *st, char *buff )
/*------------------------------------------------------------------------*/
{
  int i;

  for ( i = 0; i < 6; i++ )
    {
      if ( buff[i] == ' ' || buff[i] == '\t' ) break;
      else st->name[i] = buff[i];
    } 

  st->name[i] = '\0';

  /*fprintf(stderr,"station name %s\n", st->name );*/

  sscanf(&(buff[i+1]),"%g%g", &(st->lon), &(st->lat) );
  //cout<<st->name<<" "<<st->lon<<" "<<st->lat<<endl;
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
  void fill_one_event (EVENT *ev, char *buff,char *bufftemp )
/*------------------------------------------------------------------------*/
{ 
  char yyc[20],mmc[20],ddc[20],hc[20],mc[20],sc[20];

  sscanf(buff,"%4d%2d%2d%2d%2d%2d %g %g", &(ev->yy), &(ev->mm), &(ev->dd), &(ev->h), &(ev->m), &(ev->s), &(ev->lon), &(ev->lat));
  ev->ms = 0;
  ev->jday = jday( ev->yy, ev->mm, ev->dd );  
  ev->t0 = abs_time (ev->yy, ev->jday, ev->h, ev->m, ev->s, ev->ms );
  sprintf(ev->name,"%04d%02d%02d%02d%02d%02d\0", ev->yy, ev->mm, ev->dd, ev->h, ev->m, ev->s );

    sprintf(bufftemp,"*_%d_%d_%d_%d_%d\0",ev->yy, ev->mm, ev->dd, ev->h, ev->m );
    //sprintf(str,"mkdir -p %s\0", ev->name );
    system(str);
}

/*========================================================================*/

int change_sac_file (char *buff1)
{
SAC_HD shd;
float sig1[NPTSMAX];
float sig2[NPTSMAX];
int i=0;
float fra;
shd=sac_null;
if (read_sac(buff1,sig1,&shd,NPTSMAX))
   {
    if (shd.nzmsec!=0)
      {
      fra=1.-0.001*(float)shd.nzmsec;
      shd.nzmsec=0;
      
      shd.npts=shd.npts-1;
      if(shd.npts==0) return 0;
      
      for (i=0;i<shd.npts;i++)
        sig2[i]=sig1[i]+fra*(sig1[i+1]-sig1[i]);
 

      if (shd.nzsec==59)
       { shd.nzsec=0;
         if (shd.nzmin==59)
            {
            shd.nzmin=0;
              if (shd.nzhour==23)
                   {
                   shd.nzhour=0;
                      if ((shd.nzjday==364&&((shd.nzyear/4)*4!=shd.nzyear) )||(shd.nzjday==365))
                      {
                      shd.nzyear+=1;
                      shd.nzjday=0;
                      }

                      else
                      shd.nzjday+=1;
                   }
               else
               shd.nzhour+=1;
             }
             else
             shd.nzmin+=1;
      
        }
       else
       {shd.nzsec=shd.nzsec+1;
       }
     
       write_sac(buff1,&(sig2[0]),&shd);
       return 1;
     
       }
     else 
       {
       return 1; 
       }
     }
else
{return 0;}


}


/*========================================================================*/


SAC_DB sdb_N,sdb_E,sdb_Z;
char buff[300];

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 int main (int na, char *arg[])
/*------------------------------------------------------------------------*/
{
  if (na!=4)
    {
      cout<<"usage: sac_from_seed [SEED_path] [events.loc] [station.lst]"<<endl;
      return 0;
    }
  int ist, iev;
  char tempfn[100],bufftemp[100];
  FILE *ff,*fii,*fff;


  for ( iev = 0; iev < NEVENTS; iev++ ) 
    for ( ist = 0; ist < NSTATION; ist++ ) 
      {
	sdb_E.rec[iev][ist].n = 0;
	sdb_N.rec[iev][ist].n = 0;
	sdb_Z.rec[iev][ist].n = 0;
      }  
  printf(" NEVENTS %d  NSTATION %d\n", NEVENTS,NSTATION);

  fprintf(stderr,"initializing DB ok\n");


  if((ff = fopen(arg[3],"r"))==NULL) {
     cout<<"Can't open file "<<arg[3]<<endl;
     return 0;
  }
  cout<<"Filling stations..."<<endl;
  for ( ist = 0; ; ist++ )
    {
      if ( !fgets(buff,300,ff) ) break;

      puts(buff);

      fill_one_sta (&(sdb_E.st[ist]), buff );
      //fprintf(stderr,"station %s filled\n", sdb_E.st[ist].name );
      fill_one_sta (&(sdb_N.st[ist]), buff );
      //fprintf(stderr,"station %s filled\n", sdb_N.st[ist].name );
      fill_one_sta (&(sdb_Z.st[ist]), buff );
      //fprintf(stderr,"station %s filled\n", sdb_Z.st[ist].name );
      
    }

  sdb_E.nst = ist;
  sdb_N.nst = ist;
  sdb_Z.nst = ist;

  fclose(ff);


  if((ff = fopen(arg[2],"r"))==NULL) {
     cout<<"Can't open file "<<arg[2]<<endl;
     return 0;
  }
  FILE *ftemp;
  char eventname[30],buff2[300];
  int int_temp,int_year;
  
  for (iev=0;;)
    {  
      printf(" iev: %d\n",iev);
      
      if ( !fgets(buff,300,ff) ) break;
      sscanf(&(buff[0]),"%14s", eventname );
      sscanf(&(buff[0]),"%4d",&int_year);
      fill_one_event (&(sdb_E.ev[iev]), buff,bufftemp );
      fill_one_event (&(sdb_N.ev[iev]), buff,bufftemp );
      fill_one_event (&(sdb_Z.ev[iev]), buff,bufftemp );
      sprintf(str,"ls %s/%s* > list_seed\0", arg[1] ,bufftemp);
      system(str);
      
      fii = fopen("list_seed","r");
      if ( fscanf(fii,"%s", &buff ) != EOF )
	{
          sprintf(str,"mkdir -p %s\0", sdb_Z.ev[iev].name );
          system(str);
	  printf("\nseed file: %s\n", buff);
	  fff = fopen("from_seed","w");
	  
	  fprintf(fff,"/home/tianye/Software/rdseedv5.1/rdseed <<END\n");
	  fprintf(fff,"%s\n", buff );
	  fprintf(fff,"\n");                             /* out file */
	  fprintf(fff,"\n");                             /* volume */
	  fprintf(fff,"d\n");                            /* option */
	  fprintf(fff,"\n");                             /* summary file */
	  fprintf(fff,"\n");                             /* station list */
	  fprintf(fff,"\n");                          /* channel list */
	  fprintf(fff,"\n");                             /* network list */
	  fprintf(fff,"\n");                             /* Loc Ids */
	  fprintf(fff,"1\n");                            /* out format */
          fprintf(fff,"N\n");                            /* endtime */
	  fprintf(fff,"N\n");                            /* Output poles & zeroes */
	  fprintf(fff,"0\n");                            /* Check Reversal */
	  fprintf(fff,"\n");                             /* Select Data Type */
	  fprintf(fff,"\n");                             /* Start Time */
	  fprintf(fff,"\n");                             /* End Time */
	  fprintf(fff,"\n");                             /* Sample Buffer Length  */
	  fprintf(fff,"Y\n");                            /* Extract Responses */
	  fprintf(fff,"quit\n");
	  fprintf(fff,"END\n");
	  
	  fclose(fff);
	  system("sh from_seed");
          system("/bin/rm from_seed list_seed");
	  for ( ist = 0; ist < sdb_Z.nst; ist++ ) 
	    {
	      mk_one_rec (&sdb_N, iev, ist, buff, "BHN");
	      mk_one_rec (&sdb_E, iev, ist, buff, "BHE");
	      mk_one_rec (&sdb_Z, iev, ist, buff, "BHZ");
	    }
	  system("/bin/rm *SAC");
	  system("/bin/rm RESP*");
	  
	}
      fclose(fii);	       
      
      
    }
  



  sdb_E.nev = iev;
  sdb_N.nev = iev;
  sdb_Z.nev = iev;
  fclose(ff);


  ff = fopen("sac_db_BHN.out","wb");
  fwrite(&sdb_N, sizeof(SAC_DB), 1, ff );
  fclose(ff);
  ff = fopen("sac_db_BHE.out","wb");
  fwrite(&sdb_E, sizeof(SAC_DB), 1, ff );
  fclose(ff);
  ff = fopen("sac_db_BHZ.out","wb");
  fwrite(&sdb_Z, sizeof(SAC_DB), 1, ff );
  fclose(ff);

  sac_db_write_to_asc ( &sdb_E, "event_station_E.tbl" );
  sac_db_write_to_asc ( &sdb_N, "event_station_N.tbl" );
  sac_db_write_to_asc ( &sdb_Z, "event_station_Z.tbl" );
}
