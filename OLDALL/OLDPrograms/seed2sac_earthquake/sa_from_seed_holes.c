#define MAIN

#include <stdio.h>
#include <math.h>
#include "mysac.h"
#include "sac_db.h"


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
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, long nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
	fsac = fopen(fname, "rb");
	if ( !fsac )
	{
	 fclose (fsac);
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
	//printf("\n\n\ntobal %f \n\n\n",SHD->o);
	/*-------------  calcule de t0  ----------------*/
	{
	  int eh, em ,i;
	  float fes;
	  char koo[9];
	  
	  for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
	  koo[8] = NULL;
	  
	  SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
	    SHD->nzsec + SHD->nzmsec*.001;
	  
	  //sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);
	  
	  //SHD->o  -= (eh*3600. + em*60. + fes);
	  //printf("\n\n\ntobal %s \n\n\n",koo);
	  /*-------------------------------------------*/}
	//printf("\n\n\ntobal %f \n\n\n",SHD->o);
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


        SHD->iftype = (long)ITIME;
        SHD->leven = (long)TRUE;

        SHD->lovrok = (long)TRUE;
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
 double abs_time ( int yy, long jday, long hh, long mm, long ss, long ms )
/*--------------------------------------------------------------------------
     computes time in s relative to 1900
--------------------------------------------------------------------------*/
{
  long nyday = 0, i;

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

char fname[300][300];
double t1[300], t2[300], T1, T2;
SAC_HD sd[300], s0;
int nf;

#define NPTSMAX 1000000

float sig0[NPTSMAX], sig1[NPTSMAX];
/*-------------------------------------------------------------------------*/
 int merge_sac(char *sta, char *chan, double *t0, float *dt, long *nrec)
/*-------------------------------------------------------------------------*/
{
  FILE *fi;
  int i, n, j, N, nfirst, Nholes;

  T1 = 1.e25;
  T2 = -100.;
  
  sprintf(str,"/bin/rm *.10.LHZ.Q.SAC *.AZ.PFO..LHZ.Q.SAC RESP.*.10.LHZ RESP.AZ.PFO..LHZ \0");
  system(str);

  sprintf(str,"ls *%s*%s*SAC > list_sac\0", sta, chan);
  system(str);

  fi = fopen("list_sac","r");

  if ( !fi )
    {
      fprintf(stderr,"no list_sac\n");

      sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      system(str);
      system("/bin/rm list_sac");

      fclose(fi);
      return 0;
    }

  if ( fscanf(fi,"%s", &(fname[0]) ) == EOF )
    {
      fprintf(stderr,"void list_sac\n" );

      sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      system(str);
      system("/bin/rm list_sac");

      fclose(fi);
      return 0;
    }

  fclose(fi);


  read_sac ( "2005.334.23.33.12.5481.IU.SNZO.00.LHZ.Q.SAC", sig1, &(sd[0]), NPTSMAX);
  write_sac("gill.SAC",sig1,&sd[0]);

  fi = fopen("list_sac","r");

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


  memcpy(&s0, &(sd[nfirst]), sizeof(SAC_HD) );

  N = nint((T2-T1)/s0.delta);
  if ( N > NPTSMAX ) N = NPTSMAX;

  if ( N > 100000 ) N = 100000;

  s0.npts = N;

  *t0 = T1;

  *dt = s0.delta;

  *nrec = s0.npts;


  for ( j = 0; j < N; j++ ) sig0[j] = 1.e30;


  fi = fopen("list_sac","r");

  for ( i = 0; ; )
    {
      int nb;
      double ti;

      if ( fscanf(fi,"%s", &(fname[i]) ) == EOF ) break;

      if ( !read_sac ( fname[i], sig1, &(sd[i]), NPTSMAX-N) )
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

  if ( (float)Nholes/(float)N > 0.1 )
    {
      fprintf(stderr,"too many holes\n");

      sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
      system(str);
      system("/bin/rm list_sac");

      return 0;
    }


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

  write_sac ("merged.sac", sig0, &s0);


  sprintf(str,"/bin/rm *%s*%s*SAC\0", sta, chan);
  system(str);
  system("/bin/rm list_sac");

  return 1;
}



/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void mk_one_rec (SAC_DB *sdb, int ne, int ns, char *nseed, char *ch)
/*------------------------------------------------------------------------*/
{
  FILE *ff;
  static char resp_name[150];

  if ( sdb->rec[ne][ns].n > 0 ) return;


  ff = fopen("from_seed","w");

  fprintf(ff,"/home/weisen/rdseedv5.0/rdseed.linux <<END\n");
  fprintf(ff,"%s", nseed );
  fprintf(ff,"\n");                             /* out file */
  fprintf(ff,"\n");                             /* volume */
  fprintf(ff,"d\n");                            /* option */
  fprintf(ff,"\n");                             /* summary file */
  fprintf(ff,"%s\n", sdb->st[ns].name );        /* station list */
  fprintf(ff,"%s\n", ch );                      /* channel list */
  fprintf(ff,"\n");                             /* network list */
  fprintf(ff,"\n");                             /* Loc Ids */
  fprintf(ff,"1\n");                            /* out format */
  fprintf(ff,"N\n");                            /* endtime */
  fprintf(ff,"N\n");                            /* Output poles & zeroes */
  fprintf(ff,"0\n");                            /* Check Reversal */
  fprintf(ff,"\n");                             /* Select Data Type */
  fprintf(ff,"\n");                             /* Start Time */
  fprintf(ff,"\n");                             /* End Time */
  fprintf(ff,"\n");                             /* Sample Buffer Length  */
  fprintf(ff,"Y\n");                            /* Extract Responses */
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");


  fclose(ff);

  system("sh from_seed");

  if ( !merge_sac(sdb->st[ns].name, ch, &(sdb->rec[ne][ns].t0), &(sdb->rec[ne][ns].dt), &(sdb->rec[ne][ns].n) ) )
    {
      sdb->rec[ne][ns].n = 0;
      return;
    }

  /*---------- response file -----------*/
  sprintf(str,"ls RESP*%s*%s* > list_resp\0",  sdb->st[ns].name,  ch);
  system(str);

  ff = fopen("list_resp","r");

  if ( fscanf(ff,"%s", resp_name ) == EOF )
    {
      sdb->rec[ne][ns].n = 0;
      return;
    }

  fclose(ff);

  sprintf(sdb->rec[ne][ns].resp_fname,"%s/%s\0", sdb->ev[ne].name, resp_name);
  sprintf(str,"/bin/mv %s %s\0", resp_name, sdb->rec[ne][ns].resp_fname);
  system(str);

  system("/bin/rm list_resp");
  system("/bin/rm RESP*");


  /*------------- mooving sac file -------*/
  sprintf(str,"/bin/mv merged.sac %s/%s.%s.SAC\0", sdb->ev[ne].name, sdb->st[ns].name, ch);
  system(str);
  
  sprintf(sdb->rec[ne][ns].fname,"%s/%s.%s.SAC\0", sdb->ev[ne].name, sdb->st[ns].name, ch);
  sprintf(sdb->rec[ne][ns].ft_fname,"%s/ft_%s.%s.SAC\0", sdb->ev[ne].name, sdb->st[ns].name, ch);

  sprintf(sdb->rec[ne][ns].chan,"%s\0", ch );
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void  fill_one_sta (STATION *st, char *buff )
/*------------------------------------------------------------------------*/
{
  int i;

  for ( i = 0; i < 4; i++ )
    {
      if ( buff[i] == ' ' ) break;
      else st->name[i] = buff[i];
    } 

  st->name[i] = '\0';

  /*fprintf(stderr,"station name %s\n", st->name );*/

  sscanf(&(buff[5]),"%g%g", &(st->lon), &(st->lat) );
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
  void fill_one_event (EVENT *ev, char *buff )
/*------------------------------------------------------------------------*/
{
  sscanf(&(buff[7]),"%d", &(ev->yy) );
  sscanf(&(buff[14]),"%d", &(ev->mm) );
  sscanf(&(buff[17]),"%d", &(ev->dd) );
  sscanf(&(buff[20]),"%2d", &(ev->h) );
  sscanf(&(buff[22]),"%2d", &(ev->m) );
  sscanf(&(buff[24]),"%2d", &(ev->s) );
  sscanf(&(buff[27]),"%2d", &(ev->ms) );

  ev->ms = 10.*ev->ms;

  ev->jday = jday( ev->yy, ev->mm, ev->dd );

  ev->t0 = abs_time (ev->yy, ev->jday, ev->h, ev->m, ev->s, ev->ms );

  sprintf(ev->name,"%d_%d_%d_%d_%d_%d\0",ev->yy, ev->mm, ev->dd, ev->h, ev->m, ev->s );

  sprintf(str,"mkdir %s\0", ev->name );
  system(str);
}



/*========================================================================*/


SAC_DB sdb;
char buff[300];

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 int main (void)
/*------------------------------------------------------------------------*/
{
  int ist, iev;
  FILE *ff;

  for ( iev = 0; iev < NEVENTS; iev++ ) for ( ist = 0; ist < NSTATION; ist++ ) sdb.rec[iev][ist].n = 0;

  fprintf(stderr,"initializing DB ok\n");


  ff = fopen("station.lst","r");

  for ( ist = 0; ; ist++ )
    {
      if ( !fgets(buff,300,ff) ) break;

      puts(buff);

      fill_one_sta (&(sdb.st[ist]), buff );

      fprintf(stderr,"filling station %s\n", sdb.st[ist].name );
    }

  sdb.nst = ist;

  fclose(ff);


  /*fprintf(stderr,"stations filled\n");
    scanf("%*1d");*/



  ff = fopen("input_ev_seed","r");

  for ( iev = 0;; )
    {
      if ( !fgets(buff,300,ff) ) break;

      if ( !strncmp(" PDE", buff, 4 ) )
	{
	  fill_one_event (&(sdb.ev[iev]), buff );

	  iev++;
	}

      else
	{
	  for ( ist = 0; ist < sdb.nst; ist++ )
	    {
	      /*fprintf(stderr,"starting rdseed ev %d  st %d\n", iev-1, ist );
		scanf("%*1d");*/

	      mk_one_rec (&sdb, iev-1, ist, buff, "LHZ");

	      /*fprintf(stderr,"finishing rdseed ev %d  st %d\n", iev-1, ist );
		scanf("%*1d");*/
	    }
	}
    }

  sdb.nev = iev;

  fclose(ff);


  ff = fopen("sac_db.out","wb");
  fwrite(&sdb, sizeof(SAC_DB), 1, ff );
  fclose(ff);

  sac_db_write_to_asc ( &sdb, "event_station.tbl" );
}
