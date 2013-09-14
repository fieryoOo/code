// this program is used to build up SACDB.out
// when you have sac files for each day of one month.
// you need two files in the directory of each month,
// which are event.dat, station.lst.
// The format for these two files can be checked out
// from the example files with the same names located 
// under the directory "/home/yingjie/progs/NOISE_CODA/SAC_FROM_SEED" 


#define MAIN
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/64_mysac.h"
#include "/home/tianye/code/Programs/head/64_sac_db.h"

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
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
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
  SAC_HD shd;

  ff = fopen(fname,"w");



  for ( ie = 0; ie < sdb->nev; ie++ ) for ( is = 0; is < sdb->nst; is++ )
    {
      fprintf(ff,"%s  %s  ", sdb->ev[ie].name, sdb->st[is].name );


      if ( sdb->rec[ie][is].n <= 0 ) 
       {fprintf(ff,"NO DATA\n");
       }
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
  double abssec;
  for ( i = 1901; i < yy; i++ )
    {
      if ( 4*(i/4) == i ) nyday += 366;
      else nyday += 365;
    }
  abssec =  24.*3600.*(nyday+jday) + 3600.*hh + 60.*mm + 1.*ss + 0.001*ms;
  return abssec;
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
SAC_HD sd, s0;
int nf;


#define NPTSMAX 1000000
float sig0[NPTSMAX], sig1[NPTSMAX];

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void mk_one_rec (SAC_DB *sdb, int ne, int ns, char *nseed, char *ch)
/*------------------------------------------------------------------------*/
{
  FILE *ff;
  static char resp_name[150];
  char seisname,ch2[4];
  char filename[300]; 
  double t1;
  if ( sdb->rec[ne][ns].n > 0 ) return;

  sprintf(filename, "%s/%s.%s.SAC", sdb->ev[ne].name, sdb->st[ns].name, ch);
 
   if (ch[2]=='Z') sprintf(ch2,"BHZ");
   if (ch[2]=='N') sprintf(ch2,"BHN");
   if (ch[2]=='E') sprintf(ch2,"BHE");

  if ( access(filename,F_OK) != 0) return;

  printf("filename %s iev %d ista %d\n",filename, ne,ns);

         if ( !read_sac (filename, sig1, &sd, NPTSMAX) )
	{
	  fprintf(stderr,"file %s not found\n", filename );
	}

  printf(" sd.nzyear %d,sd.nzjday %d,sd.nzhour %d ,sd.nzmin %d,sd.nzsec %d,sd.nzmsec %d\n",
	   sd.nzyear,   sd.nzjday,   sd.nzhour,    sd.nzmin,   sd.nzsec,   sd.nzmsec);

  t1 = abs_time(sd.nzyear,sd.nzjday,sd.nzhour,sd.nzmin,sd.nzsec,sd.nzmsec);
	printf("ok here get abs\n");

  sdb->rec[ne][ns].t0 = t1;
  sdb->rec[ne][ns].dt = sd.delta;
  sdb->rec[ne][ns].n =  sd.npts;


  /***************define the path to response file HERE!!  ********************/

  //    sprintf(str,"ls %s/RESP*%s* > list_resp\0", sdb->ev[ne].name, sdb->st[ns].name);
//here we change the polezero file to RESP file
//        sprintf(str,"ls ~/work/noise/northeast/PZs/SAC_%s_%s_PZs > list_resp\0", sdb->st[ns].name,ch2);
//        sprintf(str,"ls /home/zheng/work/noise/northeast/response_CH/RESP.%s.00.%s > list_resp\0", sdb->st[ns].name,ch2);
 //	  sprintf(str,"ls /mtera/weisen/for_yong/response_CH/RESP.%s.00.%s > list_resp\0", sdb->st[ns].name,ch2);
	//  sprintf(str,"ls /Users/jiayixie/Tianshan/RESP/RESP.%s.00.%s > list_resp",sdb->st[ns].name,ch2);
 //       sprintf(str,"ls ~/work/noise/northeast/PZs/SAC_%s_%s_PZs > list_resp\0", sdb->st[ns].name,ch);
          sprintf(str,"ls /utera/tianye/data_check_amp/RESP/RESP.TA.%s..%s > list_resp",sdb->st[ns].name,ch2);
/*	printf("%s\n",str);
       system(str);

   ff = fopen("list_resp","r");

   if ( fscanf(ff,"%s", resp_name ) == EOF )
    {
      sdb->rec[ne][ns].n = 0;
      return;
      sprintf(resp_name,"no_resp");
    }
  fclose(ff);
  

  sprintf(sdb->rec[ne][ns].resp_fname,"%s",  resp_name);
*/	
  //sprintf(sdb->rec[ne][ns].ft_name,"%s/%s/%s/%s.%s..%s.%d.%d",DATAROOT,shd1.knetwk,shd1.kstnm,shd1.kstnm,shd1.knetwk,shd1.kcmpnm,shd1.nzyear,shd1.nzjday);
 /*****************  define the name format HERE!  ***********************/
	sprintf(sdb->rec[ne][ns].fname,"%s.%s.%s.SAC", sdb->ev[ne].name, sdb->st[ns].name, ch);
  sprintf(sdb->rec[ne][ns].ft_fname,"ft_%s.%s.%s.SAC", sdb->ev[ne].name, sdb->st[ns].name, ch);
  sprintf(sdb->rec[ne][ns].chan,"%s", ch );
  
  printf("make_one_rec  sdb->rec[ne][ns].t0 %f\n",sdb->rec[ne][ns].t0);
  
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 void  fill_one_sta (STATION *st, char *buff )
/*------------------------------------------------------------------------*/
{
  int i;
  for ( i = 0; i < 6; i++ )
    {
      if ( buff[i] == ' ' ) break;
      else st->name[i] = buff[i];
    } 

  st->name[i] = '\0';
  sscanf(&(buff[6]),"%g%g", &(st->lon), &(st->lat) );
  printf(" lon %f lat %f\n",  st->lon, st->lat );
}

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
  void fill_one_event (EVENT *ev, char *buff )
/*------------------------------------------------------------------------*/
{

  sscanf(&(buff[0]),"%4d", &(ev->yy) );
  sscanf(&(buff[5]),"%2d", &(ev->mm) );
  sscanf(&(buff[8]),"%2d", &(ev->dd) );
  //sscanf(&(buff[11]),"%2d", &(ev->h) );
  //sscanf(&(buff[13]),"%2d", &(ev->m) );
  //sscanf(&(buff[15]),"%2d", &(ev->s) );

  //  sscanf(&(buff[20]),"%2d", &(ev->ms) );
  //ev->ms = 10.*ev->ms;

  ev->h  = 0;
  ev->m  = 0;
  ev->ms = 0;


  printf("ev->yy %d, ev->mm %d,ev->dd %d, ev->h %d, ev->m %d, ev->s %d, ev->ms %d\n",
	  ev->yy,    ev->mm,   ev->dd,    ev->h,    ev->m,    ev->s,     ev->ms);



  ev->jday = jday( ev->yy, ev->mm, ev->dd );

  ev->t0 = abs_time (ev->yy, ev->jday, ev->h, ev->m, ev->s, ev->ms );
  sprintf(ev->name,"%d_%d_%d_%d_%d_%d",ev->yy, ev->mm, ev->dd, ev->h, ev->m, ev->s );
  printf(" ev->t0 %f\n",ev->t0);
  system(str);
}


/*========================================================================*/


SAC_DB sdb;
char buff[300];

/*////////////////////////////////////////////////////////////////////////*/
/*------------------------------------------------------------------------*/
 int main (int argc,char *argv[])
/*------------------------------------------------------------------------*/
{
  int ist, iev;
  FILE *ff;
  char channel[3],reverse_file[300];
  if( argc < 3 ) {
    printf("usage: XXX Ch file_of_unit_stnm \n proveide event.dat first\n");
	  printf("see 'HERE' to define the formate/path for 1)input seismic_sac file 2)output revered seimic_sac file\n");
    return 0;
  }

  sscanf(argv[1],"%s",channel);
	sscanf(argv[2],"%s",&reverse_file[0]);
  for ( iev = 0; iev < NEVENTS; iev++ ) for ( ist = 0; ist < NSTATION; ist++ ) sdb.rec[iev][ist].n = 0;

  fprintf(stderr,"initializing DB ok\n");


/*  if((ff = fopen("station.lst","r"))==NULL)
  {fprintf(stderr,"cannot open file station.lst\n");
	  return(0);
  }
 */
	if((ff = fopen(reverse_file,"r"))==NULL)
	{fprintf(stderr,"cannot open file %s\n",reverse_file);
		return(0);
	}
	
  for ( ist = 0; ; ist++ )
  {  fprintf(stderr,"num of ist:%d\n",ist);
      if ( !fgets(buff,300,ff) ) break;

      puts(buff);

      fill_one_sta (&(sdb.st[ist]), buff );

      fprintf(stderr,"filling station %s!!!\n", sdb.st[ist].name );
    }

  sdb.nst = ist;
  fprintf(stderr,"%d\n",sdb.nst);

  fclose(ff);
	
fprintf(stderr,"OK here1!\n");

  if((ff = fopen("event.dat","r"))==NULL)
  { fprintf(stderr,"cannot open file event.dat");
	  return(0);}
	
	fprintf(stderr,"OK here2!\n");
	
  for ( iev = 0;;iev++ )
  { fprintf(stderr,"num of iev:%d\n",iev);
      if ( !fgets(buff,300,ff) ) break;
		fprintf(stderr,"get from event.dat:\n");
		puts(buff);

	  fill_one_event (&(sdb.ev[iev]), buff );
	  fprintf(stderr,"now do each sta:::: \n");

	  for ( ist = 0; ist < sdb.nst; ist++ )
	    {
	      mk_one_rec (&sdb, iev, ist, buff, channel);
	    }
    }



  sdb.nev = iev;
  fclose(ff);

  printf(" events %d  station %d\n", iev,ist);

  ff = fopen("sac_db_unit.out","wb");
  fwrite(&sdb, sizeof(SAC_DB), 1, ff );
  printf("write out sdb\n"); 
  sac_db_write_to_asc ( &sdb, "event_station_unit.tbl" );
  fclose(ff);


}
