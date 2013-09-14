#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "/home/tianye/code/Programs/head/64_mysac.h"
#include "/home/tianye/code/Programs/head/64_sac_db.h"

//#include "/home/zheng/progs/NOISE_CODA/HEAD_NOISE/64_mysac.h"
//#include "/home/zheng/progs/NOISE_CODA/HEAD_NOISE/64_sac_db.h"



/* FUNCTION PROTOTYPES */
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax);
void write_sac (char *fname, float *sig, SAC_HD *SHD);
void one_rec_trans( SAC_DB *sd, int ne, int ns);
void one_rec_cut(SAC_DB *sd, int ne, int ns, float t1, float n );


/*--------------------------------------------------------------------------*/
SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, int nmax)
/*--------------------------------------------------------------------------*/
/* function to read sac files given the name, fname. The function outputs the time signal to the pointer sig, fills the header SHD, if the signal has fewer than nmax points */
{
  FILE *fsac;

  if((fsac = fopen(fname, "rb")) == NULL) {
    printf("could not open sac file to read %s \n", fname);
    exit(1);
  }

  if ( !fsac ) {
    /*fprintf(stderr,"file %s not found\n", fname);*/
    return NULL;
  }

  if ( !SHD ) SHD = &SAC_HEADER;

  fread(SHD,sizeof(SAC_HD),1,fsac);

  if ( SHD->npts > nmax ) {
    fprintf(stderr,"ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);
    SHD->npts = nmax;
  }

  fread(sig,sizeof(float),(int)(SHD->npts),fsac);
  fclose (fsac);

/*-------------  calcule de t0  ----------------*/
   {
     int eh, em ,i;
     float fes;
     char koo[9];

     for ( i = 0; i < 8; i++ ) {
       koo[i] = SHD->ko[i];
     }
     koo[8] = 0;

     SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
     SHD->nzsec + SHD->nzmsec*.001;

     sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

     SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

   return SHD;
}


/*--------------------------------------------------------------------------*/
void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*--------------------------------------------------------------------------*/
{
  FILE *fsac;
  int i;
  if((fsac = fopen(fname, "wb"))==NULL) {
    printf("could not open sac file to write\n");
    exit(1);
  }

  if ( !SHD ) {
    SHD = &SAC_HEADER;
  }

  SHD->iftype = (int)ITIME;
  SHD->leven = (int)TRUE;
  SHD->lovrok = (int)TRUE;
  SHD->internal4 = 6L;
  SHD->depmin = sig[0];
  SHD->depmax = sig[0];

  for ( i = 0; i < SHD->npts ; i++ ) {
    if ( SHD->depmin > sig[i] ) {
      SHD->depmin = sig[i];
    }
    if ( SHD->depmax < sig[i] ) {
      SHD->depmax = sig[i];
    }
   }

  fwrite(SHD,sizeof(SAC_HD),1,fsac);
  fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);

  fclose (fsac);
}


/*--------------------------------------------------------------------------*/
void one_rec_cut(SAC_DB *sd, int ne, int ns, float t1, float n)
/*--------------------------------------------------------------------------*/
{
  float sig1[7500000],sig2[7500000],sig3[7500000],sig4[7500000];
  double t1b, t1e, t2b, t2e, t2;
  int n1, n2,i,j,kk,npts_real,flag=0;
  char ft_name[500],ft_name_bit[500],ft_name_cut[500],str[300],tempstnm[10];
  SAC_HD shd1,shd2,shd3;
	FILE *script;


 if ( ne >= sd->nev ) return;
 if ( ns >= sd->nst  ) return;
 if ( sd->rec[ne][ns].n <= 0 ) return;

  fprintf(stderr,"t1 n %g %g\n",t1,n);
  printf("ne ns [begin from 0] %d %d\n",ne,ns);
  t2 = t1 + (n-1)*sd->rec[ne][ns].dt;
//printf("n: %g  rec.dt: %g\n",n,sd->rec[ne][ns].dt);
//printf("rec.t0: %g\n",sd->rec[ne][ns].t0);
//printf("ev.t0: %g\n",sd->ev[ne].t0);
//printf("rec.n: %g\n",sd->rec[ne][ns].n);
//printf("rec.dt: %g\n",sd->rec[ne][ns].dt);

  t1b = sd->rec[ne][ns].t0 - sd->ev[ne].t0;
  t1e = t1b + (sd->rec[ne][ns].n-1)*sd->rec[ne][ns].dt;
	
/*	sprintf(str,"/bin/cp %s s1.sac", sd->rec[ne][ns].fname ); 
	fprintf(stderr,"/bin/cp %s s1.sac\n", sd->rec[ne][ns].fname );
	system(str);	
*/
	  //strcpy(ft_name, sd->rec[ne][ns].ft_fname);
	
	/////++++++ If ONLY do the cut not need the trans, use this sentence. Otherwise, delete it
//	sprintf(str,"/bin/cp %s s1.sac", sd->rec[ne][ns].fname ); 
//	system(str);
	///////++++++
	
	strcpy(ft_name, sd->rec[ne][ns].ft_fname);
	fprintf(stderr,"t1 %lg  t2 %lg   t1b %lg  t1e %lg\n", t1, t2, t1b, t1e);
//	sprintf(ft_name_bit, "%s_bit",ft_name);
	
  fprintf(stderr,"t1 %lg  t2 %lg   t1b %lg  t1e %lg\n", t1, t2, t1b, t1e);
  

    if ( !read_sac ("s1.sac", sig1, &shd1, 7500000 ) ) {
       fprintf(stderr,"file %s not found\n", sd->rec[ne][ns].fname );
       return;
    }

    if ( !read_sac ("s1.sac", sig1, &shd2, 7500000 ) ) {
       fprintf(stderr,"file %s not found\n", sd->rec[ne][ns].fname );
       return;
    }

    if ( !read_sac ("s1.sac", sig1, &shd3, 7500000 ) ) {
       fprintf(stderr,"file %s not found\n", sd->rec[ne][ns].fname );
       return;
	}
// sig2[]contain the seimic waveform, if no signal(i.e. there is a range in sig2[] that is out of that of sig1[]), we put 1 there. sig3[], the bit file, if sig2[]has signal, sig3[]=1, else sig3[]=0
//sig4 contain only the real signals, they are attached together without any manmade sig, but the time of the sig4 is not true, the length is shorter is there was manmade sig
	//ini sig3 == 0->1 changed! ; sig2 == 1;  
	
	for(i=0;i<n;i++) sig3[i]=1.0;

  if ( (t1b>t1) || (t1e<t2) ) {
      printf("\n t1b>t1 OR t1e<t2: not all contained\n");
//	  printf("here1!!!\n");

    for(i=0;i<n;i++) sig2[i]=1;

    if( (t1b>t1) ) {
       n1 =  (int)((t1b-t1)/sd->rec[ne][ns].dt);
		kk=0;
		for (i=0; i<shd1.npts; i++) {
			if (n1+i<n) sig2[n1+i] = sig1[i];
		}
		
       for(i=0;i< n;i++) 
	   { 
//		   sig2[n1+i] =sig1[i];
		   if(fabs(sig2[i+1]-sig2[i])<0.00000001)
		   {	if(flag==0)fprintf(stderr,"==============here manmade sig!==============\n");
			   flag=flag+1;
			   sig3[i]=0.0;
			   sig3[i+1]=0.0;
		//	   continue; // in this case, sig2 and sig3 will keep the initialized values, but sig4 won't keep any
		   }

 //          sig2[n1+i] =sig1[i];
		//   sig3[n1+i] = 1;
		  if(fabs(sig3[i]-1.0)<0.1)// its still within the day we need
		  { sig4[kk]=sig2[i];
			  kk++;}
       }  // for i
		npts_real=kk;
//		abort();
    }

    else {  
//		printf("here2!!!\n");
       n1 =  (int)((t1-t1b)/sd->rec[ne][ns].dt);
		kk=0;
		for (i=n1;i<shd1.npts;i++) sig2[i-n1] = sig1[i];
		
		
       for(i=0;i< n;i++) 
	   {
	//	   fprintf(stderr,"here i =%d\n",i);
		   if(fabs(sig2[i+1]-sig2[i])<0.00000001)
		   {	if(flag==0)fprintf(stderr,"==========here manmade sig!==============\n");
			   flag= flag+1;
			   sig3[i]=0.0;
			   sig3[i+1]=0.0;
		//	   continue; // see before
		   }
	//		   sig2[i-n1] =sig1[i]; 
		//		sig3[i-n1] = 1;
		   if(fabs(sig3[i]-1.0)<0.000001)
		   {  sig4[kk]=sig2[i];
			   kk++;}
		   
       }       //for i
		npts_real=kk;
    } 
     
     shd2.npts = n;//here cut the singnal into the length we need !!!
     shd2.nzyear = sd->ev[ne].yy;
     shd2.nzjday = sd->ev[ne].jday;
     shd2.nzhour = 0;
     shd2.nzmin = 0;
     shd2.nzsec = 0;
     shd2.nzmsec = 0;
     shd2.b = t1;
	  
/*		  //get rid of the space after shd.kstnm	  
	  int mm;
	  for ( mm = 0; mm < 6; mm++ )
	  {
		  if ( shd2.kstnm[mm] == ' ' ) break;
		  else tempstnm[mm] = shd2.kstnm[mm];
	  } 
	  tempstnm[mm] = '\0';
	  
	  //mkdir for SAC file DATAROOT/NTW/STN/...
	  if((script=fopen("make_Dir.csh","w"))==NULL)
	  {fprintf(stderr,"cannot creat file make_Dir\n");
		  exit(1);}
	  fprintf(script,"if (! -e %s/%s/%s) then\n mkdir %s/%s/%s\n endif\n",DATAROOT,shd2.knetwk,tempstnm,DATAROOT,shd2.knetwk,tempstnm);
	  fclose(script);
	  
///	  system("csh make_Dir.csh");
	  fprintf(stderr,"mkdir %s/%s/%s\n",DATAROOT,shd2.knetwk,tempstnm);
	  
	  sprintf(ft_name,"%s/%s/%s/%s.%s.00.%s.%d.%d",DATAROOT,shd2.knetwk,tempstnm,tempstnm,shd2.knetwk,shd2.kcmpnm,shd2.nzyear,shd2.nzjday);
	  //fprintf(stderr,"here is the name:%s\n",ft_name);
	  fprintf(stderr,"write SAC:%s\n\n",ft_name);
//	  sprintf(ft_name_bit, "%s_bit",ft_name);
	  sprintf(ft_name_cut,"%s_cut",ft_name);
	  
	  printf("n=%g n4=%d\n",n,npts_real);
	  
///  write_sac (ft_name_cut, &(sig2[0]), &shd2 );
//  write_sac (ft_name_bit, &(sig3[0]), &shd2 );
	  write_sac(ft_name,&(sig2[0]),&shd2);
///	  shd2.npts = npts_real;
///	  write_sac(ft_name,&(sig4[0]),&shd2);
*/
	  //++++
	  write_sac (ft_name,     &(sig2[0]), &shd2 );
//	  write_sac (ft_name_bit, &(sig3[0]), &shd2 );
	  system("/bin/rm s1.sac");      

  }


  else { // the merged signal have the TOTAL signal for the day we need.
     printf("\n else: all contained  \n");
     n1 = (int)((t1-t1b)/sd->rec[ne][ns].dt);
	  kk=0;
	  for(i=0;i<n;i++)
	  {
		  sig3[i] = 1;
		  if(fabs(sig1[i+1+n1]-sig1[i+n1])<0.00000001)
		  {	if(flag==0) {fprintf(stderr,"==========here manmade sig!==============\n");}
			  flag=flag+1;
			  sig3[i]=0.0;
			  sig3[i+1]=0.0;
		//	  continue;
		  }
		//  sig3[i] = 1;
		  if(fabs(sig3[i]-1.0)<0.1)
		  {
			  sig4[kk]=sig1[i+n1];
			  kk++;
			}
		 
	  }
	  npts_real=kk;
	  
	  printf("n=%g n4=%d\n",n,npts_real);
//	  abort();

//     for(i = 0;i<n;i++)  sig3[i] = 1; //the original definition for sig3. can't figure the manmade sig in sig1
//	  printf("finish initialize sig3\n");
     shd1.npts = n;
     shd1.nzyear = sd->ev[ne].yy;
     shd1.nzjday = sd->ev[ne].jday;
     shd1.nzhour = 0;
     shd1.nzmin = 0;
     shd1.nzsec = 0;
     shd1.nzmsec = 0;
     shd1.b = t1;
	  
	  
/*	  int mm;
	  for ( mm = 0; mm < 4; mm++ )
	  {
		  if ( shd1.kstnm[mm] == ' ' ) break;
		  else tempstnm[mm] = shd1.kstnm[mm];
	  } 
	  
	  tempstnm[mm] = '\0';
	  
	
	  if((script=fopen("make_Dir.csh","w"))==NULL)
	  {fprintf(stderr,"cannot creat file make_Dir\n");
		  exit(1);}
	  fprintf(script,"if (! -e %s/%s/%s) then\n mkdir %s/%s/%s\n endif\n",DATAROOT,shd1.knetwk,tempstnm,DATAROOT,shd1.knetwk,tempstnm);
	  fclose(script);
	
	  system("csh make_Dir.csh");
	  system ("/bin/rm make_Dir.csh");
	//  fprintf(stderr,"mkdir %s/%s/%s\n",DATAROOT,shd1.knetwk,tempstnm);
	  
	  sprintf(ft_name,"%s/%s/%s/%s.%s.00.%s.%d.%d",DATAROOT,shd1.knetwk,tempstnm,tempstnm,shd1.knetwk,shd1.kcmpnm,shd1.nzyear,shd1.nzjday);
	  fprintf(stderr,"write SAC:%s\n\n",ft_name);
//	  sprintf(ft_name_bit, "%s_bit",ft_name);
	  sprintf(ft_name_cut,"%s_cut",ft_name);
	
	 
  write_sac (ft_name_cut, &(sig1[n1]), &shd1 );
//  write_sac (ft_name_bit, &(sig3[0]), &shd1 );
///	  write_sac(ft_name,&(sig1[n1]),&shd1);
	  shd1.npts = npts_real;
	  write_sac(ft_name,&(sig4[0]),&shd1);
  */
	  //++++
	  write_sac (ft_name, &(sig1[n1]), &shd1 );
//	  write_sac (ft_name_bit, &(sig3[0]), &shd1 );
	  system("/bin/rm s1.sac");

  }
  


}


/*--------------------------------------------------------------------------*/
void one_rec_trans( SAC_DB *sd, int ne, int ns)
/*--------------------------------------------------------------------------*/
{
  FILE *ff;
  float fl1, fl2, fl3, fl4;
  int n1, n2;
  char str[300];

/* ASSUME THAT THE DATA ARE WITHIN THE FOLLOWING FILTER BAND */
  fl1=1.0/170.0;		/* currently set for 5-150 s period band */
  fl2=1.0/150.0;
  fl3=1.0/5.0;
  fl4=1.0/4.0;

  if ( (fl1<=0.)||(fl2<=0.)||(fl3<=0.)||(fl4<=0.)||(fl1>=fl2)||(fl2>=fl3)||(fl3>=fl4) ) {
    fprintf(stderr,"Error with frequency limits for transfer from evalresp.\n");
    exit(1);
  }
  else {
    fprintf(stderr,"Frequency limits: %f %f %f %f\n", fl1, fl2, fl3, fl4);
  }

  if ( ne >= sd->nev ) return;
  if ( ns >= sd->nst  ) return;
  if ( sd->rec[ne][ns].n <= 0 ) return;

  sprintf(str,"/bin/cp %s s1.sac", sd->rec[ne][ns].fname ); 
//  printf("/bin/cp %s s1.sac", sd->rec[ne][ns].fname );
  system(str);

/* GENERATE TEMPORARY FILES */
  sprintf(str,"/bin/cp %s resp1", sd->rec[ne][ns].resp_fname );
//  fprintf(stderr,"/bin/cp %s resp1\n", sd->rec[ne][ns].resp_fname );
 
  system(str);

  ff = fopen("sac_bp_respcor","w");
  fprintf(ff,"/home/weisen/progs/sac/bin/sac << END\n");
  fprintf(ff,"r s1.sac\n");
  fprintf(ff,"rmean\n");
  fprintf(ff,"rtrend\n");
  //fprintf(ff,"transfer from polezero subtype resp1 to vel freqlimits %f %f %f %f\n", fl1, fl2, fl3, fl4 );
  fprintf(ff,"transfer from  EVALRESP FNAME  resp1 to vel freqlimits %f %f %f %f\n", fl1, fl2, fl3, fl4 );
  fprintf(ff,"w s1.sac\n");
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");
  fclose(ff);

  system("sh sac_bp_respcor");
  system("/bin/rm resp1");
}


/*--------------------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------------------*/
{
  FILE *ff;
  int ne, ns;
  float t1, npts;
  static SAC_DB sdb;
//	char DATAROOT[300];
	
/* CHECK INPUT ARGUMENTS */
  if ( argc < 3 ) {
    fprintf(stderr,"USAGE: cut_trans [t1] [npts]\n");
	  printf("we need the sac_db.out first, to define \n 1)the path/formate to resp \n 2) the channel going to be used AND \n 3)path/name format of the seis_SAC file, the de_resp_SAC file: sdb->rec[ne][ns].fname/ft_fname/resp_fname\n");
    exit(1);
  }
/*  t1=atof(argv[1]);
  npts=atof(argv[2]); */

  sscanf(argv[1],"%f",&t1);
  sscanf(argv[2],"%f",&npts);
///  sscanf(argv[3],"%s",&DATAROOT[0]);

fprintf(stderr,"t1-%f. npts-%f.\n", t1, npts);
fprintf(stderr,"The program assumes the results are within the 5-150 s period band.\n");

/* OPEN SAC DATABASE FILE AND READ IN TO MEMORY */
  if((ff = fopen("sac_db.out", "rb"))==NULL) {
    fprintf(stderr,"sac_db.out file not found\n");
    exit(1);
  }
  fread(&sdb, sizeof(SAC_DB), 1, ff);
  fclose(ff);
  printf("events->%d  stations-> %d\n", sdb.nev,sdb.nst);

/* REMOVE INSTRUMENT RESPONSE AND CUT TO DESIRED LENGTH */
  for ( ns = 0; ns < sdb.nst; ns++ ) for ( ne = 0; ne < sdb.nev; ne++ ) {
	fprintf(stderr,"~~~~~~~~~~~begin to tran ne:%d ns:%d~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",ne,ns);
//+++++	
	  one_rec_trans( &sdb, ne, ns); //++if not need do the trans, also change another sentence noted with +++++
   // fprintf(stderr,"back to main prog\n");
	//  fprintf(stderr,"~~~~~~~~~~~begin to cut ne:%d ns:%d~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",ne,ns);
    one_rec_cut( &sdb, ne, ns, t1, npts);
	//  one_rec_reverse();
	  fprintf(stderr,"==========end of one cut for ne:%d ns:%d========================= \n \n",ne,ns);
  }

  return 0;
}
