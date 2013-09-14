// this version does cross-correlations after  judging whether a cross-correlatins exists. 
// If a cross-correlation of one pair of stations exists,the code skips that pair. 
// Revised to retain amplitude info. Read rec files to correct for correlation time length.
//Output list file with Stack day number.

#define MAIN

//#include "/home/zheng/progs/NOISE_CODA/HEAD_NOISE/64_mysac.h"
//#include "/Users/jiayixie/progs/NOISE_CODA/HEAD_NOISE/64_mysac.h"
//#include "/home/zheng/progs/NOISE_CODA/HEAD_NOISE/64_sac_db.h"
//#include "/Users/jiayixie/progs/NOISE_CODA/HEAD_NOISE/64_sac_db.h"
#include <stdio.h>
#include <unistd.h>
#include "/home/tianye/code/Programs/head/mysac.h"
#include "/home/tianye/code/Programs/head/sac_db.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

/* Finction prorotypes */

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

/*************************************************************************/
 int check_info ( SAC_DB *sdb, int ne, int ns1, int ns2 )
/*************************************************************************/
{
  //fprintf(stderr, "ne %d ns1 %d ns2 %d\n", ne, ns1, ns2);
  if ( ne >= sdb->nev ) {
    fprintf(stderr,"cannot make correlation: too large event number\n");
    return 0;
  }
  if ( (ns1>=sdb->nst) ||(ns2>=sdb->nst)  ) {
    fprintf(stderr,"cannot make correlation: too large station number\n");
    return 0;
  }
  if ( sdb->rec[ne][ns1].n <= 0 ) {
    fprintf(stderr,"no data for station %s and event %s\n", sdb->st[ns1].name, sdb->ev[ne].name );
    return 0;
  }
  if ( sdb->rec[ne][ns2].n <= 0 ) {
    fprintf(stderr,"no data for station %s and event %s\n", sdb->st[ns2].name, sdb->ev[ne].name );
    return 0;
  }
  if ( fabs(sdb->rec[ne][ns1].dt-sdb->rec[ne][ns2].dt) > .0001 ) {
    fprintf(stderr,"incompatible DT\n");
    return 0;
  }
  return 1;
}

/*c/////////////////////////////////////////////////////////////////////////*/
 float sig[1000000];
 SAC_HD shdamp1, shdph1, shdamp2, shdph2, shd_cor,shd_temp;

/*--------------------------------------------------------------------------*/
 int do_cor( SAC_DB *sdb, int lag, int rec_flag )
/*--------------------------------------------------------------------------*/
{
  int ine, jsta1, jsta2;

  int len, ns, i, k, t, irec1, irec2, nrec1, nrec2, recB, recE; 

  int rec_b1[1000],rec_e1[1000],rec_b2[1000],rec_e2[1000],day_num[sdb->nst][sdb->nst], day_flag[sdb->nst][sdb->nst][sdb->nev];
  char filename[200], filename_temp[200], rec_name[200], amp_sac[200], phase_sac[200], filename_tp[100],filename_tp0[100], str[600];
  float amp[400000], phase[400000], cor_temp[40000], cor[40000],sactemp[40000];
  float seis_out[400000];
  double cor_pre, noise;
  FILE *ff,*fsac,*frec;

  // outermost loop over day number, then station number
      fprintf(stderr,"sdb->nev %d sdb->nst %d\n",sdb->nev,sdb->nst);

    for( jsta1 = 0; jsta1 < sdb->nst; jsta1++ )
       for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ){
          day_num[jsta1][jsta2]=0;
          for( ine = 0; ine < sdb->nev; ine++ )
             day_flag[jsta1][jsta2][ine]=0;
       }

    for( ine = 0; ine < sdb->nev; ine++ ) {

      fprintf(stderr,"sdb->nev %d\n",ine);

    // loop over "base" station number, this will be stored into common memory

    for( jsta1 = 0; jsta1 < sdb->nst; jsta1++ ) {  
      sprintf(str, "mkdir -p COR/%s",sdb->st[jsta1].name);
      system(str);

      if(sdb->rec[ine][jsta1].n > 0){
        sprintf( amp_sac, "%s.am", sdb->rec[ine][jsta1].ft_fname );
        sprintf( phase_sac, "%s.ph", sdb->rec[ine][jsta1].ft_fname );
// read amp and phase files and read into common memory
        if ( read_sac(amp_sac, amp, &shdamp1, 1000000 )==NULL ) {
  	  fprintf( stderr,"file %s did not found\n", amp_sac );
   	  goto loop1;
        }
        if ( isnan(amp[0]) ) {
	  fprintf(stderr,"%s is not a valid file\n",amp_sac);
           goto loop1;
        } 
        if ( read_sac(phase_sac, phase, &shdph1, 1000000)== NULL ) {
          fprintf( stderr,"file %s did not found\n", phase_sac );
          goto loop1;
        }

        if ( isnan(phase[0]) ) goto loop1;

	len = shdamp1.npts;

        dcommon_( &len, amp, phase ); // reads amp and phase files into common memory

        if( rec_flag ) {
           sprintf( rec_name, "%s_rec", sdb->rec[ine][jsta1].ft_fname );
           if((frec = fopen(rec_name,"r")) == NULL) {
             fprintf( stderr,"file %s did not found\n", rec_name );
             goto loop1;
           }
           for(irec1=0;;irec1++) 
              if(fscanf(frec,"%d	%d",&rec_b1[irec1],&rec_e1[irec1])!=2) break;
           nrec1=irec1;
           fclose(frec);
        }
        else {
           rec_b1[0]=0; rec_e1[0]=len-1;
           nrec1=1;
        }

/*        if ( read_sac(rec_name, rec_sig, &shdamp1, 1000000 )==NULL ) {
          fprintf( stderr,"file %s did not found\n", rec_name );
          goto loop1;
        }
        if ( isnan(rec_sig[0]) ) {
          fprintf(stderr,"%s is not a valid file\n",rec_name);
           goto loop1;
        } */

     for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ) {
         
  	  if(sdb->rec[ine][jsta2].n > 0){

	    // compute correlation
	    sprintf(amp_sac, "%s.am", sdb->rec[ine][jsta2].ft_fname);
            sprintf(phase_sac, "%s.ph", sdb->rec[ine][jsta2].ft_fname);
	  fprintf(stderr,"Correlating: %s  %s\n", sdb->rec[ine][jsta1].ft_fname,sdb->rec[ine][jsta2].ft_fname );
            // get array of floats for amp and phase of first signal

            if ( read_sac(amp_sac, amp, &shdamp2, 100000) ==NULL ) {
              fprintf(stderr,"file %s did not found\n", amp_sac );
              goto loop2;
            }
	    if( isnan(amp[0]) ) { 
	     fprintf(stderr,"%s is not a valid file\n",amp_sac);
             goto loop2;
	    } 
            if ( read_sac(phase_sac, phase, &shdph2, 100000)==NULL ) {
              fprintf(stderr,"file %s did not found\n", phase_sac );
              goto loop2;
            }
            if( isnan(phase[0]) ) goto loop2;      


	      len = shdamp2.npts;

            if( rec_flag) {
               sprintf(rec_name, "%s_rec", sdb->rec[ine][jsta2].ft_fname);
               if((frec = fopen(rec_name,"r")) == NULL) {
                 fprintf( stderr,"file %s did not found\n", rec_name );
                 goto loop2;
                }
               for(irec2=0;;irec2++)
                  if(fscanf(frec,"%d	%d",&rec_b2[irec2],&rec_e2[irec2])!=2) break;
               nrec2=irec2;
               fclose(frec);
            }
            else {
               rec_b2[0]=0; rec_e2[0]=len-1;
               nrec2=1;
            }

/*            if ( read_sac(rec_name, rec_sig2, &shdamp2, 100000) ==NULL ) {
              fprintf(stderr,"file %s did not found\n", rec_name );
              goto loop2;
            }
            if( isnan(rec_sig2[0]) ) {
              fprintf(stderr,"%s is not a valid file\n",rec_name);
             goto loop2;
            }  */


            if(!check_info(sdb, ine, jsta1, jsta2 )) {
              fprintf(stderr,"files incompatible\n");
              return 0;
            }
            else
            {

                 dmultifft_(&len, amp, phase, &lag, seis_out,&ns);
		   if(isnan(seis_out[0])) goto loop2;

/*                 cor_temp[lag]=0;
                 for(i=0;i<len;i++)
                    cor_temp[lag]+=rec_sig[i]*rec_sig2[i];
                 if(isnan(cor[lag])) goto loop2;
                 for(t=1;t<=lag;t++){
                    cor_temp[lag+t]=0; cor_temp[lag-t]=0;
                    for(i=0;i<len-t;i++){
                       cor_temp[lag+t]+=rec_sig[i]*rec_sig2[i+t];
                       cor_temp[lag-t]+=rec_sig[i+t]*rec_sig2[i];
                      }
                   }  */

                 if(isnan(cor[lag])) goto loop2;
                 for(t=0;t<=lag;t++){
                    cor_temp[lag+t]=0; cor_temp[lag-t]=0;
                    for(irec1=0;irec1<nrec1;irec1++){
                       for(irec2=0;irec2<nrec2;irec2++){
                          if(rec_b1[irec1]>=rec_e2[irec2]-t) continue;
                          if(rec_e1[irec1]<=rec_b2[irec2]-t) break;
                          recB = max(rec_b1[irec1],rec_b2[irec2]-t);
                          recE = min(rec_e1[irec1],rec_e2[irec2]-t);
                          cor_temp[lag+t] += recE - recB;
                         }
                       for(irec2=0;irec2<nrec2;irec2++){
                          if(rec_b1[irec1]>=rec_e2[irec2]+t) continue;
                          if(rec_e1[irec1]<=rec_b2[irec2]+t) break;
                          recB = max(rec_b1[irec1],rec_b2[irec2]+t);
                          recE = min(rec_e1[irec1],rec_e2[irec2]+t);
                          cor_temp[lag-t] += recE - recB;
                         }
                      }
                   }
                 cor_temp[lag] /= 2;

                 if(cor_temp[0]<15000/sdb->rec[ine][jsta1].dt || cor_temp[lag*2]<15000/sdb->rec[ine][jsta1].dt){
//                   printf("cor_time[0]: %d  cor_time[leg*2]: %d  delta: %f\n",cor_temp[0],cor_temp[lag*2],sdb->rec[ine][jsta1].dt);
                   printf("cor time less than 15000 sec. Skipped!\n");
                   goto loop2;
                  }

                 cor[lag] = seis_out[0]/cor_temp[lag];
                 for( i = 1; i< (lag+1); i++)
                 { 
     	           cor[lag-i] =  seis_out[i]/cor_temp[lag-i];
	           cor[lag+i] =  seis_out[ns-i]/cor_temp[lag+i];
	         }
               
		 noise = 0;
		 for( i = lag+2000; i< lag+2500; i++) noise += pow(cor[i],2);
		 noise = sqrt(noise/500.);
                 cor_pre = 0;
                 for( i = lag-10; i< lag+10; i++) cor_pre += pow(cor[i],2);
		 cor_pre = sqrt(cor_pre/20.);
                 if( cor_pre > noise*3) {
                    printf("Too much precursor signal around 0 time. Skipped!\n");
                    goto loop2;
                   } 

                 day_num[jsta1][jsta2]++;
                 day_flag[jsta1][jsta2][ine]++;

  	    // move and rename cor file accordingly 
	    sprintf(filename, "COR/%s/COR_%s_%s.SAC.prelim",
	      sdb->st[jsta1].name, sdb->st[jsta1].name, sdb->st[jsta2].name);

          //  sprintf(filename_temp, "COR/COR_%s_%s.SAC.%d_rec",
            //  sdb->st[jsta1].name, sdb->st[jsta2].name, ine+1);

	    if(access(filename, F_OK) == 0) { // if file alread present, do this
	      if ( !read_sac (filename, sig, &shd_cor, 1000000 ) ) {
	        fprintf(stderr,"file %s did not found\n", filename );
	        return 0;
	      }
	      // add new correlation to previous one
 	      for(k = 0; k < (2*lag+1); k++) sig[k] += cor[k];
              shd_cor.user0 = day_num[jsta1][jsta2];
	      write_sac (filename, sig, &shd_cor );
	    }
	    // if file doesn't already exist, use one of the current headers
	    // and change a few values. more may need to be added
	     else {
	      shdamp1.delta = sdb->rec[ine][jsta1].dt;
	      shdamp1.evla =  sdb->st[jsta1].lat;
	      shdamp1.evlo =  sdb->st[jsta1].lon;
	      shdamp1.stla =  sdb->st[jsta2].lat;
	      shdamp1.stlo =  sdb->st[jsta2].lon;
	      shdamp1.npts =  2*lag+1;
              shdamp1.b    = -(lag)*shdamp1.delta;
              shdamp1.user0 = day_num[jsta1][jsta2];
             // shdamp1.nzjday = shdamp2.nzjday - ine;
              shdamp1.nzjday = 1;
	      write_sac (filename, cor, &shdamp1);
           //   write_sac (filename_temp, cor_temp, &shdamp1);
	     }
 	    }   //loop over check

	 }    //loop over if jsta2
      loop2:
        ;
       }   //loop over jsta2

      }  //loop over if jsta1
    loop1:
     ;
     }  //loop over jsta1

    }  //loop over events

    sprintf(filename, "COR/Cor_dayflag.lst");
    frec = fopen(filename,"w");
    for( jsta1 = 0; jsta1 < sdb->nst-1; jsta1++ )
        for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ) {
           sprintf(filename, "%s/COR_%s_%s.SAC", sdb->st[jsta1].name, sdb->st[jsta1].name, sdb->st[jsta2].name);
           fprintf(frec,"%s   ", filename);
           for( ine = 0; ine < sdb->nev; ine++ ) 
              fprintf(frec,"%d", day_flag[jsta1][jsta2][ine]);
           fprintf(frec,"\n");
          }
    fclose(frec);
/*    for( jsta1 = 0; jsta1 < sdb->nst-1; jsta1++ ) {
        sprintf(str, "mkdir -p COR/%s",sdb->st[jsta1].name);
        system(str);
        sprintf(filename, "COR/%s/Cor_daynum.lst", sdb->st[jsta1].name);
        if((frec = fopen(filename,"w")) == NULL) {
           fprintf( stderr,"Can not open Cor_daynum list file to write\n");
           return 0;
          }
        for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ) {
           sprintf(filename, "COR_%s_%s.SAC", sdb->st[jsta1].name, sdb->st[jsta2].name);
           fprintf(frec,"%s   %d\n", filename, day_num[jsta1][jsta2]);
          }
        fclose(frec);
       }  */

  return 0;
}



SAC_DB sdb;

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
int main (int na, char *arg[])
/*--------------------------------------------------------------------------*/
{
  FILE *ff;
  int ns1 = 0, ns2 = 0,lag;
  char str[600], filename[200];

  if ( na != 3 )
    {
      fprintf(stderr,"usage: just_cor [lag] [rec_flag(0 or 1)]\n");
      exit(1);
    }

  sscanf(arg[1],"%d", &lag );


  ff = fopen("sac_db.out","rb");
  fread(&sdb, sizeof(SAC_DB), 1, ff );
  fclose(ff);

fprintf(stderr, "read info fine\n");

  // do all the work of correlations here

  do_cor(&sdb,lag,atoi(arg[2]));  

  fprintf(stderr, "finished correlations\n");

  // move COR/COR_STA1_STA2.SAC.prelim to COR/COR_STA1_STA2.SAC

  for ( ns1 = 0; ns1 < sdb.nst-1; ns1++ ) {
//      sprintf(str, "mkdir -p COR/%s",sdb.st[ns1].name);
//      system(str);
      for ( ns2 = ns1+1; ns2 < sdb.nst; ns2++ ) {
//  for ( ns2 = 1; ns2 < sdb.nst; ns2++ ) for ( ns1 = 0; ns1 < ns2; ns1++ ) {
      sprintf(filename, "COR/%s/COR_%s_%s.SAC.prelim\0", sdb.st[ns1].name, sdb.st[ns1].name, sdb.st[ns2].name);
      sprintf(str, "mv %s COR/%s/COR_%s_%s.SAC", filename, sdb.st[ns1].name, sdb.st[ns1].name, sdb.st[ns2].name);
      if(access(filename, F_OK) == 0) system(str);
      }
    }
  
  return 0;
}
