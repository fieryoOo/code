// this version does cross-correlations after  judging whether a cross-correlatins exists. 
// If a cross-correlation of one pair of stations exists,the code skips that pair. 
// Revised to retain amplitude info. Read rec files to correct for correlation time length.
//Output list file with Stack day number.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <math.h>
#include "mysac64.h"
#include "64_sac_db.h"
using namespace std;

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

/* Finction prorotypes */

extern "C"{
void dcommon_(int *len, float *amp,float *phase);

void dmultifft_(int *len,float *amp,float *phase,float *seis_out, int *ns);
}

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec);

/*************************************************************************/
 int check_info ( SAC_DB *sdb, int ne, int ns1, int ns2 )
/*************************************************************************/
{
  //fprintf(stderr, "ne %d ns1 %d ns2 %d\n", ne, ns1, ns2);
  if ( ne >= NEVENTS ) {
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

/*--------------------------------------------------------------------------*/
 int do_cor( SAC_DB *sdb, int lagtime, int rec_flag )
/*--------------------------------------------------------------------------*/
{

  int ine, jsta1, jsta2;
  int lag, lagmax = 0, len, lenmax = 0;
  int ns, i, k, t, irec1, irec2, nrec1, nrec2, recB, recE; 
  int rec_b1[1000],rec_e1[1000],rec_b2[1000],rec_e2[1000],day_num[sdb->nst][sdb->nst], day_flag[sdb->nst][sdb->nst][NEVENTS];
  char filename[200], filename_temp[200], recname[200], amp_sac[200], phase_sac[200], filename_tp[100],filename_tp0[100], str[600];
  float *sig, *amp, *phase;
  float *cor_temp = NULL, *cor = NULL, *seis_out = NULL;
  float cor_pre, noise;
  FILE *ff,*fsac,*frec;
  SAC_HD shdamp1, shdph1, shdamp2, shdph2, shd_cor,shd_temp;


   for( jsta1 = 0; jsta1 < sdb->nst; jsta1++ ) {
      sprintf(str, "mkdir -p COR/%s",sdb->st[jsta1].name);
      system(str);
      for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ){
          day_num[jsta1][jsta2]=0;
          for( ine = 0; ine < NEVENTS; ine++ ) day_flag[jsta1][jsta2][ine]=0;
      }
   }

   for( ine = 0; ine < NEVENTS; ine++ ) {
      for( jsta1 = 0; jsta1 < sdb->nst; jsta1++ ) {  
         if(sdb->rec[ine][jsta1].n <= 0) continue;
         sprintf( amp_sac, "%s.am", sdb->rec[ine][jsta1].ft_fname );
         sprintf( phase_sac, "%s.ph", sdb->rec[ine][jsta1].ft_fname );
// read amp and phase files and read into common memory
         if( read_sac(amp_sac, &amp, &shdamp1)==NULL ) {
    	    cout<<"Warning: Cannot open file "<<amp_sac<<endl;
   	    continue;
         }
         if( read_sac(phase_sac, &phase, &shdph1)== NULL ) {
            cout<<"Warning: Cannot open file "<<phase_sac<<endl;
            continue;
         }

	 len = shdamp1.npts;
         dcommon_( &len, amp, phase ); // reads amp and phase files into common memory
         free(amp);
         free(phase);

         sprintf(recname, "%s_rec", sdb->rec[ine][jsta1].ft_fname);
         if( ! read_rec(rec_flag, recname, len, rec_b1, rec_e1, &nrec1) ) {
            cout<<"Warning: Cannot open file "<<sdb->rec[ine][jsta1].ft_fname<<"_rec"<<endl;
            continue;
         }

         for( jsta2 = jsta1+1; jsta2 < sdb->nst; jsta2++ ) {
  	    if(sdb->rec[ine][jsta2].n <= 0) continue;
	    sprintf(amp_sac, "%s.am", sdb->rec[ine][jsta2].ft_fname);
            sprintf(phase_sac, "%s.ph", sdb->rec[ine][jsta2].ft_fname);
	    fprintf(stderr,"   Correlating: %s  %s\n", sdb->rec[ine][jsta1].ft_fname,sdb->rec[ine][jsta2].ft_fname );
            if ( read_sac(amp_sac, &amp, &shdamp2) ==NULL ) {
               cout<<"Warning: Cannot open file "<<amp_sac<<endl;
               continue;
            }
            if ( read_sac(phase_sac, &phase, &shdph2)==NULL ) {
               cout<<"Warning: Cannot open file "<<phase_sac<<endl;
               continue;
            }

	    len = shdamp2.npts;

            sprintf(recname, "%s_rec", sdb->rec[ine][jsta2].ft_fname);
            if( ! read_rec(rec_flag, recname, len, rec_b2, rec_e2, &nrec2) ) {
               cout<<"Warning: Cannot open file "<<sdb->rec[ine][jsta2].ft_fname<<"_rec"<<endl;
               continue;
            }

            if(!check_info(sdb, ine, jsta1, jsta2 )) {
               cout<<"Warning: incompatible am/ph record"<<endl;
               continue;
            }

            lag = (int)floor(lagtime/sdb->rec[ine][jsta2].dt+0.5);
            if(len>lenmax) {
               seis_out = (float *) realloc ( seis_out, 2*len * sizeof(float) );
               lenmax = len;
            }
            if(lag>lagmax) {
               cor_temp = (float *) realloc ( cor_temp, (2*lag+1) * sizeof(float) );
               cor = (float *) realloc ( cor, (2*lag+1) * sizeof(float) );
               lagmax = lag;
            }
            dmultifft_(&len, amp, phase, seis_out,&ns);
            free(amp);
            free(phase);
 //           if(isnan(seis_out[0])) continue;
 //           if(isnan(cor[lag])) continue;
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
               printf("cor time less than 15000 sec. Skipped!\n");
               continue;
            }

            cor[lag] = seis_out[0]/cor_temp[lag];
            for( i = 1; i< (lag+1); i++) { 
               cor[lag-i] =  seis_out[i]/cor_temp[lag-i];
               cor[lag+i] =  seis_out[ns-i]/cor_temp[lag+i];
            }
               
            noise = 0;
	    for( i = lag*9/5; i< 2*lag; i++) noise += pow(cor[i],2);
	    noise = sqrt(noise/(lag/5));
            cor_pre = 0;
            for( i = lag-10; i< lag+10; i++) cor_pre += pow(cor[i],2);
	    cor_pre = sqrt(cor_pre/20.);
            if( cor_pre > noise*5) {
               printf("       Too much precursor signal around 0 time. Skipped!\n");
               continue;
            }

            day_num[jsta1][jsta2]++;
            day_flag[jsta1][jsta2][ine]++;

  	    // move and rename cor file accordingly 
	    sprintf(filename, "COR/%s/COR_%s_%s.SAC.prelim",
	      sdb->st[jsta1].name, sdb->st[jsta1].name, sdb->st[jsta2].name);

          //  sprintf(filename_temp, "COR/COR_%s_%s.SAC.%d_rec",
            //  sdb->st[jsta1].name, sdb->st[jsta2].name, ine+1);

	    if(access(filename, F_OK) == 0) { // if file alread present, do this
	      if ( read_sac (filename, &sig, &shd_cor)==NULL ) {
	        fprintf(stderr,"file %s did not found\n", filename );
	        return 0;
	      }
	      // add new correlation to previous one
 	      for(k = 0; k < (2*lag+1); k++) sig[k] += cor[k];
              shd_cor.user0 = day_num[jsta1][jsta2];
	      write_sac (filename, sig, &shd_cor );
              free(sig);
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

       }   //loop over jsta2

     }  //loop over jsta1

    }  //loop over events
    free(seis_out);
    free(cor_temp);
    free(cor);

    sprintf(filename, "COR/Cor_dayflag.lst");
    frec = fopen(filename,"w");
    for( jsta1 = 0; jsta1 < sdb->nst-1; jsta1++ )
        for( jsta2 = (jsta1+1); jsta2 < sdb->nst; jsta2++ ) {
           sprintf(filename, "%s/COR_%s_%s.SAC", sdb->st[jsta1].name, sdb->st[jsta1].name, sdb->st[jsta2].name);
           fprintf(frec,"%s   ", filename);
           for( ine = 0; ine < NEVENTS; ine++ ) 
              fprintf(frec,"%d", day_flag[jsta1][jsta2][ine]);
           fprintf(frec,"\n");
          }
    fclose(frec);

  return 0;
}


void CrossCorr(int imo, SAC_DB *sdb, int lagtime, int tnorm_flag, int ftlen, int mintlen, int fskip) {
   int ns1, ns2;
   char str[600], filename[200];

//   int rec_flag = 0;
//   if( tnorm_flag==3 ) rec_flag = 1;
   do_cor(sdb,lagtime,ftlen);  

   fprintf(stderr, "finished correlations\n");

   for ( ns1 = 0; ns1 < sdb->nst-1; ns1++ ) {
      for ( ns2 = ns1+1; ns2 < sdb->nst; ns2++ ) {
         sprintf(filename, "COR/%s/COR_%s_%s.SAC.prelim\0", sdb->st[ns1].name, sdb->st[ns1].name, sdb->st[ns2].name);
         sprintf(str, "mv %s COR/%s/COR_%s_%s.SAC", filename, sdb->st[ns1].name, sdb->st[ns1].name, sdb->st[ns2].name);
         if(access(filename, F_OK) == 0) system(str);
      }
   }

   sprintf(str, "mv %s/COR %s/COR_old", sdb->mo[imo].name, sdb->mo[imo].name);
   system(str);
   sprintf(str, "mv COR %s", sdb->mo[imo].name);
   system(str);
  
}
