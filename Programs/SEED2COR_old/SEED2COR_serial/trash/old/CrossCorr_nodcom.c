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

// station record structure
struct starec {
   SAC_HD shd;
   float *real, *imag, dt, delta;
   int ns, nr, rec[2][1000];
};

/* Finction prorotypes */

extern "C"{
void dmultifft_(int *len, float *real1, float *imag1, float *real2, float *imag2, float *seis_out, int *ns);
}

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec);

void ReadGroup(SAC_DB *sdb, int rec_flag, int iev, int size, int ig, struct starec **stag) {
   float *am, *ph;
   char am_sac[200], ph_sac[200], recname[200];
   int i, j, ist = ig*size, nst;
   SAC_HD shdam, shdph;
   fprintf(stderr, "   Reading in records for station group %d:  ", ig);
   for(i=0, nst=0; i<size; i++, ist++) {
      stag[i]->ns = 0;
      if( ist>=sdb->nst ) continue;
      fprintf(stderr, "%s ", sdb->st[ist].name);
      if( sdb->rec[iev][ist].n<=0 ) {
         fprintf(stderr, "*** No record ***" );
         continue;
      }
      sprintf( am_sac, "%s.am", sdb->rec[iev][ist].ft_fname );
      if( read_sac(am_sac, &am, &shdam)==NULL ) {
         fprintf(stderr, "*** Warning: Cannot open file %s ***", am_sac );
         continue;
      }
      sprintf( ph_sac, "%s.ph", sdb->rec[iev][ist].ft_fname );
      if( read_sac(ph_sac, &ph, &shdph)==NULL ) {
         fprintf(stderr, "*** Warning: Cannot open file %s ***", ph_sac );
         continue;
      }
      if( shdam.npts != shdph.npts || shdam.delta != shdph.delta ) {
         fprintf(stderr, "*** Warning: ph-am mismatch ***" );
         continue;
      }
      stag[i]->ns = shdam.npts;
      stag[i]->dt = shdam.delta;
      sprintf(recname, "%s_rec", sdb->rec[iev][ist].ft_fname);
      if( ! read_rec(rec_flag, recname, sdb->rec[iev][ist].n, stag[i]->rec[0], stag[i]->rec[1], &(stag[i]->nr)) ) {
         fprintf(stderr, "*** Warning: Cannot open file %s ***", recname );
         continue;
      }
      for(j=0; j<shdam.npts; j++) {
         stag[i]->real[j] = am[j]*cos(ph[j]);
         stag[i]->imag[j] = am[j]*sin(ph[j]);
      }
      free(am);
      free(ph);
      nst++;
   }
   cout<<endl<<"      "<<nst<<" stations read in"<<endl;
}

int CalcRecCor( struct starec stag1, struct starec stag2, float *cor_rec, int lag, float dt, int mintlen) {
   int t, irec1, irec2;
   int recB, recE;
   for(t=0;t<=lag;t++){
      cor_rec[lag+t]=0; cor_rec[lag-t]=0;
      for(irec1=0;irec1<stag1.nr;irec1++){
         for(irec2=0;irec2<stag2.nr;irec2++){
            if(stag1.rec[0][irec1]>=stag2.rec[1][irec2]-t) continue;
            if(stag1.rec[1][irec1]<=stag2.rec[0][irec2]-t) break;
            recB = max(stag1.rec[0][irec1], stag2.rec[0][irec2]-t);
            recE = min(stag1.rec[1][irec1], stag2.rec[1][irec2]-t);
            cor_rec[lag+t] += recE - recB;
         }
         for(irec2=0;irec2<stag2.nr;irec2++){
            if(stag1.rec[0][irec1]>=stag2.rec[1][irec2]+t) continue;
            if(stag1.rec[1][irec1]<=stag2.rec[0][irec2]+t) break;
            recB = max(stag1.rec[0][irec1],stag2.rec[0][irec2]+t);
            recE = min(stag1.rec[1][irec1],stag2.rec[1][irec2]+t);
            cor_rec[lag-t] += recE - recB;
         }
      }
   }
   cor_rec[lag] /= 2;

   if(cor_rec[0]<mintlen/dt || cor_rec[lag*2]<mintlen/dt){
      fprintf(stderr, "*** cor time less than %d sec. Skipped! ***", mintlen);
      return 0;
   }
   return 1;
}

int CheckPrecNoise(int lag, float dt, float *cor) {
   float noise = 0., cor_pre = 0.;
   int i;
   for( i = lag*9/5; i< 2*lag; i++) noise += pow(cor[i],2);
   noise = sqrt(noise/(lag/5));
   for( i = lag-int(10/dt); i< lag+int(10/dt); i++) cor_pre += pow(cor[i],2);
   cor_pre = sqrt(cor_pre/(20/dt));
   if( cor_pre > noise*5) {
      printf("       Too much precursor signal. Skipped!\n");
      return 0;
   }
   return 1;
}

int DoCor (SAC_DB *sdb, char *dirname, int lagtime, int rec_flag, int is1, struct starec *stag1, float dt1, int is2, struct starec *stag2, float dt2, int *dnum, int *dflag, int mintlen) {
   fprintf(stderr, "%s-%s ", sdb->st[is1].name, sdb->st[is2].name);
   if(dt1!=dt2) {
      fprintf(stderr, "*** Warning: Incompatible record! ***");
      return 0;
   }
   int i, ns;
   int len = stag1->ns, lag = (int)floor(lagtime/dt1+0.5);
   float *sig, seis_out[2*len], cor_rec[2*lag+1], cor[2*lag+1];
   SAC_HD shd;

   dmultifft_(&len, stag1->real, stag1->imag, stag2->real, stag2->imag, seis_out, &ns);

   if( !CalcRecCor( *stag1, *stag2, cor_rec, lag, dt1, mintlen) ) return 0;

   for( i = 1; i< (lag+1); i++) {
      cor[lag-i] =  seis_out[i]/cor_rec[lag-i];
      cor[lag+i] =  seis_out[ns-i]/cor_rec[lag+i];
   }
   cor[lag] = seis_out[0]/cor_rec[lag];
            
   if( !CheckPrecNoise(lag, dt1, cor) ) return 0;

   *dnum = *dnum + 1;
   *dflag = 1;

   //Update SAC file
   char filename[200];
   sprintf(filename, "%s/%s/COR_%s_%s.SAC", dirname, sdb->st[is1].name, sdb->st[is1].name, sdb->st[is2].name);
   if( *dnum>1 ) {
      if ( read_sac(filename, &sig, &shd)==NULL ) {
         fprintf(stderr,"*** Warning: Cannot read file %s ***", filename );
         return 0;
      }
      for(i = 0; i < (2*lag+1); i++) sig[i] += cor[i];
      shd.user0 = *dnum;
      write_sac (filename, sig, &shd);
      free(sig);
   }
   else {
      shd.delta = dt1;
      shd.evla =  sdb->st[is1].lat;
      shd.evlo =  sdb->st[is1].lon;
      shd.stla =  sdb->st[is2].lat;
      shd.stlo =  sdb->st[is2].lon;
      shd.npts =  2*lag+1;
      shd.b    = -(lag)*dt1;
      shd.user0 = 0;
      shd.nzjday = 1;
      write_sac (filename, cor, &shd);
   }

   return 1;
}

int GroupSize(SAC_DB *sdb, int *spnpts) {
   long pnum = sysconf(_SC_PHYS_PAGES);
   long psize = sysconf(_SC_PAGE_SIZE);
   int iev, ist, nmax=0, npts;
   for(iev=0;iev<NEVENTS;iev++)
      for(ist=0;ist<sdb->nst;ist++)
         if(sdb->rec[iev][ist].n > nmax) nmax = sdb->rec[iev][ist].n;
   npts = (int)pow(2., floor(log(nmax)/log(2.)))+10;
   *spnpts = npts;
   int size = (int)floor(( 0.2*pnum*psize - 10000000 - 6*npts ) / 4. / (npts+1010));
   if(size > 0.5*sdb->nst) size = sdb->nst;
   return size;
}


void CrossCorr(int imo, SAC_DB *sdb, int lagtime, int tnorm_flag, int ftlen, int mintlen, int fskip) {
   int i1, i2, ig1, ig2, is1, is2, iev, npair;
   int dnum[sdb->nst][sdb->nst], dflag[sdb->nst][sdb->nst][NEVENTS];
   char str[600], dirname[200];

   int spnpts, size = GroupSize(sdb, &spnpts);
   //int size = sdb->nst;
   int ng = (int)ceil(sdb->nst/size);
   struct starec *stag1[size], *stag2[size];
   for(i1=0;i1<size;i1++) {
      stag1[i1] = (starec *) malloc ( sizeof(starec) );
      stag1[i1]->real = (float *) malloc ( spnpts * sizeof(float));
      stag1[i1]->imag = (float *) malloc ( spnpts * sizeof(float));
//      stag1[i1]->real = NULL; stag1[i1]->imag = NULL;
      if(ng==1) continue;
      stag2[i1] = (starec *) malloc ( sizeof(starec) );
      stag2[i1]->real = (float *) malloc ( spnpts * sizeof(float));
      stag2[i1]->imag = (float *) malloc ( spnpts * sizeof(float));
//      stag2[i1]->real = NULL; stag2[i1]->imag = NULL;
   }

   sprintf(dirname, "%s/COR", sdb->mo[imo].name);
   if( access(dirname, 0) == 0 ) {
      sprintf(str, "rm -rf %s_old", dirname);
      system(str);
      sprintf(str, "mv %s %s_old", dirname, dirname);
      system(str);
   }
   sprintf(str, "mkdir %s", dirname);
   system(str);
   for(is1=0; is1<sdb->nst; is1++) {
      sprintf(str, "mkdir -p %s/%s", dirname, sdb->st[is1].name);
      system(str);
      for(is2=0; is2<sdb->nst; is2++) {
         dnum[is1][is2]=0;
         for( iev = 0; iev < NEVENTS; iev++ ) dflag[is1][is2][iev]=0;
      }
   }

   for (iev = 0; iev < NEVENTS; iev++) {
      //if(strcmp(sdb->mo[imo].seedf[iev],"0")==0) continue;
      for ( ig1 = 0; ig1 < ng; ig1++ ) {
         ReadGroup(sdb, ftlen, iev, size, ig1, stag1);
         for ( ig2 = ig1; ig2 < ng; ig2++ ) {
            if(ig2!=ig1) ReadGroup(sdb, ftlen, iev, size, ig2, stag2);
            fprintf(stderr, "   Cross-Correlating for event %s between group %d and %d:  ", sdb->ev[iev].name, ig1, ig2);
            npair = 0;
            for(i1=0; i1<size; i1++) {
               if(stag1[i1]->ns==0) continue;
               is1 = i1+ig1*size;
               //dcommon_(&(stag1[i1]->ns), stag1[i1]->real, stag1[i1]->imag);
               for(i2=0; i2<size; i2++) {
                  is2 = i2+ig2*size;
                  if( ig1==ig2 ) {
                     if( is1>=is2 || stag1[i2]->ns==0 ) continue;
                     if( DoCor(sdb, dirname, lagtime, ftlen, is1, stag1[i1], sdb->rec[iev][is1].dt, 
                        is2, stag1[i2], sdb->rec[iev][is2].dt, &(dnum[is1][is2]), 
                        &(dflag[is1][is2][iev]), mintlen) ) npair++;
                  }
                  else {
                     if(stag2[i2]->ns==0) continue;
                     if( DoCor(sdb, dirname, lagtime, ftlen, is1, stag1[i1], sdb->rec[iev][is1].dt, 
                        is2, stag2[i2], sdb->rec[iev][is2].dt, &(dnum[is1][is2]), 
                        &(dflag[is1][is2][iev]), mintlen) ) npair++;
                  }
               }
            }
            cout<<endl<<npair<<" station pairs processed."<<endl;
         }
      }
   }

   for(is1=0;is1<size;is1++) {
      free(stag1[is1]->real); free(stag1[is1]->imag); free(stag1[is1]);
      if(ng==1) continue;
      free(stag2[is1]->real); free(stag2[is1]->imag); free(stag2[is1]);
   }

   char filename[200];
   FILE *frec;
   sprintf(filename, "%s/Cor_dayflag.lst", dirname);
   frec = fopen(filename,"w");
   for( is1 = 0; is1 < sdb->nst-1; is1++ )
      for( is2 = (is1+1); is2 < sdb->nst; is2++ ) {
         sprintf(filename, "%s/COR_%s_%s.SAC", sdb->st[is1].name, sdb->st[is1].name, sdb->st[is2].name);
         fprintf(frec,"%s   ", filename);
         for( iev = 0; iev < NEVENTS; iev++ )
            fprintf(frec,"%d", dflag[is1][is2][iev]);
         fprintf(frec,"\n");
      }
   fclose(frec);
}
