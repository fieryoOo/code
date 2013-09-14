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
   float *amp, *pha, dtrec;
   int nr, rec[2][1000];
};

/* Finction prorotypes */

extern "C"{
void dcommon_(int *len, float *amp, float *pha);

void dmultifft_(int *len, float *amp2, float *pha2, float *seis_out, int *ns);
}

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec);

int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

void ReadGroup(SAC_DB *sdb, int rec_flag, int iev, int size, int ig, struct starec **stag) {
   float *am, *ph;
   char am_sac[200], ph_sac[200], recname[200];
   int i, j, ist = ig*size, nst;
   SAC_HD shdam, shdph;
   fprintf(stderr, "### Reading in records for station group %d: ", ig);
   for(i=0, nst=0; i<size; i++, ist++) {
      stag[i]->nr = 0;
      if( ist>=sdb->nst ) continue;
      if( i%5 == 0 ) fprintf(stderr, "\n   ");
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
      stag[i]->shd = shdam;
      stag[i]->dtrec = sdb->rec[iev][ist].dt;
      if(rec_flag==2) sprintf(recname, "%s_rec2", sdb->rec[iev][ist].ft_fname);
      else sprintf(recname, "%s_rec", sdb->rec[iev][ist].ft_fname);
      if( ! read_rec(rec_flag, recname, sdb->rec[iev][ist].n, stag[i]->rec[0], stag[i]->rec[1], &(stag[i]->nr)) ) {
         fprintf(stderr, "*** Warning: Cannot open file %s ***", recname );
         continue;
      }
      for(j=0; j<shdam.npts; j++) {
         stag[i]->amp[j] = am[j];//am[j]*cos(ph[j]);
         stag[i]->pha[j] = ph[j];//am[j]*sin(ph[j]);
      }
      free(am);
      free(ph);
      nst++;
   }
   cout<<endl<<"   "<<nst<<" stations read in. ###"<<endl;
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

int CheckPrecNoise(int lag, float dt, float *cor, float dist) {
   float noise = 0., cor_pre = 0.;
   int i, ndis, nb;

   ndis = (int)floor((dist/0.8+50.)/dt+0.5);
   if( ndis > lag ) return 1;
   nb = lag*4/5;
   if( ndis > nb ) nb = ndis;
   for( i = lag+nb; i< 2*lag; i++) noise += pow(cor[i],2);
   noise = sqrt(noise/(lag-nb));

   ndis = (int)floor((dist/4.5-50.)/dt+0.5);
   if( ndis < 10./dt ) return 1;
   nb = int(100./dt);
   if(ndis < nb ) nb = ndis;
   for( i = lag-nb; i<= lag+nb; i++) cor_pre += pow(cor[i],2);
   cor_pre = sqrt(cor_pre/(2.*nb));
   if( cor_pre > noise*5) {
      fprintf(stderr, "*** Warning: Large precursor signal. Skipped! ***");
      return 0;
   }
   return 1;
}

int ComprRec(struct starec *stag1, struct starec *stag2) {
   if( stag1->dtrec!=stag2->dtrec ) return 0;
   if( stag1->shd.delta!=stag2->shd.delta ) return 0;
   if( stag1->shd.nzyear!=stag2->shd.nzyear ) return 0;
   if( stag1->shd.nzjday!=stag2->shd.nzjday ) return 0;
   if( stag1->shd.nzhour!=stag2->shd.nzhour ) return 0;
   if( stag1->shd.nzmin!=stag2->shd.nzmin ) return 0;
   if( stag1->shd.nzsec!=stag2->shd.nzsec ) return 0;
   if( abs(stag1->shd.nzmsec-stag2->shd.nzmsec)>1 ) return 0;
   return 1;
}

int DoCor (int fprcs, SAC_DB *sdb, char *dirname, int lagtime, int is1, struct starec *stag1, int is2, struct starec *stag2, short *dnum, short *dflag, int mintlen) {
   fprintf(stderr, "%s-%s ", sdb->st[is1].name, sdb->st[is2].name);
   float dt1 = stag1->dtrec;
   if( !ComprRec(stag1, stag2) ) {
      fprintf(stderr, "*** Warning: Incompatible record! Skipped. ***");
      return 0;
   }
   int i, ns;
   int len = stag1->shd.npts, lag = (int)floor(lagtime/dt1+0.5);
   float *sig, seis_out[2*len];

   dmultifft_(&len, stag2->amp, stag2->pha, seis_out, &ns);

   if(lag>ns/2) { 
      fprintf(stderr, "*** Warning: lagtime overflow, corrected back to max len ***");
      lag = ns/2;
   }
   float cor_rec[2*lag+1], cor[2*lag+1];
   if( !CalcRecCor( *stag1, *stag2, cor_rec, lag, dt1, mintlen) ) return 0;

   for( i = 1; i< (lag+1); i++) {
      cor[lag-i] =  seis_out[i]/cor_rec[lag-i];
      cor[lag+i] =  seis_out[ns-i]/cor_rec[lag+i];
   }
   cor[lag] = seis_out[0]/cor_rec[lag];
            
   //Update SAC file
   double distmp;
   char filename[200];
   SAC_HD shd;
   sprintf(filename, "%s/%s/COR_%s_%s.SAC", dirname, sdb->st[is1].name, sdb->st[is1].name, sdb->st[is2].name);
   if( *dnum>0 ) {
      if ( read_sac(filename, &sig, &shd)==NULL ) {
         fprintf(stderr,"*** Warning: Cannot read file %s ***", filename );
         return 0;
      }
      if( shd.dist == -12345.) {
         calc_dist((double)shd.evla, (double)shd.evlo, (double)shd.stla, (double)shd.stlo, &distmp);
         shd.dist = distmp;
      }
      if( fprcs ) if( !CheckPrecNoise(lag, dt1, cor, shd.dist) ) return 0;
      if( dt1!=shd.delta ) {
         fprintf(stderr,"*** Warning: Incompatible sampling rate! Skipped. ***");
         return 0;
      }
      *dnum = *dnum + 1;
      *dflag = 1;
      for(i = 0; i < (2*lag+1); i++) sig[i] += cor[i];
      shd.user0 = *dnum;
      write_sac (filename, sig, &shd);
      free(sig);
   }
   else {
      shd = stag2->shd;
      shd.delta = dt1;
      strcpy(shd.kevnm, sdb->st[is1].name);
      shd.evla =  sdb->st[is1].lat;
      shd.evlo =  sdb->st[is1].lon;
      calc_dist((double)shd.evla, (double)shd.evlo, (double)shd.stla, (double)shd.stlo, &distmp);
      shd.dist = distmp;
      if( fprcs ) if( !CheckPrecNoise(lag, dt1, cor, shd.dist) ) return 0;
      //shd.stla =  sdb->st[is2].lat;
      //shd.stlo =  sdb->st[is2].lon;
      shd.npts =  2*lag+1;
      shd.e = lag*dt1;
      shd.b    = -shd.e;
      *dnum = 1;
      *dflag = 1;
      shd.user0 = 1;
      shd.nzjday = 1;
      write_sac (filename, cor, &shd);
   }

   return 1;
}

int GroupSize(SAC_DB *sdb, float memomax, int *spnpts) {
   long pnum = sysconf(_SC_PHYS_PAGES);
   long psize = sysconf(_SC_PAGE_SIZE);
   int iev, ist, nmax=0, npts;
   for(iev=0;iev<NEVENTS;iev++)
      for(ist=0;ist<sdb->nst;ist++)
         if(sdb->rec[iev][ist].n > nmax) nmax = sdb->rec[iev][ist].n;
   npts = (int)pow(2., floor(log(nmax)/log(2.)))+10;
   *spnpts = npts;
   int size = (int)floor((memomax*pnum*psize-(10000000.+6.*npts)*sizeof(float)-(sdb->nst*sdb->nst*(NEVENTS+1)*sizeof(short)))/2./(sizeof(starec)+2.*npts*sizeof(float)));
   if(size > 0.5*sdb->nst) size = sdb->nst;
   return size;
}

int ExtractFlag(SAC_DB *sdb, char *dirname, char *str, int *is1, int *is2, char *pflagout) {
   char *name = strtok(str, " "), fname[200];
   sprintf(fname, "%s/%s", dirname, name);
   if( access( fname, R_OK) != 0 ) return -1;
   char *pflag = strtok(NULL, " ");
   strtok(name, "/");
   strtok(NULL, "_");
   char *sta1 = strtok(NULL, "_");
   char *sta2 = strtok(NULL, ".");
   int i, ii;
   for(i=0;i<sdb->nst;i++) {
      ii = *is1-i;
      if(ii>=0) if(strcmp(sdb->st[ii].name,sta1)==0) { *is1 = ii; break; }
      ii = *is1+i;
      if(ii<sdb->nst) if(strcmp(sdb->st[ii].name,sta1)==0) { *is1 = ii; break; }
   }
   if(i>=sdb->nst) return 0;
   for(i=0;i<sdb->nst;i++) {
      ii = *is2-i;
      if(ii>=0) if(strcmp(sdb->st[ii].name,sta2)==0) { *is2 = ii; break; }
      ii = *is2+i;
      if(ii<sdb->nst) if(strcmp(sdb->st[ii].name,sta2)==0) { *is2 = ii; break; }
   }
   if(i>=sdb->nst) return 0;
   sprintf(pflagout, "%s", pflag);
   return 1;
}

void PrepareCor(SAC_DB *sdb, char *dirname, int fskip, short ***dnum, short ****dflag) {
   char str[600];
   FILE *fd;
   int is1, is2, iev;
   int dex = 1;
   //short dnum[sdb->nst][sdb->nst], dflag[sdb->nst][sdb->nst][NEVENTS];
   *dnum = (short **) malloc (sdb->nst * sizeof(short *));
   *dflag = (short ***) malloc (sdb->nst * sizeof(short **));
   for(is1=0; is1<sdb->nst; is1++) {
      (*dnum)[is1] = (short *) malloc (sdb->nst * sizeof(short));
      (*dflag)[is1] = (short **) malloc (sdb->nst * sizeof(short *));
      for(is2=0; is2<sdb->nst; is2++) {
	 (*dnum)[is1][is2]=0; //mark dnum 0 for the station pairs to be processed.
	 if( sdb->st[is1].flag == sdb->st[is2].flag ) { //skip if both in group 0
	    if( sdb->st[is1].flag == 0 ) (*dnum)[is1][is2]=-1;
	 }
	 else if( sdb->st[is1].flag * sdb->st[is2].flag != 0 ) (*dnum)[is1][is2]=-1; //skip if in different groups that are both not group 0
         (*dflag)[is1][is2] = (short *) malloc (NEVENTS * sizeof(short));
         for( iev = 0; iev < NEVENTS; iev++ ) (*dflag)[is1][is2][iev]=0;
      }
   }

   if( access(dirname, 0) != 0 ) { 
      dex = 0; 
      sprintf(str, "mkdir %s", dirname);
      system(str);
      if(fskip)
         cout<<"   Cannot access "<<dirname<<". No Correlation will be skipped!"<<endl;
   }
   for(is1=0; is1<sdb->nst; is1++) {
      sprintf(str, "mkdir -p %s/%s", dirname, sdb->st[is1].name);
      system(str);
   }
   if(!dex) return;

   sprintf(str, "%s/Cor_dayflag.lst", dirname);
   if((fd=fopen(str,"r"))==NULL) {
      if(fskip)
         cout<<"   Cannot access "<<str<<". No Correlation will be skipped!"<<endl;
      return;
   }
   sprintf(str, "cp %s/Cor_dayflag.lst %s/Cor_dayflag.lst_old", dirname, dirname);
   system(str);
   char pflag[NEVENTS+1];
   if(fskip) {
      for(is1=0,is2=0; fgets(str, 300, fd)!=NULL; ) {
         if( ExtractFlag(sdb, dirname, str, &is1, &is2, pflag)<=0 ) continue;
         //for(iev=0;iev<NEVENTS;iev++) (*dflag)[is1][is2][iev] = pflag[iev]-'0';
         (*dnum)[is1][is2] = -1;
      }
      sprintf(str, "rm -f %s/Cor_dayflag.lst_tmp", dirname);
      system(str);
   }
   else {
      int sflag;
      char str2[600];
      FILE *fd2;
      fclose(fd);
      sprintf(str, "%s/Cor_dayflag.lst_old", dirname);
      fd = fopen(str, "r");
      sprintf(str, "%s/Cor_dayflag.lst_tmp", dirname);
      fd2 = fopen(str, "w");
      for(is1=0,is2=0; fgets(str, 300, fd)!=NULL; ) {
         sprintf(str2, "%s", str);
	 sflag = ExtractFlag(sdb, dirname, str2, &is1, &is2, pflag);
         if( (sflag==1 && (*dnum)[is1][is2] == -1) || sflag==0 )
            fprintf(fd2, "%s", str);
      }
      fclose(fd2);
      /*
      sprintf(str, "rm -rf %s_old", dirname);
      system(str);
      sprintf(str, "mv %s %s_old", dirname, dirname);
      system(str);
      sprintf(str, "mkdir %s", dirname);
      system(str);
      */
   }
   fclose(fd);

}

void CrossCorr(int imo, SAC_DB *sdb, int lagtime, float memomax, int ftlen, int fprcs, int mintlen, int fskip) {
   if(fskip==2) return;
   int i1, i2, ig1, ig2, is1, is2, iev, npair;
   short **dnum, ***dflag;
   char dirname[200];

   //Check memory size and allocate memories for station groups
   int spnpts, size = GroupSize(sdb, memomax, &spnpts);
   int ng = (int)ceil(sdb->nst/size);
   struct starec *stag1[size], *stag2[size];
   for(i1=0;i1<size;i1++) {
      stag1[i1] = (starec *) malloc ( sizeof(starec) );
      stag1[i1]->amp = (float *) malloc ( spnpts * sizeof(float));
      stag1[i1]->pha = (float *) malloc ( spnpts * sizeof(float));
      if(ng==1) continue;
      stag2[i1] = (starec *) malloc ( sizeof(starec) );
      stag2[i1]->amp = (float *) malloc ( spnpts * sizeof(float));
      stag2[i1]->pha = (float *) malloc ( spnpts * sizeof(float));
   }

   //Initialize day-num-flag list
   sprintf(dirname, "%s/COR", sdb->mo[imo].name);
   PrepareCor(sdb, dirname, fskip, &dnum, &dflag);

   //do cross-correlatings for each group pairs
   for (iev = 0; iev < NEVENTS; iev++) {
      if(strcmp(sdb->mo[imo].seedf[iev],"0")==0) continue;
      for ( ig1 = 0; ig1 < ng; ig1++ ) {
         ReadGroup(sdb, ftlen, iev, size, ig1, stag1);
         for ( ig2 = ig1; ig2 < ng; ig2++ ) {
            if(ig2!=ig1) ReadGroup(sdb, ftlen, iev, size, ig2, stag2);
            fprintf(stderr, "### Cross-Correlating for event %s between group %d and %d: \n   ", sdb->ev[iev].name, ig1, ig2);
            npair = 0;
            for(i1=0; i1<size; i1++) {
               if(stag1[i1]->nr==0) continue;
               is1 = i1+ig1*size;
               dcommon_(&(stag1[i1]->shd.npts), stag1[i1]->amp, stag1[i1]->pha);
               for(i2=0; i2<size; i2++) {
                  is2 = i2+ig2*size;
		  if(dnum[is1][is2]==-1) continue;
                  if( ig1==ig2 ) {
                     if( is1>=is2 || stag1[i2]->nr==0 ) continue;
                     if( DoCor(fprcs, sdb, dirname, lagtime, is1, stag1[i1],is2, stag1[i2], 
                         &(dnum[is1][is2]), &(dflag[is1][is2][iev]), mintlen) ) {
			npair++;
			if( npair%10 == 0) fprintf(stderr, "\n   ");
		     }
                  }
                  else {
                     if(stag2[i2]->nr==0) continue;
                     if( DoCor(fprcs, sdb, dirname, lagtime, is1, stag1[i1], is2, stag2[i2], 
                         &(dnum[is1][is2]), &(dflag[is1][is2][iev]), mintlen) ) {
			npair++;
			if( npair%10 == 0) fprintf(stderr, "\n   ");
		     }
                  }
               }
            }
	    if( npair==0 || npair%10 != 0 ) fprintf(stderr, "\n   ");
            cout<<npair<<" station pairs processed. ###"<<endl;
         }
      }
   }

   //clean up group memories
   for(is1=0;is1<size;is1++) {
      free(stag1[is1]->amp); free(stag1[is1]->pha); free(stag1[is1]);
      if( ng==1 ) continue;
      free(stag2[is1]->amp); free(stag2[is1]->pha); free(stag2[is1]);
   }

   //output day-num-flag list
   char filename[200], str[300];
   FILE *frec;
   sprintf(filename, "%s/Cor_dayflag.lst", dirname);
   sprintf(str, "mv %s_tmp %s", filename, filename);
   system(str);
   frec = fopen(filename,"a");
   for( is1 = 0; is1 < sdb->nst-1; is1++ ) {
      for( is2 = (is1+1); is2 < sdb->nst; is2++ ) {
         if(dnum[is1][is2]>0) {
            sprintf(filename, "%s/COR_%s_%s.SAC", sdb->st[is1].name, sdb->st[is1].name, sdb->st[is2].name);
            fprintf(frec,"%s   ", filename);
            for( iev = 0; iev < NEVENTS; iev++ ) fprintf(frec,"%d", dflag[is1][is2][iev]);
            fprintf(frec,"\n");
         }
         free(dflag[is1][is2]);
      }
      free(dnum[is1]); free(dflag[is1]);
   }
   fclose(frec);
   free(dnum); free(dflag);
}
