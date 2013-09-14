#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "mysac64.h"
#include "64_sac_db.h"
#include "Param.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////
char str[1000];
///////////////////////////////////////////////////////////////////////////

extern "C"{
void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,double *dt,int *n,float *seis_in,float *seis_out);
}

SAC_HD *read_shd (char *fname);

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (char *fname, float *sig, SAC_HD *SHD);

void UpdateTime(SAC_HD *shd);

double abs_time ( int yy, int jday, int hh, int mm, int ss, int ms );

void SetName(int ne, int ns) {
   sprintf(sdb->rec[ne][ns].fname,"%s/%s/%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].ft_fname,"%s/%s/ft_%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].chan,"%s", ch );
}

int CheckExistence(int ne, int ns) {
   SAC_HD *shd = read_shd(sdb->rec[ne][ns].fname);
   if( shd==NULL ) return 0;
   sdb->rec[ne][ns].n = shd->npts;
   sdb->rec[ne][ns].t0 = abs_time (shd->nzyear, shd->nzjday, shd->nzhour, shd->nzmin, shd->nzsec, shd->nzmsec );
   sdb->rec[ne][ns].dt = shd->delta;
   return 1;
}

int Resampling(char *sacname, float **sig2, SAC_HD *sd) {
   SAC_HD shd=sac_null;
   float *sig1=NULL;
   if (read_sac(sacname, &sig1, &shd)==NULL) return 0;
   float dt = 1./sps;
   int iinc = (int)floor(dt/shd.delta+0.5);
   //shd.delta = 1./(int)floor(1./shd.delta+0.5);
/*
   if(fabs(iinc*shd.delta-dt)>1e-7) {
      cout<<"Error: "<<sps<<" is not a factor of "<<int(1/shd.delta)<<endl;
      exit(0);
   }
*/
   if(iinc!=1) {
      int npow = 1;
      double f1 = -1., f2 = -1., f3 = sps/2.2, f4 = f3*1.2;
      double dtd = (double)shd.delta;
      filter4_(&f1,&f2,&f3,&f4,&npow,&dtd,&(shd.npts),sig1,sig1);
   }

   int i, j;
   int nptst = (int)floor((shd.npts-1)*shd.delta*sps+0.5)+10;
   float nb;
   *sig2 = (float *) malloc (nptst * sizeof(float));
   long double fra1, fra2;
   nb = ceil((shd.nzmsec*0.001+shd.b)*sps);
   i = (int)floor((nb*dt-shd.nzmsec*0.001-shd.b)/shd.delta);
   if(fabs(iinc*shd.delta-dt)<1.e-7) { //sps is a factor of 1/delta
      fra2 = (nb*dt-i*shd.delta-shd.nzmsec*0.001-shd.b)/shd.delta;
      fra1 = 1.-fra2;
      if(fra2==0)
         for(j=0;i<shd.npts;j++) {
            (*sig2)[j] = sig1[i];
            i += iinc;
         }
      else
         for(j=0;i<shd.npts-1;j++) {
            (*sig2)[j] = sig1[i]*fra1 + sig1[i+1]*fra2;
            i += iinc;
         }
   }
   else { //sps isn't a factor, slower way
      fprintf(stderr, "*** Warning: sps isn't a factor of %d, watch out for rounding error! ***", (int)floor(1/shd.delta+0.5));
      long double ti, tj;
      iinc = (int)floor(dt/shd.delta);
      ti = i*shd.delta+shd.nzmsec*0.001+shd.b;
      tj = nb*dt;
      for(j=0;i<shd.npts-1;j++) {
         fra2 = tj-ti;
         (*sig2)[j] = sig1[i] + (sig1[i+1]-sig1[i])*fra2;
         tj += dt;
         i += iinc;
         ti += iinc*shd.delta;
         if( ti+shd.delta <= tj ) { ti += shd.delta; i++; }//if(j%1000==0)cout<<i<<" "<<ti<<j<<" "<<tj<<" "<<endl;}
      }
   }
   free(sig1);
   shd.nzmsec = (int)(nb*dt*1000+0.5);
   shd.b = 0.;
   if(shd.nzmsec>=1000) UpdateTime(&shd);
   shd.delta = dt;
   shd.npts = j;
   *sd = shd;
//   sprintf(sacname,"temp.sac");
//   write_sac(sacname,*sig2,&shd);
//exit(0);
   return 1;
}

char **Seed2Sac(int ne, int ns, char *nseed, int *nfile) {

   FILE *ff = fopen("from_seed","w");
   fprintf (ff, "%s <<END\n", rdsexe);
   fprintf (ff, "%s\n", nseed);
   fprintf(ff,"\n");                             /* out file */
   fprintf(ff,"\n");                             /* volume */
   fprintf(ff,"d\n");                            /* option */
   fprintf(ff,"\n");                             /* summary file */
   fprintf(ff,"%s\n", sdb->st[ns].name );        /* station list */
   fprintf(ff,"%s\n", ch );                      /* channel list */
   fprintf(ff,"\n");                             /* network list */
   fprintf(ff,"\n");                             /* Loc Ids */
   fprintf(ff,"1\n");                            /* out format */
   fprintf(ff,"N\n");                            // new version!!!!!!!!!!
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
   system("sh from_seed >& /dev/null");

   /*---------- response file -----------*/
   
   sprintf(str,"( ls RESP.*.%s.*.%s>list_resp ) >& /dev/null", sdb->st[ns].name, ch);
   system(str);

   char resp_name[100];
   ff = fopen("list_resp","r");
   if ( fscanf(ff,"%s", resp_name ) == EOF ) {
       fclose(ff);
       //sdb->rec[ne][ns].n = 0;
       return NULL;
   }
   fclose(ff);
   sprintf(sdb->rec[ne][ns].resp_fname,"%s/%s/%s", sdb->mo[imonth].name, sdb->ev[ne].name, resp_name);
   sprintf(str,"/bin/mv %s %s", resp_name, sdb->rec[ne][ns].resp_fname);
   system(str);
 
   system("/bin/rm -f list_resp");
   sprintf(str, "/bin/rm -f RESP.*.%s.*.%s", sdb->st[ns].name, ch);
   system(str);

   /*------------ sac list -------------*/
   sprintf(str,"( ls *%s*%s*SAC>list_sac ) >& /dev/null", sdb->st[ns].name, ch);
   system(str);
   if((ff=fopen("list_sac","r"))==NULL) {
     cout<<"Cannot open file list_sac"<<endl;
     *nfile = 0;
     //sdb->rec[ne][ns].n = 0;
     return NULL;
   }
   char **filelst=NULL;
   int i;
   for(i=0;;i++) {
      filelst = (char **) realloc (filelst, (i+1)*sizeof(char*));
      filelst[i] = (char *) malloc (100*sizeof(char));
      if(fscanf(ff,"%s",filelst[i])==EOF) break;
   }
   fclose(ff);
   *nfile = i;

   if(i==0) { 
      free(filelst[0]); free(filelst); 
      //sdb->rec[ne][ns].n = 0;
      return NULL;
   }

   return filelst;
}

float av_sig (float *sig, int i, int N, int nwin ) {
   int n1, n2, j, nav = 0;
   float av = 0.;

   if ( nwin > N ) nwin = N;
   n1 = i - nwin/2;
   if ( n1 < 0 ) n1 = 0;
   n2 = n1 + nwin - 1;
   if ( n2 > N-1 ) n2 = N-1;
   n1 = n2 - nwin + 1;

   for ( j = n1; j <= n2; j++ ) 
      if ( sig[j] < 1.e29 ) {
         av += sig[j];
         nav++;
      }

   if ( nav < 1 ) av = 1.e30;
   else av = av/(float)nav;

   return av;
}

int merge_sac(float * sig[], SAC_HD *sd, int nfile, int ne, int ns)
{
   if(nfile == 1) {
      sdb->rec[ne][ns].n = sd[0].npts;
      sdb->rec[ne][ns].t0 = abs_time (sd[0].nzyear, sd[0].nzjday, sd[0].nzhour, sd[0].nzmin, sd[0].nzsec, sd[0].nzmsec );
      sdb->rec[ne][ns].dt = sd[0].delta;
      write_sac (sdb->rec[ne][ns].fname, sig[0], &sd[0]);
      system("/bin/rm -f `more list_sac`");
      return 1;
   }

   int i, nb, j, jj, N, nfirst, Nholes;
   float *sig0;
   double t1[1000], t2[1000], T1 = 1.e25, T2 = -100.;
   SAC_HD s0;

   for(i=0; i<nfile; i++) {
      t1[i] = abs_time (sd[i].nzyear, sd[i].nzjday, sd[i].nzhour, sd[i].nzmin, sd[i].nzsec, sd[i].nzmsec );
      t2[i] = t1[i] + (sd[i].npts-1)*sd[i].delta;
      if ( t1[i] < T1 ) {
         T1 = t1[i];
         nfirst = i;
      }
      if ( t2[i] > T2 ) T2 = t2[i];
   }

   memcpy(&s0, &(sd[nfirst]), sizeof(SAC_HD) );
   double dt = (int)floor(s0.delta*1e8+0.5)/1e8, tshift;
   N = (int)floor((T2-T1)/s0.delta+0.5)+1;
   s0.npts = N;
   sdb->rec[ne][ns].n = N;
   sdb->rec[ne][ns].t0 = T1;
   sdb->rec[ne][ns].dt = dt;

   sig0 = (float *) malloc (N * sizeof(float));
   for (j=0;j<N;j++) sig0[j] = 1.e30;

   for(i=0;i<nfile;i++) {
      if( fabs(sd[i].delta-s0.delta)>.0001 ) {
          fprintf(stderr, "*** Warning: sps mismatch! ***");
          continue;
      }
      nb = (int)floor((t1[i]-T1)/dt+0.5);
      tshift = fabs((sd[i].b-s0.b) + (nb*dt-(t1[i]-T1)));
      if( tshift > 1.e-3 ) {
         fprintf(stderr, "*** Warning: signal shifted by %fsec when merging! ***", tshift);
      }
      for(j=0,jj=nb;j<sd[i].npts;j++,jj++) if(sig0[jj] > 1.e29) sig0[jj] = sig[i][j];
   }

   for(Nholes=0,j=0;j<N;j++) {
      if(sig0[j]>1.e29) Nholes++;
   }
   if( (float)Nholes/(float)N > gapfrac ) {
      system("/bin/rm -f `more list_sac`");
      sdb->rec[ne][ns].n = -1;
      return 0;
   }

   int rec_b[1000], rec_e[1000];
   char recname[200];
   rec_b[0] = 0;
   for(i=1,j=0;i<N;i++) {
      if(sig0[i-1]>1.e29) { if(sig0[i]<1.e29) rec_b[j] = i; }
      else if(sig0[i]>1.e29) rec_e[j++] = i;
   }
   if(sig0[N-1]<1.e29) rec_e[j++] = N;
   sprintf(recname, "%s_rec1", sdb->rec[ne][ns].ft_fname);
   FILE *frec = fopen(recname, "w");
   for(i=0;i<j;i++) fprintf(frec, "%d %d\n", rec_b[i], rec_e[i]);
   fclose(frec);

   float av;
   int npart;
   for ( j = 0; j < N; j++ ) if ( sig0[j] > 1.e29 ) {
      for(npart=16;npart!=1;npart/=2) {
         av = av_sig (sig0, j, N, N/npart );
         if ( av < 1.e29 ) break;
      }
      if(npart==1) av=0.;
      sig0[j] = av;
   }

   write_sac (sdb->rec[ne][ns].fname, sig0, &s0);
   system("/bin/rm -f `more list_sac`");

   free(sig0);

   return 1;
}

int MakeRecord (int ne, int ns) {
  //if ( sdb->rec[ne][ns].n > 0 ) return 0;
  int i, nfile;

  char **filelst = Seed2Sac(ne, ns, sdb->mo[imonth].seedf[ne], &nfile);
  if(filelst==NULL) return 0;

  float *sigrspd[nfile];
  SAC_HD sd[nfile];
  for(i=0;i<nfile;i++) {
     Resampling(filelst[i], &sigrspd[i], &sd[i]);
  }
  if ( !merge_sac(sigrspd, sd, nfile, ne, ns) ) return 0;

  for(i=0;i<nfile;i++) { free(filelst[i]); free(sigrspd[i]); }
  free(filelst[i]); free(filelst);
  return 1;
}

void ExtractSac() {
   int ist, iev, nst, flag;

   system("/bin/mkdir -p old_sac_files");
   system("/bin/mv *.SAC old_sac_files >& /dev/null");
   system("/bin/mv RESP.* old_sac_files >& /dev/null");
   
   for(iev=0;iev<NEVENTS;iev++) {
      for( ist = 0; ist < sdb->nst; ist++ ) sdb->rec[iev][ist].n = 0;
      if(strcmp(sdb->mo[imonth].seedf[iev],"0")==0) continue;
      nst = 0;
      fprintf(stderr, "### Extracting sac files from event %s: ", sdb->ev[iev].name);
      for( ist = 0; ist < sdb->nst; ist++ ) {
         SetName(iev, ist);
         //sdb->rec[iev][ist].n = 0;
         if( fskip1==2 || fskip1==1 ) {
            flag = CheckExistence(iev, ist);
            if(fskip1==2) continue;
            else if(flag) continue;
         }
         if( MakeRecord(iev, ist) ) {
	    if( nst%20 == 0) fprintf(stderr, "\n   ");
            fprintf(stderr, "%s ", sdb->st[ist].name);
            nst++;
         }
      }
      cout<<endl<<"   "<<nst<<" stations processed. ###"<<endl;
   }

   system("/bin/rm -f list_sac list_resp >& /dev/null");
}
