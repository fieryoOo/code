#include "Param.h"

int Whiten( double f1, double f2, double f3, double f4, double dt, int n, float hlen, float *seis_in, float *seissm, float **outam, float **outph, int *nk, double *dom);

void UpdateRec(char *name, int *rec_b, int *rec_e, int nrec, int ithread) {
   int rec_b1[1000], rec_e1[1000], nrec1;
   FILE *frec;
   int irec, irec1;
   if( ! read_rec(1, name, 0, rec_b1, rec_e1, &nrec1) ) {
      reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: cannot open record file %s ***", name);
      frec = fopen(name, "w");
      for(irec=0; irec<nrec; irec++)
         fprintf(frec, "%d %d\n", rec_b[irec], rec_e[irec]);
      fclose(frec);
   }
   int recB, recE;
   char name2[100];
   sprintf(name2, "%s2", name);
   frec = fopen(name2, "w");
   for(irec=0; irec<nrec; irec++)
      for(irec1=0;irec1<nrec1;irec1++){
         if(rec_b[irec]>=rec_e1[irec1]) continue;
         if(rec_e[irec]<=rec_b1[irec1]) break;
         recB = max(rec_b[irec],rec_b1[irec1]);
         recE = min(rec_e[irec],rec_e1[irec1]);
         fprintf(frec, "%d %d\n", recB, recE);
      }
   fclose(frec);
}

void OneBit(float *sig, SAC_HD *shd) {
   int i;
   for(i=0;i<shd->npts;i++) {
      if(sig[i]>0) sig[i] = 1.;
      else sig[i] = -1.;
   }
}

void RunAvg(float *sig, SAC_HD *shd) {

   int n = shd->npts;
   float *sigw, dt = shd->delta;
   sigw = new float[n];
   double f2 = 1./Eperh, f1 = f2*0.6, f3 = 1./Eperl, f4 = f3*1.4;
//   double f1 = 1./60., f2 = 1./50., f3 = 1./10., f4 = 1./7.5;
   if( Eperl == -1 ) memcpy(sigw, sig, n*sizeof(float));
   else Filter(f1, f2, f3, f4, (double)dt, n, sig, sigw);

   int i, j, wb, we;
   int half_l = (int)floor(timehlen/dt+0.5);
   float wsum;

   // fast running average before whitening

   if(half_l*2>n-1) half_l = (n-1)/2;
   for(i=0,wsum=0.;i<=half_l;i++) wsum += fabs(sigw[i]);
   wb = 0; we = i;
   for(i=1;i<=half_l;i++,we++) {
      if(wsum!=0.) sig[i-1] /= (wsum/we);
      wsum += fabs(sigw[we]);
   }
   for(j=we;i<n-half_l;i++,wb++,we++) {
      if(wsum!=0.) sig[i-1] /= (wsum/j);
      wsum += fabs(sigw[we]) - fabs(sigw[wb]);
   }
   for(;i<n;i++,wb++) {
      if(wsum!=0.) sig[i-1] /= (wsum/(we-wb));
      wsum -= fabs(sigw[wb]);
   }
   if(wsum!=0.) sig[n-1] /= (wsum/(we-wb));
   delete [] sigw; sigw = NULL;
   //write_sac("temp1.sac", sig, shd );
   //exit(0);
}

int EqkCut(float *sig, SAC_HD *shd, char *recname, int ithread) {

   int i, n = shd->npts, ninc = n/1000, npole=0;
   for(i=0;i<n;i+=ninc) if(fabs(sig[i])<1.e-20) npole++;
   if(npole>600) {
      reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Signal time length not long enough. ***");
      return 0;
   }

   float sigw[n];
   double f2 = 1./Eperh, f1 = f2*0.8, f3 = 1./Eperl, f4 = f3*1.2;
   double dt = (double)(shd->delta);
   if( Eperl == -1 ) memcpy(sigw, sig, shd->npts*sizeof(float));
   else Filter( f1, f2, f3, f4, dt, n, sig, sigw);

   //compute noise level
   int   ii, j, is, s1k, flag;
   s1k=(int)floor(1000./dt+0.5);
   int rec_i,rec_b[1000],rec_e[1000];

   double win_max[s1k+1],win_min,window_avg,window_std;
   for(i=0;i<n;i++) sigw[i] = fabs(sigw[i]);
   for( is=0, i=0; i<= n-s1k; is++,i+=s1k){
      win_max[is]=0;
      for( ii=i; ii<i+s1k; ii++ )
         if(win_max[is]<sigw[ii]) win_max[is]=sigw[ii];
   }
   flag=0;
   if(i<n) flag=1;
   for( ii=i;ii<n;ii++ )
      if(win_max[is]<sigw[ii]) win_max[is]=sigw[ii];

   window_avg=0;ii=0;
   win_min = 1e20;
   for(i=0; i<(int)(n/s1k); i++) if(win_max[i]>1e-20 && win_min>win_max[i]) win_min=win_max[i];

   for( i =0; i< (int)(n/s1k); i++){
      if(win_max[i]>win_min*2.0 || win_max[i]<1e-20) continue;
      window_avg+=win_max[i];
      ii+=1;
   }
   if( ii < 20 ) {
      reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Time length not enough after removing earthquakes. ***");
      return 0;
   }
   window_avg=window_avg/ii;
   window_std=0;
   for( i =0; i< (int)(n/s1k); i++){
      if(win_max[i]>win_min*2.0 || win_max[i]<1e-20)continue;
      window_std+=(window_avg-win_max[i])*(window_avg-win_max[i]);
   }
   window_std=sqrt(window_std/(ii-1));

   ii=0;
   for( i =0; i < (int)(n/s1k); i++){
      if(win_max[i]>window_avg+2.0*window_std){  // || win_max[i]<1e-20){
         for ( j=0; j<s1k; j++ ) sig[i*s1k+j]=0;
         sigw[i]=0;
      }
      else sigw[i]=1;
      ii+=1;
   }
   if(win_max[i]>window_avg+2.0*window_std) {
      for ( j=i*s1k; j<n; j++ ) sig[j]=0;
      sigw[i]=0;
   }
   else sigw[i]=1;

   if( ii < 20 ) {
      reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Time length not enough after removing earthquakes. Skipped. ***");
      return 0;
   }

   rec_i=0;
   rec_b[0]=0;
   for( i=1; i<(int)(n/s1k)+flag;){
      if(sigw[i]-sigw[i-1]==1) rec_b[rec_i]=i*s1k;
      else if(sigw[i]-sigw[i-1]==-1) {
         rec_e[rec_i]=i*s1k;
         if ((rec_e[rec_i]-rec_b[rec_i])<2500/dt)
            for(ii=rec_b[rec_i];ii<rec_e[rec_i];ii++) sig[ii]=0;
         else rec_i++;
      }
      i++;
   }
   if(sigw[i-1]==1) {
      rec_e[rec_i]=n;
      if((rec_e[rec_i]-rec_b[rec_i])<1500/dt)
         for(ii=rec_b[rec_i];ii<rec_e[rec_i];ii++) sig[ii]=0;
      else rec_i++;
   }
   for(i=0;i<rec_i;i++) {
      if(rec_b[i]!=0){
         rec_b[i]+=300;
         for(ii=rec_b[i]-300;ii<rec_b[i];ii++) sig[ii]=0;
      }
      if(rec_e[i]!=n){
         rec_e[i]-=300;
         for(ii=rec_e[i]+1;ii<=rec_e[i]+300;ii++) sig[ii]=0;
      }
   }

   UpdateRec(recname, rec_b, rec_e, rec_i, ithread);

   return 1;
}

int TemperalNorm( char *fname, float **sig, SAC_HD *shd, int ithread) {

   *sig = NULL;
   if( read_sac(fname, sig, shd) == NULL ) {
      reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Cannot open file %s ***", fname);
      return 0;
   }
   char recname[300];
   sprintf(recname,"%s_rec",fname);

   if( tnorm_flag==1 ) OneBit(*sig, shd);
   else if( tnorm_flag==2 ) RunAvg(*sig, shd);
   else if( tnorm_flag==3 ) { if(!EqkCut(*sig, shd, recname, ithread)) return 0; }
   else if( tnorm_flag!=0 ) {  
      cerr<<"ERROR(TemperalNorm): Undefined normalization method!"<<endl; 
      exit(0); 
   }

/*
   char ftmp[100], *fp;
   sprintf(ftmp, "%s", fname);
   strtok(ftmp, "_");
   fp = strtok(NULL, "\n");
   sprintf(ftmp, "temp/ft_%s", fp);
   write_sac(ftmp, *sig, shd );
*/
// exit(0);

   return 1;
}

int SpectralNorm(char *fname, float *sig, SAC_HD shd, int ithread) {

   int nk, flag_whiten;
   double dom;
   double dt = (double)shd.delta;
   int n = shd.npts;
   //double f1 = 1./100., f2 = 1./80., f3 = 1./4., f4 = 1./3.;
   double f2 = 1./perh, f1 = f2*0.8, f3 = 1./perl, f4 = f3*1.2;
   float *sigw = NULL;
   //if(strcmp(fwname,"0") != 0) {
   if(frechlen==-1.) {
      SAC_HD shdw;
      if( read_sac(fwname, &sigw, &shdw) == NULL ) {
         cerr<<"ERROR(SpectralNorm): Cannot open file "<<fwname<<endl;
         free(sig);
         exit(0);
      }
      if(shdw.npts!=shd.npts || fabs(shdw.delta-shd.delta)>1.e-3) {
         cerr<<"ERROR(SpectralNorm): Smoothing spectrum is incompatibale with input signals!"<<endl;
         exit(0);
      }
   }
   //float seis_out[n];
   //int nf = (int)(log((double)n)/log(2.))+1;
   //if(nf<13) nf = 13;
   //nf = (int)pow(2,nf);
   float *outam=NULL, *outph=NULL; //[nf];
   //memset (seis_outamp,0,nf*sizeof(float));
   //memset (seis_outph,0,nf*sizeof(float));
   //pthread_mutex_lock(&fftlock);
   //whiten_(&f1,&f2,&f3,&f4,&npow,&dt,&n,&frechlen,sig,sigw,seis_out,outam,outph,&ns,&dom,&flag_whiten);
   flag_whiten = Whiten( f1, f2, f3, f4, dt, n, frechlen, sig, sigw, &outam, &outph, &nk, &dom);
   //pthread_mutex_unlock(&fftlock);
   if(flag_whiten==0) reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Skipped due to probamatic spectrum. ***");
   else {
      char nameamp[200], nameph[200];
      sprintf(nameamp, "%s.am", fname);
      sprintf(nameph, "%s.ph", fname);
      shd.npts = nk;
      shd.delta = dom;
      shd.b = 0;
      shd.iftype = IXY;
      write_sac(nameamp,outam, &shd );
      write_sac(nameph, outph,  &shd );
   }

   delete [] outam; outam = NULL;
   delete [] outph; outph = NULL;
   free(sig); free(sigw);
  
   return flag_whiten;
}

void * NormEvent (void *tid) {
   int ithread = *((int *)tid);
   int iev, ist, nst;
   char amname[200], phname[200];
   float *sig;
   SAC_HD shd;

   for(;;) {
      //get/update current event number
      pthread_mutex_lock(&cevlock);
      iev = currevn;
      currevn++;
      pthread_mutex_unlock(&cevlock);
      if(iev>=NEVENTS) break;
      //Normalizations
      if(strcmp(sdb->mo[imonth].seedf[iev],"0")==0) continue;
      nst = 0;
      reports[ithread].tail += sprintf(reports[ithread].tail, "### Am & ph records produced for event %s from thread %d: ", sdb->ev[iev].name, ithread);
      for( ist = 0; ist < sdb->nst; ist++ ) {
         sprintf(amname, "%s.am", sdb->rec[iev][ist].ft_fname);
         sprintf(phname, "%s.ph", sdb->rec[iev][ist].ft_fname);
         if( fskip3==2 || (fskip3==1 && access( amname, R_OK) != -1 && access( phname, R_OK) != -1) ) continue;
         if(sdb->rec[iev][ist].n <= 0) continue;
         if( !TemperalNorm( sdb->rec[iev][ist].ft_fname, &sig, &shd, ithread) )
            { sdb->rec[iev][ist].n = 0; continue; }
         if( !SpectralNorm( sdb->rec[iev][ist].ft_fname, sig, shd, ithread ) )
            { sdb->rec[iev][ist].n = 0; continue; }
         if( nst%20 == 0 ) reports[ithread].tail += sprintf(reports[ithread].tail, "\n   ");
         reports[ithread].tail += sprintf(reports[ithread].tail, "%s ", sdb->st[ist].name);
         nst++;
      }
      reports[ithread].tail += sprintf(reports[ithread].tail, "\n   %d stations processed. ###\n", nst);
      cout<<reports[ithread].head;
      reports[ithread].tail = reports[ithread].head;
   }

   pthread_exit(NULL);
}

void TempSpecNorm () {
   int ithread;
   //initialize report arrays
   reports = (struct NOTE *) malloc ( NTHRDS * sizeof(struct NOTE));
   for(ithread=0;ithread<NTHRDS;ithread++) {
      reports[ithread].head = (char *) malloc ( (sdb->nst+1) * 100 * sizeof(char) );
      reports[ithread].tail = reports[ithread].head;
   }

   currevn = 0;
   //Thread ids and lock
   pthread_t tid[NTHRDS];

   //Create threads to produce normalized am/ph files
   int rc, targs[NTHRDS];
   for(ithread=0;ithread<NTHRDS;ithread++) {
      targs[ithread] = ithread;
      rc = pthread_create( &tid[ithread], &attr_j, NormEvent, (void *) (&(targs[ithread])) );
      if(rc) {
         cerr<<" Thread creation failed!  ERR: "<<strerror(rc)<<endl;
         exit(0);
      }
   }
   //Wait for all threads to finish
   for(ithread=0;ithread<NTHRDS;ithread++) pthread_join(tid[ithread], NULL);

   //pthread_mutex_destroy(&rdslock);

   //free report arrays
   for(ithread=0;ithread<NTHRDS;ithread++) free(reports[ithread].head);
   free(reports);

}

