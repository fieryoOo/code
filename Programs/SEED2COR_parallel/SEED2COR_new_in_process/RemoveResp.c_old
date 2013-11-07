#include "Param.h"

//Transfer(): Remove instrument response
//CutRec(): cut the signal based on event origin time

void FDivide (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out, double *freq, double *amp, double *pha, int nf);

int CheckExistenceft(int ne, int ns) {
   SAC_HD *shd = read_shd(sdb->rec[ne][ns].ft_fname);
   if( shd==NULL ) return 0;
   sdb->rec[ne][ns].n = shd->npts;
   sdb->rec[ne][ns].dt = shd->delta;
   return 1;
}

void UpdateRecCut(char *name, int nstart, int rec_b, int rec_e) {
   int rec_b1[1000], rec_e1[1000], nrec1;
   FILE *frec;
   char recname1[200], recname[200];
   sprintf(recname1, "%s_rec1", name);
   sprintf(recname, "%s_rec", name);
   if( ! read_rec(1, recname1, 0, rec_b1, rec_e1, &nrec1) ) {
      frec = fopen(recname, "w");
      fprintf(frec, "%d	%d\n", rec_b, rec_e);
      fclose(frec);
      return;
   }
   int i;
   frec = fopen(recname, "w");
   for(i=0; i<nrec1; i++) {
      rec_b1[i] -= nstart;
      rec_e1[i] -= nstart;
      if(rec_e1[i]<=rec_b || rec_b1[i]>=rec_e) continue;
      if(rec_b1[i]<rec_b) rec_b1[i] = rec_b;
      if(rec_e1[i]>rec_e) rec_e1[i] = rec_e;
      fprintf(frec, "%d	%d\n", rec_b1[i], rec_e1[i]);
   }
   fclose(frec);
}

int CutRec(int ne, int ns, float *sig1, SAC_HD shd1) {
   int n, nstt, nend;
   float t1b, t1e, t2;

   float dt = sdb->rec[ne][ns].dt;
   n = (int)floor(tlen/dt+0.5)+1;
   t2 = t1 + (n-1)*dt;

//   if( read_sac(sdb->rec[ne][ns].ft_fname, sig1, shd1) == NULL ) {
//      //cerr<<"Cannot open file "<<sdb->rec[ne][ns].ft_fname<<endl;
//      sdb->rec[ne][ns].n = 0;
//      return 0;
//   }

   t1b = sdb->rec[ne][ns].t0 + shd1.b - sdb->ev[ne].t0;;
   t1e = t1b + (sdb->rec[ne][ns].n-1)*dt;

   nstt = (int)floor((t1-t1b)/dt+0.5);
   nend = nstt + n;
   t2 = fabs(nstt*dt - (t1-t1b));
   if( t2 > 1.e-3 )
      fprintf(stderr, "*** Warning: signal shifted by %fsec when cutting! ***", t2);

   int flag = 0;
   float *sig2=NULL;
   int rec_b, rec_e;
   rec_b = 0; rec_e = n;
   if( (nstt<0) || (nend>shd1.npts) ) {
      int i=0;
      if(nstt<0) i += -nstt;
      if(nend>shd1.npts) i += nend-shd1.npts;
      if((float)i/n > gapfrac) {
         sdb->rec[ne][ns].n = 0;
	 fRemove(sdb->rec[ne][ns].ft_fname);
         return 0;
      }
      flag = 1;
      fprintf(stderr, "*** Warning: cut range isn't fully covered. zeros padded ***");
      sig2 = (float *) malloc ( n * sizeof(float) );
      for(i=nstt;i<nend;i++) {
         if(i<0 || i>shd1.npts-1) sig2[i-nstt] = 0.;
         else sig2[i-nstt] = sig1[i];
      }
      if(nstt<0) rec_b = -nstt;
      if(nend > shd1.npts) rec_e = shd1.npts-1-nstt;
   }
   UpdateRecCut(sdb->rec[ne][ns].ft_fname, nstt, rec_b, rec_e);

   shd1.npts = n;
   sdb->rec[ne][ns].n = n;
   sdb->rec[ne][ns].dt = shd1.delta;
   //shd1.nzmsec += (int)floor((t1-t1b)*1000+0.5);
   shd1.b = 0.;
   shd1.nzhour = 0;
   shd1.nzmin = 0;
   shd1.nzsec = 0;
   shd1.nzmsec = (int)floor(t1*1000+0.5);
   UpdateTime(&shd1);
   if( flag ) {
      write_sac(sdb->rec[ne][ns].ft_fname, sig2, &shd1 );
      free(sig2);
   }
   else write_sac(sdb->rec[ne][ns].ft_fname, &(sig1[nstt]), &shd1 );
   free(sig1);

   return 1;
}
/*
int TransferSac(int ne, int ns) {

   if(ne>=NEVENTS || ns>=sdb->nst) {
      fprintf(stderr, "*** Warning: event/station # out of range! ***");
      //sdb->rec[ne][ns].n = 0;
      return 0;
   }

   FILE *ff;
   float fl2 = 1./perh*0.7, fl1 = fl2*0.8, fl3 = 1./perl*1.3, fl4 = fl3*1.2;
   //float fl1 = 1./100., fl2 = 1./80., fl3 = 1./4., fl4 = 1./3.;

   ff = fopen("sac_bp_respcor","w");
   fprintf(ff, "%s << END\n", sacexe);
   fprintf(ff, "r %s\n", sdb->rec[ne][ns].fname);
   fprintf(ff,"rmean\n");
   fprintf(ff,"rtrend\n");
   fprintf(ff,"transfer from  EVALRESP FNAME  %s to vel freqlimits %f %f %f %f\n", sdb->rec[ne][ns].resp_fname, fl1, fl2, fl3, fl4 );
   fprintf(ff,"w %s\n", sdb->rec[ne][ns].ft_fname);
   fprintf(ff,"quit\n");
   fprintf(ff,"END\n");
   fclose(ff);

   System("sh sac_bp_respcor >& /dev/null");
   if( fdel1 ) fRemove(sdb->rec[ne][ns].fname);
   return 1;
}
*/
int TransferEvr(int ne, int ns, float **sig, SAC_HD *sd) {
   // read in sac file
   *sig = NULL;
   if( (read_sac(sdb->rec[ne][ns].fname, sig, sd))==NULL ) {
      fprintf(stderr, "### Warning: cannot read sac file %s ###", sdb->rec[ne][ns].fname);
      return 0;
   }
   // running evalresp
   int nf = 100;
   char buff[300], sta[8], ch[8], net[8];
   float f2 = 1./perh*0.7, f1 = f2*0.8, f3 = 1./perl*1.3, f4 = f3*1.2;
   sscanf(sd->kstnm, "%s", sta);
   sscanf(sd->kcmpnm, "%s", ch);
   sscanf(sd->knetwk, "%s", net);
   sprintf(buff, "%s %s %s %4d %3d %f %f %d -f %s -v >& /dev/null", evrexe, sta, ch, sd->nzyear, sd->nzjday, f1, f4, nf, sdb->rec[ne][ns].resp_fname);
   system(buff);
   char nameam[50], nameph[50];
   sprintf(nameam, "AMP.%s.%s.*.%s", net, sta, ch);
   sprintf(nameph,"PHASE.%s.%s.*.%s", net, sta, ch);
   // find am file
   FILE *fam = NULL, *fph = NULL;
   int nlist;
   char *list = List(".", nameam, 0, &nlist);
   if( nlist!=1 ) {
      cerr<<"Error: "<<nlist<<" AMP file(s) found!"<<endl;
      exit(0);
   }
   sscanf(list, "%s", nameam);
   free(list);
   if( (fam = fopen(nameam, "r")) == NULL ) {
      cerr<<"Cannot open file "<<nameam<<endl;
      exit(0);
   }
   // find ph file
   list = List(".", nameph, 0, &nlist);
   if( nlist!=1 ) {
      cerr<<"Error: "<<nlist<<" PHASE file(s) found!"<<endl;
      exit(0);
   }
   sscanf(list, "%s", nameph);
   free(list);
   if( (fph = fopen(nameph, "r")) == NULL ) {
      cerr<<"Cannot open file "<<nameph<<endl;
      exit(0);
   }
   // read in am and ph data
   double pi=4*atan(1.0), pio180=pi/180.;
   double freq[nf], dtmp, amp[nf], pha[nf];
   int i = 0;
   while(i<nf) {
      if(fgets(buff, 300, fam)==NULL) break;
      sscanf(buff, "%lf %lf", &freq[i], &amp[i]);
      if(fgets(buff, 300, fph)==NULL) break;
      sscanf(buff, "%lf %lf", &dtmp, &pha[i]);
      if(dtmp!=freq[i]) {
	 cerr<<"### Warning: incompatible AMP - PHASE pair! ###"<<endl;
	 continue;
      }
      amp[i] *= 0.000000001;
      pha[i] *= pio180;
      i++;
   }
   fclose(fam); fclose(fph);
   fRemove(nameam); fRemove(nameph);
   // run rmresponse
   FDivide (f1, f2, f3, f4, (double)(sd->delta), sd->npts, *sig, *sig, freq, amp, pha, nf);
   if( fdel1 ) fRemove(sdb->rec[ne][ns].fname);
   return 1;
}

void RmRESP(){
   int ne, ns, nev, flag;
   float *sig=NULL; //pointer to signal produced by TransferEvr
   SAC_HD shd;

   for(ns=0;ns<sdb->nst;ns++) {
      nev = 0;
      fprintf(stderr, "### Removing response for the %dth station %s: ", ns+1, sdb->st[ns].name);
      for(ne=0;ne<NEVENTS;ne++) {
         if(strcmp(sdb->mo[imonth].seedf[ne],"0")==0) continue;
         if( fskip2==2 || fskip2==1 ) {
            flag = CheckExistenceft(ne, ns);
            if(fskip2==2) { 
               if( !flag ) sdb->rec[ne][ns].n = 0; 
               continue; 
            }
            else if(flag) continue;
         }
         if(sdb->rec[ne][ns].n > 0 ) if( !TransferEvr(ne, ns, &sig, &shd) ) continue; //sig allocated here
         if( !CutRec(ne, ns, sig, shd) ) continue; //sig freed here
	 if( nev %10 == 0 ) fprintf(stderr, "\n   ");
         fprintf(stderr, "%s ", sdb->ev[ne].name);
         nev++;
      }
      fprintf(stderr, "\n   %d events processed. ###\n", nev);
   }
   for(ns=0;ns<sdb->nst;ns++) for(ne=0;ne<NEVENTS;ne++) {
      if( (sdb->rec[ne][ns].resp_fname)==NULL ) continue; 
      if( fdel1 ) fRemove( sdb->rec[ne][ns].resp_fname );
      delete [] sdb->rec[ne][ns].resp_fname;
   }

   fRemove("sac_bp_respcor"); fRemove("RESP_tmp");
}
