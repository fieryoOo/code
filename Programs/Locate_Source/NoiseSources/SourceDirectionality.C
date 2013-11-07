#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <float.h>
#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_num_threads() { return 1;}
#endif

#define BLKSIZE 200

struct DATA {
   float key;
   float azi, dis, val;
};

struct RESULTS {
   float amp, ampsigma;
   float alpdec, alpsigma;
};

char finname[300];
int ibindebug = -1;
struct DATA *data;
float hwbin, alpazi, hdis, alpdis, refdis;

void Merge(struct DATA *arr, int p, int q, int r) {
   int i, j, k, n1 = q-p+1, n2 = r-q;
   struct DATA L[n1], R[n2];
   for(i=0; i<n1; i++) L[i] = arr[p+i];
   for(j=0; j<n2; j++) R[j] = arr[q+j+1];
   for(i=0,j=0,k=p;;k++)
      if(L[i].key<=R[j].key) {
         arr[k] = L[i++];
         if(i==n1) break;
      }
      else {
         arr[k] = R[j++];
         if(j==n2) break;
      }
   if(i==n1) for(k++;k<=r;k++,j++) arr[k] = R[j];
   else for(k++;k<=r;k++,i++) arr[k] = L[i];
}
void MergeSort(struct DATA *arr, int p, int r) {
   if(p<r) {
     // if(r-p>20) {
         int q=(p+r)/2;
         MergeSort(arr, p, q);
         MergeSort(arr, q+1, r);
         Merge(arr,p,q,r);
     // }
     // else InserSort(arr,p,r);
   }
}

int ReadData(char *fname, float maxdis) {
   FILE *fin;
   if( (fin=(fopen(fname, "r"))) == NULL ) {
      std::cerr<<"Error(ReadData): Cannot access file "<<fname<<std::endl;
      return -1;
   }
   int i, nblk=0;
   char buff[300];
   data = NULL;
   for(i=0;fgets(buff, 300, fin);) {
      if(nblk*BLKSIZE <= i) data = (struct DATA *) realloc (data, (++nblk)*BLKSIZE * sizeof(struct DATA));
      sscanf(buff, "%f %f %f", &(data[i].azi), &(data[i].dis), &(data[i].val));
      if( data[i].val <= 0 || data[i].dis > maxdis ) continue;
      if( data[i].dis == 0. ) data[i].dis = 1.e-10;
      //if( ibindebug > 0 ) std::cerr<<"ReadData "<<i<<": azi="<<data[i].azi<<" dis="<<data[i].dis<<" val="<<data[i].val<<std::endl;
      i++;
   }
   fclose(fin);
   return i;
}

void least_square_fit(int type, float *datx, float *daty, float *weight, int ndat, float *aout, float *sigmaaout, float *bout, float *sigmabout ) {
   if( ndat < 2 ) {
      *aout = *bout = *sigmaaout = *sigmabout = -12345.;
      return;
   }

   int i;
   float a, b, w, x, y, W=0, WX=0, WY=0, WX2=0, WY2=0, WXY=0;
   for(i=0;i<ndat;i++) {
      //w = 1./( sigma[i] * sigma[i] );
      w = weight[i];
      x = datx[i]; y = daty[i];
      W += w;
      WX += w * x; WY += w * y;
      WX2 += w * x * x; WY2 += w * y * y;
      WXY += w * x * y; 
   }

   // determine a b accordingly
   if(type == 0) {
      w = 1./(W*WX2-WX*WX);
      a = (W*WXY-WX*WY) * w;
      b = (-WX*WXY+WX2*WY) * w;
   }
   else if(type == 1) {
      w = 1./(W*WXY-WX*WY);
      a=(W*WY2-WY*WY) * w;
      b=(WY*WXY-WY2*WX) * w;
   }
   else {
      std::cout<<"Line_fit: Wrong input for indep var, stopped!"<<std::endl;
      exit(0);
   }
   *aout = a; *bout = b;

   if( ndat == 2 ) {
      *sigmaaout = *sigmabout = DBL_MAX;
      return;
   }
   // compute uncertainty in each parameter
   float S2 = 0., k; 
   float dtmp;
   // define k according to ndat from t distribution
   // assuming a 95% conf is equivalent to 2 sigma conf
   if( ndat == 3 ) k = 12.706;
   else if ( ndat == 4 ) k = 4.303;
   else if ( ndat == 5 ) k = 3.182;
   else if ( ndat == 6 ) k = 2.776;
   else k = 1.960+13.8/pow((float)ndat, 1.6); // this could be made better
   k *= 0.5; // now this is 1 sigma
   // compute uncertainties
   for(i=0;i<ndat;i++) {
      dtmp = daty[i] - a * datx[i] - b;
      S2 += dtmp * dtmp;
   }
   S2 /= ndat-2.;
   if( type == 0 ) {
      *sigmaaout = k * sqrt(S2 * W * w);
      *sigmabout = k * sqrt(S2 * WX2 * w);
   }
   else {
      S2 /= a*a;
      *sigmaaout = k * sqrt(S2 * W * w);
      *sigmabout = k * sqrt(S2 * WY2 * w);
   }
}

int FitBinDecay(float cazi, struct DATA *bindata, int npt, struct RESULTS *results) {
   if( npt < 8 ) {
      results->amp = -1.;
      return 0;
   }
   // least square fit line
   float datx[npt], daty[npt], weight[npt], dazi;
   FILE *fout;
   if( ibindebug == 1000 ) {
      char outname[300];
      sprintf(outname, "%s_%.1f", finname, cazi);
      fout = fopen(outname, "w");
   }
   for(int i=0;i<npt;i++) {
      datx[i] = bindata[i].dis;
      daty[i] = log(bindata[i].val);
      dazi = bindata[i].azi - cazi;
      weight[i] = exp(-alpazi*dazi*dazi);
      //if( cazi == ibindebug*hwbin ) std::cerr<<"FitBinDecay data "<<i<<": dis="<<bindata[i].dis<<" val="<<bindata[i].val<<" weight="<<weight[i]<<" dazi="<<dazi<<std::endl;
      if( ibindebug == 1000 ) fprintf(fout, "%f %f %f %f\n", bindata[i].dis, bindata[i].val, weight[i], dazi);
   }
   if( ibindebug == 1000 ) fclose(fout);
   float a, sigmaa, b, sigmab;
   least_square_fit(0, datx, daty, weight, npt, &a, &sigmaa, &b, &sigmab);
   //if(cazi==ibindebug*hwbin) std::cerr<<"FitBinDecay lsf results: a="<<a<<" b="<<b<<std::endl;
   // parameters and uncertainties
   results->amp = exp(b);
   results->alpdec = -a;
   results->alpsigma = sigmaa;
   results->ampsigma = results->amp * sigmab;
   return 1;
}

int CorrectDecay(struct DATA *data, int npt, struct RESULTS *results) {
   if( npt < 8 ) {
      results->amp = -1.;
      return 0;
   }
   // least square fit line
   float datx[npt], daty[npt], weight[npt], dazi;
   FILE *fout;
   int i;
   for(i=0;i<npt;i++) {
      datx[i] = data[i].dis;
      daty[i] = log(data[i].val);
      weight[i] = exp(alpdis*datx[i]*datx[i]);
   }
   float a, sigmaa, b, sigmab;
   least_square_fit(0, datx, daty, weight, npt, &a, &sigmaa, &b, &sigmab);
   // correct amplitudes for distance ( to a distance of refdis )
   for(i=0;i<npt;i++) data[i].val *= exp( a * (refdis-data[i].dis) );
   // parameters and uncertainties
   results->amp = exp(b);
   results->alpdec = -a;
   results->alpsigma = sigmaa;
   results->ampsigma = results->amp * sigmab;
   return 1;
}

int AmpAvg(struct DATA *data, int npt, int decflag, struct RESULTS *results) {
   if( npt < 3 || (decflag!=0&&decflag!=2) ) {
      results->amp = -1.;
      return 0;
   }
   int i;
   float ftmp, weight[npt], avg = 0., stdom = 0., V1 = 0., V2 = 0.;
   // average
   for(i=0;i<npt;i++) {
      ftmp = exp(alpdis*data[i].dis*data[i].dis);
      weight[i] = ftmp;
      V1 += ftmp;
      V2 += ftmp * ftmp;
   }
   if ( decflag == 2 ) for(i=0;i<npt;i++) data[i].val *= sqrt(data[i].dis/hdis);
   for(i=0;i<npt;i++) avg += data[i].val * weight[i];
   avg /= V1;
   // std of the mean
   for(i=0;i<npt;i++) {
      ftmp = avg - data[i].val;
      stdom += ftmp * ftmp * weight[i];
   }
   stdom = sqrt(stdom * V1) / (V1*V1 - V2);
   // output results
   results->amp = avg;
   results->ampsigma = stdom;
   results->alpdec = -1.;
   results->alpsigma = -1.;
   return 1;
}
 
void BinAvg(struct DATA *data, int ndat, float hwbin, int nbin, struct RESULTS *results) {
   int ipt, npt, idat, ibin;
   struct DATA *bindata;
   float azibin, azil, azih;
   float weight, weit, V2, ftmp;
   #pragma omp parallel for private(azibin, azil, azih, ipt, npt, idat, bindata, weight, weit, V2, ftmp)
   for(ibin=0; ibin<nbin; ibin++) {
      azibin = ibin * hwbin;
      azil = azibin - 3*hwbin;
      azih = azibin + 3*hwbin;
      ipt = 0; 
      bindata = new struct DATA[ndat];
      if( azil < 0. ) {
	 azil += 360.;
	 for(idat=ndat-1; idat>0 && data[idat].azi>azil; idat--) {
	    bindata[ipt] = data[idat];
	    bindata[ipt++].azi -= 360.;
	    //if(ibin==ibindebug) std::cerr<<"BinAvg data "<<ipt-1<<": idat="<<idat<<" dis="<<data[idat].dis<<std::endl;
	 }
	 azil = 0.;
      }
      if( azih > 360. ) {
	 azih -= 360.;
	 for(idat=0; idat<ndat && data[idat].azi<azih; idat++) {
	    bindata[ipt] = data[idat];
	    bindata[ipt++].azi += 360.;
	    //if(ibin==ibindebug) std::cerr<<"BinAvg data "<<ipt-1<<": idat="<<idat<<" dis="<<data[idat].dis<<std::endl;
	 }
	 azih = 360.;
      }
      for(idat=0; idat<ndat; idat++) {
	 if( data[idat].azi <= azil ) continue;
	 if( data[idat].azi >= azih ) break;
	 bindata[ipt++] = data[idat];
	 //if(ibin==ibindebug) std::cerr<<"BinAvg data "<<ipt-1<<": idat="<<idat<<" dis="<<data[idat].dis<<std::endl;
      }
      npt = ipt;
      //for(idat=0; idat<npt; idat++) bindata[idat].key = bindata[idat].dis;
      //MergeSort(bindata, 0, npt-1);
      //FitBinDecay(azibin, bindata, npt, &(results[ibin]));
      if( npt < 2 ) {
	 results[ibin].amp = -1.;
	 delete [] bindata; bindata = NULL;
	 continue;
      }
      weit = 0.; results[ibin].amp = 0.;
      for(ipt=0; ipt<npt; ipt++) {
	 weight = bindata[ipt].azi - azibin;
	 weight = exp(alpazi * weight * weight) * exp(alpdis * bindata[ipt].dis * bindata[ipt].dis);
	 weit += weight;
	 results[ibin].amp += bindata[ipt].val * weight;
      }
      if( weit < 0.5 ) {
         results[ibin].amp = -1.;
         delete [] bindata; bindata = NULL;
         continue;
      }
      results[ibin].amp /= weit;
      V2 = 0.; results[ibin].ampsigma = 0.;
      for(ipt=0; ipt<npt; ipt++) {
         weight = bindata[ipt].azi - azibin;
         weight = exp(alpazi * weight * weight) * exp(alpdis * bindata[ipt].dis * bindata[ipt].dis);
         V2 += weight * weight;
         ftmp = bindata[ipt].val - results[ibin].amp;
         results[ibin].ampsigma += ftmp*ftmp * weight;
      }
      results[ibin].ampsigma = sqrt(results[ibin].ampsigma*weit)/(weit*weit-V2);
      delete [] bindata; bindata = NULL;
   }
}

void PrintResults(struct RESULTS *results) {
   std::cout<<"RESULTS: amp = "<<results->amp<<" (+/-"<<results->ampsigma<<")  alpha = "<<results->alpdec<<" (+/-"<<results->alpsigma<<")"<<std::endl;
}

int main(int argc, char *argv[])
{
   if( argc != 4 && argc != 5 ){
      std::cerr<<"Usage: "<<argv[0]<<" [input file (azimuth distance value)] [hdist] [correct decay? (0=no, 1=fit, 2=geometric)] [optional refdis (normalize SNR to the refdis. refdis=hdist if not provided)]"<<std::endl;
      exit(-1);
   }

 // read in and sort by azimuth
   hdis = atof(argv[2]); alpdis = -0.5/(hdis*hdis);
   if( argc == 4 ) refdis = hdis;
   else refdis = atof(argv[4]);
   //std::cout<<"Normalizing SNR to an refdis of "<<refdis<<"km"<<std::endl;
   sprintf(finname, "%s", argv[1]);
   int i, ndat = ReadData(argv[1], hdis*3.);
   if( ndat <= 0 ) exit(-1);
   for(i=0; i<ndat; i++) data[i].key = data[i].azi;
   MergeSort(data, 0, ndat-1);

 // fit decay, correct amplitude for distance, and return overall amp and alp
   int decflag = atoi(argv[3]), ierr; 
   struct RESULTS results;
   if( decflag == 1 ) ierr = CorrectDecay(data, ndat, &(results));
   else if( decflag == 0 || decflag == 2 ) ierr = AmpAvg(data, ndat, decflag, &(results));
   else {
      std::cerr<<"Unknow decay flag: "<<decflag<<std::endl;
      exit(0);
   }
   if( ierr == 0 ) exit(0);
   //std::cout<<"decay "; PrintResults(&results);
   std::cout<<results.amp<<" "<<results.ampsigma<<" "<<results.alpdec<<" "<<results.alpsigma<<std::endl;

 // bin and compute alpha and alpha-sigma in each bin
   hwbin = 10.; alpazi = -0.5/(hwbin*hwbin);
   int nbin = (int)ceil(360./hwbin);
   struct RESULTS rs, binresults[nbin];
   memset(binresults, 0, nbin*sizeof(struct RESULTS));
   BinAvg(data, ndat, hwbin, nbin, binresults);

 // output azimuth-SNR-sigmaSNR
   char outname[300];
   sprintf(outname, "%s_binavg", finname);
   FILE *fout = fopen(outname, "w");
   for(i=0;i<nbin;i++) { 
      rs = binresults[i];
      if( rs.amp <= 0 ) continue;
      fprintf(fout, "%f %f %f\n", i*hwbin, rs.amp, rs.ampsigma);
   }
   fclose(fout);

   free(data);

   return 0;
}
