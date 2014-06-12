#include "SacRec.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

inline int nint( float x ) { return (int)floor(x+0.5); }
float ComputeRMS_G( const SacRec& sac, const float xc, const float hlen ) {
   const int npts = sac.shd.npts;
   const float dt = sac.shd.delta;
   const int ic = nint( ( xc - sac.shd.b ) / dt );
   const int ihlen = nint( hlen / dt );

   //std::cerr<<npts<<" "<<dt<<"   "<<hlen<<" "<<ihlen<<" "<<ic<<std::endl;
   if( ihlen < 1 ) return -12345.;
   if( ic - ihlen < 0 || ic + ihlen > npts ) return -12345.;

   int ib = ic - 3*ihlen, ie = ic + 3*ihlen + 1;
   if( ib < 0 ) ib = 0;
   if( ie > npts ) ie = npts;
   float alpha = -0.5/(hlen*hlen);
   float rms = 0., V1 = 0., V2 = 0.;
   for(int i=ib; i<ie; i++) {
      float ftmp = xc - (i*dt+sac.shd.b);
      float weight = exp(alpha*ftmp*ftmp);
      ftmp = sac.sig[i];
      rms += ftmp * ftmp * weight;
      V1 += weight; V2 += weight*weight;
   }
   //return sqrt( rms / V1 );
   return sqrt( rms * V1 / (V1*V1-V2) );
}
float ComputeITG( const SacRec& sac, const float xc, const float hlen ) {
   const int npts = sac.shd.npts;
   const float dt = sac.shd.delta;
   const int ic = nint( ( xc - sac.shd.b ) / dt );
   const int ihlen = nint( hlen / dt );

   if( ihlen < 1 ) return -12345.;
   if( ic - ihlen < 0 || ic + ihlen > npts ) return -12345.;

   int ib = ic - 3*ihlen, ie = ic + 3*ihlen + 1;
   if( ib < 0 ) ib = 0;
   if( ie > npts ) ie = npts;
   float alpha = -0.5/(hlen*hlen);
   float itg = 0.;
   for(int i=ib; i<ie; i++) {
      float ftmp = xc - (i*dt+sac.shd.b);
      float weight = exp(alpha*ftmp*ftmp);
      ftmp = sac.sig[i];
      itg += ftmp * weight;
   }
   return itg;
}
float ComputeRMS_L( const SacRec& sac, const float xb, const float xe ) {
   const int npts = sac.shd.npts;
   const float dt = sac.shd.delta;
   const int ib = nint( ( xb - sac.shd.b ) / dt );
   const int ie = nint( ( xe - sac.shd.b ) / dt );

   if( ib < 0 || ie > npts ) return -12345.;
   float rms = 0.;
   for(int i=ib; i<ie; i++) {
      float ftmp = sac.sig[i];
      rms += ftmp * ftmp;
   }
   return sqrt( rms / (ie-ib-1) );
}

int main(int argc, char* argv[]) {
   // check sacf
   if(argc != 4) {
      std::cerr<<"Usage: "<<argv[0]<<" [SAC] [tb] [te]"<<std::endl;
      exit(-1);
   }
   SacRec sac_origin(argv[1]);
   if( ! sac_origin.Load() ) {
      std::cerr<<"Error loading sacfile "<<argv[1]<<std::endl;
      exit(0);
   }

   /* ---------- SNR from freq domain --------- */
   // cut out signal and noise
   float tb = atof(argv[2]), te = atof(argv[3]);
   SacRec sac_sig, sac_noise;
   sac_origin.cut(tb, te, sac_sig);
   if( te+1000. >= sac_origin.shd.npts ) {
      std::cerr<<"out of range!"<<std::endl;
      exit(0);
   }
   sac_origin.cut(te+500., te+1000., sac_noise);

   // compute amplitude spectrum
   SacRec amp_sig, amp_noise;
   sac_sig.ToAm(amp_sig);
   sac_noise.ToAm(amp_noise);
   // normalize by time series length
   amp_sig.Mul(1./200.);
   amp_noise.Mul(1./sac_noise.shd.npts);

   // compute rms average int each period bin
   float fhlen = 0.005;
   float perl = 5., perh = 50., perstep = 0.5;
   int nper = (int)ceil( (perh-perl) / perstep );
   std::vector<float> Per(nper), SNR_F(nper), AMP_F(nper), NOI_F(nper);
   for(int iper=0; iper<nper; iper++) {
      float per = perl + iper * perstep;
      float freq = 1. / per;
      AMP_F[iper] = ComputeRMS_G( amp_sig, freq, fhlen );
      NOI_F[iper] = ComputeRMS_G( amp_noise, freq, fhlen );
      //AMP_F[iper] = ComputeRMS_L( amp_sig, freq-fhlen, freq+fhlen );
      //NOI_F[iper] = ComputeRMS_L( amp_noise, freq-fhlen, freq+fhlen );
      //AMP_F[iper] = amp_sig.sig[nint( (freq-amp_sig.shd.b)/amp_sig.shd.delta )];
      //NOI_F[iper] = amp_noise.sig[nint( (freq-amp_sig.shd.b)/amp_sig.shd.delta )];
      if( AMP_F[iper]==-12345 || NOI_F[iper]==-12345 ) {
	 std::cerr<<"invalid rms at per="<<per<<"!"<<std::endl;
	 exit(0);
      }
      Per[iper] = per;
      SNR_F[iper] = AMP_F[iper] / NOI_F[iper];
   }


   /* ---------- SNR from time domain --------- */
   std::vector<float> SNR_T(nper), AMP_T(nper), NOI_T(nper);
   for(int iper=0; iper<nper; iper++) {
      float per = perl + iper * perstep;
      float freq = 1. / per;
      SacRec sac_ft;
      sac_origin.Filter(-1., freq, fhlen, -1., sac_ft);
      float ftmp, min, max;
      sac_ft.MinMax(tb, te, ftmp, min, ftmp, max);
      AMP_T[iper] = std::max(fabs(min), fabs(max));
      NOI_T[iper] = ComputeRMS_L( sac_ft, te+500., te+1000. );
      SNR_T[iper] = AMP_T[iper] / NOI_T[iper];
   }


   /* --------- output --------- */
   std::ofstream fout("SNR_freq_time.txt");
   for(int iper=0; iper<nper; iper++) {
      fout<<Per[iper]<<"   "<<AMP_F[iper]<<" "<<NOI_F[iper]<<" "<<SNR_F[iper]<<"   "<<AMP_T[iper]<<" "<<NOI_T[iper]<<" "<<SNR_T[iper]<<"\n";
   }
   fout.close();

   return 0;
}
