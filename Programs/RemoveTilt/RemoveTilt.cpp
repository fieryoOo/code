#include "SacRec.h"
#include <iostream>

class SacPair {
public:
	SacPair( const std::string fname1, const std::string fname2 )
		: sacr(fname1), sacs(fname2) {
		sacr.Load(); sacs.Load();
		if( sacr.shd.npts != sacs.shd.npts )
			throw ErrorSR::HeaderMismatch(FuncName, "npts");
		if( sacr.shd.delta != sacs.shd.delta )
			throw ErrorSR::HeaderMismatch(FuncName, "delta");

		CalcTransferF();
	}

private:
	SacRec sacr, sacs;

	void CalcTransferF() {
		// check sigs
		if( ! (sacr.sig&&sacs.sig) )
			throw ErrorSR::EmptySig(FuncName, "");
		// FFT
		SacRec sacr_am, sacr_ph; sacr.ToAmPh(sacr_am, sacr_ph);
		SacRec sacs_am, sacs_ph; sacs.ToAmPh(sacs_am, sacs_ph);
		auto sigram = sacr_am.sig.get(), sigsam = sacs_am.sig.get();
		auto sigrph = sacr_ph.sig.get(), sigsph = sacs_ph.sig.get();
		// starting point
		//const float deltaf = sacr_am.shd.delta;
		//int ib = (int)floor(f1/deltaf);
		//int hlen = 50;
		//if(ib<hlen+1) { ib=hlen; cout<<"Warning: Starting frequency too small!"<<endl; }
		// allocate sac mems for the coherence and transfer functions
		SacRec sacC_am; sacC_am.MutateAs(sacr_am); auto sigcam = sacC_am.sig.get();
		SacRec sacT_am; sacT_am.MutateAs(sacr_am); auto sigtam = sacT_am.sig.get();
		SacRec sacT_ph; sacT_ph.MutateAs(sacr_ph); auto sigtph = sacT_ph.sig.get();
		// loop over the spectrum
		float twopi = M_PI * 2.;
		for(int i=0; i<sacr_am.shd.npts; i++) {
			//float freq = sacr_am.X(i);
			//if( freq > f4 ) break;
			float ssam = sigsam[i] * sigsam[i];
			float rram = sigram[i] * sigram[i];
			float rsam = sigsam[i] * sigram[i];
			float rsph = sigsph[i] - sigrph[i];
			if(rsph < -M_PI) rsph += twopi;
			else if(rsph >= M_PI) rsph -= twopi;
			sigcam[i] = rsam / sqrt(ssam*rram);
			sigtam[i] = rsam / ssam;
			sigtph[i] = rsph;
		}
		// output
		sacC_am.Dump("debug1.txt");
		sacT_am.Dump("debug2.txt");
		sacT_ph.Dump("debug3.txt");
	}
};


int main (int argc, char *argv[])
{
   if( argc != 8) {
      std::cout<<"Usage: "<<argv[0]<<" [SAC_target] [SAC_noise] [f1] [f2] [f3] [f4] [out_name]"<<std::endl;
      exit(-1);
   }
   
	// check input freqs
   double f1 = atof(argv[3]);
   double f2 = atof(argv[4]);
   double f3 = atof(argv[5]);
   double f4 = atof(argv[6]);
   if( f1<0 || f1>=f2 || f2>=f3 || f3>=f4 ) {
      std::cout<<"Incorrect corner frequencies"<<std::endl;
      exit(0);
   }

	SacPair sp(argv[1], argv[2]);

/*
   int i, j, ii;
   int hlen = 50, ib = (int)floor(f1/domZ);
   if(ib<hlen+1) { ib=hlen; cout<<"Warning: Starting frequency too small!"<<endl; }
   double xx=0.,yy=0.,xy_r=0.,xy_i=0.,xy1,xy2;
   double temp_ph, temp_am, ftmp, fcurrent = ib*domZ, pi = 4.0*atan(1.0);
   double freq[nkZ], phres[nkZ], amres[nkZ];
   FILE *fout;
   char nametmp[100];
   sprintf(nametmp, "%s_TF", argv[7]);
   fout = fopen(nametmp, "w");
   for(i=ib-hlen-1; i<ib+hlen; i++) {
      xx += pow(sigDam[i],2.);
      yy += pow(sigZam[i],2.);
      xy1 = sigDam[i]*sigZam[i];
      xy_r += xy1*cos(sigDph[i]-sigZph[i]);
      xy_i += xy1*sin(sigDph[i]-sigZph[i]);
   }
   for(i=ib,ftmp=f4+domZ; i<nkZ-hlen&&fcurrent<=ftmp; i++, fcurrent+=domZ) {
      xx += pow(sigDam[i+hlen],2.) - pow(sigDam[i-hlen-1],2.);
      yy += pow(sigZam[i+hlen],2.) - pow(sigZam[i-hlen-1],2.);
      xy1 = sigDam[i+hlen]*sigZam[i+hlen]; xy2 = sigDam[i-hlen-1]*sigZam[i-hlen-1];
      xy_r += xy1 * cos(sigDph[i+hlen]-sigZph[i+hlen]) - xy2 * cos(sigDph[i-hlen-1]-sigZph[i-hlen-1]);
      xy_i += xy1 * sin(sigDph[i+hlen]-sigZph[i+hlen]) - xy2 * sin(sigDph[i-hlen-1]-sigZph[i-hlen-1]);
      temp_ph=atan(xy_i/xy_r);
      if(xy_r<0) temp_ph+=pi;
      if(temp_ph>pi) temp_ph-=2*pi;
      temp_am=sqrt(xy_r*xy_r+xy_i*xy_i);
      ii = i-ib;
      freq[ii] = fcurrent;
      phres[ii] = temp_ph;
      amres[ii] = temp_am/xx;
      fprintf(fout, "%g %g %g %g\n", fcurrent, temp_am/sqrt(xx*yy), temp_am/xx, temp_ph);
   }
   fclose(fout);

   dt = shdD.delta;
   rmresponse_(&f1,&f2,&f3,&f4,sigD,freq,phres,amres);
   //rmresponse_(&f1,&f2,&f3,&f4,&(shdD.npts),&dt,sigD,freq,phres,amres);
   //write_sac ("PRED_tmp.SAC", sigD, &shdD);
   for(i=0;i<shdZ.npts;i++) sigZ[i] -= sigD[i];
   write_sac (argv[7], sigZ, &shdZ );

   delete [] sigZ;
   delete [] sigD;
*/
   return 0;
}
