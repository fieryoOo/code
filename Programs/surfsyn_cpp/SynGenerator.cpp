#include "SynGenerator.h"
#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

/* FORTRAN entrance */

#ifdef __cplusplus
extern"C" {
#endif
	void read_params_( char name_feigen[255], char name_fparam[255], int* mode, char name_fmodel[256], float* elat, float* elon, float* aM,
							 float* dt, float tm[6], int* im, int* iq, char* its, int* nd, int* npoints, int* nt, float* tmin, float* tmax, float* vmax,
							 int* nsta, char sta[2000][8], float lat[2000], float latc[2000], float lon[2000], char net[2000][8] );

	void atracer_( char name_fmodel[256], float* elon, float* elat, int* ncor, bool* applyQ, int* nsta,
						float latc[2000], float lon[2000], float cor[2000][2][500] );

	void surfread_( char name_feigen[255], char* sigR, char* sigL, char modestr[2], int* nt, int* nd, float* depth, float freq[2000],
						float cr[2000], float ur[2000], float wvr[2000], float cl[2000], float ul[2000], float wvl[2000], float v[2000][3], 
						float dvdz[2000][3], float ampr[2000], float ampl[2000], float ratio[2000], float qR[2000], float qL[2000], float I0[2000] );

	void cal_synsac_( int* ista, char* its, char* sigR, char* sigL, float cor[2000][2][500], float* tmin, float* tmax, float* vmax, 
							int* n2pow, float* fix_vel, int* iq, int* npoints, float freq[2000], float* df, int* nf, float* dt, int* nt,
							float* elatc, float* elonc, bool* key_compr, float qR[2000], float qL[2000], int* im, float* aM, float tm[6],
							float ampl[2000], float cr[2000], float ul[2000], float ur[2000], float wvl[2000], float wvr[2000],
							float v[2000][3], float dvdz[2000][3], float ratio[2000], float I0[2000], float* slatc, float* slon,
							float* sigz, float* sign, float* sige );
#ifdef __cplusplus
}
#endif


SynGenerator::SynGenerator( const fstring& name_fparam, const fstring& name_feigen, int mode, float depth ) {
	// input params
	// wave type
	char sigR = '-', sigL = '-';
	std::string wavetype("RL");
	for(int i=0; i<wavetype.size(); i++)
		switch( wavetype[i] ) {
			case 'R': sigR = '+'; break;
			case 'L': sigL = '+'; break;
			default: std::cerr<<"Warning(??): unknown wave type "<<wavetype[i]<<std::endl;
		}
	if( sigR=='-' && sigL=='-' )
		throw std::runtime_error("Error(SynGenerator::SynGenerator): invalid wave type input!");
	// mode number
	char modestr[2]; sprintf(modestr, "%d", mode);
	if( modestr[1] == '\0' ) modestr[1] = ' ';
	std::cerr<<"!!"<<modestr<<"!!"<<std::endl;
	// depth
	int idepth = depth;
	char depstr[3]; sprintf(depstr, "%d", idepth);
	// fixed velocity
	bool key_compr = false;
	float fix_vel = 2.8;
	/*
	if(argc==7) {
		key_compr = true;
		fix_vel = atof(argv[6]);
	}
	*/

	std::cout<<"Wavetype = "<<wavetype<<" mode# = "<<mode<<"  depth = "<<idepth<<std::endl;

	// call read_params
	char fstr_fmodel[256], its;
	char sta[2000][8], net[2000][8];
	int im, iq, nd, npoints, nt, nsta;
	float elat, elon, aM, dt, tm[6], tmin, tmax, vmax;
	float lat[2000], latc[2000], lon[2000];
	read_params_( name_feigen.f_str(255), name_fparam.f_str(255), &mode, fstr_fmodel, &elat, &elon, &aM, &dt, 
					  tm, &im, &iq, &its, &nd, &npoints, &nt, &tmin, &tmax, &vmax,
					  &nsta, sta, lat, latc, lon, net );
	fstring name_fmodel; name_fmodel.assignf(fstr_fmodel, 256);

	// call atracer
	int ncor;
	float cor[2000][2][500];
	bool applyQ = true;
	atracer_( fstr_fmodel, &elat, &elon, &ncor, &applyQ, &nsta, latc, lon, cor );

	// call surfread
	float freq[2000], cr[2000], ur[2000], wvr[2000], ul[2000], cl[2000], wvl[2000];
	float ampr[2000], ampl[2000], ratio[2000], qR[2000], qL[2000], I0[2000];
	float v[2000][3], dvdz[2000][3];
	surfread_( name_feigen.f_str(255), &sigR, &sigL, modestr, &nt, &nd, &depth, freq, cr, ur, wvr,
				  cl, ul, wvl, v, dvdz, ampr, ampl, ratio, qR, qL, I0 );
	std::cerr<<nt<<" input periods\n";

	// calc base size
	int n2pow, nbase = 2;
	for(n2pow=1; n2pow<13; n2pow++) {
		if( nbase > npoints ) break;
		nbase <<= 1; 
	}
	std::cerr<<"n2pow = "<<n2pow<<"  nbase = "<<nbase<<"  npoints = "<<npoints<<"\n";
	int nf = nbase / 2;
	float df = 1. / (dt*nbase), f0 = 0.;

	// main loop on stations
	float drad = M_PI / 180.;
	float elatc = elat * drad, elonc = elon * drad;
	fstring fstr_sta, fstr_net;
	for(int ista=0; ista<nsta; ista++) {
		fstr_sta.assignf(sta[ista], 8); fstr_net.assignf(net[ista], 8);
		std::cerr<<ista<<" STA = "<<fstr_sta<<" NET = "<<fstr_net<<" STLAT = "<<lat[ista]<<" STLON = "<<lon[ista]<<"\n";
		// prepare sac header
		SacRec sacz, sacn, sace;
		WriteSACHeader( sacz.shd, npoints, dt, elat, elon, fstr_sta, lat[ista], lon[ista], "NaN", fstr_net );
		sacz.MutateAs(sacz); sacn.MutateAs(sacz); sace.MutateAs(sacz);
		cal_synsac_( &ista, &its, &sigR, &sigL, cor, &tmin, &tmax, &vmax, &n2pow, &fix_vel, &iq,
						 &npoints, freq, &df, &nf, &dt, &nt, &elatc, &elonc, &key_compr, qR, qL,
						 &im, &aM, tm, ampl, cr, ul, ur, wvl, wvr, v, dvdz, ratio, I0,
						 &(latc[ista]), &(lon[ista]), sacz.sig.get(), sacn.sig.get(), sace.sig.get() );
		std::string chnprefix = dt>=1 ? "LH" : "BH", outname, chn;
		chn = chnprefix + "Z";
		outname = "sac/"+fstr_sta+"."+chn+".SAC";
		strcpy(sacz.shd.kcmpnm, chn.c_str());
		sacz.Write(outname);
	}
}

/*
void SynGenerator::Run() {
	char fparam[255] = "RUNR", feigen[255] = "R.R", fout[255] = "bred", wtype[2] = "R", mode[2] = "0";
	std::fill(fparam+4, fparam+255, ' ');
	std::fill(feigen+3, feigen+255, ' ');
	std::fill(fout+4, fout+255, ' ');
	wtype[1] = ' '; mode[1] = ' ';
	float depth = 5.6;
	surfsyn_( fparam, feigen, fout, wtype, mode, &depth );
}
*/

bool SynGenerator::Synthetic( const float lon, const float lat, const std::string& chname, 
										const float f1, const float f2, const float f3, const float f4, SacRec& sac ) {
	std::string saclist("sacS.list");
	bool found = false;
	std::ifstream fin(saclist);
	if( ! fin )
		throw std::runtime_error("cannot access file "+saclist);
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line); ss >> line;
		sac.Load( line );
		if( chname==sac.chname() && lon==sac.shd.stlo && lat==sac.shd.stla ) {
			found = true;
			break;
		}
	}
	if( ! found ) {
		sac.clear();
	} else {
		sac.Mul(-1.);	// is the synthetic upside down???
		sac.Filter(f1,f2,f3,f4);
		sac.shd.b -= minfo.t0;
		sac.Resample();
	}
	return found;
}

void SynGenerator::WriteSACHeader( SAC_HD& shd, const int npts, const float dt, const float elat, const float elon,
														const std::string& sta, const float slat, const float slon, 
														const std::string& chn, const std::string& net ) {
	shd.npts = npts;
	shd.delta = dt;
	shd.b = 0.; 
	shd.e = shd.delta*(shd.npts-1);
	shd.evdp = 0.0;
	shd.evel = 0.0;
	shd.evla = elat;
	shd.evlo = elon;
	strcpy(shd.kstnm,sta.c_str());
	shd.stdp = 0.;
	shd.stel = 0.;
	shd.stla = slat;
	shd.stlo = slon;
	shd.nzyear = 2000;
	shd.nzjday = 1;
	shd.nzhour = shd.nzmin = shd.nzsec = shd.nzmsec = 0;
	strcpy(shd.kcmpnm, chn.c_str());
	strcpy(shd.knetwk, net.c_str());
}
