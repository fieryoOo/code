#ifndef STASACS_H
#define STASACS_H

#include "SacRec.h"
#include "RDirect.h"
#include <iostream>
#include <map>
#include <tuple>
//#define DEBUG

class PointC5 : public PointC {
public:
	float z2=NaN, z3=NaN;
	PointC5(const PointC& pc)
		: PointC(pc) {}
	friend std::ostream& operator<< ( std::ostream& o, const PointC5& pt ) {
		o << pt.x << " " << pt.y << " " << pt.z << " " << pt.z2 << " " << pt.z3;
		return o;
	}
};

class StaSacs {
public:
	// orientation of horizontals are assumed to be azi_H1 (for H1) and azi_H1+90 (for H2)
	StaSacs( const std::string& fnameZ, const std::string& fnameH1, 
				const std::string& fnameH2="", const std::string& fnameD="",
				const int sactypein = 0, const float water_depth_in = -1, const float azi_H1 = 0.,
				const float Eperl = 10., const float Eperu = 40. );

	// theoretical cut-off at freq = sqrt(g/(2 pi d)) where water_depth = ingragravity_wavelength
	float fcutoffCompliance() const {
		if( water_depth < 0. ) throw ErrorSR::BadParam(FuncName, "invalid (non-positive) water depth");
		//return 1.2489 / sqrt(water_depth);
		return 1.45 / sqrt(water_depth);	// use a larger cutoff f to ensure a clean cut
	}

	float RemoveCompliance( const std::string& outinfoname = "", float tseg=2000. );

	PointC RemoveTilt( const std::string& outinfoname = "", float tseg=2000. );

	void test(float tseg, int nsm, float tb = NaN, float te = NaN);
	PointC5 RemoveTiltCompliance( const std::string& outinfoname = "", float tseg=2000. );

	RDirect RayleighDirectionality(const float dazi, const std::vector<std::pair<float, float>>& freqRangeV, 
											 float tseg=2000., float twin=1800.);

	PointC EstimateTiltDirection( const float ddeg ) const;

	// compute transfer function for each window (twin) from segments (tseg) within that window
	void ReshapeByZHPhaseDiff( float twin, float tseg );

	void Write( const std::string& foutZ, const std::string& foutZinterm="",
					const std::string foutH1="", const std::string& foutH2="", const std::string foutD="" ) const;

	void clearV() const {
		sacZamV.clear(); sacDamV.clear(); sacH1amV.clear(); sacH2amV.clear(); sacHtamV.clear(); 
		sacZphV.clear(); sacDphV.clear(); sacH1phV.clear(); sacH2phV.clear(); sacHtphV.clear();
	}

	void clear() {
		//sactype = 0; water_depth = -1.;
		trec.clear(); sacZ.clear(); 
		sacD.clear(); sacH1.clear(); sacH2.clear(); 
		sacHt.clear(); sacZinterm.clear();
		clearV();
	}

protected:
	static constexpr float NaN = SacRec::NaN;
	static constexpr float PIo2 = M_PI * 0.5; 
	static const int nfcorr_min = 10;
	static constexpr float fb = 0.01;		// signal below fb is discarded
	static constexpr float cohmin = 0.5;
	static constexpr float fmax_tilt = 0.105;
	static constexpr float fb_avg = 0.01;	// cohavg is computed between fb_avg and fe_avg Hz
	static constexpr float fe_avg = 0.09;	
	//static constexpr float cohmin = 0.2;
	//static constexpr float fmax_tilt = 0.4;
	//static constexpr float fb_avg = 0.2;	// cohavg is computed between fb_avg and fe_avg Hz
	//static constexpr float fe_avg = 0.4;	

private:
	const int sactype;
	const float water_depth, fmax_comp;
	const float _Eperl, _Eperu;
	SacRec sacZinterm;	// stores intermediate Z after removing the stronger of compliance/tilt
	SacRec sacZ, sacH1, sacH2, sacD, sacHt;
	mutable float _tseg = -1.;
	mutable std::vector<std::pair<float, float>> trec; 
	mutable std::vector<SacRec> sacZamV, sacH1amV, sacH2amV, sacDamV, sacHtamV;
	mutable std::vector<SacRec> sacZphV, sacH1phV, sacH2phV, sacDphV, sacHtphV;
	// use a map of pointers to make things cleaner
	const std::map<char, std::tuple<SacRec*,std::vector<SacRec>*,std::vector<SacRec>*>> sacsM{
		{'Z', std::make_tuple(&sacZ, &sacZamV, &sacZphV)},
		{'D', std::make_tuple(&sacD, &sacDamV, &sacDphV)},
		{'1', std::make_tuple(&sacH1, &sacH1amV, &sacH1phV)},
		{'2', std::make_tuple(&sacH2, &sacH2amV, &sacH2phV)},
		{'t', std::make_tuple(&sacHt, &sacHtamV, &sacHtphV)},
	};

	int nint( float val ) const { return (int)floor(val+0.5); }

	// compute valid time segments (removing earthquakes)
	SacRec DetectNoiseWindows( float Eperl, float Eperu, bool ataper = false ) const;

	// segmentize a single channel
	void Segmentize(const SacRec& sac, std::vector<SacRec>& sacamV, std::vector<SacRec>& sacphV) const;
	std::pair<int, int> GetSegRange(float &tb, float &te, std::vector<SacRec>& sacV) const;

	void FitIntoParabola( SacRec& sac1, const SacRec& sacsigma, const float fmin, const float fmax, const int nparab = 1 ) const;

	void CalcTransferF( const SacRec& sac1, std::vector<SacRec>& sac1amV, std::vector<SacRec>& sac1phV,
							  const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
							  SacRec& Coh, SacRec& Adm, SacRec& Pha, const float tb=NaN, const float te=NaN ) const;
	void CalcTransferF( const char c1, const char c2, SacRec& Coh, SacRec& Adm, SacRec& Pha, 
							  float tb=NaN, float te=NaN ) const;
	void CalcTransferF2( const SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha ) const;
	float CohAvg(const SacRec& Coh, const SacRec& Pha, float fb, float fe, float pha0) const;

	float EstimateCoherence( const char c1, const char c2, const float fmax, 
									 SacRec& Coh, SacRec& Adm, SacRec& Pha, bool autoflip = false );
	PointC EstimateTilt( SacRec& Coh, SacRec& Adm, SacRec& Pha );

	void ApplyCorrection( const char c, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha, 
								 float fmax, const float pha0 = 0. );

};


/* --------------------- Public Methods --------------------- */

StaSacs::StaSacs( const std::string& fnameZ, const std::string& fnameH1, 
						const std::string& fnameH2, const std::string& fnameD,
						const int sactypein, const float water_depth_in, const float azi_H1,
						const float Eperl, const float Eperu )
						: sacZ(fnameZ), sacH1(fnameH1), sacH2(fnameH2), sacD(fnameD), sactype(sactypein)
						, water_depth(water_depth_in), fmax_comp(fcutoffCompliance()), _Eperl(Eperl), _Eperu(Eperu) {
	// lambda function to load sac and convert to displacement
	auto LoadSAC = [&]( SacRec& sac ) {
		sac.Load(); //sac.cut(60000., 85000.);
		switch(sactype) {
			case 2:	// acc: integrate twice
				sac.Integrate();
			case 1:	// vel: integrate once
				sac.Integrate();
			case 0:	// dis: do nothing
				break;
			default:
				throw ErrorSR::BadParam( FuncName, "unknown sactype("+std::to_string(sactype)+")." );
		}
		#ifdef DEBUG
		auto chn = sac.chname();
		sac.Write("debug_"+chn+".SAC"); 
		#endif
	};
	LoadSAC(sacZ);
	try {	// load horizontal channels
		LoadSAC(sacH1); LoadSAC(sacH2); SACRotate(sacH1, sacH2, -azi_H1);
	} catch( ErrorSR::Base& e ) {
		std::cerr<<"Warning(StaSacs): empty/non-accessable sacH1/H2 ("<<e.what()<<")"<<std::endl;
	}
	try { 
		// load pressure gauge data and convert to MPa
		sacD.Load(); sacD.Mul(1.0e-6); sacD.ZoomToEvent(sacZ.shd);
		//sacD.cut(60000., 85000.);
	} catch( ErrorSR::Base& e ) {
		std::cerr<<"Warning(StaSacs): empty/non-accessable sacD ("<<e.what()<<")"<<std::endl;
	}
	if( sacH1.sig && (sacZ.shd.npts!=sacH1.shd.npts || sacZ.shd.delta!=sacH1.shd.delta ) )
		throw ErrorSR::HeaderMismatch(FuncName, sacH1.fname);
	if( sacH2.sig && (sacZ.shd.npts!=sacH2.shd.npts || sacZ.shd.delta!=sacH2.shd.delta ) )
		throw ErrorSR::HeaderMismatch(FuncName, sacH2.fname);
	if( sacD.sig && (sacZ.shd.npts!=sacD.shd.npts || sacZ.shd.delta!=sacD.shd.delta ) )
		throw ErrorSR::HeaderMismatch(FuncName, sacD.fname);
}

float StaSacs::RemoveCompliance( const std::string& outinfoname, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(_Eperl, _Eperu);

	// estimate compliance noise
	SacRec Coh_c, Adm_c, Pha_c;
	float cohavg = EstimateCoherence( 'Z', 'D', fmax_comp, Coh_c, Adm_c, Pha_c, true );	// flip sacD when pha is close to PI

	// output 1
	std::vector<SacRec> sacV; // use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	if( output_coh ) sacV = {Coh_c, Adm_c, Pha_c};

	// compliance correction on Z (Z -= coh_c*D)
	ApplyCorrection( 'Z', sacD, Coh_c, Adm_c, Pha_c, fmax_comp );

	// output 2
	if(output_coh) { 
		sacV.push_back(std::move(Adm_c)); sacV.push_back(std::move(Pha_c));
		DumpSACs( sacV, outinfoname );
	}

	return cohavg;
}

float StaSacs::EstimateCoherence( const char c1, const char c2, const float fmax, 
											 SacRec& Coh, SacRec& Adm, SacRec& Pha, bool autoflip ) {
	// compute transfer F for compliance noise
	CalcTransferF(c1, c2, Coh, Adm, Pha);
	// compute average coherence
   float fe = fe_avg<fmax?fe_avg:fmax, flen = 0.5 * (fe-fb_avg);
   float coh1 = flen<0.02 ? CohAvg( Coh, Pha, fb_avg, fe, 0. ) 
					 : std::max( CohAvg( Coh, Pha, fb_avg, fb_avg+flen, 0. ),
									 CohAvg( Coh, Pha, fb_avg+flen, fe, 0. ) );
	if( autoflip ) {
		float coh2 = flen<0.02 ? CohAvg( Coh, Pha, fb_avg, fe, M_PI ) 
						 : std::max( CohAvg( Coh, Pha, fb_avg, fb_avg+flen, M_PI ),
										 CohAvg( Coh, Pha, fb_avg+flen, fe, M_PI ) );
		if( coh2 > coh1 ) {
			float adder = Pha.shd.user1<0 ? M_PI : -M_PI;
			Pha.Transform( [&](float& val){val += adder;} ); 
			coh1 = coh2; Pha.shd.user1 += adder; 
			std::get<0>(sacsM.at(c2))->Mul(-1.);
			std::get<1>(sacsM.at(c2))->clear();
		}
	}
	return coh1;
}

PointC StaSacs::RemoveTilt( const std::string& outinfoname, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(_Eperl, _Eperu);

	// estimate tilt noise
	SacRec Coh_t, Adm_t, Pha_t;
	PointC Pres = EstimateTilt( Coh_t, Adm_t, Pha_t );

	// output 1
	std::vector<SacRec> sacV; // use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	if( output_coh ) sacV = {Coh_t, Adm_t, Pha_t};

	// tilt correction on Z (Z -= coh_t*H1)
	if( Pres.y>cohmin ) ApplyCorrection( 'Z', sacHt, Coh_t, Adm_t, Pha_t, fmax_tilt );

	// output 2
	if(output_coh) {
		sacV.push_back(std::move(Adm_t)); sacV.push_back(std::move(Pha_t));
		DumpSACs( sacV, outinfoname );
	}

	return Pres;
}

PointC StaSacs::EstimateTilt( SacRec& Coh, SacRec& Adm, SacRec& Pha ) {
	// estimate tilt direction (and the corresponding coherence)
	// where Pres.x = direction_t, Pres.y = coh_t
	PointC Pres = EstimateTiltDirection( 1.0 );
	sacHt = SACProject(sacH1, sacH2, Pres.x); sacHtamV.clear();
	// compute transfer F for tilt noise
	Pres.y = EstimateCoherence( 'Z', 't', fmax_tilt, Coh, Adm, Pha );
	return Pres;
}

void StaSacs::test(float tseg, int nsm, float tb, float te) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ modified
	sacZ = DetectNoiseWindows(11., 20., true);
//sacZ.Write("debugZ.SAC");
//sacD.Mul(-1.);
	std::vector<SacRec> sacV; sacV.resize(3);
	SacRec Coh, Adm, Pha;
	CalcTransferF('Z', 'D', Coh, Adm, Pha, tb, te);
	sacV[0] = std::move(Coh); sacV[1] = std::move(Adm); sacV[2] = std::move(Pha);
	DumpSACs( sacV, "testinfo1.txt" );
	if( tb!=NaN && te!=NaN ) { sacZ.cut(tb, te); sacD.cut(tb, te); }
	CalcTransferF2(sacZ, sacD, Coh, Adm, Pha);
	sacV[0] = std::move(Coh); sacV[1] = std::move(Adm); sacV[2] = std::move(Pha);
	DumpSACs( sacV, "testinfo2.txt" );
}

PointC5 StaSacs::RemoveTiltCompliance( const std::string& outinfoname, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(_Eperl, _Eperu);

	// estimate tilt noise
	SacRec Coh_t, Adm_t, Pha_t;
	PointC5 Pres = EstimateTilt( Coh_t, Adm_t, Pha_t );

	// estimate comp noise
	SacRec Coh_c, Adm_c, Pha_c;
	Pres.z = !sacD.sig ? 0. :
				EstimateCoherence( 'Z', 'D', fmax_comp, Coh_c, Adm_c, Pha_c, true );

	// output 1
	std::vector<SacRec> sacV; // use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	if( output_coh ) {
		sacV = {Coh_t, Adm_t, Pha_t};
		sacV.resize(sacD.sig?15:5);
		if( sacD.sig ) { sacV[5] = Coh_c; sacV[6] = Adm_c; sacV[7] = Pha_c; }
	}

	// and apply to the one with larger average coherence first
	auto &coh_t = Pres.y, &coh_c = Pres.z;
	if( coh_t<cohmin && coh_c<cohmin ) {
		//output_coh = false;
	} else if( coh_t > coh_c ) {
		// tilt correction on Z (Z -= coh_t*H1)
		ApplyCorrection( 'Z', sacHt, Coh_t, Adm_t, Pha_t, fmax_tilt );
		sacZinterm = sacZ;
		if(output_coh) { sacV[3] = std::move(Adm_t); sacV[4] = std::move(Pha_t); }
		if( sacD.sig ) { // stop if sacD is empty
			// tilt correction on D (D -= coh_t*H1)
			Pres.z2 = EstimateCoherence( 'D', 't', fmax_tilt, Coh_t, Adm_t, Pha_t, true );
			if(Pres.z2>cohmin) ApplyCorrection( 'D', sacHt, Coh_t, Adm_t, Pha_t, fmax_tilt );
			// compliance correction on Z (Z -= coh_c*D)
			if(output_coh) { sacV[8] = std::move(Adm_c); sacV[9] = std::move(Pha_c); }
			Pres.z3 = EstimateCoherence( 'Z', 'D', fmax_comp, Coh_c, Adm_c, Pha_c, true );
			if(output_coh) { sacV[10] = Coh_c; sacV[11] = Adm_c; sacV[12] = Pha_c; }
			if(Pres.z3>cohmin) ApplyCorrection( 'Z', sacD, Coh_c, Adm_c, Pha_c, fmax_comp );
			if(output_coh) { sacV[13] = std::move(Adm_c); sacV[14] = std::move(Pha_c); }
		}
	} else {
		// compliance correction on Z (Z -= coh_c*D)
		ApplyCorrection( 'Z', sacD, Coh_c, Adm_c, Pha_c, fmax_comp );
		sacZinterm = sacZ;
		if(output_coh) { sacV[8] = std::move(Adm_c); sacV[9] = std::move(Pha_c); }
		// compliance correction on H1 (H1 -= coh_c*D)
		float coh_cH1 = EstimateCoherence( '1', 'D', fmax_comp, Coh_c, Adm_c, Pha_c );
		if(coh_cH1>cohmin) ApplyCorrection( '1', sacD, Coh_c, Adm_c, Pha_c, fmax_comp );
		// compliance correction on H2 (H2 -= coh_c*D)
		float coh_cH2 = EstimateCoherence( '2', 'D', fmax_comp, Coh_c, Adm_c, Pha_c );
		if(coh_cH2>cohmin) ApplyCorrection( '2', sacD, Coh_c, Adm_c, Pha_c, fmax_comp );
		Pres.z2 = sqrt(coh_cH1*coh_cH1 + coh_cH2*coh_cH2);
		if(output_coh) { sacV[3] = std::move(Adm_t); sacV[4] = std::move(Pha_t); }
		// tilt correction on Z (Z -= coh_t*H1)
		PointC Pres2 = EstimateTilt( Coh_t, Adm_t, Pha_t ); 
		Pres.x = Pres2.x; Pres.z3 = Pres2.y;	// update tilt orientation and coh
		if(output_coh) { sacV[10] = Coh_t; sacV[11] = Adm_t; sacV[12] = Pha_t; }
		if(Pres.z3>cohmin) ApplyCorrection( 'Z', sacHt, Coh_t, Adm_t, Pha_t, fmax_tilt );
		if(output_coh) { sacV[13] = std::move(Adm_t); sacV[14] = std::move(Pha_t); }
	}

	// output all three sets of Coh-Adm-Pha when requested
	if(output_coh) DumpSACs( sacV, outinfoname );

	return Pres;
}

RDirect StaSacs::RayleighDirectionality(const float dazi, const std::vector<std::pair<float, float>>& freqRangeV, 
													 float tseg, float twin) {
	_tseg = tseg;
	float tb = sacZ.shd.b, te = sacZ.shd.e; //, twin = 1800.;
	//DetectNoiseWindows(_Eperl, _Eperu);
	trec = { {tb, te} }; clearV();

	//sacZ.RunAvg( 40., _Eperl, _Eperu ); //sacZ.OneBit();
	std::vector<SacRec> sacV{sacZ, sacH1, sacH2};
	//RunAvg(40., 15., 30., sacV);
	RunAvg(40., 12., 20., sacV, true);
	SacRec &sacZt(sacV[0]), &sacH1t(sacV[1]), &sacH2t(sacV[2]); 
	int nrotate = ceil(90./dazi), ncycle = nrotate*4, ntwin = (int)(te-tb-twin)/(0.5*twin) + 1;
	RDirect rdirect( freqRangeV );
	//SacRec sacH1t(sacH1), sacH2t(sacH2);
	for(int irotate=0; irotate<nrotate; irotate++) {
		float azi1 = irotate * dazi, azi2 = azi1 + 90.;
		float azi3 = azi2 + 90., azi4 = azi3 + 90.;
		//const SacRec &sacH1t(sacH1), &sacH2t(sacH2);
		SacRec Coh, Adm, Pha;
		auto computeCohs = [&](const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
									  float tb, float te, float azipos, float azineg) {
			CalcTransferF(sacZt, sacZamV, sacZphV, sac2, sac2amV, sac2phV, Coh, Adm, Pha, tb, te);
/*
if( tb==26300 && te==28100. && azipos==70. ) {
	std::cerr<<"found!"<<std::endl;
	SacRec Coht; Coh.Smooth(0.01, Coht);
	Coht.Write("Coh_6500sec_120deg.SAC");
	SacRec Admt; Adm.Smooth(0.01, Admt);
	Admt.Write("Adm_6500sec_120deg.SAC");
	SacRec Phat(Pha); Phat.Wrap();
	Phat.Write("Pha_6500sec_120deg.SAC");
}
*/
			auto& donut = rdirect[std::make_pair(tb, te)];
			auto &resVpos = donut[azipos], &resVneg = donut[azineg];
			for( const auto& pair : freqRangeV ) {
				//cohAvg( pair.first, pair.second, resVpos[i], resVneg[i] ); i++;
				resVpos.push_back( CohAvg( Coh, Pha, pair.first, pair.second, -PIo2 ) );
				resVneg.push_back( CohAvg( Coh, Pha, pair.first, pair.second, PIo2 ) );
			}
		};
		for(float twinb=tb,twine=tb+twin; twine<=te; twinb+=0.5*twin, twine+=0.5*twin) {
			try {
				computeCohs(sacH1t, sacH1amV, sacH1phV, twinb, twine, azi1, azi3);
				computeCohs(sacH2t, sacH2amV, sacH2phV, twinb, twine, azi2, azi4); 
			} catch(ErrorSR::Base &e) {}
		}
		sacH1amV.clear();	sacH2amV.clear();
		SACRotate(sacH1t, sacH2t, dazi); 
	}
	sacZamV.clear();
	return rdirect;
}

PointC StaSacs::EstimateTiltDirection( const float ddeg ) const {
	// calc the average_coh - azimuth curve (by rotating the two horizontal channels)
	//std::cerr<<"zeropha: "<<zeropha<<std::endl;
	int nrotate = ceil(90./ddeg); 
	std::vector<PointC> cohV(nrotate*4);
	SacRec sacH1t(sacH1), sacH2t(sacH2);
	float fe_tilt = fe_avg<fmax_tilt ? fe_avg : fmax_tilt;
	for(int irotate=0; ; irotate++) {
		auto &p1 = cohV[irotate], &p2 = cohV[irotate+nrotate];
		auto &p3 = cohV[irotate+nrotate*2], &p4 = cohV[irotate+nrotate*3];
		p1.x = irotate*ddeg; p2.x = p1.x + 90.;
		p3.x = p1.x + 180.; p4.x = p1.x + 270.;
		SacRec Coh, Adm, Pha;
		//p1.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1t, sacH1amV, sacH1phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		//p2.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2t, sacH2amV, sacH2phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH1t, sacH1amV, sacH1phV, Coh, Adm, Pha);
		p1.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, 0. );
		p3.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, M_PI );
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH2t, sacH2amV, sacH2phV, Coh, Adm, Pha);
		p2.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, 0. );
		p4.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, M_PI );
		sacH1amV.clear(); sacH2amV.clear();
		//std::cerr<<p1<<"\n"<<p2<<"\n"<<p3<<"\n"<<p4<<std::endl;
		if(irotate>=nrotate-1) break;
		SACRotate(sacH1t, sacH2t, ddeg); 
	}
	// find the 3 points near the peak
	std::nth_element( cohV.begin(), cohV.begin()+2, cohV.end(), 
			[](const PointC& p1, const PointC& p2){ return p1.y>p2.y; } );
	// find the maximum
	PointC &pmax = *( std::max_element( cohV.begin(), cohV.begin()+3, 
				[](const PointC& p1, const PointC& p2){ return p1.y<p2.y; } ) );
	// return if coh-avg-peak < cohmin
	//if( pmax.y < cohmin )	return pmax;
	PointC &p1 = cohV[0], &p2 = cohV[1], &p3 = cohV[2];
	// assume a 360 deg periodicity
	p2.x -= 360. * nint((p2.x-p1.x)/360.);
	p3.x -= 360. * nint((p3.x-p1.x)/360.);
	// check if p1 p2 p3 are continguous
	if( fabs(p2.x-p1.x) + fabs(p3.x-p1.x) + fabs(p3.x-p2.x) > ddeg*4.1 ) return pmax;
	// find tilt direction (with peak coh)
	return Parabola(p1, p2, p3).Vertex();
}

// compute transfer function for each window (twin) from segments (tseg) within that window
void StaSacs::ReshapeByZHPhaseDiff( float twin, float tseg ) {
	if( twin < tseg ) throw ErrorSR::BadParam( FuncName, "twin has to be >= tseg" );
	_tseg = tseg;
	// zero-out earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(_Eperl, _Eperu);
	// loop through segments
	SacRec Coh1, Adm1, Pha1, Coh2, Adm2, Pha2;
	float twinb_e, twine_e = sacZ.shd.b;
	for( float twinb=sacZ.shd.b; twinb<sacZ.shd.e; twinb+=twin ) {
		bool succ = true;
		try {
			CalcTransferF('Z', '1', Coh1, Adm1, Pha1, twinb, twinb+twin );
			CalcTransferF('Z', '2', Coh2, Adm2, Pha2, twinb, twinb+twin );
		} catch( ErrorSR::InsufData& e ) { succ = false; }
		// zero out [twine_e_old, twinb_e)
		sacZ.Transform( [&](float& val){val = 0.;}, sacZ.Index(twine_e), sacZ.Index(Coh1.shd.user2) );
		twinb_e = Coh1.shd.user2; twine_e = Coh1.shd.user3;
		if( succ ) {
			// FFT on sacZ for reshaping
			SacRec sacZw, sacZwam, sacZwph;
			sacZ.cut(twinb_e, twine_e, sacZw); sacZw.ToAmPh(sacZwam, sacZwph);
			auto sigZwam = sacZwam.sig.get();
			// spectrum weights from s1 ( defined as (pi/2 - dist(pha, +/-pi/2)) * coh )
			auto sigPha = Pha1.sig.get();
			Coh1.Transform2i( [&](const int i, float& val) { 
				val *= ( PIo2 - fabs(fabs(fmod(sigPha[i], M_PI)) - PIo2) );	} );
				//val *= ( fabs(fabs(fmod(sigPha[i], M_PI)) - PIo2) );	} );
			Coh1.Smooth(0.01, false, fb); 
			// spectrum weights from s2
			sigPha = Pha2.sig.get();
			Coh2.Transform2i( [&](const int i, float& val) { 
				val *= ( PIo2 - fabs(fabs(fmod(sigPha[i], M_PI)) - PIo2) );	} );
				//val *= ( fabs(fabs(fmod(sigPha[i], M_PI)) - PIo2) );	} );
			Coh2.Smooth(0.01, false, fb);
			// use whichever is greater
			Coh1.PullUpTo(Coh2);
			// apply the weighting
			auto sigCoh = Coh1.sig.get();
			float dfast = sacZwam.shd.delta, dslow = Coh1.shd.delta;
			float ffast = 0, fslow = dslow;
			float wlast = sigCoh[0];
			for(int ifast=0,islow=1; islow<Coh1.shd.npts; islow++,fslow+=dslow) {
				float w = sigCoh[islow];
				//float w = sigCoh[islow] * (fslow<0.105&&fslow>=0.025 ? 50 : 1);
				for(; ffast<=fslow; ffast+=dfast,ifast++)
					sigZwam[ifast] *= wlast + (w-wlast) * (ffast-fslow+dslow) / dslow;
				wlast = w;
			}
			sacZw.FromAmPh(sacZwam, sacZwph);
			sacZ.merge(sacZw);
		} else {
			sacZ.Transform( [&](float& val){val *= 0.1;}, sacZ.Index(twinb_e), sacZ.Index(twine_e) );
		}
	}
	sacZ.Transform( [&](float& val){val = 0.;}, sacZ.Index(twine_e), sacZ.shd.npts );
}

void StaSacs::Write( const std::string& foutZ, const std::string& foutZinterm,
							const std::string foutH1, const std::string& foutH2, const std::string foutD ) const {
	// lambda function to convert sac back to original type and write
	SacRec sact;
	auto WriteSAC = [&]( const SacRec& sac, const std::string& outname ) {
		if( outname.empty() || !sac.sig ) return;
		sact = sac;
		switch(sactype) {
			case 2:	// acc: differentiate twice
				sact.Differentiate();
			case 1:	// vel: differentiate once
				sact.Differentiate();
			case 0:	// dis: do nothing
				break;
			default:
				throw ErrorSR::BadParam( FuncName, "unknown sactype("+std::to_string(sactype)+")." );
		}
		sact.Write(outname);
	};
	WriteSAC(sacZ, foutZ); WriteSAC(sacZinterm, foutZinterm);
	WriteSAC(sacH1, foutH1); WriteSAC(sacH2, foutH2); 
	if( ! foutD.empty() ) { sact = sacD; sact.Write(foutD); }
}


/* --------------------- Private Methods --------------------- */

// compute valid time segments (removing earthquakes)
SacRec StaSacs::DetectNoiseWindows( float Eperl, float Eperu, bool ataper ) const {
	clearV();	// clear all old sac-seg vectors (which were based on the old rec windows)
	SacRec sacZcut; std::vector<int> recb, rece;
	sacZ.EqkCut(sacZcut, recb, rece, Eperl, Eperu, 2.5, ataper); // apply taper
	for(int i=0; i<recb.size(); i++) {
		float tb = sacZ.X(recb[i]), te = sacZ.X(rece[i]), tmid = 0.5*(tb+te);
		trec.push_back({tb, te});
		//if( docut ) sacZ.cosTaperL(tb, tmid, false); sacZ.cosTaperR(tmid, te, false);
		//std::cerr<<"   rec window("<<sacZ.fname<<"): "<<sacZ.X(recb[i])<<" "<<sacZ.X(rece[i])<<std::endl; 
	}
	return sacZcut;
}

// segmentize a single channel
void StaSacs::Segmentize(const SacRec& sac, std::vector<SacRec>& sacamV, std::vector<SacRec>& sacphV) const {
	if( !sacamV.empty() && !sacphV.empty() && sacamV.size()==sacphV.size() ) return;
	if( ! sac.sig ) throw ErrorSR::EmptySig(FuncName, sac.fname);
	if( _tseg <= 0 ) throw ErrorSR::BadParam(FuncName, "Invalid/Uninitialized _tseg("+std::to_string(_tseg)+")");
	// clear old segments
	sacamV.clear(); sacphV.clear();
	float ttaper = std::min(_tseg*0.4, 200.);
	for(int i=0; i<trec.size(); i++) {
		float twin_b = trec[i].first, twin_e = trec[i].second, tb;
		//std::cout<<twin_b<<" "<<twin_e<<std::endl;
		for(tb=twin_b; tb<=twin_e-_tseg; tb+=_tseg) {
			float te = tb+_tseg;
			if( tb<sac.shd.b || te>sac.shd.e ) continue;
			SacRec sac_tmp, sac_am, sac_ph;
			// cut and FFT
			sac.cut(tb,te, sac_tmp);
			// Tukey window to reduce broadband spectral leakage
			sac_tmp.cosTaperL(tb, tb+ttaper);
			sac_tmp.cosTaperR(te-ttaper, te);
			sac_tmp.ToAmPh(sac_am, sac_ph);
			sac_am.Mul(1.0e-5);	// to prevent floating-point overflow
			auto& shd = sac_am.shd; shd.user2 = tb; shd.user3 = te;
			sacamV.push_back(std::move(sac_am)); 
			sacphV.push_back(std::move(sac_ph));
		}
	}
	//std::cerr<<sac.fname<<" "<<sacamV.size()<<" "<<sacphV.size()<<std::endl;
}

std::pair<int, int> StaSacs::GetSegRange(float &tb, float &te, std::vector<SacRec>& sacV) const {
	if( sacV.empty() ) throw ErrorSR::EmptySig(FuncName, "empty sac-seg vector");
	int nseg = sacV.size(), isegb = 0, isege = nseg;
	if( tb!=NaN && te!=NaN ) {
		isegb = -1;	float teeff;
		for(int iseg=0; iseg<nseg; iseg++) {
			const auto& shd = sacV[iseg].shd;
			const float tmid = 0.5 * (shd.user2+shd.user3);
			if( isegb<0 && tmid>=tb ) {
				isegb = iseg; tb = shd.user2;
			}
			if( tmid >= te ) { isege = iseg; break; }
			teeff = shd.user3;
		}
		if( isegb==-1 ) { tb = teeff; isegb=isege; }
		te = teeff;
	} else {
		tb = sacV[0].shd.user2; te = sacV.back().shd.user3;
	}
	return std::make_pair(isegb, isege);
}

void StaSacs::FitIntoParabola( SacRec& sac1, const SacRec& sacsigma, const float fmin, const float fmax, const int nparab ) const {
	// square sigma to get weights
	SacRec sacw(sacsigma); sacw.Mulf(sacsigma);
	// dump data to vector of PointCs
	std::vector<PointC> dataV;
	sac1.Dump( dataV, fmin, fmax );
	// check header
	if( dataV.size() != sacw.Index(fmax) - sacw.Index(fmin) ) 
		throw std::runtime_error(std::string("Error(")+FuncName+"): npts mismatch");
	auto sigwsac = sacw.sig.get();
	// assign weight into dataV
	// and set logscale to x axis
	int isig = sacw.Index(fmin); 
	for( auto &p : dataV ) { 
		p.x = std::log(p.x); 
		p.z = sigwsac[isig++]; 
	}
	// divide [fmin,fmax) into nparab segments such that each segment has nvalid/nparab valid points
	// count valid points (with weight>0.)
	std::vector<int> validC(dataV.size());
	validC[0] = dataV[0].z>0.01;
	for(int i=1; i<validC.size(); i++) 
		validC[i] = validC[i-1] + (int)(dataV[i].z>0.01);
	// fit parabola(s)
	float cstep = validC.back() / (float)nparab;
	for(int i=0; i<nparab; i++) {
		auto Ib = dataV.begin() + ( std::lower_bound(validC.begin(), validC.end(), (int)(cstep*i+0.5)) - validC.begin() );
		auto Ie = dataV.begin() + ( std::upper_bound(validC.begin(), validC.end(), (int)(cstep*(i+1)+0.5)) - validC.begin() );
		Parabola parab; float rmsw; parab.Fit(Ib, Ie, rmsw);
		// write parabola to sac1
		float f1 = i==0 ? fmin*0.9 : std::exp(Ib->x);
		float f2 = i==nparab-1 ? fmax*1.1 : std::exp(Ie->x);
		sac1.Transform2( [&](const float x, float& y) {
				y = parab[std::log(x)];
		},	sac1.Index(f1), sac1.Index(f2) );
	}
	//sac1.Write("debug_Adm_BHZBH1.SAC_pb");
}

void StaSacs::CalcTransferF2( const SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha ) const {
	SacRec sac1_am, sac1_ph; sac1.ToAmPh(sac1_am, sac1_ph);
	SacRec sac2_am, sac2_ph; sac2.ToAmPh(sac2_am, sac2_ph);
	//sac1_am.Mul(1.0e-5);	sac2_am.Mul(1.0e-5); // to prevent floating-point overflow
	// calculate autospectral density functions Gss, Grr and
	// the one-sided cross-spectral density function Grs
	SacRec Grr(sac1_am);  Grr.Mulf(sac1_am);
	SacRec Gss(sac2_am);  Gss.Mulf(sac2_am);
	SacRec GrsA(sac2_am); GrsA.Mulf(sac1_am);
	SacRec GrsP(sac2_ph); GrsP.Subf(sac1_ph);
	// convert GrsAmp&GrsPha to GrsR&GrsI in place
   AmPhToReIm( GrsA, GrsP );
	// freqency domain smoothing for statistical stability
	int nsm = 20; float fhlen = nsm * GrsA.shd.delta;
	Grr.Smooth(fhlen, true, fb); Gss.Smooth(fhlen, true, fb);
	GrsA.Smooth(fhlen, false, fb); GrsP.Smooth(fhlen, false, fb);
	// convert GrsR&GrsI back to GrsAmp&GrsPha
	ReImToAmPh( GrsA, GrsP );
	// construct coherence, admittance, and phase from the spectral density functions
std::cerr<<" fhlen = "<<fhlen<<std::endl;
GrsA.Write("debugRSA.SAC"); GrsP.Write("debugRSP.SAC");
Grr.Write("debugRR.SAC"); Gss.Write("debugSS.SAC");
	Pha = std::move(GrsP); Adm = GrsA; Adm.Divf( Gss );
	Coh = Adm; Coh.Mulf( GrsA ); Coh.Divf( Grr ); Coh.sqrt();
Adm.Write("debugAdm.SAC"); Coh.Write("debugCoh.SAC");
	// smooth Coh
	Coh.Smooth(0.002, false, fb); 
}

void StaSacs::CalcTransferF( const char c1, const char c2, SacRec& Coh, SacRec& Adm, SacRec& Pha, 
									  const float tb, const float te ) const {
	auto sacs1 = sacsM.at(c1); const SacRec& sac1 = *(std::get<0>(sacs1));
	auto &sac1amV = *(std::get<1>(sacs1)), &sac1phV = *(std::get<2>(sacs1));
	auto sacs2 = sacsM.at(c2); const SacRec& sac2 = *(std::get<0>(sacs2));
	auto &sac2amV = *(std::get<1>(sacs2)), &sac2phV = *(std::get<2>(sacs2));
	CalcTransferF( sac1, sac1amV, sac1phV, sac2, sac2amV, sac2phV, Coh, Adm, Pha, tb, te );
}
void StaSacs::CalcTransferF( const SacRec& sac1, std::vector<SacRec>& sac1amV, std::vector<SacRec>& sac1phV,
										const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
										SacRec& Coh, SacRec& Adm, SacRec& Pha, float tb, float te ) const {
	Segmentize( sac1, sac1amV, sac1phV );
   Segmentize( sac2, sac2amV, sac2phV );
	// check consistency of segment Numbers
	int nseg = sac1amV.size();
	if( nseg!=sac1phV.size() || nseg!=sac2amV.size() || nseg!=sac2phV.size() ) {
		std::string info = std::to_string(nseg)+" "+std::to_string(sac1phV.size())+" "+std::to_string(sac2amV.size())+" "+std::to_string(sac2phV.size());
		throw std::runtime_error(std::string("Error(")+FuncName+"): No. segments mismatch ("+info+")");
	}

	// only use segments with (shd.b+shd.e)/2 between [tb, te)
	int isegb, isege; std::tie(isegb, isege) = GetSegRange(tb, te, sac1amV);
	// minimum No. segments requirement
	if( isege-isegb < 5 ) {
		Coh.shd.user2 = tb; Coh.shd.user3 = te;
		throw ErrorSR::InsufData(FuncName, "No. segments = "+std::to_string(isege-isegb));
	}

	// calculate autospectral density functions Gss, Grr and
	// the one-sided cross-spectral density function Grs,
	// summing up all segments
	SacRec Gss, Grr, GrsR, GrsI;
	const SacRec& sac0 = sac2amV[0];
	Gss.MutateAs(sac0); Grr.MutateAs(sac0); 
	GrsR.MutateAs(sac0); GrsI.MutateAs(sac0);
	for(int iseg=isegb; iseg<isege; iseg++) {
		//const auto &shd = sac1amV[iseg].shd; std::cerr<<"   seg window: "<<shd.user2<<" "<<shd.user3<<std::endl;
		// compute Gss, Grr, Grs_am, Grs_ph for the current seg
		SacRec sacrr_am(sac1amV[iseg]); sacrr_am.Mulf(sac1amV[iseg]);
		SacRec sacss_am(sac2amV[iseg]); sacss_am.Mulf(sac2amV[iseg]);
		SacRec sacrs_am(sac2amV[iseg]); sacrs_am.Mulf(sac1amV[iseg]);
		SacRec sacrs_ph(sac2phV[iseg]); sacrs_ph.Subf(sac1phV[iseg]);
		// Add Gss and Grr to the sum
		Grr.Addf( sacrr_am ); Gss.Addf( sacss_am );
		// convert am&ph to re&im
		AmPhToReIm( sacrs_am, sacrs_ph );
		// and add
		GrsR.Addf( sacrs_am ); GrsI.Addf( sacrs_ph );
	}
	// convert GrsR&GrsI to GrsAmp&GrsPha in place
	ReImToAmPh( GrsR, GrsI );
	// construct coherence, admittance, and phase from the spectral density functions
	Pha = std::move(GrsI); Adm = GrsR; Adm.Divf( Gss );
	Coh = Adm; Coh.Mulf( GrsR ); Coh.Divf( Grr ); Coh.sqrt();
	Coh.shd.user2 = tb; Coh.shd.user3 = te;
	// smooth Coh
	//Coh.Smooth(0.002, false, fb);
	// re-centralize phase
	float phac = Pha.MeanPha(fb_avg, fe_avg);
	Pha.shd.user1 = phac;
	const float twopi = M_PI*2.;
	Pha.Transform( [&](float& val){
		if(val<phac-M_PI) val+=twopi;
		else if(val>=phac+M_PI) val-=twopi;
		} );
}

float StaSacs::CohAvg(const SacRec& Coh, const SacRec& Pha, float fb, float fe, float pha0) const {
	int ib = Coh.Index(fb), ie = Coh.Index(fe);
	if( ib >= ie ) throw ErrorSR::BadParam(FuncName, "fb >= fe");
	float mean = 0.;
	auto sigcoh = Coh.sig.get(), sigpha = Pha.sig.get();
	for(int i=ib; i<ie; i++)
		mean += sigcoh[i] * std::max(cos(sigpha[i]-pha0),0.);
	return mean / (ie-ib);
}

void StaSacs::ApplyCorrection( const char c, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha, 
										 float fmax, const float pha0 ) {
	// get reference to sac1 and clean up seg vectors
	auto sacs1 = sacsM.at(c); SacRec& sac1 = *(std::get<0>(sacs1));
	auto &sac1amV = *(std::get<1>(sacs1)), &sac1phV = *(std::get<2>(sacs1));
	sac1amV.clear(); sac1phV.clear(); // seg vectors for sac1 need to be updated
	#ifdef DEBUG
	auto chn1 = sac1.chname();
	auto chn2 = sac2.chname();
	Adm.Write("debug_Adm_"+chn1+chn2+".SAC");
	Pha.Write("debug_Pha_"+chn1+chn2+".SAC");
	Coh.Write("debug_Coh_"+chn1+chn2+".SAC");
	#endif

	// smooth/modify Coh as sigma for parabola fitting
	const auto sigpha = Pha.sig.get(); 
	Coh.Transform2i( [&](const int i, float& val){ 
		val *= cos(sigpha[i]-pha0); 
	}, Coh.Index(fb), Coh.Index(fmax) );
	Coh.Smooth(0.002, false, fb);

	// invalidate fpoints with coh<cohmin-0.2
	// and shrink [fmin, fmax]
	int ib = Coh.Index(fb), ie = Coh.Index(fmax), nvalid = 0;
	Coh.Transform( [&](float& val){ 
		if(val<cohmin-0.2) val = 0.; 
		else nvalid++;
	}, ib, ie );
	float coreperc = 0.05;	
	auto sigcoh = Coh.sig.get();
	for(int ivalid=nvalid*coreperc; ib<ie&&ivalid!=0; ib++ ) ivalid -= sigcoh[ib]>0.;
	float fmin = Coh.X(ib);
	for(int ivalid=nvalid*coreperc; ie>ib&&ivalid!=0; ie-- ) ivalid -= sigcoh[ie]>0.;
	fmax = Coh.X(ie);
	float fex = (fmax-fmin)*1.2*coreperc/(1.-2*coreperc); 
	fmin -= fex; fmax += fex;

	// fit phase with a single parabola
	FitIntoParabola( Pha, Coh, fmin, fmax, 1 );
	// and admittance with two parabolas
	FitIntoParabola( Adm, Coh, fmin, fmax, (fmax-fmin>0.04 ? 2 : 1) );
	//Adm.Smooth(0.002, false, fmin);
	#ifdef DEBUG
	Coh.Write("debug_CohP_"+chn1+chn2+".SAC");
	Pha.Write("debug_PhaP_"+chn1+chn2+".SAC");
	Adm.Write("debug_AdmP_"+chn1+chn2+".SAC");
	#endif

	// apply transfer function in between fmin and fmax Hz
	SacRec sac2_am, sac2_ph; sac2.ToAmPh(sac2_am, sac2_ph);
	// predicted spectrum
	int npts_ratio = (int)floor(sac2_am.shd.npts/(float)Adm.shd.npts+0.5);
	SacRec AdmInt, PhaInt;
	Adm.Interpolate(npts_ratio, AdmInt); 
	Pha.Interpolate(npts_ratio, PhaInt);
	SacRec sacP_ph(sac2_ph); sacP_ph.Subf(PhaInt); sacP_ph.Wrap(); //sacP_ph.Mul(-1.);
	SacRec sacP_am(sac2_am); sacP_am.Mulf(AdmInt);
	// Tukey window
	sacP_am.cosTaperL(fmin*0.9, fmin);
	sacP_am.cosTaperR(fmax, fmax*1.1);
	SacRec sacP; sacP.shd = sac1.shd;
	sacP.FromAmPh(sacP_am, sacP_ph);
	sac1.Subf(sacP);
	#ifdef DEBUG
	sacP.Write("debug_Pred_"+chn1+chn2+".SAC");
	sac1.Write("debug_Crct_"+chn1+chn2+".SAC");
	#endif
}

#endif
