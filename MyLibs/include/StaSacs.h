#ifndef STASACS_H
#define STASACS_H

#include "SacRec.h"
#include <iostream>
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
				const int sactypein = 0, const float water_depth_in = -1, const float azi_H1 = 0. );

	// cut-off at freq = sqrt(g/(2 pi d)) where water_depth = ingragravity_wavelength
	float fcutoffCompliance() const {
		if( water_depth < 0. ) throw ErrorSR::BadParam(FuncName, "invalid (non-positive) water depth");
		return 1.2489 / sqrt(water_depth);
	}

	float RemoveCompliance( const std::string& outinfoname = "", bool zeropha = true, 
									float Eperl=10., float Eperu=40., float tseg=2500. );

	PointC RemoveTilt( const std::string& outinfoname = "", bool zeropha = true,
							 float Eperl=10., float Eperu=40., float tseg=2500. );

	void test(float tseg, int nsm, float tb = NaN, float te = NaN);
	PointC5 RemoveTiltCompliance( const std::string& outinfoname = "", float Eperl=10., float Eperu=40., float tseg=2500. );

	std::vector<std::vector<float>> RayleighDirectionality(const float dazi, 
			const std::vector<std::pair<float, float>>& freqRangeV) const;

	PointC EstimateTiltDirection( const bool zeropha, const float ddeg ) const;

	// compute transfer function for each window (twin) from segments (tseg) within that window
	void ReshapeByZHPhaseDiff( float twin, float tseg, float Eperl=10., float Eperu=40. );

	//std::vector<int> RemoveTilts( float Eperl=10., float Eperu=40., float tseg=2500. );

	void Write( const std::string& foutZ, const std::string& foutZinterm="",
					const std::string foutH1="", const std::string& foutH2="", const std::string foutD="" );

	void clear() {
		//sactype = 0; water_depth = -1.;
		trec.clear();
		sacZ.clear(); sacZinterm.clear(); sacH1.clear(); sacH2.clear(); sacD.clear();
		sacZamV.clear(); sacH1amV.clear(); sacH2amV.clear(); sacDamV.clear();
		sacZphV.clear(); sacH1phV.clear(); sacH2phV.clear(); sacDphV.clear();
	}

protected:
	static constexpr float NaN = SacRec::NaN;
	static constexpr float PIo2 = M_PI * 0.5; 
	static const int nfcorr_min = 10;
	static constexpr float fb = 0.01;		// signal below fb is discarded
	static constexpr float cohmin = 0.5;
	static constexpr float fmax_tilt = 0.105;
	static constexpr float fb_avg = 0.01;	// cohavg is computed between fb_avg and fe_avg Hz
	static constexpr float fe_avg = 0.05;	
	//static constexpr float cohmin = 0.2;
	//static constexpr float fmax_tilt = 0.4;
	//static constexpr float fb_avg = 0.2;	// cohavg is computed between fb_avg and fe_avg Hz
	//static constexpr float fe_avg = 0.4;	

private:
	SacRec sacZ, sacH1, sacH2, sacD;
	SacRec sacZinterm;	// stores intermediate Z after removing the stronger of compliance/tilt
	const int sactype;
	const float water_depth;
	std::vector<std::pair<float, float>> trec; 
	float _tseg = -1.;
	mutable std::vector<SacRec> sacZamV, sacH1amV, sacH2amV, sacDamV;
	mutable std::vector<SacRec> sacZphV, sacH1phV, sacH2phV, sacDphV;

	int nint( float val ) const { return (int)floor(val+0.5); }

	// segmentize a single channel
	void Segmentize(const SacRec& sac, std::vector<SacRec>& sacamV, std::vector<SacRec>& sacphV) const;

	// compute valid time segments (removing earthquakes)
	void DetectNoiseWindows( float Eperl=10., float Eperu=40., bool docut = false ) {
		SacRec sacZcut; std::vector<int> recb, rece;
		sacZ.EqkCut(sacZcut, recb, rece, Eperl, Eperu, 2.5, docut);	// apply taper only when docut==true
		if( docut ) sacZ = std::move(sacZcut);
		for(int i=0; i<recb.size(); i++) {
			float tb = sacZ.X(recb[i]), te = sacZ.X(rece[i]), tmid = 0.5*(tb+te);
			trec.push_back({tb, te});
			//if( docut ) sacZ.cosTaperL(tb, tmid, false); sacZ.cosTaperR(tmid, te, false);
std::cerr<<"   rec window: "<<sacZ.X(recb[i])<<" "<<sacZ.X(rece[i])<<std::endl; 
		}
	}

	void FitIntoParabola( SacRec& sac1, const SacRec& sacsigma, const float fmax, const int nparab = 1 ) const;

	void CalcTransferF( const SacRec& sac1, std::vector<SacRec>& sac1amV, std::vector<SacRec>& sac1phV,
								const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
								SacRec& Coh, SacRec& Adm, SacRec& Pha, const float tb=NaN, const float te=NaN ) const;
	void CalcTransferF2( const SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha ) const;
	float CohAvg(SacRec& Coh, SacRec& Pha, float fb, float fe, float pha0) const;

	void ApplyCorrection( SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha, const float fmax );

};


/* --------------------- Public Methods --------------------- */

StaSacs::StaSacs( const std::string& fnameZ, const std::string& fnameH1, 
						const std::string& fnameH2, const std::string& fnameD,
						const int sactypein, const float water_depth_in, const float azi_H1 )
						: sacZ(fnameZ), sacH1(fnameH1), sacH2(fnameH2), sacD(fnameD)
						, sactype(sactypein), water_depth(water_depth_in) {
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
		throw ErrorSR::HeaderMismatch(FuncName, "H1");
	if( sacH2.sig && (sacZ.shd.npts!=sacH2.shd.npts || sacZ.shd.delta!=sacH2.shd.delta ) )
		throw ErrorSR::HeaderMismatch(FuncName, "H2");
	if( sacD.sig && (sacZ.shd.npts!=sacD.shd.npts || sacZ.shd.delta!=sacD.shd.delta ) )
		throw ErrorSR::HeaderMismatch(FuncName, "D");
}

float StaSacs::RemoveCompliance( const std::string& outinfoname, bool zeropha, 
											float Eperl, float Eperu, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(Eperl, Eperu);

	// if output_coh is on, use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	std::vector<SacRec> sacV; 
	if(output_coh) sacV.resize(5);

	// compute transfer F for compliance noise
	SacRec Coh_c, Adm_c, Pha_c;
	float fmax_comp = 0.;
	// cut-off at freq = sqrt(g/(2 pi d)) where water_depth = ingragravity_wavelength
	fmax_comp = fcutoffCompliance();
	//float cohavg = CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c, fmax_comp, zeropha);
	CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
	float fe = fe_avg<fmax_comp?fe_avg:fmax_comp;
	float cohavg = CohAvg( Coh_c, Pha_c, fb_avg, fe, 0. );
	if(output_coh) { sacV[0] = Coh_c; sacV[1] = Adm_c; sacV[2] = Pha_c; }

	// compliance correction on Z (Z -= coh_c*D)
	ApplyCorrection( sacZ, sacD, Coh_c, Adm_c, Pha_c, fmax_comp ); sacZamV.clear();
	if(output_coh) { 
		sacV[3] = std::move(Adm_c); sacV[4] = std::move(Pha_c); 
		DumpSACs( sacV, outinfoname );
	}

	return cohavg;
}

PointC StaSacs::RemoveTilt( const std::string& outinfoname, bool zeropha,
									 float Eperl, float Eperu, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(Eperl, Eperu);

	// if output_coh is on, use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	std::vector<SacRec> sacV; 
	if(output_coh) sacV.resize(5);

	// compute transfer F for tilt noise
	// estimate tilt direction (and the corresponding coherence)
	// where Pres.x = direction_t, Pres.y = coh_t
	PointC Pres = EstimateTiltDirection( zeropha, 1.0 );
	SacRec Coh_t, Adm_t, Pha_t;
	if( Pres.y>cohmin || output_coh ) {	// re-calc transfer F in the tilt direction
		SACRotate(sacH1, sacH2, Pres.x); sacH1amV.clear(); sacH2amV.clear();
		//Pres.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t, fmax_tilt, zeropha);
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t, fmax_tilt);
		float fe = fe_avg<fmax_tilt?fe_avg:fmax_tilt;
		Pres.y = CohAvg( Coh_t, Pha_t, fb_avg, fe, 0. );
	}
	if(output_coh) { sacV[0] = Coh_t; sacV[1] = Adm_t; sacV[2] = Pha_t; }

	// tilt correction on Z (Z -= coh_t*H1)
	if( Pres.y>cohmin ) ApplyCorrection( sacZ, sacH1, Coh_t, Adm_t, Pha_t, fmax_tilt ); sacZamV.clear();

	// rotate H1 and H2 back to their initial azimuth
	if( Pres.y>cohmin || output_coh ) { SACRotate(sacH1, sacH2, -Pres.x); sacH1amV.clear(); sacH2amV.clear(); }

	// output Coh-Adm-Pha when requested
	if(output_coh) {
		sacV[3] = std::move(Adm_t); sacV[4] = std::move(Pha_t);
		DumpSACs( sacV, outinfoname );
	}

	return Pres;
}

void StaSacs::test(float tseg, int nsm, float tb, float te) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ modified
	DetectNoiseWindows(11., 20., true);
//sacZ.Write("debugZ.SAC");
//sacD.Mul(-1.);
	if( tb!=NaN && te!=NaN ) { sacZ.cut(tb, te); sacD.cut(tb, te); }
	std::vector<SacRec> sacV; sacV.resize(3);
	SacRec Coh, Adm, Pha;
	CalcTransferF2(sacZ, sacD, Coh, Adm, Pha);
	sacV[0] = std::move(Coh); sacV[1] = std::move(Adm); sacV[2] = std::move(Pha);
	DumpSACs( sacV, "testinfo2.txt" );
	CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh, Adm, Pha);
	sacV[0] = std::move(Coh); sacV[1] = std::move(Adm); sacV[2] = std::move(Pha);
	DumpSACs( sacV, "testinfo1.txt" );
}

PointC5 StaSacs::RemoveTiltCompliance( const std::string& outinfoname, float Eperl, float Eperu, float tseg ) {
	_tseg = tseg;
	// search for earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(Eperl, Eperu);

	// if output_coh is on, use a vector to store sacs to be written later
	bool output_coh = !outinfoname.empty();
	std::vector<SacRec> sacV; 
	if(output_coh) sacV.resize(sacD.sig?15:5);

	// f range for computing coh_avg
	float fmax_comp = fcutoffCompliance();
	float fe_tilt = fe_avg<fmax_tilt ? fe_avg : fmax_tilt;
	float fe_comp = fe_avg<fmax_comp ? fe_avg : fmax_comp;

	// compute transfer F for tilt noise
	// estimate tilt direction (and the corresponding coherence)
	// where Pres.x = direction_t, Pres.y = coh_t
	PointC5 Pres = EstimateTiltDirection( true, 1.0 );
	SACRotate(sacH1, sacH2, Pres.x); sacH1amV.clear(); sacH2amV.clear();
	SacRec Coh_t, Adm_t, Pha_t;
	//if( Pres.y>cohmin || output_coh ) {	// re-calc transfer F in the tilt direction
	//Pres.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t, fmax_tilt);
	CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
	Pres.y = CohAvg( Coh_t, Pha_t, fb_avg, fe_tilt, 0. );
	//}
	if(output_coh) { sacV[0] = Coh_t; sacV[1] = Adm_t; sacV[2] = Pha_t; }
	auto &coh_t = Pres.y;

	// compute transfer F for compliance noise
	SacRec Coh_c, Adm_c, Pha_c;
	if( sacD.sig ) {
		// cut-off at freq = sqrt(g/(2 pi d)) where water_depth = infragravity_wavelength
		//Pres.z = CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c, fmax_comp);
		CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
		Pres.z = CohAvg( Coh_c, Pha_c, fb_avg, fe_comp, 0. );
		if(output_coh) { sacV[5] = Coh_c; sacV[6] = Adm_c; sacV[7] = Pha_c; }
	} else {
		Pres.z = 0.;
	}
	auto &coh_c = Pres.z;

	//std::cerr<<" Initial coh_t & coh_c = "<<coh_t<<" "<<coh_c<<std::endl;

	// and apply to the one with larger average coherence first
	if( coh_t<cohmin && coh_c<cohmin ) {
		//output_coh = false;
	} else if( coh_t > coh_c ) {
		// tilt correction on Z (Z -= coh_t*H1)
		ApplyCorrection( sacZ, sacH1, Coh_t, Adm_t, Pha_t, fmax_tilt ); sacZamV.clear();
		sacZinterm = sacZ;
		if(output_coh) { sacV[3] = std::move(Adm_t); sacV[4] = std::move(Pha_t); }
		if( sacD.sig ) { // stop if sacD is empty
			// tilt correction on D (D -= coh_t*H1)
			//Pres.z2 = CalcTransferF(sacD, sacDamV, sacDphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t, fmax_tilt);
			CalcTransferF(sacD, sacDamV, sacDphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
			Pres.z2 = CohAvg( Coh_t, Pha_t, fb_avg, fe_tilt, 0. );
			if(Pres.z2>cohmin) { ApplyCorrection( sacD, sacH1, Coh_t, Adm_t, Pha_t, fmax_tilt ); sacDamV.clear(); }
			// compliance correction on Z (Z -= coh_c*D)
			if(output_coh) { sacV[8] = std::move(Adm_c); sacV[9] = std::move(Pha_c); }
			//Pres.z3 = CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c, fmax_comp);
			CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
			Pres.z3 = CohAvg( Coh_c, Pha_c, fb_avg, fe_comp, 0. );
			if(output_coh) { sacV[10] = Coh_c; sacV[11] = Adm_c; sacV[12] = Pha_c; }
			if(Pres.z3>cohmin) { ApplyCorrection( sacZ, sacD, Coh_c, Adm_c, Pha_c, fmax_comp ); sacZamV.clear(); }
			if(output_coh) { sacV[13] = std::move(Adm_c); sacV[14] = std::move(Pha_c); }
		}
	} else {
		// compliance correction on Z (Z -= coh_c*D)
		ApplyCorrection( sacZ, sacD, Coh_c, Adm_c, Pha_c, fmax_comp ); sacZamV.clear();
		sacZinterm = sacZ;
		if(output_coh) { sacV[8] = std::move(Adm_c); sacV[9] = std::move(Pha_c); }
		// compliance correction on H1 (H1 -= coh_c*D)
		//float coh_cH1 = CalcTransferF(sacH1, sacH1amV, sacH1phV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c, fmax_comp);
		CalcTransferF(sacH1, sacH1amV, sacH1phV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
		float coh_cH1 = CohAvg( Coh_c, Pha_c, fb_avg, fe_comp, 0. );
		if(coh_cH1>cohmin) { ApplyCorrection( sacH1, sacD, Coh_c, Adm_c, Pha_c, fmax_comp ); sacH1amV.clear(); }
		// compliance correction on H2 (H2 -= coh_c*D)
		//float coh_cH2 = CalcTransferF(sacH2, sacH2amV, sacH2phV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c, fmax_comp);
		CalcTransferF(sacH2, sacH2amV, sacH2phV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
		float coh_cH2 = CohAvg( Coh_c, Pha_c, fb_avg, fe_comp, 0. );
		if(coh_cH2>cohmin) { ApplyCorrection( sacH2, sacD, Coh_c, Adm_c, Pha_c, fmax_comp ); sacH2amV.clear(); }
		Pres.z2 = sqrt(coh_cH1*coh_cH1 + coh_cH2*coh_cH2);
		// tilt correction on Z (Z -= coh_t*H1)
		if(output_coh) { sacV[3] = std::move(Adm_t); sacV[4] = std::move(Pha_t); }
		PointC Pres2 = EstimateTiltDirection( true, 1.0 );
		SACRotate(sacH1, sacH2, Pres2.x); sacH1amV.clear(); sacH2amV.clear();
		Pres.x += Pres2.x;	// this is the second rotation on H1&H2
		//Pres.z3 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t, fmax_tilt);
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
		Pres.z3 = CohAvg( Coh_t, Pha_t, fb_avg, fe_tilt, 0. );
		if(output_coh) { sacV[10] = Coh_t; sacV[11] = Adm_t; sacV[12] = Pha_t; }
		if(Pres.z3>cohmin) { ApplyCorrection( sacZ, sacH1, Coh_t, Adm_t, Pha_t, fmax_tilt ); sacZamV.clear(); }
		if(output_coh) { sacV[13] = std::move(Adm_t); sacV[14] = std::move(Pha_t); }
	}

	// rotate H1 and H2 back to their initial azimuth
	SACRotate(sacH1, sacH2, -Pres.x); sacH1amV.clear(); sacH2amV.clear();

	// output all three sets of Coh-Adm-Pha when requested
	if(output_coh) DumpSACs( sacV, outinfoname );

	return Pres;
}

std::vector<std::vector<float>> StaSacs::RayleighDirectionality(const float dazi, 
		const std::vector<std::pair<float, float>>& freqRangeV) const {
	int nrotate = ceil(90./dazi); 
	std::vector<std::vector<float>> resVV(nrotate*4, std::vector<float>(freqRangeV.size()+1));
	SacRec sacH1t(sacH1), sacH2t(sacH2);
	for(int irotate=0; ; irotate++) {
		float azi = irotate * dazi;
		auto &resV1 = resVV[irotate]; resV1[0] = azi;
		auto &resV2 = resVV[irotate+nrotate]; resV2[0] = azi+90.;
		auto &resV3 = resVV[irotate+2*nrotate]; resV3[0] = azi+180.;
		auto &resV4 = resVV[irotate+3*nrotate]; resV4[0] = azi+270.;
		SacRec Coh, Adm, Pha;
		//lambda: compute average coh between fb and fe Hz
		auto CohAvg = [&](float fb, float fe, float &coh_pos, float &coh_neg) {
			coh_pos = coh_neg = 0.;	int npos = 0, nneg = 0;
			auto sigcoh = Coh.sig.get(), sigpha = Pha.sig.get();
			for(int i=Coh.Index(fb); i<Coh.Index(fe); i++) {
				if( sigpha[i] >= 0 ) {
					float weight = 1. - fabs(sigpha[i]/PIo2 - 1.);
					coh_neg += sigcoh[i]*weight; nneg++;
				} else {
					float weight = 1. - fabs(-sigpha[i]/PIo2 - 1.);
					coh_pos += sigcoh[i]*weight; npos++;
				}
			}
			coh_pos = npos>0.1*nneg ? coh_pos/npos : 0.;
			coh_neg = nneg>0.1*npos ? coh_neg/nneg : 0.;
			//coh_pos /= npos+nneg;
			//coh_neg /= npos+nneg;
		};
		auto computeCohs = [&](const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
									  std::vector<float> &resVpos, std::vector<float> &resVneg) {
			CalcTransferF(sacZ, sacZamV, sacZphV, sac2, sac2amV, sac2phV, Coh, Adm, Pha);
			int i=1; for( const auto& pair : freqRangeV ) {
				CohAvg( pair.first, pair.second, resVpos[i], resVneg[i] ); i++;
			}
		};
		computeCohs(sacH1t, sacH1amV, sacH1phV, resV1, resV3); sacH1amV.clear();
		computeCohs(sacH2t, sacH2amV, sacH2phV, resV2, resV4); sacH2amV.clear();
		if(irotate>=nrotate-1) break;
		SACRotate(sacH1t, sacH2t, dazi); 
	}
	return resVV;
}

PointC StaSacs::EstimateTiltDirection( const bool zeropha, const float ddeg ) const {
	// calc the average_coh - azimuth curve (by rotating the two horizontal channels)
	//std::cerr<<"zeropha: "<<zeropha<<std::endl;
	int nrotate = ceil(90./ddeg); 
	std::vector<PointC> cohV(nrotate*2);
	SacRec sacH1t(sacH1), sacH2t(sacH2);
	float fe_tilt = fe_avg<fmax_tilt ? fe_avg : fmax_tilt;
	for(int irotate=0; ; irotate++) {
		auto &p1 = cohV[irotate], &p2 = cohV[irotate+nrotate];
		p1.x = irotate*ddeg; p2.x = p1.x + 90.;
		SacRec Coh, Adm, Pha;
		//p1.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1t, sacH1amV, sacH1phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		//p2.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2t, sacH2amV, sacH2phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH1t, sacH1amV, sacH1phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		p1.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, 0. );
		CalcTransferF(sacZ, sacZamV, sacZphV, sacH2t, sacH2amV, sacH2phV, Coh, Adm, Pha, fmax_tilt, zeropha);
		p2.y = CohAvg( Coh, Pha, fb_avg, fe_tilt, 0. );
		sacH1amV.clear(); sacH2amV.clear();
		//std::cerr<<p1<<"\n"<<p2<<std::endl;
		if(irotate>=nrotate-1) break;
		SACRotate(sacH1t, sacH2t, ddeg); 
	}
	// find the 3 points near the peak
	std::nth_element( cohV.begin(), cohV.begin()+2, cohV.end(), 
			[](const PointC& p1, const PointC& p2){ return p1.y>p2.y; } );
	// return if coh-avg-peak < cohmin
	PointC &pmax = *( std::max_element( cohV.begin(), cohV.begin()+3, 
				[](const PointC& p1, const PointC& p2){ return p1.y<p2.y; } ) );
	//if( pmax.y < cohmin )	return pmax;
	// assume a 180 deg periodicity
	PointC &p1 = cohV[0], &p2 = cohV[1], &p3 = cohV[2];
	p2.x -= 180. * nint((p2.x-p1.x)/180.);
	p3.x -= 180. * nint((p3.x-p1.x)/180.);
	// check if p1 p2 p3 are continguous
	if( fabs(p2.x-p1.x) + fabs(p3.x-p1.x) + fabs(p3.x-p2.x) > ddeg*4.1 ) {
		return pmax;
		//std::stringstream ss;
		//ss<<"non continguous peak points ("<<p1<<" - "<<p2<<" - "<<p3<<")";
		//throw ErrorSR::BadParam(FuncName, ss.str());
	}
	// find the tilt direction (with peak coh)
	return Parabola(p1, p2, p3).Vertex();
}

// compute transfer function for each window (twin) from segments (tseg) within that window
void StaSacs::ReshapeByZHPhaseDiff( float twin, float tseg, float Eperl, float Eperu ) {
	if( twin < tseg ) throw ErrorSR::BadParam( FuncName, "twin has to be >= tseg" );
	_tseg = tseg;
	// zero-out earthquakes and save noise windows, sacZ not modified
	DetectNoiseWindows(Eperl, Eperu);
	// loop through segments
	SacRec Coh1, Adm1, Pha1, Coh2, Adm2, Pha2;
	float twinb_e, twine_e = sacZ.shd.b;
	for( float twinb=sacZ.shd.b; twinb<sacZ.shd.e; twinb+=twin ) {
		bool succ = true;
		try {
			CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh1, Adm1, Pha1, twinb, twinb+twin );
			CalcTransferF(sacZ, sacZamV, sacZphV, sacH2, sacH2amV, sacH2phV, Coh2, Adm2, Pha2, twinb, twinb+twin );
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
							const std::string foutH1, const std::string& foutH2, const std::string foutD ) {
	// lambda function to convert sac back to original type and write
	auto WriteSAC = [&]( SacRec& sac, const std::string& outname ) {
		if( outname.empty() || !sac.sig ) return;
		switch(sactype) {
			case 2:	// acc: differentiate twice
				sac.Differentiate();
			case 1:	// vel: differentiate once
				sac.Differentiate();
			case 0:	// dis: do nothing
				break;
			default:
				throw ErrorSR::BadParam( FuncName, "unknown sactype("+std::to_string(sactype)+")." );
		}
		sac.Write(outname);
	};
	WriteSAC(sacZ, foutZ); WriteSAC(sacZinterm, foutZinterm);
	WriteSAC(sacH1, foutH1); WriteSAC(sacH2, foutH2); WriteSAC(sacD, foutD);
}


/* --------------------- Private Methods --------------------- */

// segmentize a single channel
void StaSacs::Segmentize(const SacRec& sac, std::vector<SacRec>& sacamV, std::vector<SacRec>& sacphV) const {
	if( !sacamV.empty() && !sacphV.empty() && sacamV.size()==sacphV.size() ) return;
	if( ! sac.sig ) throw ErrorSR::EmptySig(FuncName, sac.fname);
	if( _tseg <= 0 ) throw ErrorSR::BadParam(FuncName, "Invalid/Uninitialized _tseg("+std::to_string(_tseg)+")");
	// clear old segments
	sacamV.clear(); sacphV.clear();
	float ttaper = std::min(_tseg*0.25, 200.);
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

void StaSacs::FitIntoParabola( SacRec& sac1, const SacRec& sacsigma, const float fmax, const int nparab ) const {
	// square sigma to get weights
	SacRec sacw(sacsigma); sacw.Mulf(sacsigma);
	// dump data to vector of PointCs
	std::vector<PointC> dataV;
	sac1.Dump( dataV, fb, fmax );
	// check header
	if( dataV.size() != sacw.Index(fmax) - sacw.Index(fb) ) 
		throw std::runtime_error(std::string("Error(")+FuncName+"): npts mismatch");
	auto sigwsac = sacw.sig.get();
	// assign weight into dataV
	// and set logscale to x axis
	int isig = sacw.Index(fb); 
	float wmin = cohmin * cohmin;
	for( auto &p : dataV ) { 
		p.x = std::log(p.x); 
		p.z = sigwsac[isig++]; 
		if( p.z < wmin ) p.z = 0.;
	}
	// divide [fb,fmax) into nparab segments such that each segment has nvalid/nparab valid points
	// count valid points (with weight>0.3)
	std::vector<int> validC(dataV.size());
	validC[0] = dataV[0].z>wmin;
	for(int i=1; i<validC.size(); i++) 
		validC[i] = validC[i-1] + (int)(dataV[i].z>wmin);
	// fit parabola(s)
	float cstep = validC.back() / (float)nparab;
	for(int i=0; i<nparab; i++) {
		auto Ib = dataV.begin() + ( std::lower_bound(validC.begin(), validC.end(), (int)(cstep*i+0.5)) - validC.begin() );
		auto Ie = dataV.begin() + ( std::upper_bound(validC.begin(), validC.end(), (int)(cstep*(i+1)+0.5)) - validC.begin() );
		Parabola parab; float rmsw; parab.Fit(Ib, Ie, rmsw);
		// write parabola to sac1
		float f1 = i==0 ? fb*0.9 : std::exp(Ib->x);
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

void StaSacs::CalcTransferF( const SacRec& sac1, std::vector<SacRec>& sac1amV, std::vector<SacRec>& sac1phV,
										const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
										SacRec& Coh, SacRec& Adm, SacRec& Pha, const float tb, const float te ) const {
	Segmentize( sac1, sac1amV, sac1phV );
	Segmentize( sac2, sac2amV, sac2phV );
	// check consistency of segment Numbers
	int nseg = sac1amV.size();
	if( nseg!=sac1phV.size() || nseg!=sac2amV.size() || nseg!=sac2phV.size() )
		throw std::runtime_error(std::string("Error(")+FuncName+"): No. segments mismatch");
	// only use segments with (shd.b+shd.e)/2 between [tb, te)
	int isegb = 0, isege = nseg;
	float tbeff = sac1.shd.b, teeff = sac1.shd.e;
	if( tb!=NaN && te!=NaN ) {
		isegb = -1;
		for(int iseg=0; iseg<nseg; iseg++) {
			const auto& shd = sac1amV[iseg].shd;
			const float tmid = 0.5 * (shd.user2+shd.user3);
			if( isegb<0 && tmid>=tb ) {
				isegb = iseg; tbeff = shd.user2;
			}
			if( tmid >= te ) { isege = iseg; break; }
			teeff = shd.user3;
		}
		if( isegb==-1 ) { tbeff = teeff; isegb=isege; }
	}
	// minimum No. segments requirement
	if( isege-isegb < 5 ) {
		Coh.shd.user2 = tbeff; Coh.shd.user3 = teeff;
		throw ErrorSR::InsufData(FuncName, "No. segments = "+std::to_string(isege-isegb));
	}
	//nseg -= 1;	// last element is reserved for the original sac

	// calculate autospectral density functions Gss, Grr and
	// the one-sided cross-spectral density function Grs,
	// summing up all segments
	SacRec Gss, Grr, GrsR, GrsI;
	const SacRec& sac0 = sac2amV[0];
	Gss.MutateAs(sac0); Grr.MutateAs(sac0); 
	GrsR.MutateAs(sac0); GrsI.MutateAs(sac0);
	for(int iseg=isegb; iseg<isege; iseg++) {
const auto &shd = sac1amV[iseg].shd;
std::cerr<<"   seg window: "<<shd.user2<<" "<<shd.user3<<std::endl;
		// compute Gss, Grr, Grs_am, Grs_ph for the current seg
		SacRec sacrr_am(sac1amV[iseg]); sacrr_am.Mulf(sac1amV[iseg]);
		SacRec sacss_am(sac2amV[iseg]); sacss_am.Mulf(sac2amV[iseg]);
		SacRec sacrs_am(sac2amV[iseg]); sacrs_am.Mulf(sac1amV[iseg]);
		SacRec sacrs_ph(sac2phV[iseg]); sacrs_ph.Subf(sac1phV[iseg]);
		// Add Gss and Grr to the sum
		Grr.Addf( sacrr_am ); Gss.Addf( sacss_am );
		// project sacrs_am&sacrs_ph onto phase==0.
		/*
		float pha0 = 0.;
		auto sigrsam = sacrs_am.sig.get(), sigrsph = sacrs_ph.sig.get();
		for(int i=0; i<sacrs_am.shd.npts; i++) {
			sigrsam[i] *= cos(sigrsph[i]-pha0);
			sigrsph[i] = pha0;
		}
		*/
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
	Coh.shd.user2 = tbeff; Coh.shd.user3 = teeff;
	// smooth Coh
	Coh.Smooth(0.002, false, fb);
	// re-centralize phase
	float phac = Pha.MeanPha();
	Pha.shd.user1 = phac;
	const float twopi = M_PI*2.;
	Pha.Transform( [&](float& val){
		if(val<phac-M_PI) val+=twopi;
		else if(val>=phac+M_PI) val-=twopi;
		} );
}

float StaSacs::CohAvg(SacRec& Coh, SacRec& Pha, float fb, float fe, float pha0) const {
	int ib = Coh.Index(fb), ie = Coh.Index(fe);
	if( ib >= ie ) throw ErrorSR::BadParam(FuncName, "fb >= fe");
	float mean = 0.;
	auto sigcoh = Coh.sig.get(), sigpha = Pha.sig.get();
	for(int i=ib; i<ie; i++)
		mean += sigcoh[i] * std::max(cos(sigpha[i]-pha0),0.);
	return mean / (ie-ib);
}

void StaSacs::ApplyCorrection( SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha, const float fmax ) {
	#ifdef DEBUG
	auto chn1 = sac1.chname();
	auto chn2 = sac2.chname();
	Adm.Write("debug_Adm_"+chn1+chn2+".SAC");
	Pha.Write("debug_Pha_"+chn1+chn2+".SAC");
	Coh.Write("debug_Coh_"+chn1+chn2+".SAC");
	#endif
	// discard points with phase in the wrong quad
	auto sigpha = Pha.sig.get(), sigcoh = Coh.sig.get();
	const float phac = Pha.shd.user1, pio4 = M_PI*0.25;
	for(int i=0; i<Pha.Index(fmax, true); i++) {
		float phadiff = fabs(sigpha[i]-phac);
		if( phadiff > pio4 ) sigcoh[i] = 0.;
	}
	// fit phase with a single parabola
	FitIntoParabola( Pha, Coh, fmax, 1 );
	// and admittance with two parabolas
	FitIntoParabola( Adm, Coh, fmax, 2 );
	//Adm.Smooth(0.01);
	#ifdef DEBUG
	Coh.Write("debug_CohP_"+chn1+chn2+".SAC");
	Pha.Write("debug_PhaP_"+chn1+chn2+".SAC");
	Adm.Write("debug_AdmP_"+chn1+chn2+".SAC");
	#endif
	// apply transfer function in between 0 and 0.1 Hz
	//const SacRec &sac2_am = sac2amV.back(), &sac2_ph = sac2phV.back();
	SacRec sac2_am, sac2_ph; sac2.ToAmPh(sac2_am, sac2_ph);
	// predicted spectrum
	int npts_ratio = (int)floor(sac2_am.shd.npts/(float)Adm.shd.npts+0.5);
	SacRec AdmInt, PhaInt;
	Adm.Interpolate(npts_ratio, AdmInt); 
	Pha.Interpolate(npts_ratio, PhaInt);
	SacRec sacP_ph(sac2_ph); sacP_ph.Subf(PhaInt); sacP_ph.Wrap(); //sacP_ph.Mul(-1.);
	SacRec sacP_am(sac2_am); sacP_am.Mulf(AdmInt);
	// Tukey window
	sacP_am.cosTaperL(fb*0.9, fb*1.1);
	sacP_am.cosTaperR(fmax*0.9, fmax*1.1);
	SacRec sacP; sacP.shd = sac1.shd;
	sacP.FromAmPh(sacP_am, sacP_ph);
	sac1.Subf(sacP);
	#ifdef DEBUG
	sacP.Write("debug_Pred_"+chn1+chn2+".SAC");
	sac1.Write("debug_Crct_"+chn1+chn2+".SAC");
	#endif
	// clear seg vectors for sac1 (which has been corrected)
	// sac1amV.clear(); sac1phV.clear();
}

#endif
