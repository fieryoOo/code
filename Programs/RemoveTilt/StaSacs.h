#ifndef STASACS_H
#define STASACS_H

#include "SacRec.h"
#include <iostream>
#define DEBUG

class StaSacs {
public:
	StaSacs( const std::string& fnameZ, const std::string& fnameH1, 
				const std::string& fnameH2, const std::string& fnameD="", const int sactypein = 0 )
		: sacZ(fnameZ), sacH1(fnameH1), sacH2(fnameH2), sacD(fnameD), sactype(sactypein) {
		// lambda function to load sac and convert to displacement
		auto LoadSAC = [&]( SacRec& sac ) {
			sac.Load();
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
		LoadSAC(sacZ); LoadSAC(sacH1); LoadSAC(sacH2);
		try { // load pressure gauge data and convert to MPa
			sacD.Load(); sacD.Mul(1.0e-6);
		} catch( ErrorSR::Base& e ) {
			std::cerr<<"Warning(StaSacs): empty/non-accessable sacD "<<sacD.fname<<std::endl;
		}
		if( sacZ.shd.npts!=sacH1.shd.npts || sacZ.shd.npts!=sacH2.shd.npts )
			throw ErrorSR::HeaderMismatch(FuncName, "npts");
		if( sacZ.shd.delta!=sacH1.shd.delta || sacZ.shd.delta!=sacH2.shd.delta )
			throw ErrorSR::HeaderMismatch(FuncName, "delta");
	}

	PointC RemoveTiltCompliance( float Eperl=10., float Eperu=40., float tseg=2500. ) {
		_tseg = tseg;
		// search for earthquakes and save noise windows
		DetectNoiseWindows(Eperl, Eperu);

		// compute transfer F for tilt noise
		// estimate tilt direction (and the corresponding coherence)
		// where Pres.x = direction_t, Pres.y = coh_t
		PointC Pres = EstimateTiltDirection( 1.0 );
		SacRec Coh_t, Adm_t, Pha_t;
		if( Pres.y > cohmin ) {	// re-calc transfer F in the tilt direction
			SACRotate(sacH1, sacH2, Pres.x); sacH1amV.clear();
			Pres.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
		}
		auto &coh_t = Pres.y;

		// compute transfer F for compliance noise
		SacRec Coh_c, Adm_c, Pha_c;
		Pres.z = (!sacD.sig) ? 0. :
					CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
		auto &coh_c = Pres.z;

		std::cerr<<" Initial coh_t & coh_c = "<<coh_t<<" "<<coh_c<<std::endl;

		// and apply to the one with larger average coherence first
		if( coh_t<cohmin && coh_c<cohmin ) {
		} else if( coh_t > coh_c ) {
			// tilt correction on Z (Z -= coh_t*H1)
			ApplyCorrection( sacZ, sacH1, Coh_t, Adm_t, Pha_t ); sacZamV.clear();
			if( sacD.sig ) { // stop if no sacD is empty
				// tilt correction on D (D -= coh_t*H1)
				float coh_tD = CalcTransferF(sacD, sacDamV, sacDphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
				if(coh_tD>cohmin) { ApplyCorrection( sacD, sacH1, Coh_t, Adm_t, Pha_t ); sacDamV.clear(); }
				// compliance correction on Z (Z -= coh_c*D)
				coh_c = CalcTransferF(sacZ, sacZamV, sacZphV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
				if(coh_c>cohmin) { ApplyCorrection( sacZ, sacD, Coh_c, Adm_c, Pha_c ); sacZamV.clear(); }
			}
		} else {
			// compliance correction on Z (Z -= coh_c*D)
			ApplyCorrection( sacZ, sacD, Coh_c, Adm_c, Pha_c ); sacZamV.clear();
			// compliance correction on H1 (H1 -= coh_c*D)
			float coh_cH = CalcTransferF(sacH1, sacH1amV, sacH1phV, sacD, sacDamV, sacDphV, Coh_c, Adm_c, Pha_c);
			if(coh_cH>cohmin) { ApplyCorrection( sacH1, sacD, Coh_c, Adm_c, Pha_c ); sacH1amV.clear(); }
			// tilt correction on Z (Z -= coh_t*H1)
			coh_t = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh_t, Adm_t, Pha_t);
			if(coh_t>cohmin) { ApplyCorrection( sacZ, sacH1, Coh_t, Adm_t, Pha_t ); sacZamV.clear(); }
		}

		return Pres;
	}

	PointC EstimateTiltDirection( const float ddeg ) const {
		// calc the average_coh - azimuth curve (by rotating the two horizontal channels)
		int nrotate = ceil(90./ddeg);
		std::vector<PointC> cohV(nrotate*2);
		SacRec sacH1t(sacH1), sacH2t(sacH2);
		for(int irotate=0; ; irotate++) {
			auto &p1 = cohV[irotate], &p2 = cohV[irotate+nrotate];
			p1.x = irotate*ddeg; p2.x = p1.x + 90.;
			SacRec Coh, Adm, Pha;
			p1.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1t, sacH1amV, sacH1phV, Coh, Adm, Pha);
			p2.y = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2t, sacH2amV, sacH2phV, Coh, Adm, Pha);
			if(irotate>=nrotate-1) break;
			SACRotate(sacH1t, sacH2t, ddeg); 
			sacH1amV.clear(); sacH2amV.clear();
		}
		// find the 3 points near the peak
		std::nth_element( cohV.begin(), cohV.begin()+2, cohV.end(), 
								[](const PointC& p1, const PointC& p2){ return p1.y>p2.y; } );
		// return if coh-avg-peak < cohmin
		PointC &pmax = *( std::max_element( cohV.begin(), cohV.begin()+3, 
								[](const PointC& p1, const PointC& p2){ return p1.y<p2.y; } ) );
		if( pmax.y < cohmin )	return pmax;
		// assume a 180 deg periodicity
		PointC &p1 = cohV[0], &p2 = cohV[1], &p3 = cohV[2];
		p2.x -= 180. * nint((p2.x-p1.x)/180.);
		p3.x -= 180. * nint((p3.x-p1.x)/180.);
		// check if p1 p2 p3 are continguous
		if( fabs(p2.x-p1.x) + fabs(p3.x-p1.x) + fabs(p3.x-p2.x) > ddeg*4.1 ) {
			std::stringstream ss;
			ss<<"non continguous peak points ("<<p1<<" - "<<p2<<" - "<<p3<<")";
			throw ErrorSR::BadParam(FuncName, ss.str());
		}
		// find the tilt direction (with peak coh)
	   return Parabola(p1, p2, p3).Vertex();
	}

	/*
	std::vector<int> RemoveTilts( float Eperl=10., float Eperu=40., float tseg=2500. ) {
		_tseg = tseg;
		DetectNoiseWindows(Eperl, Eperu);
		// a vector to store (a total of 4) nfcorrs (which indicates the level of coherence between each channel)
		// to be returned
		std::vector<int> nfcorrV(4);
		// compute transfer function on both horizontals
		SacRec Coh1, Adm1, Pha1;
		int ncorr1 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh1, Adm1, Pha1);
		nfcorrV[0] = ncorr1;
		SacRec Coh2, Adm2, Pha2;
		int ncorr2 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2, sacH2amV, sacH2phV, Coh2, Adm2, Pha2);
		nfcorrV[1] = ncorr2;
		// and apply to the one with larger correlation coef first
		if( ncorr1 > ncorr2 ) {
			Coh2.clear(); Adm2.clear(); Pha2.clear();
			if(ncorr1>=nfcorr_min) { ApplyCorrection( sacZ, sacH1, Coh1, Adm1, Pha1 ); sacZamV.clear(); }
			ncorr1 = CalcTransferF(sacH2, sacH2amV, sacH2phV, sacH1, sacH1amV, sacH1phV, Coh1, Adm1, Pha1);
			nfcorrV[2] = ncorr1;
			if(ncorr1>=nfcorr_min) { ApplyCorrection( sacH2, sacH1, Coh1, Adm1, Pha1 ); sacH2amV.clear(); }
			ncorr1 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2, sacH2amV, sacH2phV, Coh1, Adm1, Pha1);
			nfcorrV[3] = ncorr1;
			if(ncorr1>=nfcorr_min) { ApplyCorrection( sacZ, sacH2, Coh1, Adm1, Pha1 ); sacZamV.clear(); }
		} else {
			Coh1.clear(); Adm1.clear(); Pha1.clear();
			if(ncorr2>=nfcorr_min) { ApplyCorrection( sacZ, sacH2, Coh2, Adm2, Pha2 ); sacZamV.clear(); }
			ncorr2 = CalcTransferF(sacH1, sacH1amV, sacH1phV, sacH2, sacH2amV, sacH2phV, Coh2, Adm2, Pha2);
			nfcorrV[2] = ncorr2;
			if(ncorr2>=nfcorr_min) { ApplyCorrection( sacH1, sacH2, Coh2, Adm2, Pha2 ); sacH1amV.clear(); }
			ncorr2 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh2, Adm2, Pha2);
			nfcorrV[3] = ncorr2;
			if(ncorr2>=nfcorr_min) { ApplyCorrection( sacZ, sacH1, Coh2, Adm2, Pha2 ); sacZamV.clear(); }
		}
		return nfcorrV;
	}
	*/

	void Write( const std::string& foutZ, const std::string foutH1="",
					const std::string& foutH2="", const std::string foutD="" ) {
		// lambda function to convert sac back to original type and write
		auto WriteSAC = [&]( SacRec& sac, const std::string& outname ) {
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
		if( ! foutZ.empty() ) WriteSAC(sacZ, foutZ);
		if( ! foutH1.empty() ) WriteSAC(sacH1, foutH1);
		if( ! foutH2.empty() ) WriteSAC(sacH2, foutH2);
		if( ! foutD.empty() ) WriteSAC(sacD, foutD);
	}

protected:
	static const int nfcorr_min = 10;
	static constexpr float cohmin = 0.5;
	static constexpr float fmax = 0.1;

private:
	SacRec sacZ, sacH1, sacH2, sacD;
	const int sactype;
	std::vector<int> recb, rece; 
	float _tseg = -1.;
	mutable std::vector<SacRec> sacZamV, sacH1amV, sacH2amV, sacDamV;
	mutable std::vector<SacRec> sacZphV, sacH1phV, sacH2phV, sacDphV;

	int nint( float val ) const { return (int)floor(val+0.5); }
	// segmentize a single channel
	void Segmentize(const SacRec& sac, std::vector<SacRec>& sacamV, std::vector<SacRec>& sacphV) const {
		if( !sacamV.empty() && !sacphV.empty() && sacamV.size()==sacphV.size() ) return;
		// clear old segments
		sacamV.clear(); sacphV.clear();
		for(int i=0; i<recb.size(); i++) {
			float twin_b = sac.X(recb[i]), twin_e = sac.X(rece[i]);
			//std::cout<<twin_b<<" "<<twin_e<<std::endl;
			for(float tb=twin_b; tb<twin_e-_tseg; tb+=_tseg) {
				float te = tb+_tseg;
				SacRec sac_tmp, sac_am, sac_ph;
				// cut and FFT
				sac.cut(tb,te, sac_tmp); 
				sac_tmp.ToAmPh(sac_am, sac_ph);
				sacamV.push_back(std::move(sac_am)); 
				sacphV.push_back(std::move(sac_ph));
			}
		}
	}

	void DetectNoiseWindows( float Eperl=10., float Eperu=40. ) {
		// compute valid time segments (removing earthquakes)
		SacRec sacZcut;
		sacZ.EqkCut(sacZcut, recb, rece, Eperl, Eperu, false);	// no taper
		sacZcut.clear();
	}

	void FitIntoParabola( SacRec& sac1, const SacRec& sacsigma, const float fmax = 0.1 ) const {
		// square sigma to get weights
		SacRec sacw(sacsigma); sacw.Mulf(sacsigma);
		// dump data to vector of PointCs
		std::vector<PointC> dataV;
		sac1.Dump( dataV );
		// assign weight into dataV
		if( dataV.size() != sacw.shd.npts ) 
			throw std::runtime_error(std::string("Error(")+FuncName+"): npts mismatch");
		auto sigwsac = sacw.sig.get();
	   int i;
		for(i=0; i<sacw.Index(fmax); i++)	dataV[i].z = sigwsac[i];
	   for(; i<dataV.size(); i++) dataV[i].z = 0.;
		// fit parabola
	   Parabola parab; float rmsw;
		parab.Fit(dataV, rmsw);
	   // write parabola to sac1
		auto sig1sac = sac1.sig.get();
	   sac1.Transform2( [&](const float x, float& y) {
		   y = parab[x];
	   } );
		//sac1.Write("debug_Adm_BHZBH1.SAC_pb");
	}

	float CalcTransferF( const SacRec& sac1, std::vector<SacRec>& sac1amV, std::vector<SacRec>& sac1phV,
							 const SacRec& sac2, std::vector<SacRec>& sac2amV, std::vector<SacRec>& sac2phV,
							 SacRec& Coh, SacRec& Adm, SacRec& Pha ) const {
		Segmentize( sac1, sac1amV, sac1phV );
		Segmentize( sac2, sac2amV, sac2phV );
		// check No. segments
		int nseg = sac1amV.size();
		if( nseg < 15 )
			throw std::runtime_error(std::string("Error(")+FuncName+","+sac1.fname+"): No. segments("+
											 std::to_string(nseg)+") not enough for computing transfer function");
		if( nseg!=sac1phV.size() || nseg!=sac2amV.size() || nseg!=sac2phV.size() )
			throw std::runtime_error(std::string("Error(")+FuncName+"): No. segments mismatch");
		//nseg -= 1;	// last element is reserved for the original sac

		// calculate autospectral density functions Gss, Grr and
		// the one-sided cross-spectral density function Grs,
		// summing up all segments
		SacRec Gss, Grr, GrsR, GrsI;
		const SacRec& sac0 = sac2amV[0];
		Gss.MutateAs(sac0); Grr.MutateAs(sac0); 
		GrsR.MutateAs(sac0); GrsI.MutateAs(sac0);
		auto sigre = GrsR.sig.get(), sigim = GrsI.sig.get();
		for(int iseg=0; iseg<nseg; iseg++) {
			//sac1amV[iseg]
			// compute Gss, Grr, Grs_am, Grs_ph for the current seg
			SacRec sacrr_am(sac1amV[iseg]); sacrr_am.Mulf(sac1amV[iseg]);
			SacRec sacss_am(sac2amV[iseg]); sacss_am.Mulf(sac2amV[iseg]);
			SacRec sacrs_am(sac2amV[iseg]); sacrs_am.Mulf(sac1amV[iseg]);
			SacRec sacrs_ph(sac2phV[iseg]); sacrs_ph.Subf(sac1phV[iseg]);
			// Add Gss and Grr to the sum
			Grr.Addf( sacrr_am ); Gss.Addf( sacss_am );
			// Add Grs to the sum (converting from am&ph to re&im)
			auto sigam = sacrs_am.sig.get(), sigph = sacrs_ph.sig.get();
			for(int i=0; i<sacrs_ph.shd.npts; i++) {
				float &amp = sigam[i], &pha = sigph[i];
				sigre[i] += amp * cos(pha);
				sigim[i] += amp * sin(pha);
			}
		}
		// convert GrsR&GrsI to GrsAmp&GrsPha in place
		ReImToAmPh( GrsR, GrsI );
		// construct coherence, admittance, and phase from the spectral density functions
		Pha = std::move(GrsI); Adm = GrsR; Adm.Divf( Gss );
		Coh = Adm; Coh.Mulf( GrsR ); Coh.Divf( Grr ); Coh.sqrt();
		// construct smooth transfer function over the whole freq range:
		// smooth Coh
		SacRec Cohs; Coh.Smooth(0.002, Cohs);
		Coh.clear();
		//SacRec& Cohs = Coh;
		// count for points near 0 or pi (with coh>cohmin)
		float pio4 = M_PI*0.25, pi3o4 = M_PI*0.75, twopi = M_PI*2.;
		auto sigpha = Pha.sig.get(), sigcoh = Cohs.sig.get();
		int count0 = 0, countpi = 0;
		for(int i=0; i<Pha.Index(fmax); i++) {
			if( sigcoh[i]<cohmin ) continue;
			float phabs = fabs(sigpha[i]);
			if( phabs<pio4 ) count0++;
			else if(phabs>pi3o4) countpi++;
		}
		int nvalid = std::max(count0, countpi);
		// compute average coh between 0.01 and 0.05 Hz
		float cohavg = 0.;
		int i, ib = Pha.Index(0.01);
		for(i=ib; i<Pha.Index(0.05); i++) cohavg += sigcoh[i];
		cohavg /= (i-ib);
		// return if no enough valid points
		//if( nvalid < nfcorr_min ) return nvalid;
		// theoretically, pha should be either 0 or pi:
		float phac = count0>countpi ? 0. : M_PI;
		// correct 2pi
		if( phac > M_PI*0.5 ) Pha.Transform( [&](float& val){if(val<0.)val+=twopi;} );
		// discard points with phase in the wrong quad
		float phamin = phac-pio4, phamax = phac+pio4;
		for(int i=0; i<Pha.Index(fmax); i++) {
			float phadiff = fabs(sigpha[i]-phac);
			if( phadiff > pio4 ) sigcoh[i] = 0.;
		}
		Coh = std::move(Cohs);
		return cohavg;
	}

	void ApplyCorrection( SacRec& sac1, const SacRec& sac2, SacRec& Coh, SacRec& Adm, SacRec& Pha ) {
		#ifdef DEBUG
		auto chn1 = sac1.chname();
		auto chn2 = sac2.chname();
		Adm.Write("debug_Adm_"+chn1+chn2+".SAC");
		Pha.Write("debug_Pha_"+chn1+chn2+".SAC");
		Coh.Write("debug_Coh_"+chn1+chn2+".SAC");
		#endif
		// fit phase with parabola
		FitIntoParabola( Pha, Coh, fmax );
		// and admittance
		FitIntoParabola( Adm, Coh, fmax );
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
		Adm.Interpolate(npts_ratio); 
		Pha.Interpolate(npts_ratio);
		SacRec sacP_ph(sac2_ph); sacP_ph.Subf(Pha); sacP_ph.Wrap(); //sacP_ph.Mul(-1.);
		SacRec sacP_am(sac2_am); sacP_am.Mulf(Adm);
		sacP_am.cosTaperR(fmax-0.005, fmax+0.005);
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

};

#endif
