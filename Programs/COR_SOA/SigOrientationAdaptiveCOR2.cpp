#include "SacRec.h"
#include "RDirect.h"
#include "DisAzi.h"
#include "MyOMP.h"
#include <vector>
#include <string>
#include <fstream>

// check if frange matches
typedef std::vector<RDirect::Range> RV;
bool IsSameRange(const RV &rV0, const RV &rV1) {
	if( rV0.size() != rV1.size() ) return false;
	for( int i=0; i<rV0.size(); ++i ) {
		auto &r0 = rV0[i]; auto &r1 = rV1[i];
		if( r0.first!=r1.first || r0.second!=r1.second ) return false;
	}
	return true;
};

int main(int argc, char *argv[]) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [saclst(sac1 sac2 frdCoh1 frdAdm1 frdCoh2 frdAdm2)]" <<std::endl;
		return -1;
	}

	// load from saclist
	std::ifstream fin(argv[1]);
	if( ! fin ) {
		std::cerr<<"Error(main): I/O failed on file "<<argv[1]<<std::endl;
		return -2;
	}
	std::vector<std::string> flistV;
	for(std::string line; std::getline(fin, line); )
		flistV.push_back(std::move(line));

	// loop over flist
	SacRec sacCC1, sacCC2;
	SacRec sacSW1, sacSW2;
	//std::string lag("pos"); float azi0 = 0.;
	std::string lag("neg"); float azi0 = 180.;
	std::string sta1, sta2;
	#pragma omp parallel for schedule (dynamic, 1)
	for(int iline=0; iline<flistV.size(); iline++) {
		const auto& line = flistV[iline];
	//for( const auto& line : flistV ) {
		std::stringstream ss(line);
		std::string sac1name, sac2name, rdCoh1name, rdAdm1name, rdCoh2name, rdAdm2name;
		if( ! (ss >> sac1name >> sac2name >> rdCoh1name >> rdAdm1name >> rdCoh2name >> rdAdm2name) ) continue;

		#pragma omp critical
		std::cout<<sac1name<<" - "<<sac2name<<std::endl;

		// load in sac and rdirect files
		SacRec sac1(sac1name); sac1.Load(); sta1 = sac1.stname();
		SacRec sac2(sac2name); sac2.Load(); sta2 = sac2.stname();
		RDirect rdCoh1, rdAdm1, rdCoh2, rdAdm2; 
		try {
			rdCoh1.Load(rdCoh1name); rdAdm1.Load(rdAdm1name);
			rdCoh2.Load(rdCoh2name); rdAdm2.Load(rdAdm2name);
		} catch( std::exception &e ) {
			std::cerr<<"Warning(main): failed to load rd file(s): "<<e.what()<<std::endl;
			continue;
		}

		// check franges
		const auto &frangeV = rdCoh1.fRanges();
		if( !IsSameRange(frangeV, rdAdm1.fRanges()) || !IsSameRange(frangeV, rdCoh2.fRanges()) || !IsSameRange(frangeV, rdAdm2.fRanges()) )
			throw std::runtime_error("Error(main): franges mismatch between rdCoh1&rdAdm1&rdCoh2&rdAdm2!");

		// normalize sac1 and sac2 before cross-correlation
		sac1.RunAvg(50., 12., 30.); sac1.Whiten(0.015, 0.5);
		sac2.RunAvg(50., 12., 30.); sac2.Whiten(0.015, 0.5);

		// compute azimuths
		Path<float> path(sac1.shd.stlo, sac1.shd.stla, sac2.shd.stlo, sac2.shd.stla);
		float azi1 = path.Azi1()+azi0, azi2 = path.Azi2()+azi0;

		// loop over time windows in rdirect
		SacRec sacSWd1, sacSWd2;
		for( const auto& tdPair : rdCoh1 ) {
			auto twin = tdPair.first;
			if( rdAdm1.find(twin) == rdAdm1.end() ) continue;
			if( rdCoh2.find(twin) == rdCoh2.end() ) continue;
			if( rdAdm2.find(twin) == rdAdm2.end() ) continue;
			//std::cout<<twin.first<<" "<<twin.second<<std::endl;
			// lambda: cut, FFT, and apply weightings
			// only time-freq windows satisfying (Adm_al*Adm_bl^coh < adm < Adm_au*Adm_bu^coh) are used!
			auto FFTWeighted = [&]( const SacRec& sac, const float azi, const RDirect& rdCoh, const RDirect& rdAdm, 
											const float Adm_al, const float Adm_bl, const float Adm_au, const float Adm_bu,
											SacRec& sac_am, SacRec& sac_ph, SacRec& sacS, SacRec& sacSWd ) {
				// get spectrum weights
				auto Cohs = rdCoh.AvgsAlong(twin, azi);
				auto Adms = rdAdm.AvgsAlong(twin, azi);
//std::cerr<<twin.first<<" "<<twin.second<<"\nCohs: ";
//int ifr=0; for( auto coh : Cohs ) std::cerr<<frangeV[ifr++].first<<" "<<coh<<"   "; std::cerr<<sac.stname()<<"\nAdms: ";
//ifr = 0;   for( auto adm : Adms ) std::cerr<<frangeV[ifr++].first<<" "<<adm<<"   "; std::cerr<<sac.stname()<<std::endl;
				// cut out associated sac segments
				float ttaper = std::min(200., (twin.second-twin.first)*0.2);
				sac.cut(twin.first, twin.second, sacS);
				sacS.cosTaperL(twin.first, twin.first+ttaper);
				sacS.cosTaperR(twin.second-ttaper, twin.second);
				// FFT
				sacS.ToAmPh(sac_am, sac_ph);
//SacRec sac_sw = sac_am;
				// weighting amplitude spectrums
				float flast = 0., wlast = 0.; int ilast = 0;
				int ifreq=0; for( auto fwin : frangeV ) {
					//std::cout<<fwin.first<<" "<<fwin.second<<" "<<SWeightV[ifreq++]<<std::endl;
					// get coh and adm
					float coh = Cohs[ifreq], adm = Adms[ifreq++];
					// define weight
					//static const float alpha = -0.5 / (0.13*0.13);
					//float wcur = 0.75-coh; wcur = wcur>0. ? exp(alpha*wcur*wcur) : 1.;
					static const float alpha = -0.5 / (0.3*0.3);
					float wcur = 1.-coh; wcur = exp(alpha*wcur*wcur);
					// only time-freq windows satisfying (Adm_al*Adm_bl^coh < adm < Adm_au*Adm_bu^coh) are used!
//					if( fabs(fwin.first-0.05)<0.001 && fabs(fwin.second-0.0707)<0.001 ) {
						if( adm<Adm_al*pow(Adm_bl,coh) || adm>Adm_au*pow(Adm_bu,coh) ) wcur = 0.;
//					}
					// central freq
					float fcur = (fwin.first+fwin.second) * 0.5;
					int icur = sac_am.Index( fcur ); fcur = sac_am.X(icur);
					if( ifreq == 1 ) {
						sac_am.Transform( [&](float &val){	val *= wcur; }, ilast, icur );
//sac_sw.Transform( [&](float &val){	val = wcur; }, ilast, icur );
					} else {
						float slope = (wcur-wlast)/(fcur-flast);
						sac_am.Transform2( [&](const float f, float &val){
								val *= wlast + slope * (f-flast);
						}, ilast, icur );
//sac_sw.Transform2( [&](const float f, float &val){
//	val = wlast + slope * (f-flast);
//}, ilast, icur );
					}
					ilast = icur; flast = fcur; wlast = wcur;
				}
				sac_am.Transform( [&](float &val){ val *= wlast; }, ilast );
//sac_sw.Transform( [&](float &val){ val = wlast; }, ilast );
//sacSWd.Addf(sac_sw);
				//sacS.FromAmPh(sac_am, sac_ph);
				//sacS.Write("debug_"+sac.stname()+".SAC");
			};
			// CC with weightings
			SacRec sac1_am, sac1_ph, sac1S;
			//FFTWeighted( sac1, azi1, rdCoh1, rdAdm1, 0.1, 12., 1.5, 5., sac1_am, sac1_ph, sac1S, sacSWd1 );
			//FFTWeighted( sac1, azi1, rdCoh1, rdAdm1, 0.005, 8., 0.08, 3., sac1_am, sac1_ph, sac1S, sacSWd1 );
			FFTWeighted( sac1, azi1, rdCoh1, rdAdm1, 0.1, 1., 6., 1., sac1_am, sac1_ph, sac1S, sacSWd1 );
			SacRec sac2_am, sac2_ph, sac2S, sac2_sw;
			FFTWeighted( sac2, azi2, rdCoh2, rdAdm2, 0.1, 1., 6., 1., sac2_am, sac2_ph, sac2S, sacSWd2 );
			// do the cross-correlation
			SacRec sacCC2seg = CrossCorrelateSACs( sac1_am, sac1_ph, sac2_am, sac2_ph, sac1S.shd, sac2S.shd );
			#pragma omp critical
			sacCC2.Addf(sacCC2seg);
			// CC without weightings
			SacRec sacCC1seg = sac1S.CrossCorrelate(sac2S);
			#pragma omp critical
			sacCC1.Addf(sacCC1seg);
		}
/*
#pragma omp critical
{
sacSWd1.Write(sac1.fname+"_sw"); sacSW1.Addf(sacSWd1);
sacSWd2.Write(sac2.fname+"_sw"); sacSW2.Addf(sacSWd2);
}
*/

	}
	std::string oname("COR_"+sta1+"_"+sta2);
	sacCC1.Write(oname+"_seg.SAC"); 
	sacCC2.Whiten(0.02, 0.4); 
	sacCC2.Write(oname+"_aziNorm.SAC_"+lag); 

	return 0;
}
