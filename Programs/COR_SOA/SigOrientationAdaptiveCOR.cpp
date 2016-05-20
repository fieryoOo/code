#include "SacRec.h"
#include "RDirect.h"
#include "DisAzi.h"
#include "MyOMP.h"
#include <vector>
#include <string>
#include <fstream>

int main(int argc, char *argv[]) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [saclst(sac1 sac2 fRDirect1 fRDirect2)]" <<std::endl;
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
	#pragma omp parallel for schedule (dynamic, 1)
	for(int iline=0; iline<flistV.size(); iline++) {
		const auto& line = flistV[iline];
	//for( const auto& line : flistV ) {
		std::stringstream ss(line);
		std::string sac1name, sac2name, rd1name, rd2name;
		if( ! (ss >> sac1name >> sac2name >> rd1name >> rd2name) ) continue;

		#pragma omp critical
		std::cerr<<sac1name<<" - "<<rd1name<<"   "<<sac2name<<" - "<<rd2name<<std::endl;

		// load in sac and rdirect files
		SacRec sac1(sac1name); sac1.Load();
		SacRec sac2(sac2name); sac2.Load();
		RDirect rdirect1; rdirect1.Load(rd1name);
		RDirect rdirect2; rdirect2.Load(rd2name);

		// normalize sac1 and sac2 before cross-correlation
		sac1.RunAvg(200., 12., 20.); sac1.Whiten(0.015, 0.5);
		sac2.RunAvg(200., 12., 20.); sac2.Whiten(0.015, 0.5);

		// compute azimuths
		Path<float> path(sac1.shd.stlo, sac1.shd.stla, sac2.shd.stlo, sac2.shd.stla);
		float azi1 = path.Azi1()+azi0, azi2 = path.Azi2()+azi0;

		// loop over time windows in rdirect
		const auto &frangeV = rdirect1.fRanges();
		SacRec sacSWd1, sacSWd2;
		for( const auto& tdPair : rdirect1 ) {
			auto twin = tdPair.first;
			if( rdirect2.find(twin) == rdirect2.end() ) continue;
			//std::cout<<twin.first<<" "<<twin.second<<std::endl;
			// lambda: cut, FFT, and apply weightings
			auto FFTWeighted = [&]( const SacRec& sac, const float azi, const RDirect& rdirect, 
					SacRec& sac_am, SacRec& sac_ph, SacRec& sacS, SacRec& sacSWd ) {
				// get spectrum weights
				auto SWeightV = rdirect.SpecWeight(twin, azi);
				// cut out associated sac segments
				float ttaper = std::min(200., (twin.second-twin.first)*0.25);
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
					float fcur = (fwin.first+fwin.second) * 0.5;
					int icur = sac_am.Index( fcur ); fcur = sac_am.X(icur);
					float wcur = SWeightV[ifreq++];
					if( ifreq == 0 ) {
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
			FFTWeighted( sac1, azi1, rdirect1, sac1_am, sac1_ph, sac1S, sacSWd1 );
			SacRec sac2_am, sac2_ph, sac2S, sac2_sw;
			FFTWeighted( sac2, azi2, rdirect2, sac2_am, sac2_ph, sac2S, sacSWd2 );
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
	sacCC1.Write("COR_I03D_I05D_seg.SAC"); 
	std::string oname("COR_I03D_I05D_aziNorm.SAC_"+lag);
	//sacCC2.Whiten(0.015, 0.5); 
	sacCC2.Write(oname); 
	//sacSW1.Write(oname+"_sw1"); sacSW2.Write(oname+"_sw2");

	return 0;
}
