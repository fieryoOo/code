#include "SacRec.h"
#include "RDirect.h"
#include "DisAzi.h"

int main(int argc, char *argv[]) {
	if( argc != 5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [fRDirect1] [fRDirect2]" <<std::endl;
		return -1;
	}


	// load in sac and rdirect files
	SacRec sac1(argv[1]); sac1.Load();
	SacRec sac2(argv[2]); sac2.Load();
	RDirect rdirect1; rdirect1.Load(argv[3]);
	RDirect rdirect2; rdirect2.Load(argv[4]);

	// compute azimuths
	Path<float> path(sac1.shd.stlo, sac1.shd.stla, sac2.shd.stlo, sac2.shd.stla);
	float azi1 = path.Azi1(), azi2 = path.Azi2();

	// loop over time windows in rdirect
	SacRec sacCC1, sacCC2;
	const auto &frangeV = rdirect1.fRanges();
	for( const auto& tdPair : rdirect1 ) {
		auto twin = tdPair.first;
		if( rdirect2.find(twin) == rdirect2.end() ) continue;
		std::cout<<twin.first<<" "<<twin.second<<std::endl;
		// lambda: cut, FFT, and apply weightings
		auto FFTWeighted = [&]( const SacRec& sac, const float azi, const RDirect& rdirect, 
										SacRec& sac_am, SacRec& sac_ph, SacRec& sacS ) {
			// get spectrum weights
			auto SWeightV = rdirect.SpecWeight(twin, azi);
			// cut out associated sac segments
			sac.cut(twin.first, twin.second, sacS);
			// FFT
			sacS.ToAmPh(sac_am, sac_ph);
			// weighting amplitude spectrums
			float flast = 0., wlast = 0.; int ilast = 0;
			int ifreq=0; for( auto fwin : frangeV ) {
				std::cout<<fwin.first<<" "<<fwin.second<<" "<<SWeightV[ifreq++]<<std::endl;
				float fcur = (fwin.first+fwin.second) * 0.5;
				int icur = sac_am.Index( fcur ); fcur = sac_am.X(icur);
				float wcur = SWeightV[ifreq++];
				if( ifreq == 0 ) {
					sac_am.Transform( [&](float &val){ val *= wcur; }, ilast, icur );
				} else {
					float slope = (wcur-wlast)/(fcur-flast);
					sac_am.Transform2( [&](const float f, float &val){
						val *= wlast + slope * (f-flast);
					}, ilast, icur );
				}
				ilast = icur; flast = fcur; wlast = wcur;
			}
			sac_am.Transform( [&](float &val){ val *= wlast; }, ilast );
//sacS.FromAmPh(sac_am, sac_ph);
//sacS.Write("debug_"+sac.stname()+".SAC");
		};
		// CC with weightings
		SacRec sac1_am, sac1_ph, sac1S;
		FFTWeighted( sac1, azi1, rdirect1, sac1_am, sac1_ph, sac1S );
		SacRec sac2_am, sac2_ph, sac2S;
		FFTWeighted( sac2, azi2, rdirect2, sac2_am, sac2_ph, sac2S );
		// do the cross-correlation
		SacRec sacCCseg2 = CrossCorrelateSACs( sac1_am, sac1_ph, sac2_am, sac2_ph, sac1S.shd, sac2S.shd );
		sacCC2.Addf(sacCCseg2);
		// CC without weightings
		SacRec sacCCseg1 = sac1S.CrossCorrelate(sac2S);
		sacCC1.Addf(sacCCseg1);
	}
	sacCC1.Write("debug_CC1.SAC");
	sacCC2.Write("debug_CC2.SAC");

	return 0;
}
