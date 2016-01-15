#include "SacRec.h"
#include "Curve.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc!=3 && argc!=4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sacname] [DISPname] [hfactor (optional, default to 1)]"<<std::endl;
		exit(-1);
	}

	// load sac file
	std::string sacname(argv[1]);
	SacRec sac;
	try {
		sac.Load(sacname);
		if( sac.shd.dist < 0 ) throw std::runtime_error("Error(main): undefined distance in header");
	} catch(const std::exception& e) {
		std::cerr<<"Error(main): "<<e.what()<<std::endl;
		return -2;
	}

	// load dispersion info
	const float dis = sac.shd.dist;
	/*
	std::map<float, float> Tdisp;
	std::ifstream fin(argv[2]);
	if( ! fin )
		throw std::runtime_error(std::string("Error(main): IO failed on ")+argv[2]);
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line);
		float per, grv, ftmp;
		if( ! (ss >> ftmp >> ftmp >> per >> grv) ) continue;
		Tdisp[per] = dis/grv;
	}
	fin.close();
	*/
	std::vector<PointC> dispV;
	std::ifstream fin(argv[2]);
	if( ! fin )
		throw std::runtime_error(std::string("Error(main): IO failed on ")+argv[2]);
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line);
		float per, grv, ftmp;
		if( ! (ss >> ftmp >> ftmp >> per >> grv) ) continue;
		PointC pc; 
		pc.x = per; pc.y = grv; pc.z = dis/grv;
		dispV.push_back(pc);
	}
	fin.close();
	Curve<PointC> Tdisp( std::move(dispV) );

	// main loop for spectrum amplitudes
	const float hfactor = argc==3 ? 1 : atof(argv[3]);
	std::string outname = sacname + "_amps";
	std::ofstream fout(outname);
	//const float hlen = 50., hhlen = 75.;
	const float e = sac.shd.e;
	float amp1 = 0., amp2 = 0., amp3 = 0.;
	for( const auto disp : Tdisp ) {
		//float per = disp.first, Tarr = disp.second;
		float per = disp.x, U = disp.y, Tarr = disp.z;
		// empirical window length (not well tested)
		const float hlen = hfactor * (20. + per*1.5) + dis*0.065;
		const float hhlen = hlen * 1.5;
		float f = 1./per;
		if( Tarr>=0 && Tarr<e ) {
			/* cut in time domain and measure in freq domain */
			SacRec sac2(sac);
			//sac2.cosTaperL( Tarr-hhfactor*per, Tarr-hfactor*per );
			//sac2.cosTaperR( Tarr+hfactor*per, Tarr+hhfactor*per );
			sac2.cosTaperL( Tarr-hhlen, Tarr-hlen );
			sac2.cosTaperR( Tarr+hlen, Tarr+hhlen );
			sac2.ToAm(); amp1 = sac2.Sig(f);

			/* filt in freq domain and measure in time domain */
			sac.GaussianFilt(f, 0.008, sac2);
			//float fc = 6.3/(per*sqrt(dis)), fl = f-fc, fh = f+fc;
			//sac.BandpassBTWFilt(fl, fh, 3, sac2, true);	// zero phase filter
			sac2.Envelope(); 
			amp2 = sac2.Apeak( Tarr-per, Tarr+per );
			// correct amplitude
			//amp3 = amp2 * (per*sqrt( Tdisp.Deriv(per) * dis )) / (2*U) ;
		}
		fout<<per<<" "<<amp1<<" "<<amp2<<" "<<Tarr<<" "<<dis<<"\n";
	}

	return 0;
}
