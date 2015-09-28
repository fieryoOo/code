#include "SeedRec.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

bool isNumber( const std::string& s ) {
	return !s.empty() && std::find_if( s.begin(), s.end(), [](char c) {
		return !std::isdigit(c);
	} ) == s.end();
}

int main( int argc, char* argv[] ) {
	// check inputs
	if( argc != 7 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [seed name] [sta list (or name)] [chanel list (or name)] [event name (14 digits No)] [sps out] [sac outtype (0=displacement, 1=velocity, 2=acceleration)]"<<std::endl;
		exit(-1);
	}

	try {
		// input params
		const std::string seedname(argv[1]);
		const std::string fname_sta(argv[2]);
		const std::string fname_cha(argv[3]);
		const std::string evname(argv[4]);
		if( !isNumber(evname) || evname.length()!=14 )
			throw std::runtime_error( std::string("Invalid input: evname = ") + argv[4] );
		const int sps = atoi(argv[5]);
		if( sps<=0 ) 
			throw std::runtime_error( std::string("Invalid input: sps = ") + argv[5] );
		const int sactype = atoi(argv[6]);

		// other params
		const std::string netname("*");	// *: extract for any network found
		const std::string rec_outname; // empty: do not output

		// lambda: load file
		auto loadF = [&]( const std::string& fname, std::vector<std::string>& V ) {
			std::ifstream fin(fname);
			if( fin ) {
				for( std::string line; std::getline(fin, line); ) {
					std::stringstream ss(line); ss >> line;
					V.push_back( line );
				}
			} else {
				V.push_back( fname );
			}
		};
		// load station list
		std::vector<std::string> staV;
		loadF( fname_sta, staV );
		// load channel list
		std::vector<std::string> chaV;
		loadF( fname_cha, chaV );

		// main loop
		SeedRec seedrec(seedname);
		for( const auto& staname : staV )
			for( const auto& chname : chaV ) {
				// out name
				std::string outname( staname + "." + chname + ".SAC" );
				if( netname != "*" ) outname = netname + "." + outname;
				const std::string resp_outname( outname + "_RESP" );
				// extract
				SacRec sac;
				float gapfrac;
				if( seedrec.ExtractSac(staname, netname, chname, sps, rec_outname, resp_outname, gapfrac, sac) ) {
					sac.Write( outname );
					// remove resp
					float perl = sps * 2.2, perh = 1000;
					sac.RmRESP( resp_outname, perl, perh, sactype );
					sac.ZoomToEvent( evname, -12345., -12345., 0., 5000. );
					sac.Write( "ft_" + outname );
				}
			}
	} catch( std::exception& e ) {
		std::cerr<<e.what()<<std::endl;
	}

	return 0;
}
