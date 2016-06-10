#include "SacRec.h"
#include "SacList.h"
#include "Threads.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <sys/stat.h>

bool SysExists(const std::string& filename) {
   struct stat info; return stat(filename.c_str(), &info)==0;
}

int main(int argc, char *argv[]) {
	if( argc != 3 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [saclist (one sacname per line)] [outdir (has to exist)]"<<std::endl;
		exit(-1);
	}

	std::string outdir(argv[2]);
	if( ! SysExists(outdir) ) {
		std::cerr<<"### main: cannot access output dir "<<outdir<<" ###"<<std::endl;
		exit(-2);
	}

	std::cout<<"### main: grouping sac files by header station names... ###"<<std::endl;
	SacList saclst(argv[1]);
	
	auto corlist = saclst.corList();
	std::cout<<"### main: "<<corlist.size()<<" sac pairs ready to be cross-correlated ###"<<std::endl;

	// lambda function for pre-processings
	auto normalize = [&](SacRec& sac) {
		sac.RunAvg(40., 12., 25.);
		float ttaper = std::min(500., (sac.shd.e-sac.shd.b)*0.25);
		sac.cosTaperL(sac.shd.b, sac.shd.b+ttaper);
		sac.cosTaperR(sac.shd.e-ttaper, sac.shd.e);
		sac.Whiten(0.01, 0.5);
	};

	SacPool sacpool(10.);	// wait for 10 sec before output/clearup (after all sac for a single sacS in the SacPool are Consumed)

	// lambda for CC between two sac files
	auto CrossCorr = [&](const std::string &sacname1, const std::string &sacname2) {
		// load in sac headers
		SacRec sac1(sacname1); sac1.LoadHD();
		SacRec sac2(sacname2); sac2.LoadHD();

		// skip if COR file already exists
		std::string fcorname = outdir+"/COR_"+sac1.stname()+"_"+sac2.stname()+".SAC";
		if( SysExists(fcorname) ) return;

		// notify SacPool of the newly started CC
		sacpool.WaitForSac(fcorname);

		// load in sacs
		sac1.Load(); sac2.Load();

		try {
			// pre-processing
			normalize(sac1); normalize(sac2);

			// do CC
			SacRec sacCCday = sac1.CrossCorrelate(sac2); sacCCday.cut(-3000., 3000.);

			// stack sacCCday
			sacpool.ConsumeSac(fcorname, sacCCday);
		} catch( const ErrorSR::Base &e) {
			std::cerr<<"Warning(main): CC failed on "<<sac1.fname<<" "<<sac2.fname<<" ("<<e.what()<<"). Skipped."<<std::endl;
		}
	};

	// main loop
	Threads threads(12);
	for( const auto &corpair : corlist ) {
		threads.Call(CrossCorr, corpair.first, corpair.second);
	}
	sacpool.Stop();	// notify SacPool theres no more sac coming
	//threads.join();	// wait for all threads to be done (called in dstr)

	return 0;
}
