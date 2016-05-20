#include "SacRec.h"
#include "MyOMP.h"
#include <vector>
#include <string>
#include <fstream>

int main(int argc, char *argv[]) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [saclst(sac1 sac2)]" <<std::endl;
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
	SacRec sacCC;
	std::string sta1, sta2;
	#pragma omp parallel for schedule (dynamic, 1)
	for(int iline=0; iline<flistV.size(); iline++) {
		const auto& line = flistV[iline];
		std::stringstream ss(line);
		std::string sac1name, sac2name;
		if( ! (ss >> sac1name >> sac2name) ) continue;

		#pragma omp critical
		std::cerr<<sac1name<<"  "<<sac2name<<std::endl;

		// load in sac and rdirect files
		SacRec sac1(sac1name); sac1.Load();
		SacRec sac2(sac2name); sac2.Load();

		// check station names
		#pragma omp critical
		{
		if( sta1.empty() ) sta1 = sac1.stname();
		if( sta2.empty() ) sta2 = sac2.stname();
		}
		if( sta1!=sac1.stname() || sta2!=sac2.stname() )
			throw ErrorSR::HeaderMismatch("main", "inconsistent station name(s)");
		
		// normalize sac1 and sac2 before cross-correlation
		auto normalize = [&](SacRec& sac) {
			sac.RunAvg(40., 15., 30.);
			float ttaper = std::min(500., (sac.shd.e-sac.shd.b)*0.25);
			//sac.cosTaperL(sac.shd.b, sac.shd.b+ttaper);
			//sac.cosTaperR(sac.shd.e-ttaper, sac.shd.e);
			sac.Whiten(0.015, 0.5);
		};
		normalize(sac1); normalize(sac2);

		// do CC
		SacRec sacCCday = sac1.CrossCorrelate(sac2); sacCCday.cut(-3000., 3000.);
		#pragma omp critical
		sacCC.Addf(sacCCday);
	}
	sacCC.Write("COR_"+sta1+"_"+sta2+".SAC"); 

	return 0;
}
