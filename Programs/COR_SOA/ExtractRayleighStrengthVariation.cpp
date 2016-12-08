// given a list of RDirect (and RAdmittance?), extract the maximum Coherence and associated Admittance info over the given time period

#include "SacRec.h"
#include "RDirect.h"
#include "DisAzi.h"
#include "MyOMP.h"
#include <array>
#include <vector>
#include <string>
#include <fstream>

int Jday ( int y, int m, int d ) {
	int i, jd = 0;
	for( i = 1; i < m; i++ ) {
		if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
		else if (i==2) {
			if ( (y%400==0) || (y%100!=0&&y%4==0) ) jd += 29;
			else jd += 28;
		}
		else jd += 30;
	}
	return jd + d;
}

int AbsDay(int year0, int y, int m, int d) {
	if( year0 > y )
		throw std::runtime_error("Error(AbsDay): year0("+std::to_string(year0)+") > year("+std::to_string(y)+")");
   int nyday = 0;
   for( int i=year0; i<y; i++ ) {
      if ( (i%400==0) || (i%100!=0&&i%4==0) ) nyday += 366;
      else nyday += 365;
   }
	return nyday + Jday(y, m, d);
}

/*
double DayTime(int h, int m, int s, int ms = 0.) { return 3600.*h + 60.*m + s + 0.001*ms; }

double AbsTime(int year0, int y, int mo, int d, int h, int mi, int s, int ms = 0.) {
   return 24.*3600.*AbsDay(year0, y, mo, d) + DayTime(h, mi, s, ms);
}
*/


int main(int argc, char *argv[]) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [infile_list (yyyy mo dd RDirect RAdmittance)]" <<std::endl;
		return -1;
	}

	// load file list
	std::ifstream fin(argv[1]);
	if( ! fin ) {
		std::cerr<<"Error(main): I/O failed on file "<<argv[1]<<std::endl;
		return -2;
	}
	std::vector<std::string> flistV;
	for(std::string line; std::getline(fin, line); )
		flistV.push_back(std::move(line));

	// open fout
	std::string outname(argv[1]); outname += "_MaxCohs";
	std::ofstream fout(outname);

	// loop over flist
	int year0 = 2011;
	bool isFirst = true; std::vector<RDirect::Range> frangeV0;
	for(int iline=0; iline<flistV.size(); iline++) {
		const auto& line = flistV[iline];
		std::stringstream ss(line);
		int y, m, d;
		std::string Cohname, Admname;
		if( ! (ss >> y >> m >> d >> Cohname >> Admname) ) continue;
		std::cout<<"Working on "<<line<<std::endl;

		// time
		double time0 = 24.*3600.*AbsDay(year0, y, m, d);

		// load in rdirect files
		RDirect rdCoh, rdAdm;
		try {
			rdCoh.Load(Cohname); rdAdm.Load(Admname);
		} catch( std::exception &e ) {
			std::cout<<"Warning(main): Skipped. Failed to load RDirect file "<<Cohname<<" "<<Admname<<std::endl;
			continue;
		}

		// check if frange matches
		typedef std::vector<RDirect::Range> RV;
		auto IsSameRange = [](const RV &rV0, const RV &rV1) -> bool {
			if( rV0.size() != rV1.size() ) return false;
			for( int i=0; i<rV0.size(); ++i ) {
				auto &r0 = rV0[i]; auto &r1 = rV1[i];
				if( r0.first!=r1.first || r0.second!=r1.second ) return false;
			}
			return true;
		};
		// between rdCoh and rdAdm
		if( ! IsSameRange(rdCoh.fRanges(), rdAdm.fRanges()) )
			throw std::runtime_error("Error(main): frange mismatch between rdCoh&rdAdm!");
		const auto &frangeV = rdCoh.fRanges();
		// and between rdCoh0 and rdCoh
		if( isFirst ) {
			frangeV0 = frangeV; fout<<"tl tu   \\";
			for( const auto &fr : frangeV0 ) 
				fout<<"frange= "<<fr.first<<" "<<fr.second<<"   ";
			fout<<"\n";	isFirst = false;
		} else {
			if( ! IsSameRange(frangeV0, frangeV) )
				throw std::runtime_error("Error(main): frangeV mismatch between rdCoh&rdCoh0!");
		}

		// loop over time windows in rdirect
		for( const auto &tdPair : rdCoh ) {
			// current time window
			auto twin = tdPair.first;
			double tl = time0 + twin.first, tu = time0 + twin.second;
			// find the associated tdPair of rdAdm
			auto IrdAdm = rdAdm.find(twin);
			if( IrdAdm == rdAdm.end() ) 
				throw std::runtime_error("Error(main): twin mismatch between rdCoh&rdAdm!");
			// store azi-max_coh-adms of all franges in a vector
			std::vector<std::array<float, 3>> cohmaxV(frangeV.size(), {-1., -1., -1.});
			// loop over azi (for cohV) within the current donut
			auto &donutCoh = tdPair.second;
			auto &donutAdm = IrdAdm->second;
			for( const auto &aVPair : donutCoh ) {
				float azi = aVPair.first;
				auto &cohV = aVPair.second;
				auto IadmV = donutAdm.find(azi>=180.?azi-180:azi);
				if( IadmV == donutAdm.end() ) 
					throw std::runtime_error("Error(main): azi mismatch between rdCoh&rdAdm!");
				int icoh=0; for( auto &coh : cohV ) {
					auto &adm = (IadmV->second)[icoh];
					auto &cmaxA = cohmaxV[icoh++];
					if( coh > cmaxA[1] ) cmaxA = {azi, coh, adm};
				}
			}
			fout<<tl<<" "<<tu<<"   ";
			for( const auto &acaA : cohmaxV )
				fout<<acaA[0]<<" "<<acaA[1]<<" "<<acaA[2]<<"   ";
			fout<<y<<" "<<m<<" "<<d<<" "<<twin.first<<" - "<<twin.second<<std::endl;
		}
	}

	std::cout<<outname<<std::endl;

	return 0;
}
