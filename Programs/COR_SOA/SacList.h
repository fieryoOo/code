#ifndef SACLIST_H
#define SACLIST_H

#include "SacRec.h"
#include "MyOMP.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <chrono>

class SacList {
public:
	SacList(const std::string& flist) {
		// load in all sac names and group by station name
		std::ifstream fin(flist);
		for(std::string line; std::getline(fin, line);) {
			std::stringstream ss(line); 
			if( ! (ss >> line) ) continue;
			try {
				SacRec sac(line); sac.LoadHD();
				int absJday = sac.AbsDay(2000); std::string name = sac.stname();
				std::transform(name.begin(),name.end(),name.begin(),::toupper);
				_saclistH[name].emplace( absJday, line );
			} catch( ErrorSR::Base &e ) {
				std::cerr<<"Warning(main) sac header error in "<<line<<" ("<<e.what()<<"). skipped"<<std::endl;
			}
		}

		// list of station-names
		int nsacs = 0;
		for(auto &sVpair : _saclistH) {
			_stalistV.push_back(sVpair.first);
			nsacs += sVpair.second.size();
		}
		std::cout<<"SacList::SacList: "<<nsacs<<" sacs from "<<_saclistH.size()<<" stations loaded"<<std::endl;
		std::sort(_stalistV.begin(), _stalistV.end());
	}

	const std::vector<std::string>& staList() const { return _stalistV; }

	//void corList(const std::string &staname1, const std::string &staname2) {
	std::vector<std::pair<std::string,std::string>> corList(const std::string &staname1, const std::string &staname2) {
		auto IH1 = _saclistH.find(staname1), IH2 = _saclistH.find(staname2);
		if( IH1==_saclistH.end() || IH2==_saclistH.end() )
			throw std::runtime_error(std::string("SacList::corList: station ")+staname1+" or "+staname2+" not found");
		std::vector<std::pair<std::string,std::string>> corlist;
		auto IM1 = IH1->second.begin(), IM2 = IH2->second.begin();
		while( IM1!=IH1->second.end() && IM2!=IH2->second.end() ) {
			auto ajday1 = IM1->first, ajday2 = IM2->first;
			if( ajday1 < ajday2 ) {
				IM1++;
			} else if( ajday1 > ajday2 ) {
				IM2++;
			} else {
				corlist.push_back({IM1->second, IM2->second});
				IM1++; IM2++;
			}
		}
		return corlist;
	}

	std::vector<std::pair<std::string,std::string>> corList(const std::string &staname1) {
		std::vector<std::pair<std::string,std::string>> corlist;
		auto Isl = std::lower_bound(_stalistV.begin(), _stalistV.end(), staname1);
		for( ; Isl!=_stalistV.end(); Isl++ ) {
			auto curlst = corList(staname1, *Isl);
			corlist.insert(corlist.end(), std::make_move_iterator(curlst.begin()), std::make_move_iterator(curlst.end()));
		}
		return corlist;
	}

	std::vector<std::pair<std::string,std::string>> corList() {
		std::vector<std::pair<std::string,std::string>> corlist;
		for( const auto &staname : _stalistV ) {
			auto curlst = corList(staname);
			corlist.insert(corlist.end(), std::make_move_iterator(curlst.begin()), std::make_move_iterator(curlst.end()));
		}
		return corlist;
	}

private:
	// hash table that stores one <jday,sacname> map for each station
	std::unordered_map<std::string, std::map<int, std::string>> _saclistH;
	std::vector<std::string> _stalistV;
};



class SacPool {
public:
	// waittime: time to wait after all sac for a single sacS are Consumed and before output/clearup
	SacPool(const float waittime = 10.) 
		: t(&SacPool::SACManager, this), tcheck1(waittime*1500), tcheck2(waittime*100) {}

	~SacPool() { t.join(); }

	void Stop() { isWaiting = false; }

	void WaitForSac( const std::string &sacname ) {
		std::lock_guard<std::mutex> lock(mhash);
		auto ISC = _SacCountH.find(sacname); 
		if( ISC == _SacCountH.end() )	{
			_SacCountH[sacname].shd.user5 = 1;
			_SacCountH[sacname].fname = sacname;
		} else { (ISC->second).shd.user5 += 1; }
		//std::cerr<<"WaitForSac: nwait("<<sacname<<") = "<<_SacCountH[sacname].shd.user5<<std::endl;
	}

	void ConsumeSac( const std::string& sacname, const SacRec& sac ) {
		std::lock_guard<std::mutex> lock(mhash);
		auto ISC = _SacCountH.find(sacname); 
		if( ISC == _SacCountH.end() )
			throw std::runtime_error(std::string("Error(")+__FUNCTION__+"): "+sacname+" not exist in pool");
		if(sac.sig) {
			auto &sacS = ISC->second;
			int count = sacS.shd.user5;
			sacS.Addf(sac);
			sacS.shd.user5 = count; sacS.fname = sacname;
		}
		(ISC->second).shd.user5-=1; 
		//std::cerr<<"Consumed sac: nwait("<<sacname<<") = "<<(ISC->second).shd.user5<<std::endl;
	}

private:
	std::thread t;
	std::mutex mhash;
	bool isWaiting = true;
	std::unordered_map<std::string, SacRec> _SacCountH;
	int tcheck1, tcheck2;	// in milliseconds

	void SACManager() {
		while( isWaiting || !_SacCountH.empty() ) {	// check _SacCountH every tcheck1 seconds
			std::this_thread::sleep_for(std::chrono::milliseconds(tcheck1));
			//std::cout<<" ... "<<_SacCountH.size()<<std::endl;
			std::lock_guard<std::mutex> lock(mhash);
			for(const auto SC : _SacCountH) {
				int count = (SC.second).shd.user5;
				if( count < 0 ) {
					throw std::runtime_error(std::string("Error(")+__FUNCTION__+"): negative count on sac "+SC.first);
				} else if( count == 0 ) {
					std::thread(&SacPool::CheckNOutput, this, SC.first).detach();
				}
			}
		}
	}

	void CheckNOutput( const std::string& sacname ) {
		// repeat the check ten times (tcheck2) before clear and Output
		for(int i=0; i<10; i++) {
			std::this_thread::sleep_for(std::chrono::milliseconds(tcheck2));
			std::lock_guard<std::mutex> lock(mhash);
			auto ISC = _SacCountH.find(sacname);
			if(ISC==_SacCountH.end() || (ISC->second).shd.user5!=0) return;
		}
		// test passed. clear and output
		std::lock_guard<std::mutex> lock(mhash);
		//std::cerr<<"CheckNOutput 2: "<<sacname<<std::endl;
		auto ISC = _SacCountH.find(sacname);
		if( (ISC->second).sig ) (ISC->second).Write();
		_SacCountH.erase(ISC);
		std::cout<<"*** Stacked sac "<<sacname<<" output and cleared. ***"<<std::endl;
	}
};

/*
int main() {
	Rand rand(1,10);
	SacPool sacpool(10.0);
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i=0; i<30; i++) {
		SacRec sac("/lustre/janus_scratch/yeti4009/ASN_OBS/SAC/COR_J23A_I03D.SAC");
		auto sacname = "COR_"+std::to_string(rand.UniformI())+".SAC";
		std::cout<<i<<" ";
		sacpool.WaitForSac(sacname);
		sac.Load(); 
		std::this_thread::sleep_for(std::chrono::milliseconds((int)(rand.UniformI()*1000)));
		sacpool.ConsumeSac(sacname, sac);
	}
	sacpool.Stop();

	return 0;
}
*/

#endif
