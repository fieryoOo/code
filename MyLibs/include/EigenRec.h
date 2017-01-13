#ifndef EIGENREC_H
#define EIGENREC_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <functional>

// single layer (for a single period)
struct DepData{
	float dep, eig, deig;

	DepData( const float dep ) : dep(dep) {}
	DepData( const std::string& line ) {
		if( (sscanf(line.c_str(), "%f %f %f", &dep, &eig, &deig)) != 3 )
			throw std::runtime_error("DepData::DepData: format error: "+line);
	}

	friend std::ostream& operator<< ( std::ostream& o, const struct DepData& dd ) {
		o << dd.dep << " " << dd.eig << " " << dd.deig;
		return o;
	}

	friend bool operator<(const DepData &dp1, const DepData &dp2) { return dp1.dep < dp2.dep;	}
};

// single period all layers
struct PerData {
	float per, phv, grv, wvn, ac;
	float depmin = -1, depmax = -1, ddep = -1;
	int ifline_1, ifline_2, ifline_3;	// _1=H_start, _2=H_end, _2+1=V_start, 3=V_end
	std::string info1, info2;
	std::vector<DepData> depDHV, depDVV;

	friend std::ostream& operator<< ( std::ostream& o, const struct PerData& pd ) {
		o << "per = " << pd.per << " info1 = " << pd.info1 << " info2 = " << pd.info2;
		return o;
	}

	bool sort() {
		// loaded data should already be sorted
		//std::sort(depDHV.begin(), depDHV.end() );
		//std::sort(depDVV.begin(), depDVV.end() );
		int size = depDHV.size();
		if( size < 2 ) return false;
		depmin = depDHV[0].dep; depmax = depDHV.back().dep; ddep = depDHV[1].dep - depmin;
		// check data format
		for(int i=2; i<size; i++)
			if( depDHV[i].dep - depDHV[i-1].dep != ddep )
				throw std::runtime_error("PerData::sort: format error. inconsistent ddep");
		if( depDVV.empty() ) return true;
		if( size != depDVV.size() )
				throw std::runtime_error("PerData::sort: format error. depDHV - depDVV sizes mismatch");
		for(int i=1; i<size; i++)
			if( depDVV[i].dep - depDVV[i-1].dep != ddep )
				throw std::runtime_error("PerData::sort: format error. inconsistent ddep");
	}

	void WriteEigen( const std::string& outname ) {
		const auto &depDV = depDVV.empty() ? depDHV : depDVV;
		if( depDV.empty() ) return;
		std::ofstream fout(outname);
		if( ! fout )
			throw std::runtime_error("PerData::WriteEigen: IO failed on " + outname);
		for( const auto& dd : depDV ) fout << dd << "\n";
	}

	/*
	bool AtDepH( const float dep, float &eig, float &deig ) {
		DepData dd0{dep};
		auto Iu = std::lower_bound( depDHV.begin(), depDHV.end(), dd0 );
		if( Iu == depDHV.end() ) return false;
		if( Iu->dep == dep ) { eig = Iu->eig; deig = Iu->deig; return true; }
		if( Iu == depDHV.begin() ) return false;
		auto Il = Iu - 1;
		float dfactor = (dep-(Il->dep)) / ((Iu->dep)-(Il->dep));
		eig = (Il->eig) + dfactor * ((Iu->eig)-(Il->eig));
		deig = (Il->deig) + dfactor * ((Iu->deig)-(Il->deig));
		return true;
	}
	bool AtDepV( const float dep, float &eig, float &deig ) {
		DepData dd0{dep};
		auto Iu = std::lower_bound( depDVV.begin(), depDVV.end(), dd0 );
		if( Iu == depDVV.end() ) return false;
		if( Iu->dep == dep ) { eig = Iu->eig; deig = Iu->deig; return true; }
		if( Iu == depDVV.begin() ) return false;
		auto Il = Iu - 1;
		float dfactor = (dep-(Il->dep)) / ((Iu->dep)-(Il->dep));
		eig = (Il->eig) + dfactor * ((Iu->eig)-(Il->eig));
		deig = (Il->deig) + dfactor * ((Iu->deig)-(Il->deig));
		return true;
	}
	*/

	inline int nint(const float &val) const { return (int)floor(val+0.5); }
	bool AtDep( float dep, float &eigH, float &deigH ) const {
		if( dep < 0. )
			for(int i=0; i<depDHV.size(); ++i)
				if( depDHV[i].eig!=0 ) {	dep = depDHV[i].dep; break; }
		if( depDHV.empty() || dep<depmin-toler || dep>depmax+toler ) return false;
		int i1dep = nint( (dep-depmin)/ddep ); const auto &dd1 = depDHV[i1dep];
		if( fabs(dd1.dep-dep) < toler ) {
			eigH = dd1.eig; deigH = dd1.deig;
		} else {
			const auto &dd2 = depDHV[dep<dd1.dep ? i1dep-1 : i1dep+1];
			float dfactor = (dep-dd1.dep) / (dd2.dep-dd1.dep);
			eigH = dd1.eig + dfactor * (dd2.eig-dd1.eig);
			deigH = dd1.deig + dfactor * (dd2.deig-dd1.deig);
		}
		return true;
	}

	bool AtDep( float dep, float &eigH, float &deigH, float &eigV, float &deigV ) const {
		if( dep < 0. )
			for(int i=0; i<depDVV.size(); ++i)
				if( depDVV[i].eig!=0 ) {	dep = depDVV[i].dep; break; }
		if( depDHV.empty() || depDVV.empty() || dep<depmin-toler || dep>depmax+toler ) return false;
		int i1dep = nint( (dep-depmin)/ddep ); 
		const auto &dd1H = depDHV[i1dep], &dd1V = depDVV[i1dep];
		if( fabs(dd1H.dep-dep) < toler ) {
			eigH = dd1H.eig; deigH = dd1H.deig;
			eigV = dd1V.eig; deigV = dd1V.deig;
		} else {
			int i2dep = dep<dd1H.dep ? i1dep-1 : i1dep+1;
			const auto &dd2H = depDHV[i2dep], &dd2V = depDVV[i2dep];
			float dfactor = (dep-dd1H.dep) / (dd2H.dep-dd1H.dep);
			eigH = dd1H.eig + dfactor * (dd2H.eig-dd1H.eig);
			deigH = dd1H.deig + dfactor * (dd2H.deig-dd1H.deig);
			eigV = dd1V.eig + dfactor * (dd2V.eig-dd1V.eig);
			deigV = dd1V.deig + dfactor * (dd2V.deig-dd1V.deig);
		}
		return true;
	}

protected:
	static constexpr float toler = 0.001;
};

// data needed for computing source effects at a single source-depth
// formatted to feed fortran subroutines
//rad_pattern_r_( &(stk), &(dip), &(rak), &(er.nper), &(er.dper), er.per, er.eigH, er.deigH, er.eigV, er.deigV, er.wvnR, er.campR,
struct SourceData {
	int nper; float dper;
	std::vector<float> per, wvn, ac;	// period, wavenumber, amp-factor
	std::vector<float> eigH, deigH, eigV, deigV;
};

// single wavetype single mode all periods all layers
class EigenRec {
public:
	std::vector<PerData> perDV;
	SourceData sd;

	//EigenRec( const std::string& fname, const int imod = 1, const bool loadeig = false, std::ostream& sout = std::cout )
	//	: fname(fname), _imod(imod), _sout(sout) {
	EigenRec( const std::string& fname = "", const int imod = 1, const bool loadeig = false )
		: fname(fname), _imod(imod) {
		if(! fname.empty()) preLoad( loadeig );
	}

	void reLoad( const std::string& fnamein, const int imod = 1, const bool loadeig = false ) {
		fname = fnamein; _imod = imod;
		preLoad( loadeig );
	}

	const PerData& GetPD( const float per, const bool load = true ) {
		std::vector<PerData>::iterator found;
		if( ! findPer( per, found ) )
			throw std::runtime_error("EigenRec::GetPD: requested period not found " + std::to_string(per));
		if(load && found->depDHV.empty()) loadPD( *found );
		return *found;
	}

	void FillSD() {
		// Rayl line 1: read(linetmp, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
		// Love line 1: read(linetmp, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
		sd.nper = perDV.size(); sd.dper = perDV[1].per-perDV[0].per;
		sd.per.clear(); sd.per.reserve(sd.nper);
		sd.wvn.clear(); sd.wvn.reserve(sd.nper);
		sd.ac.clear(); sd.ac.reserve(sd.nper);
		for( const auto &pd : perDV ) {
			sd.per.push_back(pd.per);
			sd.wvn.push_back(pd.wvn);
			sd.ac.push_back(pd.ac);
		}
	}

	void FillSDAtDep(const float dep=-1.) {
		sd.nper = perDV.size();
		sd.eigH.assign( sd.nper, 0. ); sd.deigH.assign( sd.nper, 0. ); 
		if( type == "Love" ) { 
			for( int i=0; i<sd.nper; i++ ) perDV[i].AtDep( dep, sd.eigH[i], sd.deigH[i] );
		} else { 
			sd.eigV.assign( sd.nper, 0. ); sd.deigV.assign( sd.nper, 0. ); 
			for( int i=0; i<sd.nper; i++ ) perDV[i].AtDep( dep, sd.eigH[i], sd.deigH[i], sd.eigV[i], sd.deigV[i] );
		}
	}

private:
	std::string fname, type;	// type = Rayleigh/Love
	int _imod = 1, ifline = -1;	// keep track of file line#
	//std::ostream& _sout = std::cout;
	//std::vector<PerData>::iterator iPD = perDV.begin();

	/* Eigen file format:
	line0: Rayleigh/Love mode 1
	line1: @@@...
	line2: per, cl(k),ul(k),wvl(k),ampl(k),qL(k)
	line3: I0, ?, ?, ?
	line4 - line@: dep va(eigen) vd(derivative of eigen)
	*/
	void preLoad( const bool loadeig ) {
		std::ifstream fin(fname); ifline = 0;
		if( ! fin )
			throw std::runtime_error("EigenRec::PreLoad: IO failed on " + fname);

		// find the "Rayleigh/Love mod #mod" line to start with
		for( std::string line; Egetline(fin, line); ) {
			int imod; 
			if( ! isStartLine( line, type, imod ) ) continue;
			if( imod != _imod ) continue;
			break;	// found
		}
		if( fin.eof() )	// starting line not found!
			throw std::runtime_error("EigenRec::PreLoad: starting line for mod#="+std::to_string(_imod)+" not found!");

		// skip/check the first line
		std::string line; Egetline(fin, line);	// skip the first "@@@" line
		if( !isAtLine(line) )
			throw std::runtime_error("EigenRec::PreLoad: format error!");

		// start loop on individual periods
		PerData pd;
		while( preloadPD(fin, pd, loadeig) ) {
			perDV.push_back(pd);
			//_sout<<pd<<std::endl;
		}
		fin.close(); ifline = 0;

		// sort by period
		std::sort(perDV.begin(), perDV.end(), [&](const PerData& pd1, const PerData& pd2) {
			return pd1.per < pd2.per;
		} );
		//_sout<<"### EigenRec::PreLoad: "<<perDV.size()<<" periods pre-loaded for "<<type<<" at mod# "<<_imod<<" ###"<<std::endl;

	}

	bool Egetline( std::ifstream& f, std::string& line ) {
		ifline++;
		return (bool)std::getline(f, line);
	}

	bool findPer( const float per, std::vector<PerData>::iterator& found ) {
		found = std::lower_bound( perDV.begin(), perDV.end(), per, [](const PerData& pd, const float val) {
			return pd.per < val;
		} );
		return ( found!=perDV.end() && per==found->per );
	}

	bool isStartLine( const std::string& line ) const {
		int itmp; std::string stmp;
		return isStartLine( line, stmp, itmp );
	}
	bool isStartLine( const std::string& line, std::string& wtype, int& imod ) const {
		std::stringstream ss(line);
		std::string stmp;
		bool succ = true;
		succ &= !(ss >> wtype >> stmp >> imod).fail();
		succ &= (wtype=="Rayleigh" || wtype=="Love");
		return succ;
	}
	bool isAtLine( const std::string& line ) {
		return line.find("@@@") != std::string::npos;
	}
	bool isDollarLine( const std::string& line ) {
		return line.find("$$$") != std::string::npos;
	}
	bool preloadPD( std::ifstream& fin, struct PerData& pd, const bool loadeig ) {
		auto pdepDV = &(pd.depDHV);
		auto loadline = loadeig ?
							 [&](const std::string &line){ pdepDV->push_back(DepData(line)); } :
							 std::function<void(const std::string&)>( [](const std::string &line){} );
		pd.depDHV.clear(); pd.depDVV.clear();
		pd.ifline_1 = pd.ifline_2 = ifline+3;
		for( std::string line; Egetline(fin, line); ) {
			if( isAtLine(line) || isStartLine(line) ) break;		// end of current period
			if( ifline-pd.ifline_1 == -2 ) pd.info1 = line;			// store first line
			else if( ifline-pd.ifline_1 == -1 ) pd.info2 = line;	// store second line
			else if( isDollarLine(line) ) {								// dollar line found (Rayl)
				pd.ifline_2 = ifline;		// Rayl format: horizontal $$$ vertical
				pdepDV = &(pd.depDVV);		// done with horizontal eigen, start extracting vertical
			} else {
				loadline( line );
			}
		}
		if( pd.ifline_1 == pd.ifline_2 ) pd.ifline_2 = ifline;
		pd.ifline_3 = ifline; pd.sort();
		//std::cerr<<pd.depDHV.size()<<" "<<pd.depDVV.size()<<std::endl;
		bool succ = (ifline - pd.ifline_1 >= 0);
		if( succ ) { std::stringstream ss(pd.info1); ss>>pd.per>>pd.phv>>pd.grv>>pd.wvn>>pd.ac; }
		return succ;
	}
	void loadPD( PerData& pd ) {
		int b = pd.ifline_1, e = pd.ifline_2;
		if( b<0 || e<0 || b>e )
			throw std::runtime_error("EigenRec::loadPD: invalid params: " + std::to_string(b) + " " + std::to_string(e));
		std::ifstream fin(fname); ifline = 0;
		if( ! fin )
			throw std::runtime_error("EigenRec::LoadPD: IO failed on " + fname);

		// get to the b-1 th line
		for( std::string line; Egetline(fin, line); ) {
			if( ifline+1 == b ) break;
		}
		// read Horizontal
		pd.depDHV.clear();
		for( std::string line; Egetline(fin, line)&&ifline<e; ) 
			pd.depDHV.push_back( DepData(line) );
		// read Vertical
		pd.depDVV.clear(); ++ifline;
		for( std::string line; Egetline(fin, line)&&ifline<pd.ifline_3; )
			pd.depDVV.push_back( DepData(line) );
		pd.sort();
		//_sout<<"### EigenRec::loadPD: "<<pd.depDV.size()<<" depths("<<b<<"-"<<e<<") loaded at "<<pd.per<<" sec period. ###"<<std::endl;
	}
};

#endif
