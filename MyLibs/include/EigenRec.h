#ifndef EIGENREC_H
#define EIGENREC_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

struct DepData{
	float dep, eig, deig;
	DepData( const std::string& line ) {
		if( (sscanf(line.c_str(), "%f %f %f", &dep, &eig, &deig)) != 3 )
			throw std::runtime_error("DepData::DepData: format error: "+line);
	}
	friend std::ostream& operator<< ( std::ostream& o, const struct DepData& dd ) {
		o << dd.dep << " " << dd.eig << " " << dd.deig;
		return o;
	}
};
struct PerData {
	float per;
	int ifline_b, ifline_e;
	std::string info1, info2;
	std::vector<DepData> depDV;

	void WriteEigen( const std::string& outname ) {
		if( depDV.size() == 0 ) return;
		std::ofstream fout(outname);
		if( ! fout )
			throw std::runtime_error("EigenRec::WriteEigen: IO failed on " + outname);
		for( const auto& dd : depDV ) fout << dd << "\n";
	}
	friend std::ostream& operator<< ( std::ostream& o, const struct PerData& pd ) {
		o << "per = " << pd.per << " info1 = " << pd.info1 << " info2 = " << pd.info2;
		return o;
	}
};

class EigenRec {
public:
	EigenRec( const std::string& fname, const int imod = 1, std::ostream& sout = std::cout )
		: fname(fname), _imod(imod), _sout(sout) {
		preLoad();
	}

	const PerData& GetPD( const float per, const bool load = true ) {
		std::vector<PerData>::iterator found;
		if( ! findPer( per, found ) )
			throw std::runtime_error("EigenRec::GetPD: requested period not found " + std::to_string(per));
		if(load) loadPD( *found );
		return *found;
	}

protected:

private:
	const std::string fname;
	std::string type;	// type = Rayleigh/Love
	const int _imod = 1;
	std::ostream& _sout = std::cout;
	int ifline = -1;	// keep track of file line#
	std::vector<PerData> perDV;
	//std::vector<PerData>::iterator iPD = perDV.begin();

	/* Eigen file format:
	line0: Rayleigh/Love mode 1
	line1: @@@...
	line2: per, cl(k),ul(k),wvl(k),ampl(k),qL(k)
	line3: I0, ?, ?, ?
	line4 - line@: dep va(eigen) vd(derivative of eigen)
	*/
	void preLoad() {
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
		while( preloadPD(fin, pd) ) {
			perDV.push_back(pd);
			//_sout<<pd<<std::endl;
		}
		fin.close(); ifline = 0;

		// sort by period
		std::sort(perDV.begin(), perDV.end(), [&](const PerData& pd1, const PerData& pd2) {
			return pd1.per < pd2.per;
		} );
		_sout<<"### EigenRec::PreLoad: "<<perDV.size()<<" periods pre-loaded for "<<type<<" at mod# "<<_imod<<" ###"<<std::endl;

	}

	bool Egetline( std::ifstream& f, std::string& line ) {
		ifline++;
		return std::getline(f, line);
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
	bool preloadPD( std::ifstream& fin, struct PerData& pd ) {
		pd.ifline_b = pd.ifline_e = ifline+3;
		for( std::string line; Egetline(fin, line); ) {
			if( isAtLine(line) || isStartLine(line) ) break;
			if( ifline-pd.ifline_b == -2 ) pd.info1 = line;
			else if( ifline-pd.ifline_b == -1 ) pd.info2 = line;
			else if( isDollarLine(line) ) pd.ifline_b = ifline+1;
		}
		pd.ifline_e = ifline;
		bool succ = (pd.ifline_e - pd.ifline_b >= 0);
		if( succ ) sscanf( pd.info1.c_str(), "%f", &(pd.per) );
		return succ;
	}
	void loadPD( PerData& pd ) {
		int b = pd.ifline_b, e = pd.ifline_e;
		if( b<0 || e<0 || b>e )
			throw std::runtime_error("EigenRec::loadPD: invalid params: " + std::to_string(b) + " " + std::to_string(e));
		std::ifstream fin(fname); ifline = 0;
		if( ! fin )
			throw std::runtime_error("EigenRec::LoadPD: IO failed on " + fname);

		// get to the b-1 th line
		for( std::string line; Egetline(fin, line); ) {
			if( ifline+1 == b ) break;
		}
		pd.depDV.clear();
		for( std::string line; Egetline(fin, line); ) {
			if( ifline == e ) break;
			pd.depDV.push_back( DepData(line) );
		}
		_sout<<"### EigenRec::loadPD: "<<pd.depDV.size()<<" depths("<<b<<"-"<<e<<") loaded at "<<pd.per<<" sec period. ###"<<std::endl;
	}
};

#endif
