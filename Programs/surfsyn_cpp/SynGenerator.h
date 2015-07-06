#ifndef SYNGENERATOR_H
#define SYNGENERATOR_H

#include "SacRec.h"
#include "ModelInfo.h"
#include <string>

class fstring : public std::string {
public:
	fstring() : std::string() {}
	fstring( const char* strin ) : std::string(strin) {}
	fstring( const std::string& strin ) : std::string(strin) {}
	~fstring() { clearfstring(); }

	char* f_str() const { return f_str( size()+1 ); }
	char* f_str( const size_t fsize ) const {
		if( fsize < size()+1 )
			throw std::runtime_error("Error(fstring::f_str): input size too small");
		clearfstring();
		pfstring = new char[fsize];
		sprintf(pfstring, "%s", c_str());
		std::fill(pfstring+size(), pfstring+fsize, ' ');
		return pfstring;
	}

	void assignf( const char* strin, const size_t strsize ) {
		int i=strsize-1;
		for(; i>=0; i--)
			if( strin[i] != ' ' ) break;
		this->assign( strin, i+1 );
	}

private:
	mutable char* pfstring = nullptr;
	void clearfstring() const { if(pfstring!=nullptr) delete pfstring; pfstring=nullptr; }
};


class SynGenerator {
public:
	SynGenerator( const fstring& name_fparam, const fstring& name_feigen, const int mode, const float depth );

	bool Synthetic( const float lon, const float lat, const std::string& chname,
						 const float f1, const float f2, const float f3, const float f4, SacRec& sac );

private:
	std::string fmodel;
	ModelInfo minfo;
	static void WriteSACHeader( SAC_HD& shd, const int npts, const float dt, const float elat, const float elon,
										 const std::string& sta, const float slat, const float slon, 
										 const std::string& chn, const std::string& net );
};


#endif
