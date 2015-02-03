#include "SysTools.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

//#include <sys/types.h>
#include <sys/stat.h>

class FileList {
public:
	FileList( const std::string path, const std::string prefix = "", 
				 const std::string pattern = "*" ) {
		struct stat st;
		if( stat(path.c_str(), &st) == 0 ) {
			if( st.st_mode & S_IFDIR ) {	//is a directory
				ConstructList( path, pattern );
				AnalyzeList( _list1, _list2, prefix );
				return;
			} else if( st.st_mode & S_IFREG ) {	//is a file
				ReadList(path);
				return;
			}
		}
		throw std::runtime_error("Initialization failed!");
	}
	bool NextFile() { ifile++; }
	bool Rewind() { ifile = 0; }
	bool GetFilename1( std::string& filename ) { return GetFilename(filename, _list1); }
	bool GetFilename2( std::string& filename ) { return GetFilename(filename, _list2); }
	bool GetFilename( std::string& filename, const std::vector<std::string>& list ) {
		if( ifile >= list.size() ) {
			filename.clear();
			return false;
		} else {
			filename = list[ifile];
			return true;
		}
	}

protected:
	static constexpr int epimax = 6000;	//(hope it's big enough..)
private:
	int ifile = 0;
	std::vector<std::string> _list1, _list2;
	std::vector<int> _fno;

	void ConstructList( const std::string& dir, const std::string pattern ) {
		if( ! List( dir.c_str(), pattern.c_str(), 0, _list1 ) )
			throw std::runtime_error("Listing failed!");
	}

	void AnalyzeList( const std::vector<std::string>& list, std::vector<std::string>& listout, const std::string prefix = "" ) {
		listout.resize( list.size() );
		int iname = 0;
		for( auto& fname : list ) {
			// guess episode number
			int fnumber = NoFromFile( fname );
			if( fnumber > 0 ) {
				// get file suffix
				std::size_t lastdot = fname.find_last_of(".");
				std::string suffix;
				if( lastdot != std::string::npos )
					suffix = fname.substr(lastdot);
				listout[iname] = prefix + std::to_string(fnumber) + suffix;
			} else {
				listout[iname] = fname;
			}
			iname++;
		}
	}

	void ReadList( const std::string& fname ) {
		std::ifstream fin(fname);
		if( ! fin )
			throw std::runtime_error("Cannot read from file "+fname);
		for(std::string line; std::getline(fin, line); ) {
			std::stringstream sin(line);
			std::string name1, name2;
			sin >> name1 >> name2;
			if( name1.empty() || name2.empty() ) continue;
			_list1.push_back( name1 );
			_list2.push_back( name2 );
		}
	}

	int NoFromFile( const std::string& name ) {
		//std::string name("asd8990b.~@9x897.12397lasdn98311.23asdno132.rmvb");
		// locate and record all integers into a vector
		std::vector<int> numberV;
		std::vector<std::size_t> nobV, noeV;
		std::size_t nob = 0, noe = 0;
		while( true ) {
			nob = name.find_first_of("0123456789", noe);
			if( nob == std::string::npos ) break;
			noe = name.find_first_not_of("0123456789", nob);
			if( noe == std::string::npos ) noe == name.size();
			numberV.push_back( atoi( name.substr(nob, noe-nob).c_str() ) );
			nobV.push_back(nob); 
			noeV.push_back(noe);
		}
		if(numberV.size() == 0 ) return -12345;
		// invalidate any number pairs that are connected by 
		//	either '.'(fractional) or 'x'(resolution)
		std::vector<bool> valid(numberV.size(), true);
		for(int ino=0; ino<numberV.size()-1; ino++ ) {
			if( noeV[ino]+1 != nobV[ino+1] ) continue;
			int ic = noeV[ino];
			if( name[ic] == '.' || name[ic] == 'x' || name[ic] == 'X' ) {
				valid[ino] = valid[ino+1] = false;
			}
		}
		// return the first remaining number that's 1~epimax 
		for(int ino=0; ino<numberV.size(); ino++ ) {
			if(! valid[ino]) continue;
			int& number = numberV[ino];
			if( number>0 && number<epimax )
				return number;
		}
		return -12345; // no number found
	}
};
