/* ------------------------------------------------------------------------
 * cpp version of grep
 * takes either standard input or a file (with -f filein)
 * -w, --whole: match every single character of the whole line to the input string (same as ^string1$ ^string2$ ...)
 * -a, --adjacent: find a sequence of lines in the input contents, such that each of the line in the sequence, in the original order, corresponds to one of the input string
 * -v, --invert: invert the sense of matching, to select non-matching lines
 * '^' be the 1st character of the searching string: have to match from the begining
 * '$' be the last character of the searching string: have to match till the end
 * ------------------------------------------------------------------------ */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>


class Greper {
public:
	/*
	Greper( const std::istream& inputsin,
			  std::vector< std::string >&& triggersin,
			  const bool MatchAllin, const bool MatchAdjin )
		: inputs(inputsin), triggers(std::move(triggersin)),
		  MatchAll(MatchAllin), MatchAdj(MatchAdjin) {
	}
	*/
	Greper( std::vector<std::string>&& contentsin )
		: contents( std::move(contentsin) ) {
		// remove coloring characters
		for(auto itc=contents.begin(); itc<contents.end(); ) {
			auto line = *itc;
			std::size_t cbeg = 0, cend;
			while( true ) {
				cbeg = line.find("\e[");
				cend = line.find("m", cbeg+2);
				if( cbeg==std::string::npos || cend==std::string::npos ) break;
				cbeg = line.erase( line.begin()+cbeg, line.begin()+cend+1) - line.begin();
			}
			if( line.empty() ) {
				itc = contents.erase(itc);
			} else {
				contents_m.push_back(line);
				itc++;
			}
		}
	}

	void Grep( const std::vector<std::string>::iterator ibeg_trig,
				  const std::vector<std::string>::iterator iend_trig,
				  std::vector<std::string>& outputs,
				  const bool MatchAll = false, const bool MatchAdj = false, const bool InvMatch = false) const {
		if( contents_m.empty() || ibeg_trig>=iend_trig ) return;
		if( MatchAdj )
			grepAdjacent(ibeg_trig, iend_trig, outputs, MatchAll, InvMatch);
		else
			grepDisjoint(ibeg_trig, iend_trig, outputs, MatchAll, InvMatch);
	}
	
private:
	std::vector< std::string > contents;
	std::vector< std::string > contents_m;

	// check if 'line' is the same as 'trigger'
	static bool lineIs(std::string const& line, std::string trigger) {
		if( !trigger.empty() && trigger[0]=='^' ) trigger.erase(trigger.begin());
		auto itlast = trigger.end()-1;
		if( !trigger.empty() && *(itlast)=='$' ) trigger.erase(itlast);
		//std::cerr<<"checking: "<<line<<" "<<trigger<<std::endl;
		//if( line==trigger ) std::cerr<<"true!"<<std::endl;
		//else std::cerr<<"false!"<<std::endl;
		return line == trigger;
	}

	// check if 'line' contains 'trigger'
	static bool lineContains(std::string const& line, std::string trigger) {
		if ( trigger.empty() ) return true;
		if ('^' == trigger[0]) {
			trigger.erase(trigger.begin());
			if (!trigger.empty() && '$' == trigger[trigger.size() - 1])
				return lineIs(line, trigger);
			return 0 == line.find(trigger);
		}
		if ('$' == trigger[trigger.size() - 1]) {
			trigger.erase(trigger.end() - 1);
			std::string::size_type const pos = line.rfind(trigger);
			return std::string::npos != pos && trigger.size() + pos == line.size();
		}
		return std::string::npos != line.find(trigger);
	}

	// search for triggers all together 
	void grepAdjacent(const std::vector<std::string>::iterator ibeg_trig,
							const std::vector<std::string>::iterator iend_trig,
							std::vector<std::string>& outputs, 
							bool matchWholeLines, bool InvertMatch) const {
		std::size_t trigsize = iend_trig - ibeg_trig;
		auto itc = contents_m.begin();	// start from the begining of contents_m
		while( itc < contents_m.end() ) {// and stop at the end
			// search for the next match of the whole trig sequence
			auto itsearch = std::search(itc, contents_m.end(), ibeg_trig, iend_trig, matchWholeLines?lineIs:lineContains);
			if( InvertMatch ) {	// if invert matching
				// output all contents before the matching point
				for(; itc<itsearch; itc++)
					outputs.push_back( *(GetIcontent(itc)) );
				itc += trigsize;	// and skip the matches
			} else {	// otherwise
				// output all matches
				if(itsearch<contents_m.end()) {
					for(itc=itsearch; itc<itsearch+trigsize; itc++)
						outputs.push_back( *(GetIcontent(itc)) );
				} else {
					itc = itsearch;
				}
			}
		}
	}

	// search for each trigger independentlly
	void grepDisjoint(const std::vector<std::string>::iterator ibeg_trig,
                     const std::vector<std::string>::iterator iend_trig,
                     std::vector<std::string>& outputs, 
							bool matchWholeLines, bool InvertMatch) const {
		auto itc = contents_m.begin();	// start from the begining of contents_m
		while( itc < contents_m.end() ) {// and stop at the end
			// search for the next match of any trigger
			auto itsearch = std::find_first_of(itc, contents_m.end(), ibeg_trig, iend_trig, matchWholeLines?lineIs:lineContains);
			if( InvertMatch ) {	// if invert matching
				// output all contents before the match
            for(; itc<itsearch; itc++)
               outputs.push_back( *(GetIcontent(itc)) );
            itc++;  // and skip the match
			} else { //	otherwise
				if(itsearch<contents_m.end()) outputs.push_back( *(GetIcontent(itsearch)) );
				itc = itsearch + 1;
			}
		}
	}

	// return an iterator to the corresponding content of the current content_m
	std::vector<std::string>::const_iterator GetIcontent(const std::vector<std::string>::const_iterator icm) const {
		std::vector<std::string>::const_iterator ic = contents.end();
		if( icm>=contents_m.begin() && icm<contents_m.end() )
			ic = contents.begin() + ( icm - contents_m.begin() );
		return ic;
	}
};


int main( int argc, char const* argv[] ) {
	// prompt for inputs
	if( argc < 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [options(-w -a -v -f)] [strings]"<<std::endl;
		exit(-1);
	}
	// record options/fetch strings to search
	bool matchWholeLines = false; 
	bool useAdjacentMatches = false;
	bool InvertMatch = false;
	std::string filename;
	std::vector<std::string> triggers;
	int iarg=1;
	while( iarg<argc ) {
		std::string argstr(argv[iarg]);
		if( argstr == "-w" ) {
			matchWholeLines = true;
		} else if( argstr == "-a" ) {
			useAdjacentMatches = true;
		} else if( argstr == "-v" ) {
			InvertMatch = true;
		} else if( argstr == "-f" ) {
			filename = argv[++iarg];
		} else {
			triggers.push_back(argstr);
		}
		iarg++;
	}
	// set input stream to either fin (option -f) or std::cin
	std::ifstream fin;
	if( ! filename.empty() ) {
		fin.open(filename);
		if( ! fin ) {
			std::cerr<<"\e[1;31mError(main): cannot read from file "<<filename<<"\e[0m"<<std::endl;
			exit(-2);
		}
	}
	std::istream& contents_in = filename.empty() ? std::cin : fin;
	// fetch original contents from standard_input or file_in
	std::vector<std::string> contents;
	for( std::string line; std::getline(contents_in, line); ) {
		contents.push_back(line);
	}
	// initialize the greper with the contents vector
	Greper greper( std::move(contents) );
	// grep and store outputs into a new vector
	std::vector<std::string> outputs;
	greper.Grep( triggers.begin(), triggers.end(), outputs, 
					 matchWholeLines, useAdjacentMatches, InvertMatch );
	// output
	for(auto line : outputs) std::cout<<line<<std::endl;
	return 0;
}
