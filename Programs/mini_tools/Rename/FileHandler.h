#include "SysTools.h"
#include <string>

class FileHandler : std::string {
public:
	FileHandler() : std::string() {}
	FileHandler( const std::string& stringin ) : std::string(stringin) {}
	bool Rename( const std::string targetname, std::string& actualnewname ) {
		actualnewname = targetname;
		Move( (*this).c_str(), actualnewname.c_str() );
	}
};

