#include <iostream>
#include <sstream>

int main() {
	std::stringstream sin("rdsexe\t/home/tian\t 3.3 aa bb cc dd ee");
	std::string str1, str2, str3;
	float ftmp;

	std::cerr<<"1: "<<sin.str()<<std::endl;
	bool succeed = sin >> str1 >> str2 >> ftmp;

	std::getline(sin, str3);
	std::cerr<<succeed<<"   "<<str1<<" "<<str2<<" "<<ftmp<<std::endl;
	std::cerr<<str3<<std::endl;

	std::cerr<<"2: "<<sin.str()<<std::endl;

	return 0;
}
