#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

int main() {
	std::ifstream fin("temp.txt");
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss;
		std::string s1, s2, s3;
		float f1, f2, f3;
		ss.str(line); bool suc1 = (ss >> s1 >> f2 >> f3); ss.clear();
		std::cerr<<"trial 1: "<<line<<" "<<s1<<" "<<f2<<" "<<f3<<" "<<suc1<<"\n";
		ss.str(line); bool suc2 = (ss >> s1 >> s2 >> s3); ss.clear();
		std::cerr<<"trial 2: "<<line<<" "<<s1<<" "<<s2<<" "<<s3<<" "<<suc2<<"\n";
	}

	return 0;
}
