#include <iostream>
#include <string>
#include <sstream>

int main() {
	int i1 = -1, i2 = -2;
	//std::string str("12345-332"), tok;
	std::string str, tok;
	std::stringstream ss(str);
	if( std::getline(ss,tok,'-') ) i1 = std::stoi(tok);
	if( std::getline(ss,tok,'-') ) i2 = std::stoi(tok);
	std::cout<<str<<"   "<<i1<<" "<<i2<<std::endl;
}
