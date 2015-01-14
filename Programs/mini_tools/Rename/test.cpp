#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

int main() {
	std::stringstream sin("asdacv.asd");
	std::string name1, name2;
	sin >> name1 >> name2;
	std::cerr<<name1.empty()<<" "<<name2.empty()<<std::endl;
	return 0;
}
