#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

int main( int argc, char* argv[] ) {
   std::string buff, buff2;
   /* using ifstream directly */
   std::ifstream fin("temp.txt");
   while( fin >> buff >> buff2 ) {
      std::cerr<<buff<<" "<<buff2<<std::endl;
      fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
   }   
   fin.close();

   /* read into a vector */
   fin.open("temp.txt", std::fstream::in);
   //std::vector<std::string> filevec( (std::istream_iterator<std::string>(fin)), (std::istream_iterator<std::string>()) );
   //std::vector<std::string> filevec;
   //copy( (std::istream_iterator<std::string>(fin)), (std::istream_iterator<std::string>()), std::back_inserter< std::vector<std::string> >(filevec) );
   std::vector<std::string> filevec;
   //std::string::iterator striter;
   for( std::string line; std::getline(fin, line); ) {
      filevec.push_back( line.substr(0, line.find_first_of("\t ")) );
   }
   fin.close();
   std::vector<std::string>::iterator fveciter;
   for(fveciter=filevec.begin(); fveciter<filevec.end(); fveciter++) std::cerr<<*fveciter<<std::endl;

   typedef std::numeric_limits<int> limit_int;
   typedef std::numeric_limits<double> limit_double;
   std::cerr<<limit_double::min()<<" "<<limit_double::max()<<" "<<limit_double::infinity()<<std::endl;
   return 0;

}

