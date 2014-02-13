#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>

int main() {
   char *fnames = new char[100];
   fnames = "aaa bbb ccc ddd";
   //std::stringstream list(fnames);
   std::stringstream list(fnames, std::ios_base::app | std::ios_base::in | std::ios_base::out );
   list << "fieryoOo";
   std::string fname("abc/");
   std::string stmp;

   while( list >> stmp ) {
      fname.append(stmp);
      std::cout<<fname<<std::endl;
   }

   int itmp = 2012319813;
   stmp = "Working_" + std::to_string(itmp);
   std::cerr<<stmp<<std::endl;

   return 0;
}
