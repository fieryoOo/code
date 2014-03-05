#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

int main(int argc, char *argv[]) {
   if( argc != 4 ) {
      std::cerr<<"usage: File_Combiner [in_file_lst] [column # (0 for entire row; -1 for non-comp join)] [out_file]"<<std::endl;
      exit(-1);
   }

   // check colunm number
   if( atoi(argv[2]) != atof(argv[2]) || atoi(argv[2]) < -1 ) {
      std::cerr<<"Integer(>=-1) expected for column #"<<std::endl;
      exit(0);
   }
   int kcol = atoi(argv[2]);

   // read in file list
   std::ifstream fin(argv[1]);
   if( ! fin ) {
      std::cerr<<"Cannot read from file "<<argv[1]<<std::endl;
      exit(0);
   }
   std::vector<std::string> filelist;
   for( std::string line; std::getline(fin, line); ) {
      std::stringstream sin(line);
      std::string fname;
      sin >> fname;
      filelist.push_back(fname);
   }
   fin.close();

   // read in from each file
   std::vector< std::vector<std::string> > stringlists;
   stringlists.resize( filelist.size() );
   for( int ifile=0; ifile<filelist.size(); ifile++ ) {
      fin.clear(); fin.open( filelist.at(ifile).c_str() );
      if( ! fin ) {
	 std::cerr<<"Warning: Cannot read from file "<<filelist.at(ifile)<<". Skipped!"<<std::endl;
	 continue;
      }
      for(std::string line; std::getline(fin, line); ) {
	 stringlists.at(ifile).push_back(line); 
      }
      fin.close();
   }

   // iterate through each slist and search for matches 
   for( int ilist=0; ilist<stringlists.size(); ilist++ ) {		// for each list
      std::vector<std::string>& slist = stringlists.at(ilist);
      for( int is=0; is<slist.size(); is++ ) {				// for each string
	 std::string str1;
	 // string 1 for comparison
	 if( kcol == 0 ) str1 = slist.at(is);
	 else {
	    std::stringstream sin1( slist.at(is) );
	    int itok;
	    for( itok=0; itok<kcol; itok++ ) { if( ! (sin1 >> str1) ) break; }
	    if( itok != kcol ) continue;
	 }
	 for( int jlist=ilist+1; jlist<stringlists.size(); jlist++ ) {	// for each of the other list
	    std::vector<std::string>& slist2 = stringlists.at(jlist);
	    for( int is2=0; is2<slist2.size(); is2++ ) {		// for each string2 in list2
	       // string 2 for comparison
	       std::string str2;
	       if( kcol == 0 ) str2 = slist2.at(is2);
	       else {
		  std::stringstream sin2( slist2.at(is2) );
		  int itok;
		  for(itok=0; itok<kcol; itok++ ) { if( ! (sin2 >> str2) ) break; }
		  if( itok != kcol ) continue;
	       }
	       // compare string1 and string2
	       if( str1 == str2 ) {
		  slist.at(is) += "\t" + slist2.at(is2);
		  slist2.erase( slist2.begin() + is2 );
	       }
	    }
	 }
      }
   }

   // output
   std::ofstream fout(argv[3]);
   for( int ilist=0; ilist<stringlists.size(); ilist++ ) {
      std::vector<std::string>& slist = stringlists.at(ilist);
      for( int is=0; is<slist.size(); is++ ) {
	 fout<<slist.at(is)<<std::endl;
      }
   }
   fout.close();
   
}
