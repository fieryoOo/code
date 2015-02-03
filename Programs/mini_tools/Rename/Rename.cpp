#include "FileList.h"
#include "FileHandler.h"
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>


int main() {
	std::string prefix = "OnePiece_";
	// take rename_list or directory name from standard input
	std::string linein;
	std::cout<<"Input either rename.list or directory to be analyzed "
				<<"(leave blank for current directory): ";
	std::getline(std::cin, linein);
	std::string pathin = linein.empty() ? "." : linein;
	FileList flst(pathin, prefix);

	// output the rename list (auto analyzed when constructing FileList)
	std::string listname("backup_Rename_list.txt");
	std::ofstream fout(listname);
	if(! fout) {
		std::cerr<<"Error(main): cannot write to file "<<listname<<std::endl;
		exit(-2);
	}
	std::string fname1, fname2;
	while( flst.GetFilename1(fname1) && flst.GetFilename2(fname2) ) {
		if( fname1 != fname2 )
			fout<<fname1<<"\t"<<fname2<<std::endl;
		flst.NextFile();
	}
	fout.close();

	// show results of the analysis and prompt for corrections;
	std::cout<<"Please check the produced rename list ("<<listname<<") and make corrections\n";
	usleep(2e6);
	linein = "vi " + listname;
	system(linein.c_str());

	// prompt to continue;
	std::cout<<"Continue moving (Save any corrections before proceed)?";
	std::getline(std::cin, linein);
	char c1 = linein.empty() ? 'N' : linein[0];
	if( c1 != 'y' && c1 != 'Y' ) return -1;

	// reload the rename list to reflect any corrections made by the user
	FileList flst2(listname);
	std::string newname;
	while( flst2.GetFilename1(fname1) && flst2.GetFilename2(fname2) ) {
		if( fname1 != fname2 ) {
			FileHandler fileh(fname1);
			fileh.Rename(fname2, newname);
			std::cout<<fname1<<" -> "<<newname<<std::endl;
		}
		flst2.NextFile();
	}


	return 0;
}
