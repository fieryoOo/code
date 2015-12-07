/* convert sac files into freq domain and compute the ratio of amplitude spectrums
 * output = sac1.am / sac2.am
 */
#include "SacRec.h"

int main( int argc, char* argv[] ) {
	if( argc!=4 && argc!=6 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [sac_outname (cout=print to screen)] [tmin(optional)] [tmax(optional)]"<<std::endl;
		exit(-1);
	}

	try {
		// load sac
		SacRec sac1( argv[1] ), sac2( argv[2] );
		sac1.Load(); sac2.Load();

		if( argc == 6 ) {
			float tmin = atof(argv[4]), tmax = atof(argv[5]);
			sac1.cut(tmin, tmax); sac2.cut(tmin, tmax);
		}

		// FFT
		SacRec sac1_am, sac2_am;
		sac1.ToAm(sac1_am); sac2.ToAm(sac2_am);
		sac1.clear(); sac2.clear();

		sac1_am.Divf(sac2_am);

		// output
		std::string outname(argv[3]);
		if( outname == "cout" ) {
			sac1_am.Dump();
		} else {
			sac1_am.Write( outname );
		}

	} catch( std::exception& e ) {
		std::cerr<<e.what()<<std::endl;
		return -2;
	} catch(...) {
		std::cerr<<"Error(main): exception detected!"<<std::endl;
		return -3;
	}

	return 0;
}

