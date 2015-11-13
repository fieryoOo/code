#ifndef SEARCHER_H
#define SEARCHER_H

#include "MyOMP.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <random>
#include <algorithm>
#include <functional>

namespace Searcher {
	static float _poc = -1.;	// percentage-of-completion of the current search

	// Interfaces. Required by the searcher!!
	template < class MI >
	class IModelSpace {	
		virtual void Perturb( MI& ) const = 0;
		virtual void SetMState( const MI& ) = 0;
		// also: ModelSpace has to be assignable to MI
	};

	template < class MI >
	class IDataHandler {
		virtual void Energy( const MI&, float& E, int& Ndata ) const = 0;
	};

	// empirical formula for alpha
	inline float Alpha( const int nsearch, const float Tfactor ) {
		if( Tfactor <= 0 ) return 0.;
		return std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided alpha
	}

	// to save search information
	template < class MI >
	struct SearchInfo {
		int isearch;	// search#
		int ithread;	// thread#
		float T;			// current temperature
		MI info;			// model info
		int Ndata;		// #data
		float E;			// energy of current model
		int accepted;	// accepted in search?

		SearchInfo() {}
		SearchInfo( int isearch, float T, MI info, int Ndata, float E, int accepted = 0 )
			: isearch(isearch), T(T), info(info), Ndata(Ndata), E(E), accepted(accepted), ithread(omp_get_thread_num()) {}

		friend bool operator< ( const SearchInfo& s1, const SearchInfo& s2 ) { return s1.isearch<s2.isearch; }
		friend std::ostream& operator<<( std::ostream& o, const SearchInfo<MI>& mi ) {
			o << "search# " << std::setw(3) << mi.isearch << " (T=" << mi.T << ", tid="<<std::setw(2)<<mi.ithread
			  <<"):\tminfo = (" << mi.info << ")\tN = " << mi.Ndata << "\tE = " << mi.E;
			if( mi.accepted == 1 ) o << " (accepted)";
			else if( mi.accepted == 0 ) o << " (rejected)";
			else if( mi.accepted == -1 ) o << " (skipped)";
			else o << " (best)";
			return o;
		}
	};


   // accept (likelihood function)
   static bool Accept(const float E, const float Enew, const double T, const float probability) {
      bool accept;
      if( Enew<E || probability<exp((E-Enew)/T) ) {
         accept = true;
      } else {
         accept = false;
      }
      //std::cerr<<"E="<<E<<" Enew="<<Enew<<"  accept="<<accept<<" "<<ftmp<<" "<<exp((E-Enew)/T)<<std::endl;
      return accept;
   }


	// Simulated Allealing. Three classes required for 1) modelinfo, 2) modelspace, and 3) datahandler
	// a vector of search info is returned upon request
	// Tinit: initial temperature ( controls initial temperature. called when Tinit is real )
	// Tfactor: T_initial = E_initial * Tfactor ( controls initial temperature. called when Tfactor is int )
	// alpha: T_current = T_last * alpha ( controls rate of temperature decay )
	// outSI: -1=no, 0=all, 1=accepted
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch, float alpha, const int Tfactor,
																	  std::ostream& sout = std::cout, short outSI = -1, bool saveV = false,
																	  const int istart = 0 ) {
		if( alpha < 0. )
			alpha = std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided
		return SimulatedAnnealing<MI>( ms, dh, nsearch, alpha, -(float)Tfactor, sout, outSI, saveV, istart );
	}

	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch, const float alpha, const float Tinit,
																	  std::ostream& sout = std::cout, short outSI = -1, bool saveV = false, 
																	  const int istart = 0 ) {
		// initialize random number generator
		_poc = 0.; float pocinc = 1./(nsearch+2);
		std::default_random_engine generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() );
		std::uniform_real_distribution<float> d_uniform(0., 1.);
		//std::normal_distribution<float> d_normal(0., 1.);
		auto rnd = std::bind( d_uniform, generator1 );

		// force MS, and DH to be derived from the provided interfaces at compile time
		const IModelSpace<MI>& ims = ms;
		const IDataHandler<MI>& idh = dh;

		bool outacc = outSI>=0;
		bool outrej = outSI==0;

		// search starts
		if( ! (alpha==1 && Tinit==2.0) )	// not MonteCarlo
			sout<<"### Starting simulated annealing search (#search="<<nsearch<<" alpha="<<alpha<<" Tinit="<<Tinit<<") "<<std::endl;
			sout<<"### with model state =\n"<<ms<<std::endl;

		// initial energy (force MS to be assignable to MI at compile time)
		int Ndata0, Ndata;
		float E; dh.Energy( ms, E, Ndata0 );
		float Ebest = E;
		// initial temperature
		float T = Tinit>0. ? Tinit : std::fabs(E*Tinit);
		// SearchInfo vector
		std::vector< SearchInfo<MI> > VSinfo; 
		VSinfo.reserve( nsearch/2*(outacc+outrej) );
		// save initial
		SearchInfo<MI> sibest( istart, T, ms, Ndata0, E, true );
		if( saveV ) VSinfo.push_back( sibest );
		sout<<sibest<<"\n";
		_poc += pocinc;
		// main loop
		#pragma omp parallel for schedule(dynamic, 1) private(Ndata)
		for( int i=istart; i<istart+nsearch; i++ ) {
			float Enew; int isaccepted = 0;	// 0 for rejected
			MI minew;
			try {
				ms.Perturb( minew );
				dh.Energy( minew, Enew, Ndata );
				Enew *= ((float)Ndata0 / Ndata);
			} catch (const std::exception& e) {
				std::cerr<<e.what()<<"   skipped!\n";
				isaccepted = -1;	// -1 for skipped
			}
			SearchInfo<MI> si( i+1, T, minew, Ndata, Enew, isaccepted );
			#pragma omp critical
			{ // critical begins
			//if( Enew<E || rnd()<exp((E-Enew)/T) ) {
			if( isaccepted==0 && Accept(E, Enew, T, rnd()) ) {
				// update Energy and model state
				ms.SetMState( minew ); E = Enew;
				// mark as accepted
				si.accepted = 1;
				// update (if is) the best
				if( Enew < Ebest ) {
					Ebest = Enew;
					sibest = si;
				}
				if( outacc ) {
					// output search info
					if(saveV) VSinfo.push_back( si );
					sout<<si<<"\n";
				}
			} else if( outrej ) {
				if(saveV) VSinfo.push_back( si );
				sout<<si<<"\n";
			}
			// temperature decrease
			T *= alpha;
			_poc += pocinc;	// update perc-of-completion
			if( i % 100 == 0 ) sout.flush();
			} // critical ends
		}
		// output final result
		sibest.accepted = 2;	// 2 for best fitting model
		sout<<sibest<<std::endl;
		std::sort( VSinfo.begin(), VSinfo.end() );
		// set model state to the best fitting model
		ms.SetMState( sibest.info );
		_poc = -1.;
		return VSinfo;
	}


	// Monte Carlo search: simulated annealing with 1) temperature fixed at 2 and 2) all search info returned
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > MonteCarlo( MS& ms, DH& dh, const int nsearch, std::ostream& sout = std::cout ) {
		std::cout<<"### Monte Carlo search started. (#search = "<<nsearch<<") ###"<<std::endl;
		return SimulatedAnnealing<MI>( ms, dh, nsearch, 1., 2.0f, sout, 0 );		// save all search info
	}

	template < class MI, class MS, class DH >
	void MonteCarlo( MS& ms, DH& dh, const int nsearch, const std::string& outname ) {
		std::cout<<"### Monte Carlo search in process... (#search="<<nsearch<<") ###\n"
					<<"### results being streamed into file "<<outname<<" ###"<<std::endl;
		std::ofstream fout( outname, std::ofstream::app );
		omp_set_nested(true);
		#pragma omp parallel sections
		{	// parallel S
			#pragma omp section
			SimulatedAnnealing<MI>( ms, dh, nsearch, 1., 2.0f, fout, -1 );	// do not save search info
			#pragma omp section
			{	// section S
			std::this_thread::sleep_for( std::chrono::seconds(1) );
			while( _poc >= 0. ) {
				std::cout<<"*** In process... "<<std::setprecision(1)<<std::setw(4)<<std::fixed<<_poc*100<<"\% completed... ***"<<std::endl<<"\x1b[A";
				std::this_thread::sleep_for( std::chrono::seconds(20) );
			}
			std::cout<<"### 100.0\% completed ###\t\t\t\n";
			}	// section E
		}	// parallel E
	}
}

#endif
