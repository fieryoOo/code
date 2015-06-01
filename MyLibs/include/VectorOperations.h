#ifndef VOPERATIONS_H
#define VOPERATIONS_H

#include "StackTrace.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

namespace VO {

	/* exceptions */
	class Base : public std::runtime_error {
		public:
			Base(const std::string message)
				: runtime_error(message) {
					PrintStacktrace();
				}
	};

	class EmptyData : public Base {
		public:
			EmptyData(const std::string funcname, const std::string info = "")
				: Base("Error("+funcname+"): Empty data input ("+info+").") {}
	};

	class BadData : public Base {
		public:
			BadData(const std::string funcname, const std::string info = "")
				: Base("Error("+funcname+"): Bad data ("+info+").") {}
	};

   class BadFile : public Base {
   public:
      BadFile(const std::string funcname, const std::string info = "")
         : Base("Warning("+funcname+"): Unable to access file ("+info+").") {}
   };

   class BadParam : public Base {
   public:
      BadParam(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Bad parameters ("+info+").") {}
   };

   class SizeMismatch : public Base {
   public:
      SizeMismatch(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Incompatible sizes ("+info+").") {}
   };


	/* operations */
	template < class T >
	void Output( const std::vector<T>& dataV, const std::string& fname, const bool app = false ) {
		std::ofstream fout;
		if( app ) fout.open( fname, std::ofstream::app );
		else fout.open( fname );
		if( ! fout )
			throw BadFile(FuncName, "open "+fname);

		if( app )fout<<"\n\n";
		for( const auto& t : dataV )
			fout << t << "\n";
	}

	template < class T >
	bool isSorted( const std::vector<T>& dataV ) {
		int dsize = dataV.size();
		if( dsize > 1 ) {
			for( int i=0; i<dsize-1; i++ )
				if( dataV[i+1] < dataV[i] ) return false;
		}
		return true;
	}

	template < class T >
   void PeriodicExtension( const std::vector<T>& datain, const float binwidth, std::vector<T>& dataout ) {
      // check if datain is sorted 
		if( ! isSorted(datain) )
			throw BadData(FuncName, "not sorted");

      // clear the new vector (of Azidata) dataout
      dataout.clear();
      // prepend measurements with azi>360-binhwidth*2 to the beginning
      T bound; bound.azi = 360-binwidth;
      auto Ibound = std::upper_bound( datain.begin(), datain.end(), bound );
      for(auto iter=Ibound; iter<datain.end(); iter++) {
         auto Ttmp = *iter;
         Ttmp.azi -= 360.;
         dataout.push_back( Ttmp );
      }
      // copy datain in the middle
      size_t presize = dataout.size();
      dataout.resize( presize + datain.size() );
      std::copy( datain.begin(), datain.end(), dataout.begin() + presize );
      // append measurements with azi<binwidth*2 to the end
      bound.azi = binwidth;
      Ibound = std::lower_bound( datain.begin(), datain.end(), bound );
      for(auto iter=datain.begin(); iter<Ibound; iter++) {
         auto Ttmp = *iter;
         Ttmp.azi += 360;
         dataout.push_back( Ttmp );
      }
   }

	// computes the mean and L1 norm of a vector of user-defined class
	// implemented arithmetic operations and .sqrt functions are required for the class
	template < class T >
   bool MeanL1( const typename std::vector<T>::iterator id_lbound, 
					 const typename std::vector<T>::iterator id_ubound,	// data range
					  T& mean, T& L1, const bool stdofmean ) {
		//size_t dsize = data.size();
		size_t dsize = id_ubound - id_lbound;
      if( dsize == 0 ) return false;
      // compute mean
      mean = T{0.};
      for(auto iter=id_lbound; iter<id_ubound; iter++)
         mean += *iter;
      mean /= dsize;
      // stop if only one or two points in data vector
      if( dsize < 3 ) {
         L1 = T{};
         return true;
      }
      // compute L1 norm
      float V2 = 0.;
      L1 = T{0.};
      for(auto iter=id_lbound; iter<id_ubound; iter++)
         L1 += fabs( *iter - mean );
		float denom = dsize-1;
		if( stdofmean ) denom *= sqrt((double)dsize);
      L1 /= denom;
      return true;
	}

	// computes the mean and std of a vector of user-defined class
	// implemented arithmetic operations and .sqrt functions are required for the class
	template < class T >
   bool MeanSTD( const typename std::vector<T>::iterator id_lbound, 
					  const typename std::vector<T>::iterator id_ubound,	// data range
					  const std::vector<float>::iterator iw_lbound, 
					  const std::vector<float>::iterator iw_ubound,			// weight range
					  T& mean, T& std, const bool stdofmean ) {
		//size_t dsize = data.size();
		size_t dsize = id_ubound - id_lbound;
      if( dsize == 0 ) return false;
      if( dsize != iw_ubound - iw_lbound ) 
			throw SizeMismatch(FuncName, "data - weight");
      // compute mean
      float V1 = 0.;
      mean = T{0.};
      for(auto iter=id_lbound; iter<id_ubound; iter++) {
			float weight = *( iw_lbound + (iter-id_lbound) );
         mean += (*iter) * weight;
         V1 += weight;
      }
      mean /= V1;
      // stop if only one or two points in data vector
      if( dsize < 3 ) {
         std = T{};
         return true;
      }
      // compute std
      float V2 = 0.;
      std = T{0.};
      for(auto iter=id_lbound; iter<id_ubound; iter++) {
			float weight = *( iw_lbound + (iter-id_lbound) );
         T Ttmp = (*iter) - mean;
         std += Ttmp * Ttmp * weight;
         V2 += weight * weight;
      }
      float denom = V1*V1-V2;
      //if( frdm <= 0. ) {
      // mean = NaN;
      // return false;
      //}
		if( stdofmean ) denom *= dsize;
      std = sqrt( std * V1 / denom );
      return true;
   }
	template < class T >
   bool MeanSTD( const typename std::vector<T>::iterator id_lbound, 
					  const typename std::vector<T>::iterator id_ubound,		// data range
					  T& mean, T& std, const bool stdofmean ) {
      std::vector<float> weit( id_ubound-id_lbound, 1. );
      return MeanSTD( id_lbound, id_ubound, weit.begin(), weit.end(), mean, std, stdofmean );
   }

   // sub-BinAverage: takes the extended data vector and computes the mean and std in each bin
	// aziavg: output mid (false) or mean (true) azimuth
	// stdofmean: output std-dev (false) or std-dev of the mean (true)
	template < class T >
   void BinAvg( std::vector<T>& datain, std::vector<T>& meanoutV, std::vector<T>& stdoutV,
					 const float binstep, const float binhwidth, const float MIN_BIN_SIZE = 1,
					 const size_t norm_order = 2, const bool aziavg = false, const bool stdofmean = false ) {
      // check if datain is sorted 
		if( ! isSorted(datain) )
			throw BadData(FuncName, "not sorted");
		// check norm order
      // compute average in each bin
      int nbin = (int)ceil(360 / binstep);
      meanoutV.clear(); meanoutV.resize(nbin);
      stdoutV.clear(); stdoutV.resize(nbin);
      auto IT_lbound = datain.begin(), IT_ubound = datain.begin();
		bool (*MeanFunc)(const typename std::vector<T>::iterator, const typename std::vector<T>::iterator, T&, T&, const bool);
		//MeanFunc = norm_order==1 ? &MeanL1 : &MeanSTD;
		if( norm_order == 1 ) {
			MeanFunc = &MeanL1;
		} else if( norm_order == 2 ) {
			MeanFunc = &MeanSTD;
		} else {
			throw BadParam(FuncName, "undefined norm order");
		}
      for(float iazi=0; iazi<nbin; iazi++) {
         float azi = iazi * binstep;
         // search for dataV window boundaries
         T T_lbound, T_ubound;
         T_lbound.azi = azi-binhwidth; T_ubound.azi = azi+binhwidth;
         IT_lbound = std::lower_bound( IT_lbound, datain.end(), T_lbound );
         IT_ubound = std::upper_bound( IT_ubound, datain.end(), T_ubound );
         // compute mean&variance using operator+&* of class T
         // bins with insufficient data are invalidated
         if( IT_ubound - IT_lbound < MIN_BIN_SIZE ) continue;
         MeanFunc( IT_lbound, IT_ubound, meanoutV[iazi], stdoutV[iazi], stdofmean );
			auto& outazi = meanoutV[iazi].azi;
         if( aziavg ) {
            if( outazi >= 360. ) {
               outazi -= 360.;
            } else if( outazi < 0. ) {
               outazi += 360.;
            }
         } else {
            outazi = azi;
         }
      }
   }

	// exclude points based on the input bin-avg data (discard if out of 3-sigma)
	// a method named isWithin( T& lowerbound, T& upperbound ) is required
	template < class T, class T2 >
   bool SelectData( const std::vector<T>& datain, std::vector<T>& dataout, 
						  const std::vector<T2>& binmeanV, const std::vector<T2>& binstdV,
						  const float exfactor ) {
      // check if datain is sorted 
		if( ! ( isSorted(datain) && isSorted(binmeanV) ) )
			throw BadData(FuncName, "not sorted");
		// check sizes
		size_t binsize = binmeanV.size();
		if( binsize != binstdV.size() )
			throw SizeMismatch(FuncName, "binmeanV - binstdV");
		// check address
		if( &datain == &dataout )
			throw BadParam(FuncName, "datain and dataout cannot be identical");
		// search 
      auto IT_lbound = datain.begin(), IT_ubound = datain.begin();
      const float flarge = 0.9 * std::numeric_limits<float>::max();
      // bin ranges provided by binmeanV and binstdV
		dataout.clear();
		for( int i=0; i<binsize; i++ ) {
			auto& binmean = binmeanV[i];
			auto& binstd = binstdV[i];
         // data window (in datain) for the current azi range (from binavgV)
			T2 T2lbound = binstd * exfactor;
			T2 T2ubound = binmean + T2lbound;
			T2lbound = binmean - T2lbound;
         IT_lbound = std::lower_bound( IT_lbound, datain.end(), T2lbound );
         IT_ubound = std::upper_bound( IT_ubound, datain.end(), T2ubound );
			// check all Ts in the range
			for( auto iter=IT_lbound; iter<IT_ubound; iter++ ) {
				// store into dataout if within defined boundaries
				if( (*iter).isWithin( T2lbound, T2ubound) )
					dataout.push_back( (*iter) );
			}
      }
   }

};


#endif
