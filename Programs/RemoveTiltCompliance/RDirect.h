#ifndef RDIRECT_H
#define RDIRECT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

class RDirect {
public:
	typedef std::pair<float, float> Range;
	class rangeCmp {
	public:
		bool operator()(const Range& r1, const Range& r2) const { 
			return (r1.first+r1.second)<(r2.first+r2.second); 
		}
	};
	typedef std::map< float, std::vector<float> > Donut;	// Donut: azi-cohV pair
	typedef std::map< Range, Donut, rangeCmp > TDMap;		// TDMap: twin-Donut pair

public:
	RDirect( const std::vector<Range>& frangeV = {} ) : _frangeV(frangeV) {}
	
	void Write(const std::string& oname) const;

	void Load(const std::string& iname);

	size_t size() const { return _data.size(); }

	size_t sizedonut() const {
		if( _data.empty() ) return 0;
		else return _data.begin()->second.size();
	}

	typename TDMap::iterator		 begin()			{ return _data.begin(); }
	typename TDMap::const_iterator begin() const { return _data.begin(); }
	typename TDMap::const_iterator cbegin() const { return _data.begin(); }
	typename TDMap::iterator		 end()		 { return _data.end(); }
	typename TDMap::const_iterator end() const { return _data.end(); }
	typename TDMap::const_iterator cend() const { return _data.end(); }

	Donut& operator[] (const Range& r) { return _data[r]; }

	TDMap::const_iterator find( const Range& r ) const { return _data.find(r); }

	const std::vector<Range>& fRanges() const { return _frangeV; }

	float aziWeight( float dazi, float &wsum ) const;
	std::vector<float> AvgsAlong( const Range &twin, const float aziTarget ) const;

protected:
	static constexpr float PIoDEG = M_PI / 180.;

private:
	TDMap _data;
	std::vector<Range> _frangeV;
};


void RDirect::Write(const std::string& oname) const {
	// open file
	std::ofstream fout(oname);
	if( ! fout ) throw std::runtime_error("RDirect::Write: IO failed");
	// output freq ranges
	for( const auto& fwin : _frangeV ) fout<<fwin.first<<" "<<fwin.second<<"   ";
	fout<<"\n";
	// output a row of coherences (different frequency) for each (twin, azi)
	for( const auto& tdPair : _data ) {	// by time
		Range twin = tdPair.first;
		for( const auto& acPair : tdPair.second ) {	// by azimuth
			float azi = acPair.first;
			fout<<twin.first<<" "<<twin.second<<" "<<azi<<" ";
			// by freq
			for( const auto& val : acPair.second ) fout<<val<<" ";
			fout<<"\n";
		}
	}
}

void RDirect::Load(const std::string& iname) {
	std::ifstream fin(iname);
	if( ! fin ) throw std::runtime_error("RDirect::Load: IO failed");
	// load in freq ranges from the first line
	std::string line; std::getline(fin, line);
	std::stringstream ss(line);
	_frangeV.clear();
	for(float f1,f2; ss>>f1>>f2; ) _frangeV.push_back(std::make_pair(f1, f2));
	// load coherences (each line is a single (twin, azi) )
	_data.clear();
	auto lsize = _frangeV.size();
	while( std::getline(fin, line) ) {
		std::stringstream ss(line);
		float t1, t2, azi;
		if( ! (ss >> t1 >> t2 >> azi) )
			throw std::runtime_error("RDirect::Load: format error in line "+line);
		if( azi < 0. ) azi += 360.;
		else if( azi >= 360.) azi -= 360.;
		auto &cohV = _data[std::make_pair(t1,t2)][azi];
		for(float coh; ss>>coh; ) cohV.push_back(coh);
		if( cohV.size() != lsize )
			throw std::runtime_error("RDirect::Load: num of freq points does not match");
	}
}

float RDirect::aziWeight( float dazi, float &wsum ) const {
	// weighting function 1: cos based, weight(azi=0) = 1 and weight(azi=+/-pi/2) = -0.5
	/*
	float weight = cos( 2.*PIoDEG * (dazi) );
	if( weight < 0. ) weight *= 0.5;
	*/
	// weighting function 2: gaussian based
	// gaussian coefs
	static const float h1 = 10., alpha1 = - 0.5 / (h1*h1), bound1 = 30.;
	static const float h2 = 25., alpha2 = - 0.5 / (h2*h2), bound2 = 120.;
	//static const float h1 = 30., alpha1 = - 0.5 / (h1*h1), bound1 = 90.;
	//static const float h2 = 20., alpha2 = - 0.5 / (h2*h2), bound2 = 90.;
	dazi = fabs(dazi);
	if( dazi > 180. ) dazi = 360. - dazi;
	float weight = 0.;
	if( dazi < bound1 ) {	// close to target direction:
		weight = exp(alpha1*dazi*dazi);	// positive weight
		wsum += weight;
	} else if( dazi < bound2 ) {	// close to perpendicular:
		//weight = -0.0 * exp(alpha2*(dazi-90.)*(dazi-90.));	// negative weight
		weight = -0.2 * exp(alpha2*(dazi-90.)*(dazi-90.));	// negative weight
	}
	return weight;
}

// get averages (for all freqs) at a given twin&azi
std::vector<float> RDirect::AvgsAlong( const Range &twin, const float aziTarget ) const {
	const auto &donut = _data.at(twin);
	if( donut.size() < 5 ) 
		throw std::runtime_error("Error(RDirect::AvgsAlong): incomplete donut with size(No. of Azi) = "
																								+std::to_string(donut.size()));
	// for each freq range, compute a weighted average of coh over all azimuths
	std::vector<float> SWeightV( donut.begin()->second.size() );
	//std::vector<float> cohmaxV( donut.begin()->second.size() );
	float wsum = 0.;
	for( const auto &acPair : donut ) {
		float dazi = aziTarget-acPair.first;
		float weight = aziWeight( dazi, wsum );
		if( weight <= 0. ) continue;
		const auto &cohV = acPair.second;
		for(int i=0; i<cohV.size(); i++) {
			SWeightV[i] += cohV[i] * weight;
			/*
			if( cohV[i] < cohmaxV[i] ) continue;
			cohmaxV[i] = cohV[i];
			//SWeightV[i] = cohV[i] * weight;
			SWeightV[i] = cohV[i] * cos(dazi*M_PI/180.);
			*/
			//SWeightV[i] = std::max(SWeightV[i], cohV[i]*(float)cos(dazi*M_PI/180.));
		}
	}
	//for( auto& sweight : SWeightV ) sweight = std::max(0.1f, sweight/wsum);
	for( auto& sweight : SWeightV ) sweight /= wsum;
	return SWeightV;
}

#endif
