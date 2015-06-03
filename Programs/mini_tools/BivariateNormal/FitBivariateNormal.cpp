#include "VectorOperations.h"
#include "Point.h"
#include "DisAzi.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

#define FuncName __FUNCTION__

using namespace Eigen;


struct CSData { 
	CSData( const float probin, const float chiSin )
		: prob(probin), chiS(chiSin) {}
	float prob, chiS;

	friend bool operator< (const CSData& cd1, const CSData& cd2) { return cd1.prob<cd2.prob; }
};
class DataHandler {
public:
	DataHandler( const std::string fname = "" ) {
		if( ! fname.empty() ) Load(fname);
	}

	void Load( const std::string& fname ) {
		std::ifstream fin(fname);
		if( !fin ) throw std::runtime_error( std::string("Error(") + FuncName + "): cannot read from file " + fname );
		for( std::string line; std::getline(fin, line); ) {
			Point<float> pt;
			if( ! pt.LoadLine(line) ) continue;
			if( pt.lon<0. ) pt.lon += 360.;
			_dataGeo.push_back(pt);
			if( _lonmin > pt.lon ) _lonmin = pt.lon;
			if( _lonmax < pt.lon ) _lonmax = pt.lon;
			if( _latmin > pt.lat ) _latmin = pt.lat;
			if( _latmax < pt.lat ) _latmax = pt.lat;
		}
		// convert longitude/latitude to km
		GeoToXY();
	}

	void OutputXY( const std::string fname ) {
		if( _data.size() != _dataGeo.size() )
			throw std::runtime_error( std::string("Error(") + FuncName + "): Geo-XY data size mismatch " );
			
		std::ofstream fout(fname);
		for(int i=0; i<_data.size(); i++) {
			fout<<_data[i]<<" "<<_dataGeo[i]<<"\n";
		}
	}

	void Boundaries( float& lonmin, float& lonmax, float& latmin, float& latmax ) {
		lonmin = _lonmin; lonmax = _lonmax;
		latmin = _latmin; latmax = _latmax;
	}

	void ErrorEllipse( const float prob, float& x0, float& y0, float& a, float& b, float& theta ) {
		if( _lambda1==NaN || _lambda2==NaN || _theta==NaN )
			ComputeEigen();
		x0 = _meanlon;
		y0 = _meanlat;
		float chi_S = chiS( prob );
		a = 2. * sqrt(_lambda1 * chi_S);
		b = 2. * sqrt(_lambda2 * chi_S);
		theta = _theta;
	}

	float PDF( const float lon, const float lat ) {
		Point<float> Pxy = PgeoToPxy(Point<float>(_meanlon,_meanlat), Point<float>(lon,lat));
		float ftmpx = ( Pxy.lon - meanX() ) / stdX();
		float ftmpy = ( Pxy.lat - meanY() ) / stdY();
		return CPDF1() * exp( CPDF2() * (ftmpx*ftmpx + ftmpy*ftmpy - 2*CC()*ftmpx*ftmpy) );
	}

	Point<float> ToPxy( const float lon, const float lat ) {
		return PgeoToPxy(Point<float>(_meanlon,_meanlat), Point<float>(lon,lat));
	}

	float CPDF1() {
		if( _cpdf1 == NaN ) ComputeCPDF();
		return _cpdf1;
	}

	float CPDF2() {
		if( _cpdf2 == NaN ) ComputeCPDF();
		return _cpdf2;
	}

	float CC() {
		if( _cc == NaN ) ComputeCC();
		return _cc;
	}

	float meanX() {
		if( _meanx == NaN ) ComputeMeanSTD();
		return _meanx;
	}

	float meanY() {
		if( _meany == NaN ) ComputeMeanSTD();
		return _meany;
	}

	float stdX() {
		if( _stdx == NaN ) ComputeMeanSTD();
		return _stdx;
	}

	float stdY() {
		if( _stdy == NaN ) ComputeMeanSTD();
		return _stdy;
	}

	size_t size() { return _data.size(); }

protected:
	static constexpr float NaN = -12345.;
	static const std::vector<CSData>& chiSV;

private:
	std::vector< Point<float> > _dataGeo, _data;
	float _meanlon = NaN, _meanlat = NaN;
	float _meanx = NaN, _meany = NaN;
	float _stdx = NaN, _stdy = NaN;
	float _cov = NaN, _cc = NaN;
	float _cpdf1 = NaN, _cpdf2 = NaN;
	float _lambda1 = NaN, _lambda2 = NaN, _theta = NaN;
	float _lonmin = 360., _lonmax = 0.;
	float _latmin = 90., _latmax = -90.;

	float chiS( const float prob ) {
		if( prob<0.00001 || prob>0.99999 )
			throw std::runtime_error( std::string("Error(") + FuncName + "): invalid probability (shoud be 0.00001-0.99999)" );
		// search in chiSV
		auto iup = std::lower_bound( chiSV.begin(), chiSV.end(), CSData(prob, 0.) );
		if( iup == chiSV.begin() ) return iup->chiS;
		if( iup == chiSV.end() ) return (iup-1)->chiS;
		auto ilo = iup - 1;
		// interpolate
		float chiS_prob = ilo->chiS + (iup->chiS - ilo->chiS) * (prob - ilo->prob) / (iup->prob - ilo->prob);
		return chiS_prob;
	}

	template <class T>
	Point<T> PgeoToPxy( const Point<T> Pcenter, const Point<T> Pgeo ) {
		auto path = Path<T>(Pcenter, Pgeo);
		T dis = path.Dist(), azi = path.Azi1()*M_PI/180.;
		return Point<T>( dis*sin(azi), dis*cos(azi) );
	}
	void GeoToXY() {
		ComputeMeanGeo();
		Point<float> Pcenter(_meanlon, _meanlat);
		_data.clear();
		for( const auto& p_geo : _dataGeo ) {
			_data.push_back( PgeoToPxy(Pcenter, p_geo) );
			//std::cerr<<p_geo<<" "<<_data.back()<<"\n";
		}
	}

	void ComputeMeanGeo() {
		Point<float> pmean, pstd;
		if( VO::MeanSTD( _dataGeo.begin(), _dataGeo.end(), pmean, pstd, false ) ) {
			_meanlon = pmean.lon; _meanlat = pmean.lat;
			//_stdlon = pstd.lon; _stdlat = pstd.lat;
		}
	}

	void ComputeMeanSTD() {
		Point<float> pmean, pstd;
		if( VO::MeanSTD( _data.begin(), _data.end(), pmean, pstd, false ) ) {
			_meanx = pmean.lon; _meany = pmean.lat;
			_stdx = pstd.lon; _stdy = pstd.lat;
		}
	}

	void ComputeCov() {
		_cov = 0.;
		for( const auto& p : _data )
			_cov += (p.lon-meanX()) * (p.lat-meanY());
		_cov /= (_data.size() - 1);
	}
	void ComputeCC() {
		if( _cov == NaN ) ComputeCov();
		_cc = _cov / ( stdX() * stdY() );
	}

	void ComputeCPDF() {
		_cpdf2 = - 0.5 / ( sqrt(1. - CC()*CC()) );
		_cpdf1 = 0.5 / ( M_PI * stdX() * stdY() * sqrt(1. - CC()*CC()) );
	}

	void ComputeEigen() {
		// covariance matrix
		//size_t Nptsmo = _data.size() - 1;
		float varX = stdX() * stdX();
		float varY = stdY() * stdY();
		if( _cov == NaN ) ComputeCov();
		MatrixXd covM(2,2);
		covM << varX, _cov,
				 _cov, varY;
		// eigen values
		EigenSolver<MatrixXd> es(covM);
		VectorXcd eigvals = es.eigenvalues();
		if( eigvals(0).imag()!=0. || eigvals(1).imag()!=0. ) {
			std::cerr<<"Exception: complex eigen value(s)!"<<std::endl;
			exit(-3);
		}
		_lambda1 = eigvals(0).real();
		_lambda2 = eigvals(1).real();
		// eigen vectors
		MatrixXd eigvecs = es.eigenvectors().real();
		VectorXd eigvec(2);
		if( _lambda1 < _lambda2 ) {
			std::swap(_lambda1, _lambda2);
			eigvec = eigvecs.col(1);
		} else {
			eigvec = eigvecs.col(0);
		}
		// rotation angle
		_theta = atan2(eigvec(1), eigvec(0)) * 180./M_PI;
		//float A = 1. / (stdX()*stdX());
		//float B = -2. * CC() / (stdX()*stdY());
		//float C = 1. / (stdY()*stdY());
		//_theta = 0.5 * atan2(B, A-C) * 180./M_PI;
	}

};

static const std::vector<CSData> chiSV_data = 
	{ CSData(0.00001, 0.000020), CSData(0.10000, 0.210721), CSData(0.20000, 0.446287),
	  CSData(0.30000, 0.713350), CSData(0.40000, 1.021651), CSData(0.50000, 1.386294),
	  CSData(0.60000, 1.832581), CSData(0.70000, 2.407946), CSData(0.80000, 3.218876),
	  CSData(0.89000, 4.414550), CSData(0.95000, 5.991465), CSData(0.98000, 7.824046),
	  CSData(0.99500, 10.59663), CSData(0.99950, 15.20180), CSData(0.99999, 23.02585) };
const std::vector<CSData>& DataHandler::chiSV = chiSV_data;

int main( int argc, char *argv[] ) {
	if( argc != 3 && argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [infile(>=2 columns)] [prob confidence] [prob outfile (optional)]"<<std::endl;
		return -1;
	}

	DataHandler dh( argv[1] );

	//std::string fname( argv[1] ); fname += "_XY";
	//dh.OutputXY( fname );
	//std::cerr<<dh.size()<<" data points loaded"<<std::endl;
	//std::cerr<<"correlation_coef = "<<dh.CC()<<std::endl;

	float lon0, lat0, a, b, theta;
	dh.ErrorEllipse( atof(argv[2]), lon0, lat0, a, b, theta );
	std::cout<<lon0<<" "<<lat0<<" "<<theta<<" "<<a<<" "<<b<<std::endl;

	if( argc == 3 ) return 0;

	// compute pdfs if requested
	// define output boundaries
	float lonmin, lonmax, latmin, latmax;
	dh.Boundaries( lonmin, lonmax, latmin, latmax );
	float lonspan = lonmax - lonmin, latspan = latmax - latmin;
	float mperc = 0.3;
	lonmin -= lonspan*mperc; lonmax += lonspan*mperc;
	latmin -= latspan*mperc; latmax += latspan*mperc;
	float step=0.005*(lonmax-lonmin), Agrid = step*step;
	// compute
	std::ofstream fout( argv[3] );
	//float cdf = 0.;
	for( float lon=lonmin; lon<=lonmax; lon+=step ) {
		for( float lat=latmin; lat<latmax; lat+=step ) {
			float pdf = dh.PDF(lon,lat);
			//cdf += pdf * Agrid;
			fout<<lon<<" "<<lat<<" "<<pdf<<" "<<dh.ToPxy(lon, lat)<<"\n";
		}
		fout<<"\n";
	}
	//std::cerr<<" cdf = "<<cdf<<std::endl;

	return 0;
}
