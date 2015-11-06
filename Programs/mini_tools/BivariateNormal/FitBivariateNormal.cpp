#include "VectorOperations.h"
#include "Point.h"
#include "DisAzi.h"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

#define FuncName __FUNCTION__

using namespace Eigen;


struct CSData { 
	CSData( const double probin, const double chiSin )
		: prob(probin), chiS(chiSin) {}
	double prob, chiS;

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
			Point<double> pt;
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

	void Boundaries( double& lonmin, double& lonmax, double& latmin, double& latmax ) {
		lonmin = _lonmin; lonmax = _lonmax;
		latmin = _latmin; latmax = _latmax;
	}

	void ErrorEllipse( const double prob, double& x0, double& y0, double& a, double& b, double& theta ) {
		if( _lambda1==NaN || _lambda2==NaN || _theta==NaN )
			ComputeEigen();
		x0 = _meanlon;
		y0 = _meanlat;
		double chi_S;
		if( prob < 0 )	{	// input assumed to be multiple of standard-deviation
			chi_S = prob * prob;
		} else {
			chi_S = chiS( prob );
		}
		a = 2. * sqrt(_lambda1 * chi_S);
		b = 2. * sqrt(_lambda2 * chi_S);
		theta = _theta;
	}

	double PDF( const double lon, const double lat ) {
		Point<double> Pxy = PgeoToPxy(Point<double>(_meanlon,_meanlat), Point<double>(lon,lat));
		double ftmpx = ( Pxy.lon - meanX() ) / stdX();
		double ftmpy = ( Pxy.lat - meanY() ) / stdY();
		return CPDF1() * exp( CPDF2() * (ftmpx*ftmpx + ftmpy*ftmpy - 2*CC()*ftmpx*ftmpy) );
	}

	Point<double> ToPxy( const double lon, const double lat ) {
		return PgeoToPxy(Point<double>(_meanlon,_meanlat), Point<double>(lon,lat));
	}

	double CPDF1() {
		if( _cpdf1 == NaN ) ComputeCPDF();
		return _cpdf1;
	}

	double CPDF2() {
		if( _cpdf2 == NaN ) ComputeCPDF();
		return _cpdf2;
	}

	double CC() {
		if( _cc == NaN ) ComputeCC();
		return _cc;
	}

	double meanX() {
		if( _meanx == NaN ) ComputeMeanSTD();
		return _meanx;
	}

	double meanY() {
		if( _meany == NaN ) ComputeMeanSTD();
		return _meany;
	}

	double stdX() {
		if( _stdx == NaN ) ComputeMeanSTD();
		return _stdx;
	}

	double stdY() {
		if( _stdy == NaN ) ComputeMeanSTD();
		return _stdy;
	}

	size_t size() { return _data.size(); }

protected:
	static constexpr double NaN = -12345.;
	static const std::vector<CSData>& chiSV;

private:
	std::vector< Point<double> > _dataGeo, _data;
	double _meanlon = NaN, _meanlat = NaN;
	double _meanx = NaN, _meany = NaN;
	double _stdx = NaN, _stdy = NaN;
	double _cov = NaN, _cc = NaN;
	double _cpdf1 = NaN, _cpdf2 = NaN;
	double _lambda1 = NaN, _lambda2 = NaN, _theta = NaN;
	double _lonmin = 360., _lonmax = 0.;
	double _latmin = 90., _latmax = -90.;

	double chiS( const double prob ) {
		if( prob<0.00001 || prob>0.99999 )
			throw std::runtime_error( std::string("Error(") + FuncName + "): invalid probability (shoud be 0.00001-0.99999)" );
		// search in chiSV
		auto iup = std::lower_bound( chiSV.begin(), chiSV.end(), CSData(prob, 0.) );
		if( iup == chiSV.begin() ) return iup->chiS;
		if( iup == chiSV.end() ) return (iup-1)->chiS;
		auto ilo = iup - 1;
		// interpolate
		double chiS_prob = ilo->chiS + (iup->chiS - ilo->chiS) * (prob - ilo->prob) / (iup->prob - ilo->prob);
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
		Point<double> Pcenter(_meanlon, _meanlat);
		_data.clear();
		for( const auto& p_geo : _dataGeo ) {
			_data.push_back( PgeoToPxy(Pcenter, p_geo) );
			//std::cerr<<p_geo<<" "<<_data.back()<<"\n";
		}
	}

	void ComputeMeanGeo() {
		Point<double> pmean, pstd;
		if( VO::MeanSTD( _dataGeo.begin(), _dataGeo.end(), pmean, pstd, false ) ) {
			_meanlon = pmean.lon; _meanlat = pmean.lat;
			//_stdlon = pstd.lon; _stdlat = pstd.lat;
		}
	}

	void ComputeMeanSTD() {
		Point<double> pmean, pstd;
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
		double varX = stdX() * stdX();
		double varY = stdY() * stdY();
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
		_lambda1 = eigvals(0).real();	// variance #1
		_lambda2 = eigvals(1).real();	// variance #2
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
		//double A = 1. / (stdX()*stdX());
		//double B = -2. * CC() / (stdX()*stdY());
		//double C = 1. / (stdY()*stdY());
		//_theta = 0.5 * atan2(B, A-C) * 180./M_PI;
	}

};

// the degree two chi-square probability
static const std::vector<CSData> chiSV_data = 
	{ CSData(0.00001, 0.000020), CSData(0.10000, 0.210721), CSData(0.20000, 0.446287),
	  CSData(0.30000, 0.713350), CSData(0.40000, 1.021651), CSData(0.50000, 1.386294),
	  CSData(0.60000, 1.832581), CSData(0.70000, 2.407946), CSData(0.80000, 3.218876),
	  CSData(0.89000, 4.414550), CSData(0.95000, 5.991465), CSData(0.98000, 7.824046),
	  CSData(0.99500, 10.59663), CSData(0.99950, 15.20180), CSData(0.99999, 23.02585) };
const std::vector<CSData>& DataHandler::chiSV = chiSV_data;

int main( int argc, char *argv[] ) {
	if( argc != 3 && argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [infile(>=2 columns)] [prob confidence(0.-1.) or #sigma(-int)] [prob outfile (optional)]"<<std::endl;
		return -1;
	}

	DataHandler dh( argv[1] );

	//std::string fname( argv[1] ); fname += "_XY";
	//dh.OutputXY( fname );
	//std::cerr<<dh.size()<<" data points loaded"<<std::endl;
	//std::cerr<<"correlation_coef = "<<dh.CC()<<std::endl;

	double lon0, lat0, a, b, theta;
	dh.ErrorEllipse( atof(argv[2]), lon0, lat0, a, b, theta );
	std::cout<<std::setprecision(9);
	std::cout<<lon0<<" "<<lat0<<" "<<theta<<" "<<a<<" "<<b<<std::endl;

	if( argc == 3 ) return 0;

	// compute pdfs if requested
	// define output boundaries
	double lonmin, lonmax, latmin, latmax;
	dh.Boundaries( lonmin, lonmax, latmin, latmax );
	double lonspan = lonmax - lonmin, latspan = latmax - latmin;
	double mperc = 0.3;
	lonmin -= lonspan*mperc; lonmax += lonspan*mperc;
	latmin -= latspan*mperc; latmax += latspan*mperc;
	double step=0.005*(lonmax-lonmin), Agrid = step*step;
	// compute
	std::ofstream fout( argv[3] );
	//double cdf = 0.;
	fout<<std::setprecision(9);
	for( double lon=lonmin; lon<=lonmax; lon+=step ) {
		for( double lat=latmin; lat<latmax; lat+=step ) {
			double pdf = dh.PDF(lon,lat);
			//cdf += pdf * Agrid;
			fout<<lon<<" "<<lat<<" "<<pdf<<" "<<dh.ToPxy(lon, lat)<<"\n";
		}
		fout<<"\n";
	}
	//std::cerr<<" cdf = "<<cdf<<std::endl;

	return 0;
}
