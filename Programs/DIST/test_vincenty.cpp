#include "DisAzi.h"
#include <iostream>

int main(int argc, char* argv[]) {
	/*
	if (argc!=5) {
      std::cout<<"Usage: "<<argv[0]<<" [lat1] [lon1] [lat2] [lon2]"<<endl;
      return -1;
	}
	*/
	double lon1 = -120., lat1 = 43., dist = 1., azi = 166.;
	Path<double> p1( Point<double>(lon1, lat1), dist, azi );
	std::cout<<p1.P2()<<std::endl;
	std::cout<<p1<<std::endl;;

	Path<double> p2( Point<double>(lon1, lat1), p1.P2() );
	std::cout<<p2.Dist()<<" "<<p2.Azi1()<<std::endl;
	std::cout<<p2<<std::endl;

	// p1.P2() = -119.9970 42.9913
	//std::cout<< Path<double>(Point<double>(60.003, -42.9913), p1.P2()).DistF() <<std::endl;
	//std::cout<< Path<double>(Point<double>(60.003, -42.9913), p1.P2()).Dist() <<std::endl;
	//return -2;
	for(double ilon=40.0; ilon<=80.0; ilon+=0.1) {
		for(double ilat=-63.0; ilat<=-23.0; ilat+=0.1) {
			auto p = Path<double>(Point<double>(ilon, ilat), p1.P2());
			double dis1 = p.DistF();
			double dis2 = p.Dist();
			std::cout<<ilon<<" "<<ilat<<" "<<dis1-dis2<<" "<<dis1<<" "<<dis2<<" "<<p.iter<<std::endl;
		}
	}
	return 0;
}
